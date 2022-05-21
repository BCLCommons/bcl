// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
// include header for this class
#include "scorestat/bcl_scorestat_sspred_agreement.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_roc_curve.h"
#include "score/bcl_score_environment_predictions.h"
#include "score/bcl_score_sse_predictions.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_methods.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSPredAgreement::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSPredAgreement())
    );

    //! @brief default constructor
    SSPredAgreement::SSPredAgreement() :
      m_Method( sspred::GetMethods().e_PSIPRED),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSPredAgreement *SSPredAgreement::Clone() const
    {
      return new SSPredAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSPredAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSPredAgreement::GetOutFilePostfix() const
    {
      static const std::string s_suffix( "sspred_agreement.histograms");
      return s_suffix;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSPredAgreement::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSPredAgreement::GetAlias() const
    {
      static std::string s_name( "SSPredAgreement");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSPredAgreement::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // declare and initialize containers for holding sspred agreement statistics
      // histograms for secondary structure confidences
      // stores for each sspred method, each native sstype, the probability histograms for all sstypes
      storage::Map< sspred::Method, storage::Map< biol::SSType, storage::List< storage::Pair< double, double> > > >
        sspred_native_roc;
      storage::Map< sspred::Method, storage::Map< biol::EnvironmentType, storage::List< storage::Pair< double, double> > > >
        env_native_roc;

      const storage::Vector< biol::EnvironmentType> &env_types( biol::GetEnvironmentTypes().GetReducedTypes());
      const storage::Vector< biol::SSType> &ss_types( biol::GetSSTypes().GetReducedTypes());

      // initialize these containers
      for
      (
        auto method_itr( m_Method.Begin()), method_itr_end( m_Method.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        auto &map_o_lists( sspred_native_roc[ *method_itr]);
        auto &map_o_env_lists( env_native_roc[ *method_itr]);

        // for all sstypes, put in seed values so that even if we don't see this state we can just assume that the PPV is =
        // the predictions
        for
        (
          auto ss_type_itr_a( ss_types.Begin()), ss_type_itr_a_end( ss_types.End());
          ss_type_itr_a != ss_type_itr_a_end;
          ++ss_type_itr_a
        )
        {
          map_o_lists[ *ss_type_itr_a].PushBack( storage::Pair< double, double>( 0.0, 0.0));
          map_o_lists[ *ss_type_itr_a].PushBack( storage::Pair< double, double>( 1.0, 1.0));
        }

        for
        (
          auto env_type_itr( env_types.Begin()), env_type_itr_end( env_types.End());
          env_type_itr != env_type_itr_end;
          ++env_type_itr
        )
        {
          map_o_env_lists[ *env_type_itr].PushBack( storage::Pair< double, double>( 0.0, 0.0));
          map_o_env_lists[ *env_type_itr].PushBack( storage::Pair< double, double>( 1.0, 1.0));
        }
      }

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get a modifiable copy of current protein model
        assemble::ProteinModel &current_protein_model( *( *protein_model_itr)->HardCopy());
        sspred::PDB::SetEnvironmentTypes( current_protein_model, true);
        sspred::CIPhiPsi().Calculate( current_protein_model, true);

        // get the membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get model name
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        const std::string final_pdb_path( io::File::SplitToPathAndFileName( model_name).First());
        const std::string pdb_id( io::File::RemoveLastExtension( io::File::RemovePath( model_name)));

        // get all chains in current protein model
        util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageStd( "Skip chain " + util::Format()( ( *chain_itr)->GetChainID()) + " in " + model_name);
            continue;
          }

          // create reference on current chain
          biol::AASequence &current_chain( *( *chain_itr)->GetSequence());

          // determine sspred methods for which prediction files are available
          const storage::Set< sspred::Method> available_methods
            ( sspred::MethodHandler::AvailablePredictionFiles( m_Method, ( *chain_itr)->GetChainID(), pdb_id, final_pdb_path));

          // print out message if there are prediction files missing
          if( m_Method.GetSize() != available_methods.GetSize())
          {
            BCL_MessageCrt
            (
              "not all requested SSPRED Methods have prediction files available for pdb and chain " +
              model_name + " " + ( *chain_itr)->GetChainID()
            );
          }

          // try to read in all available sspred files
          if( !sspred::MethodHandler::ReadPredictionsForAASequence( available_methods, current_chain, pdb_id, final_pdb_path))
          {
            BCL_MessageCrt
            (
              "can't read in prediction file for all requested SSPRED Methods for pdb and chain " +
              model_name + " " + ( *chain_itr)->GetChainID()
            );
          }
          else // if all files were successfully read in
          {
            // iterate over the amino acid residues
            for
            (
              biol::AASequence::const_iterator
                aa_itr( current_chain.Begin()), aa_itr_end( current_chain.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              auto true_ss_env( ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi));
              if( !true_ss_env.IsDefined())
              {
                true_ss_env = ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_PDB);
              }
              if( !true_ss_env.IsDefined())
              {
                continue;
              }
              const biol::SSType native_ss_type( true_ss_env->GetOneStateSSPrediction());
              const biol::EnvironmentType native_env_type( true_ss_env->GetOneStateTMPrediction()->GetReducedType());
              for
              (
                auto method_itr( available_methods.Begin()), method_itr_end( available_methods.End());
                method_itr != method_itr_end;
                ++method_itr
              )
              {
                const util::SiPtr< const sspred::MethodInterface> &preds( ( *aa_itr)->GetSSPrediction( *method_itr));
                if( !preds.IsDefined())
                {
                  continue;
                }
                if( native_ss_type.IsDefined())
                {
                  auto sspreds( preds->GetThreeStatePrediction());
                  auto itr_sspreds( sspreds.Begin());
                  for( auto ss_itr( ss_types.Begin()), ss_itr_end( ss_types.End()); ss_itr != ss_itr_end; ++ss_itr, ++itr_sspreds)
                  {
                    sspred_native_roc[ *method_itr][ *ss_itr].PushBack( storage::Pair< double, double>( *itr_sspreds, *ss_itr == native_ss_type ? 1.0 : 0.0));
                  }
                }
                if( native_env_type.IsDefined())
                {
                  auto envs( preds->GetThreeStateTMPrediction());
                  auto itr_env( envs.Begin());
                  for
                  (
                    auto env_itr( env_types.Begin()), env_itr_end( env_types.End());
                    env_itr != env_itr_end;
                    ++env_itr, ++itr_env
                  )
                  {
                    env_native_roc[ *method_itr][ *env_itr].PushBack( storage::Pair< double, double>( *itr_env, *env_itr == native_env_type ? 1.0 : 0.0));
                  }
                }
              }
            }
          } // end of iteration over all chains
        } // end of iteration over current protein model
      } // end of iteration over protein ensemble

      storage::Map< sspred::Method, storage::Map< biol::SSType, math::PiecewiseFunction> >
        sspred_native_localppv;
      storage::Map< sspred::Method, storage::Map< biol::EnvironmentType, math::PiecewiseFunction> >
        env_native_localppv;

      // initialize these containers
      for
      (
        auto method_itr( m_Method.Begin()), method_itr_end( m_Method.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        auto &map_o_ss_lists( sspred_native_roc[ *method_itr]);
        auto &map_o_env_lists( env_native_roc[ *method_itr]);

        for( auto ss_itr( ss_types.Begin()), ss_itr_end( ss_types.End()); ss_itr != ss_itr_end; ++ss_itr)
        {
          sspred_native_localppv[ *method_itr][ *ss_itr] = math::ROCCurve( map_o_ss_lists[ *ss_itr], 0.5, true).GetLocalPPVCurve();
        }
        for( auto env_itr( env_types.Begin()), env_itr_end( env_types.End()); env_itr != env_itr_end; ++env_itr)
        {
          env_native_localppv[ *method_itr][ *env_itr] = math::ROCCurve( map_o_env_lists[ *env_itr], 0.5, true).GetLocalPPVCurve();
        }
      }

      // write statistics
      std::stringstream stream;

      // output file stream
      io::OFStream write;

      // write mean and sd over n residues
      io::File::MustOpenOFStream( write, score::SSEPredictions::GetDefaultHistogramFilename());
      write << sspred_native_localppv;
      io::File::CloseClearFStream( write);
      io::File::MustOpenOFStream( write, score::EnvironmentPredictions::GetDefaultHistogramFilename());
      write << env_native_localppv;
      io::File::CloseClearFStream( write);

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSPredAgreement::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics for secondary structure prediction agreement."
      );

      parameters.AddInitializer
      (
        "sspred_methods",
        "secondary structure prediction methods to be analyzed",
        io::Serialization::GetAgent( &m_Method),
        "(PSIPRED)"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      return parameters;
    }

  } // namespace scorestat
} // namespace bcl

