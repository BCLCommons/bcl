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

// include header for this class
#include "scorestat/bcl_scorestat_phipsi.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PhiPsi::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new PhiPsi())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &PhiPsi::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram2D",
        GetStaticClassName< PhiPsi::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &PhiPsi::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "phi_psi_angles_by_sstype.histogram2D",
        GetStaticClassName< PhiPsi::OutputOption>()
      };
      return s_output_filename_extensions[ size_t(OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsi::PhiPsi() :
      m_OutputOption( e_Histogram2D),
      m_NumberOfBins( 12),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    PhiPsi *PhiPsi::Clone() const
    {
      return new PhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &PhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &PhiPsi::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the number of bins of Phi/Psi dihedral historgams
    //! @return the number of bins of Phi/Psi dihedral histograms
    const size_t &PhiPsi::GetNumberOfBins() const
    {
      return m_NumberOfBins;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &PhiPsi::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &PhiPsi::GetAlias() const
    {
      static std::string s_name( "PhiPsi");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string PhiPsi::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics of phi/psi dihedrals
      storage::Vector< storage::Vector< storage::Vector< math::Histogram2D> > > phi_psi_sstypes_histograms
      (
        2,
        storage::Vector< storage::Vector< math::Histogram2D> >
        (
          biol::GetSSTypes().GetEnumCount(),
          storage::Vector< math::Histogram2D>
          (
            biol::GetAATypes().GetEnumCount(),
            math::Histogram2D
            (
              storage::VectorND< 2, double>( -math::g_Pi, -math::g_Pi),
              storage::VectorND< 2, double>( 2 * math::g_Pi / m_NumberOfBins, 2 * math::g_Pi / m_NumberOfBins),
              storage::VectorND< 2, size_t>( m_NumberOfBins, m_NumberOfBins)
            )
          )
        )
      );

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_model( *protein_model_itr);
        sspred::CIPhiPsi().Calculate( *sp_model, true);

        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains in current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageStd( "Skip undesired chain: " + util::Format()( ( *chain_itr)->GetChainID()));
            continue;
          }

          // get all sses in current chain
          const util::SiPtrVector< const assemble::SSE> &all_sses( ( *chain_itr)->GetSSEs());

          // iterate over all sse in current chain
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip undefined sses
            if( !( *sse_itr)->GetType().IsDefined())
            {
              BCL_MessageStd( "Skip undefined sse: " + ( *sse_itr)->GetType().GetName());
              continue;
            }

            // get all amino acids in current sse
            const util::ShPtrVector< biol::AABase> &all_amino_acids( ( *sse_itr)->GetData());

            // iterate over all amino acid in current sse
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                aa_itr( all_amino_acids.Begin()), aa_itr_end( all_amino_acids.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // skip unnatural amino acids or those that are first or last in the sse
              if
              (
                !( *aa_itr)->GetType()->IsNaturalAminoAcid() ||
                aa_itr == all_amino_acids.Begin() ||
                aa_itr == aa_itr_end - 1
              )
              {
                continue;
              }

              // calculate phi/psi using C atom from preceding residue and N atom from following residue
              storage::VectorND< 2, double> phi_psi_vector
              (
                ( *aa_itr)->CalculatePhi( ( *( aa_itr - 1))->GetAtom( biol::GetAtomTypes().C)),
                ( *aa_itr)->CalculatePsi( ( *( aa_itr + 1))->GetAtom( biol::GetAtomTypes().N))
              );

              if( m_OutputOption == e_Histogram2D)
              {
                auto sspre( ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi));
                biol::EnvironmentType env_type
                (
                  sp_membrane.IsDefined() && sspre.IsDefined()
                  ? sspre->GetOneStateTMPrediction()->GetReducedType()
                  : biol::GetEnvironmentTypes().e_Solution
                );
                phi_psi_sstypes_histograms( env_type == biol::GetEnvironmentTypes().e_MembraneCore ? 1 : 0)
                ( ( *sse_itr)->GetType())
                ( ( *aa_itr)->GetType()).PushBack( phi_psi_vector);
              }
            } // end of iteration over amino acid residues
          } // end of iteration over sses
        } // end of iteration over chains
      } // end of iteration over protein models

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_Histogram2D)
      {
        // iterate over sse types
        for( size_t in_membrane( 0), im_max( 2); in_membrane < im_max; ++in_membrane)
        {
          stream << ( in_membrane ? "MEMBRANE" : "SOLUTION") << '\n';
          for
          (
            biol::SSTypes::const_iterator
              sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
            sse_type_itr != sse_type_itr_end;
            ++sse_type_itr
          )
          {
            // write sse type
            stream << *sse_type_itr << '\n';

            // iterate over aa types
            for
            (
              biol::AATypes::const_iterator
                aa_itr( biol::GetAATypes().Begin()),
                aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // write one letter code of aa
              stream << ( *aa_itr)->GetOneLetterCode() << '\n';
              // write histograms
              stream << phi_psi_sstypes_histograms( in_membrane)( *sse_type_itr)( *aa_itr);
            }
          }
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PhiPsi::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes phi/psi dihedral statistics."
      );

      parameters.AddInitializer
      (
        "resolution",
        "number of bins of the phi/psi histogram",
        io::Serialization::GetAgent( &m_NumberOfBins),
        "12"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram2D"
      );

      return parameters;
    }

  } // namespace scorestat
} // namespace bcl

