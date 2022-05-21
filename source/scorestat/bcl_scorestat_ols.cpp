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
#include "scorestat/bcl_scorestat_ols.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_histogram.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> OLS::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new OLS( false))
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> OLS::s_InstanceEnv
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new OLS( true))
    );

    //! @brief ChainOption as string
    //! @param CHAIN_OPTION the ChainOption
    //! @return the string for the ChainOption
    const std::string &OLS::GetChainOptionName( const ChainOption &CHAIN_OPTION)
    {
      static std::string s_names[] =
      {
        "One",
        "All",
        GetStaticClassName< OLS::ChainOption>()
      };
      return s_names[ size_t( CHAIN_OPTION)];
    }

    //! @brief Output filename as string
    //! @param ChainOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &OLS::GetOutputFileName( const ChainOption &CHAIN_OPTION, const bool SPLIT_ENVIRONMENT)
    {
      static std::string s_names_default[] =
      {
        "ols_one_chain.wo_env.histograms",
        "ols_all_chain.wo_env.histograms",
        GetStaticClassName< OLS::ChainOption>()
      };

      static std::string s_names_env[] =
      {
        "ols_one_chain.with_env.histograms",
        "ols_all_chain.with_env.histograms",
        GetStaticClassName< OLS::ChainOption>()
      };

      return SPLIT_ENVIRONMENT ? s_names_env[ size_t( CHAIN_OPTION)] : s_names_default[ size_t( CHAIN_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OLS::OLS()
    {
      // nothing else to do
    }

    //! @brief explicit constructor
    OLS::OLS( const bool SPLIT_ENVIRONMENT) :
        m_ChainOption( e_OneChain),
        m_ChainIds( ""),
        m_SequenceExclusion( 2),
        m_Radius( 4.75),
        m_SplitEnvironment( SPLIT_ENVIRONMENT)
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    OLS *OLS::Clone() const
    {
      return new OLS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &OLS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &OLS::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_ChainOption, m_SplitEnvironment);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &OLS::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the sequence exclusion considered when creating neighbor list
    //! @return the sequence exclusion considered when creating neighbor list
    const size_t &OLS::GetSequenceExclusion() const
    {
      return m_SequenceExclusion;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &OLS::GetAlias() const
    {
      static std::string s_name_one( "OLS"), s_name_two( "OLSByEnvironment");
      return m_SplitEnvironment ? s_name_two : s_name_one;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string OLS::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for neighbor count one chain, categorized by environment types
      storage::Vector< storage::Vector< math::Histogram> > desired_histograms
      (
        biol::GetEnvironmentTypes().GetEnumCount(),
        storage::Vector< math::Histogram>( biol::GetAATypes().GetEnumCount(), math::Histogram( 0, 0.01, 100))
      );

      // initialize an ols object
      assemble::AASasaOLS ols;
      ols.SetMinimalSequenceSeparation( m_SequenceExclusion);
      ols.SetThresholdRange( math::Range< double>( 0, m_Radius));

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get the membrane
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all amino acids and cb coordinates
        util::SiPtrVector< const biol::AABase> all_chain_amino_acids;

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all amino acids
          all_chain_amino_acids.Append( ( *chain_itr)->GetAminoAcids());
        }

        // get neighbor lists from all chains for computing neighbor vector statistics
        const assemble::AANeighborListContainer all_chain_neighbor_list_container
        (
          all_chain_amino_acids,
          ols.GetDistanceCutoff(),
          ols.GetMinimalSequenceSeparation(),
          true
        );

        // iterate through all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // ship undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get amino acids in current chain
          const util::SiPtrVector< const biol::AABase> amino_acids_current_chain( ( *chain_itr)->GetAminoAcids());

          // get neighbor lists from current chain for computing neighbor vector statistics
          const assemble::AANeighborListContainer current_chain_neighbor_list_container
          (
            amino_acids_current_chain,
            ols.GetDistanceCutoff(),
            ols.GetMinimalSequenceSeparation(),
            false
          );

          const assemble::AANeighborListContainer &desired_neighbor_list_container
          (
            m_ChainOption == e_OneChain
            ? current_chain_neighbor_list_container
            : all_chain_neighbor_list_container
          );

          // iterate through current chain
          for
          (
            util::SiPtrVector< const biol::AABase>::const_iterator
              aa_itr( amino_acids_current_chain.Begin()), aa_itr_end( amino_acids_current_chain.End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // proceed only if there are coordinates given and the aa is a natural aa
            if
            (
              !( *aa_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // get current environment type
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() && m_SplitEnvironment
              ? sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates())
              : biol::GetEnvironmentTypes().e_Solution
            );

            // calculate desired ols
            const double &desired_ols
            (
              ols( desired_neighbor_list_container.Find( ( **aa_itr))->second)
            );

            desired_histograms( current_environment_type)( ( *aa_itr)->GetType()).PushBack( desired_ols);

          } // end of iterating through all aas in current chain
        } // end of iterating through all chains
      } // end of iterating through all protein models

      // write output
      std::stringstream stream;

      stream << ols.GetThresholdRange() << '\n';
      stream << ols.GetMinimalSequenceSeparation() << '\n';

      // iterate over environment types
      for
      (
        biol::EnvironmentTypes::const_iterator
          ent_itr( biol::GetEnvironmentTypes().Begin()), ent_itr_end( biol::GetEnvironmentTypes().End());
        ent_itr != ent_itr_end;
        ++ent_itr
      )
      {
        // only interested in soluble environment type
        if( !m_SplitEnvironment)
        {
          if( ( *ent_itr) != biol::GetEnvironmentTypes().e_Solution)
          {
            continue;
          }
        }
        else
        {
          stream << *ent_itr << '\n';
        }

        // for current environment type iterate over all amino acid types
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\n';
          stream << desired_histograms( *ent_itr)( *aa_itr);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer OLS::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes neighbor count statistics."
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
        "sequence_exclusion",
        "minimal sequence separation to be considered in the calculation",
        io::Serialization::GetAgent( &m_SequenceExclusion),
        "2"
      );

      parameters.AddInitializer
      (
        "radius",
        "radius for the overlapping sphere algorithm",
        io::Serialization::GetAgent( &m_Radius),
        "4.75"
      );

      parameters.AddInitializer
      (
        "chain_option",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_ChainOption),
        GetChainOptionName( e_OneChain)
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl

