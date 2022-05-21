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
#include "scorestat/bcl_scorestat_neighbor_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborCount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborCount( false))
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborCount::s_InstanceEnv
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborCount( true))
    );

    //! @brief ChainOption as string
    //! @param CHAIN_OPTION the ChainOption
    //! @return the string for the ChainOption
    const std::string &NeighborCount::GetChainOptionName( const ChainOption &CHAIN_OPTION)
    {
      static std::string s_names[] =
      {
        "One",
        "All",
        GetStaticClassName< NeighborCount::ChainOption>()
      };
      return s_names[ size_t( CHAIN_OPTION)];
    }

    //! @brief Output filename as string
    //! @param ChainOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &NeighborCount::GetOutputFileName( const ChainOption &CHAIN_OPTION, const bool SPLIT_ENVIRONMENT)
    {
      static std::string s_names_default[] =
      {
        "nc_one_chain.wo_env.histograms",
        "nc_all_chain.wo_env.histograms",
        GetStaticClassName< NeighborCount::ChainOption>()
      };

      static std::string s_names_env[] =
      {
        "nc_one_chain.with_env.histograms",
        "nc_all_chain.with_env.histograms",
        GetStaticClassName< NeighborCount::ChainOption>()
      };

      return SPLIT_ENVIRONMENT ? s_names_env[ size_t( CHAIN_OPTION)] : s_names_default[ size_t( CHAIN_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NeighborCount::NeighborCount()
    {
      // nothing else to do
    }

    //! @brief explicit constructor
    NeighborCount::NeighborCount( const bool SPLIT_ENVIRONMENT) :
        m_ChainOption( e_OneChain),
        m_ChainIds( ""),
        m_SequenceExclusion( 2),
        m_NCLowerBound( 4.0),
        m_NCUpperBound( 11.4),
        m_SplitEnvironment( SPLIT_ENVIRONMENT)
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    NeighborCount *NeighborCount::Clone() const
    {
      return new NeighborCount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &NeighborCount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &NeighborCount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_ChainOption, m_SplitEnvironment);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &NeighborCount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the sequence exclusion considered when creating neighbor list
    //! @return the sequence exclusion considered when creating neighbor list
    const size_t &NeighborCount::GetSequenceExclusion() const
    {
      return m_SequenceExclusion;
    }

    //! @brief returns the lower bound for neighbor count
    //! @return the lower bound for neighbor count
    const double &NeighborCount::GetNCLowerBound() const
    {
      return m_NCLowerBound;
    }

    //! @brief returns the upper bound for neighbor count
    //! @return the upper bound for neighbor count
    const double &NeighborCount::GetNCUpperBound() const
    {
      return m_NCUpperBound;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &NeighborCount::GetAlias() const
    {
      static std::string s_name_one( "NeighborCount"), s_name_two( "NeighborCountByEnvironment");
      return m_SplitEnvironment ? s_name_two : s_name_one;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string NeighborCount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for neighbor count one chain, categorized by environment types
      storage::Vector< storage::Vector< math::Histogram> > desired_histograms
      (
        biol::GetEnvironmentTypes().GetEnumCount(),
        storage::Vector< math::Histogram>( biol::GetAATypes().GetEnumCount(), math::Histogram( 0, 1.0, 50))
      );

      // initialize a neighbor count object
      assemble::AANeighborCount neighbor_count;
      neighbor_count.SetMinimalSequenceSeparation( m_SequenceExclusion);
      neighbor_count.SetThresholdRange( math::Range< double>( m_NCLowerBound, m_NCUpperBound));

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
        const assemble::AANeighborListContainer all_chain_neighbor_list_nc
        (
          all_chain_amino_acids,
          neighbor_count.GetDistanceCutoff(),
          neighbor_count.GetMinimalSequenceSeparation(),
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
          const assemble::AANeighborListContainer current_chain_neighbor_list_nc
          (
            amino_acids_current_chain,
            neighbor_count.GetDistanceCutoff(),
            m_SequenceExclusion,
            false
          );

          const assemble::AANeighborListContainer &desired_neighbor_list
          (
            m_ChainOption == e_OneChain
            ? current_chain_neighbor_list_nc
            : all_chain_neighbor_list_nc
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
              sp_membrane.IsDefined() && m_SplitEnvironment ?
                  sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates()) :
                  biol::GetEnvironmentTypes().e_Solution
            );

            // calculate desired neighbor count
            const double &desired_neighbor_count
            (
              neighbor_count( desired_neighbor_list.Find( ( **aa_itr))->second)
            );

            desired_histograms( current_environment_type)( ( *aa_itr)->GetType()).PushBack( desired_neighbor_count);

          } // end of iterating through all aas in current chain
        } // end of iterating through all chains
      } // end of iterating through all protein models

      // write output
      std::stringstream stream;

      stream << neighbor_count.GetThresholdRange() << '\n';
      stream << neighbor_count.GetMinimalSequenceSeparation() << '\n';

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
    io::Serializer NeighborCount::GetSerializer() const
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
        "lower_bound",
        "lower bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NCLowerBound),
        "4.0"
      );

      parameters.AddInitializer
      (
        "upper_bound",
        "upper bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NCUpperBound),
        "11.4"
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

