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
#include "scorestat/bcl_scorestat_aa_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AACount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AACount())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &AACount::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        GetStaticClassName< AACount::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &AACount::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "aa_count.tbl",
        GetStaticClassName< AACount::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AACount::AACount():
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    AACount *AACount::Clone() const
    {
      return new AACount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &AACount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AACount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &AACount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AACount::GetAlias() const
    {
      static std::string s_name( "AACount");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AACount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // aa count distribution according to environmental types, such as membrane core, transition, soluble
      storage::Vector< storage::Vector< size_t> > aa_region_count
      (
        biol::GetAATypes().GetEnumCount(),
        storage::Vector< size_t>( biol::GetEnvironmentTypes().GetEnumCount(), size_t( 0))
      );

      // aa count distribution according to secondary structure elements and environmental types
      storage::Vector< storage::Vector< storage::Vector< size_t> > > aa_sse_region_count
      (
        biol::GetSSTypes().COIL.GetIndex() + 1,
        storage::Vector< storage::Vector< size_t> >
        (
          biol::GetAATypes().GetEnumCount(),
          storage::Vector< size_t>( biol::GetEnvironmentTypes().GetEnumCount(), size_t( 0))
        )
      );

      // iterate through all protein models
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? sp_model_filename->GetData() : "");

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        if( !sp_membrane.IsDefined())
        {
          BCL_MessageDbg( util::Format()( model_filename) + " does not have an associated membrane, assuming it's a soluble protein.");
        }

        // iterate over all chains
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
            continue;
          }

          // get all sses in current chain
          const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
            &all_sses( ( *chain_itr)->GetData());

          // iterate over all sses in current chain
          for
          (
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && ( *sse_itr)->GetType() == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !( *sse_itr)->GetType().IsDefined())
            {
              continue;
            }

            // get the type of current sse
            const biol::SSType &current_sse_type( ( *sse_itr)->GetType());

            // get all amino acids in current sse
            const util::ShPtrVector< biol::AABase> &all_amino_acids( ( *sse_itr)->GetData());

            // iterate over all amino acids in current sse
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                aa_itr( all_amino_acids.Begin()), aa_itr_end( all_amino_acids.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // get environment type in which the current amino acid can be found
              const biol::EnvironmentType current_environment_type
              (
                sp_membrane.IsDefined() ?
                    sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates()) :
                      biol::GetEnvironmentTypes().e_Solution
              );

              // skip undefined environment
              if( !current_environment_type.IsDefined())
              {
                continue;
              }

              ++aa_region_count( ( *aa_itr)->GetType())( current_environment_type);
              ++aa_sse_region_count( current_sse_type)( ( *aa_itr)->GetType())( current_environment_type);
            } // end of iterating over all amino acids
          } // end of iterating over all sses in current chain
        } // end of iterating over all chains in current protein model
      } // end of iterating over all protein models in the ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
        stream << "# aa count in each environment type" << '\n';
        stream << '\n';
        stream << "environment_type" << '\t';

        for
        (
          biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\t';
        }
        stream << '\n';

        for
        (
          biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().Begin()),
            env_itr_end( biol::GetEnvironmentTypes().End());
          env_itr != env_itr_end;
          ++env_itr
        )
        {
          stream << *env_itr << '\t';

          // all amino acid types
          for
          (
            biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            stream << aa_region_count( *aa_itr)( *env_itr) << '\t';
          }
          stream << '\n';
        }

        stream << '\n';
        stream << "# aa count in each environment type categorized by sse type" << '\n';

        for
        (
          biol::SSTypes::const_iterator sse_itr( biol::GetSSTypes().Begin()),
            sse_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          stream << '\n';
          stream << ( *sse_itr)->GetName() << '\n';
          stream << "environment_type" << '\t';

          for
          (
            biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            stream << ( *aa_itr)->GetOneLetterCode() << '\t';
          }
          stream << '\n';

          for
          (
            biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().Begin()),
              env_itr_end( biol::GetEnvironmentTypes().End());
            env_itr != env_itr_end;
            ++env_itr
          )
          {
            stream << *env_itr << '\t';

            // all amino acid types
            for
            (
              biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
                aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              stream << aa_sse_region_count( *sse_itr)( *aa_itr)( *env_itr) << '\t';
            }
            stream << '\n';
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
    io::Serializer AACount::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes amino acid occurrence statistics."
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
        "Table"
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl
