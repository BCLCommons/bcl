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
#include "scorestat/bcl_scorestat_sse_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSECount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSECount())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSECount::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        GetStaticClassName< SSECount::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSECount::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extentions[] =
      {
        "sse_count.tbl",
        GetStaticClassName< SSECount::OutputOption>()
      };
      return s_output_filename_extentions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSECount::SSECount() :
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSECount *SSECount::Clone() const
    {
      return new SSECount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSECount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSECount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSECount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSECount::GetAlias() const
    {
      static std::string s_Name( "SSECount");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSECount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // create table for holding sse count statistics
      storage::Table< double> sse_count_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "HELIX", "STRAND", "COIL"
          )
        )
      );

      // statistics of sse count
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, size_t> > sse_region_count;

      // iterate over all protein models in the ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( ( **protein_model_itr));

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> &sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        if( !sp_membrane.IsDefined())
        {
          BCL_MessageDbg( "Skip: " + util::Format()( model_filename) + " does not have an associated membrane");
          continue;
        }

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains of current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains, if m_ChainIds is empty, analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( "Skip undesired chain: " + util::Format()( ( *chain_itr)->GetChainID()));
            continue;
          }

          // get all sse in current chain
          const util::SiPtrVector< const assemble::SSE> &all_sses( ( *chain_itr)->GetSSEs());

          // iterate over all sses in current chain
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // get current sse type
            const biol::SSType &current_sse_type( ( *sse_itr)->GetType());

            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && current_sse_type == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !current_sse_type.IsDefined())
            {
              continue;
            }

            // get environment type of current sse
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() ?
                  sp_membrane->DetermineEnvironmentType( ( *sse_itr)->GetCenter()) :
                    biol::GetEnvironmentTypes().e_Solution
            );

            ++sse_region_count[ ( *sse_itr)->GetType()][ current_environment_type];
          } // end of iteration over all sses in current chain
        } // end of iteration over all chains in current protein model
      } // end of iteration over protein ensemble

      // store sse statistics to a table
      if( m_OutputOption == e_Table)
      {
        // iterate over environment types
        for
        (
          biol::EnvironmentTypes::const_iterator
            env_itr( biol::GetEnvironmentTypes().Begin()), env_itr_end( biol::GetEnvironmentTypes().End());
          env_itr != env_itr_end;
          ++env_itr
        )
        {
          sse_count_table.InsertRow
          (
            ( *env_itr)->GetName(),
            storage::Vector< double>::Create
            (
              sse_region_count[ biol::GetSSTypes().HELIX][ *env_itr],
              sse_region_count[ biol::GetSSTypes().STRAND][ *env_itr],
              sse_region_count[ biol::GetSSTypes().COIL][ *env_itr]
            )
          );
        }

      }

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
         sse_count_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSECount::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics of occurrence of secondary structure elements."
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

