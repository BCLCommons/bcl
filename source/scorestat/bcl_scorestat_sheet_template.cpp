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
#include "scorestat/bcl_scorestat_sheet_template.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_collector_topology_interface.h"
#include "assemble/bcl_assemble_collector_topology_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SheetTemplate::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SheetTemplate())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SheetTemplate::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "List",
        GetStaticClassName< SheetTemplate::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SheetTemplate::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "sheet_template.list",
        GetStaticClassName< SheetTemplate::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SheetTemplate::SheetTemplate() :
      m_OutputOption( e_List),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SheetTemplate *SheetTemplate::Clone() const
    {
      return new SheetTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SheetTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SheetTemplate::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SheetTemplate::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SheetTemplate::GetAlias() const
    {
      static std::string s_name( "SheetTemplate");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SheetTemplate::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // vector to hold sheet templates
      storage::Vector< assemble::FoldTemplate> sheet_templates;

      // initialize collector
      const util::ShPtr< assemble::CollectorTopologyInterface> sp_collector
      (
        new assemble::CollectorTopologySheet()
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
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get pdb id
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        std::string pdb_id( io::File::RemovePath( model_name));
        std::transform( pdb_id.begin(), pdb_id.end(), pdb_id.begin(), toupper);

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

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

          // create a chain model from current chain
          const assemble::ProteinModel &chain_model( util::ShPtr< assemble::Chain>( ( *chain_itr)->HardCopy()));

          // collect sheets
          const util::ShPtrVector< assemble::Domain> &sheets( assemble::CollectorSheet().Collect( chain_model));

          BCL_MessageDbg( "number sheets found: " + util::Format()( sheets.GetSize()));
          size_t sheet_ctr( 0);

          // iterate over sheets
          for
          (
            util::ShPtrVector< assemble::Domain>::const_iterator
              sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
            sheet_itr != sheet_itr_end;
            ++sheet_itr
          )
          {
            ++sheet_ctr;

            BCL_MessageDbg
            (
              "Sheet #" + util::Format()( sheet_ctr) + "\n"
              + ( *sheet_itr)->GetTopology()->GetOrderedIdentification()
            );

            // if less than 2 sses, skip
            if( ( *sheet_itr)->GetNumberSSEs() < 2)
            {
              continue;
            }

            // get the elements vector
            const util::SiPtrVector< const assemble::SSEGeometryInterface> elements_vector
            (
              ( *sheet_itr)->GetTopology()->GetElements()
            );

            // initialize ShPtrVector for geometries
            util::ShPtrVector< assemble::SSEGeometryPhiPsi> geometry_vector;

            // iterate over sses in current elements_vector
            for
            (
              util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
                element_itr( elements_vector.Begin()), element_itr_end( elements_vector.End());
              element_itr != element_itr_end;
              ++element_itr
            )
            {
              // get sses in current sheet
              const util::SiPtrVector< const assemble::SSE> &sheet_sses( ( *sheet_itr)->GetSSEs());

              // iterate over all sses in the sheet
              for
              (
                util::SiPtrVector< const assemble::SSE>::const_iterator
                  sheet_sse_itr( sheet_sses.Begin()), sheet_sse_itr_end( sheet_sses.End());
                sheet_sse_itr != sheet_sse_itr_end;
                ++sheet_sse_itr
              )
              {
                // if geometries are equal
                if( **element_itr == **sheet_sse_itr)
                {
                  geometry_vector.PushBack
                  (
                    util::ShPtr< assemble::SSEGeometryInterface>( new assemble::SSEGeometryPhiPsi( **sheet_sse_itr))
                  );
                }
              }
            }

            // create the fold template
            const assemble::FoldTemplate sheet_template
            (
              geometry_vector,
              sp_collector,
              pdb_id,
              false
            );

            // insert the template if it is defined
            if( sheet_template.HasDefinedGeometries())
            {
              sheet_templates.PushBack( sheet_template);
            }
          }
        } // end of iteration over all chains
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_List)
      {
        stream << sheet_templates.GetSize() << '\n';
        // iterate over sheet_templates
        for
        (
          storage::Vector< assemble::FoldTemplate>::const_iterator
            sheet_template_itr( sheet_templates.Begin()), sheet_template_itr_end( sheet_templates.End());
          sheet_template_itr != sheet_template_itr_end;
          ++sheet_template_itr
        )
        {
          sheet_template_itr->WriteCompact( stream);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SheetTemplate::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes sheet template statistics."
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
        "List"
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl

