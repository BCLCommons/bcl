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
#include "assemble/bcl_assemble_fold_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_collector_topology_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "scorestat/bcl_scorestat_fold_template.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FoldTemplate::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new FoldTemplate())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &FoldTemplate::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
          "List",
          GetStaticClassName< FoldTemplate::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &FoldTemplate::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_file_extensions[] =
      {
          "fold_template.list",
          GetStaticClassName< FoldTemplate::OutputOption>()
      };
      return s_output_file_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FoldTemplate::FoldTemplate() :
        m_OutputOption( e_List),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    FoldTemplate *FoldTemplate::Clone() const
    {
      return new FoldTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &FoldTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &FoldTemplate::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &FoldTemplate::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &FoldTemplate::GetAlias() const
    {
      static std::string s_name( "FoldTemplate");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string FoldTemplate::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the proein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "The protein ensemble is empty.");

      // create variables for storing output information
      storage::Vector< assemble::FoldTemplate> fold_templates;
      storage::Map
      <
        storage::Pair< size_t, size_t>,
        storage::Vector< storage::Pair< std::string, double> >
      > sorted_pdbs;

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &protein_model( **protein_model_itr);

        // get pdb id
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        std::string pdb_id( io::File::RemovePath( model_name));
        std::transform( pdb_id.begin(), pdb_id.end(), pdb_id.begin(), toupper);

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( protein_model.GetChains().Begin()), chain_itr_end( protein_model.GetChains().End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip chains that are not desired
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // create a protein model from this chain only
          assemble::ProteinModel chain_model( util::ShPtr< assemble::Chain>( ( *chain_itr)->HardCopy()));

          // create the fold template
          const util::ShPtr< assemble::CollectorTopologyInterface> sp_collector
          (
            new assemble::CollectorTopologyCombined( false)
          );
          chain_model.Transform( math::Inverse( chain_model.GetOrientation()));
          const assemble::FoldTemplate fold_template( chain_model, sp_collector, pdb_id);

          // push back the fold template information if it has at least one geometry, all defined geometries, and an appropriate Rg
          if
          (
            !fold_template.GetGeometries().IsEmpty() &&
            fold_template.HasDefinedGeometries() &&
            fold_template.GetRadiusOfGyration() < 3.0 * double( fold_template.GetGeometries().GetSize())
          )
          {
            fold_templates.PushBack( fold_template);
            sorted_pdbs
            [
               storage::Pair< size_t, size_t>
              (
                fold_template.GetHelicalGeometries().GetSize(), fold_template.GetStrandGeometries().GetSize()
              )
            ].PushBack
            (
              storage::Pair< std::string, double>( fold_template.GetPDBID(), fold_template.GetRadiusOfGyration())
            );
          }
        } // end of iterating over all chains
      } // end of iterating over the current protein model

      // write list of pdbs with fold template information
      io::OFStream write;
      io::File::MustOpenOFStream( write, "fold_template_pdbs.list");
      write << sorted_pdbs;
      io::File::CloseClearFStream( write);

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_List)
      {
        stream << fold_templates.GetSize() << '\n';
        // iterate over fold templates
        for
        (
          storage::Vector< assemble::FoldTemplate>::const_iterator
            fold_template_itr( fold_templates.Begin()), fold_template_itr_end( fold_templates.End());
          fold_template_itr != fold_template_itr_end;
          ++fold_template_itr
        )
        {
          fold_template_itr->WriteCompact( stream);
        }
      }
      return stream.str();
    } // end of operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FoldTemplate::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes fold template statistics."
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

