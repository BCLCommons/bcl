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
#include "scorestat/bcl_scorestat_contact_order.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_contact_order.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContactOrder::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ContactOrder())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &ContactOrder::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        "Histogram",
        GetStaticClassName< ContactOrder::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &ContactOrder::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "contact_order.tbl",
        "contact_order.histograms",
        GetStaticClassName< ContactOrder::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

    //! @brief default constructor
    ContactOrder::ContactOrder() :
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    ContactOrder *ContactOrder::Clone() const
    {
      return new ContactOrder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &ContactOrder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ContactOrder::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &ContactOrder::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ContactOrder::GetAlias() const
    {
      static std::string s_name( "ContactOrder");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string ContactOrder::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // table for holding contact order statistics
      storage::Table< double> contact_order_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "nr_aas", "nr_aas_sses", "co_chain_sse", "co_chain_seq", "co_seq", "co_chain_sse_sqr", "co_chain_seq_sqr", "co_seq_sqr"
          )
        )
      );

      // create all types of contact orders
      const contact::Order co_relative_sse( contact::Order::e_RelativeAAsUsed, "co_relative_sse", false);
      const contact::Order co_relative_seq( contact::Order::e_RelativeSequenceLength, "co_relative_length", false);
      const contact::Order co_relative_sse_sqr( contact::Order::e_RelativeSqrAAsUsed, "co_relative_sse_sqr", false);
      const contact::Order co_relative_seq_sqr( contact::Order::e_RelativeSqrSequenceLength, "co_relative_seq_sqr", false);

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
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get current protein model name
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_name( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( sp_model_name.IsDefined() ? sp_model_name->GetData() : "");

        // get all chains
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all sequences
        const util::SiPtrVector< const biol::AASequence> &all_sequences( current_protein_model.GetSequences());

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
            continue;
          }

          // insert rows
          if( !contact_order_table.HasRow( model_name + ( *chain_itr)->GetChainID()))
          {
            contact_order_table.InsertRow( model_name + ( *chain_itr)->GetChainID());
          }

          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "nr_aas_sses"] = ( *chain_itr)->GetNumberAAs();
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_sse"] = co_relative_sse.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_seq"] = co_relative_seq.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_sse_sqr"] = co_relative_sse_sqr.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_seq_sqr"] = co_relative_seq_sqr.ContactOrder( **chain_itr);
        } // end of iteration over all chains

        // iterate over all sequences
        for
        (
          util::SiPtrVector< const biol::AASequence>::const_iterator
            aaseq_itr( all_sequences.Begin()), aaseq_itr_end( all_sequences.End());
          aaseq_itr != aaseq_itr_end;
          ++aaseq_itr
        )
        {
          if( !contact_order_table.HasRow( model_name + ( *aaseq_itr)->GetChainID()))
          {
            contact_order_table.InsertRow( model_name + ( *aaseq_itr)->GetChainID());
          }

          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "nr_aas"] = ( *aaseq_itr)->GetSize();
          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "co_seq"] = co_relative_seq.ContactOrder( **aaseq_itr);
          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "co_seq_sqr"] = co_relative_seq_sqr.ContactOrder( **aaseq_itr);
        } // end of iteration over sequences
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
        // contact order summary
        stream << "contact order summary" << '\n';
        contact_order_table.WriteFormatted( stream);
        stream << '\n';
      }
      else if( m_OutputOption == e_Histogram)
      {
        // contact order histograms
        storage::Map< std::string, math::Histogram> co_histograms( score::ContactOrder::HistogramsFromColumns( contact_order_table));

        // iterate over histograms
        for
        (
          storage::Map< std::string, math::Histogram>::const_iterator
            histogram_itr( co_histograms.Begin()), histogram_itr_end( co_histograms.End());
          histogram_itr != histogram_itr_end;
          ++histogram_itr
        )
        {
          // skip the first two columns
          if( histogram_itr->first == "nr_aas" || histogram_itr->first == "nr_aas_sses")
          {
            continue;
          }

          stream << histogram_itr->first + " histogram" << '\n';
          stream << histogram_itr->second << '\n';
        }
      }

      return stream.str();
    } // namespace scorestat

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContactOrder::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes contact order statistics."
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
        "Histogram"
      );

      return parameters;
    }

  } // namespace scorestat
} // namespace bcl

