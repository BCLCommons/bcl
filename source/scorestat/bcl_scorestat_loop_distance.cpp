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
BCL_StaticInitializationFiascoFinder

// include header for this class
#include "scorestat/bcl_scorestat_loop_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_loop.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopDistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopDistance())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopDistance::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_Names[] =
      {
        "Table",
        "Histogram",
        "LogHistogram",
        GetStaticClassName< LoopDistance::OutputOption>()
      };
      return s_Names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopDistance::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_file_extensions[] =
      {
        "loop_distance.tbl",
        "loop_distance_raw.histograms",
        "loop_distance_log.histograms",
        GetStaticClassName< OutputOption>()
      };

      return s_output_file_extensions[ OUTPUT_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopDistance::LoopDistance() :
      m_OutputOption( e_Table),
      m_OutFilePostFix( ".loopdistance"),
      m_BinSize( 1.0),
      m_ChainIds( "")
    {
    }

    //! @brief virtual copy constructor
    LoopDistance *LoopDistance::Clone() const
    {
      return new LoopDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopDistance::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the binsize for the histogram
    //! @return the binsize for the histogram
    const double &LoopDistance::GetBinSize() const
    {
      return m_BinSize;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &LoopDistance::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopDistance::GetAlias() const
    {
      static const std::string s_Name( "LoopDistance");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string LoopDistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // three output formats
      storage::Vector< math::Histogram> loop_dist_histogram
      (
        500, math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize))
      );

      storage::Table< double> loop_dist_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sequence_dist", "euclidean_dist", "de/log(ds)", "consecutive_sses"
          )
        )
      );

      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

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
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? io::File::RemovePath( model_name_ptr->GetData()) : "");

        // iterate over all chains of the protein model
        for
        (
           util::ShPtrVector< assemble::Chain>::const_iterator
             chain_itr( protein_model.GetChains().Begin()), chain_itr_end( protein_model.GetChains().End());
           chain_itr != chain_itr_end;
           ++chain_itr
        )
        {
          // get current chain
          util::ShPtr< assemble::Chain> current_chain( *chain_itr);

          BCL_MessageDbg
          (
            "Calculating loop distance statistics for chain, chain.sse.size="
            + util::Format()( current_chain->GetData().GetSize())
          );

          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( current_chain->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg
            (
              std::string( "Skip chain: ") + current_chain->GetChainID() + std::string( ", chain ids to use: ") + m_ChainIds
            );

            continue;
          }

          // iterate over all sses in the current chain
          size_t sse_a_number( 0);
          for
          (
             storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
               sse_a_itr( current_chain->GetData().Begin()), sse_itr_end( current_chain->GetData().End());
              sse_a_itr != sse_itr_end;
              ++sse_a_itr
          )
          {

            // get the first sse of the current chain
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_b_itr( sse_a_itr);

            // skip undefined sse
            if( !( *sse_a_itr)->IsDefined())
            {
              continue;
            }

            size_t sse_b_number( sse_a_number);
            bool is_sse_a_b_consecutive( false);
            if( sse_b_itr != sse_itr_end)
            {
              ++sse_b_itr;
              ++sse_b_number;
              is_sse_a_b_consecutive = true;
            }

            // iterate over the rest sses
            for( ; sse_b_itr != sse_itr_end; ++sse_b_itr)
            {
              // skip undefined sse
              if( !( *sse_b_itr)->IsDefined())
              {
                BCL_MessageDbg( "Skip: sse_b is undefined");
                continue;
              }

              BCL_MessageDbg
              (
                "Calculate loop distance statistics for sse pair " +
                util::Format()( ( *sse_a_itr)->GetIdentification()) +
                " " + util::Format()( ( *sse_b_itr)->GetIdentification())
              );

              // calculate sequence distance and euclidean distance
              const storage::Pair< size_t, double> seq_euc_distance
              (
                score::Loop::SequenceAndEuclideanDistance( **sse_a_itr, **sse_b_itr)
              );

              // store sequence distance and euclidean distance in histogram
              if( seq_euc_distance.First() >= loop_dist_histogram.GetSize())
              {
                if( m_OutputOption == e_Histogram)
                {
                  loop_dist_histogram.LastElement().PushBack( seq_euc_distance.Second());
                }
                else if( m_OutputOption == e_LogHistogram)
                {
                  loop_dist_histogram.LastElement().PushBack( score::Loop::NormalizeDistance( seq_euc_distance));
                }
              }
              else
              {
                if( m_OutputOption == e_Histogram)
                {
                  loop_dist_histogram( seq_euc_distance.First()).PushBack( seq_euc_distance.Second());
                }
                else if( m_OutputOption == e_LogHistogram)
                {
                  loop_dist_histogram( seq_euc_distance.First()).PushBack( score::Loop::NormalizeDistance( seq_euc_distance));
                }
              }

              // add loop distance statistics to table
              if
              (
                  ( **sse_a_itr).GetType() != biol::GetSSTypes().COIL
                  && ( **sse_b_itr).GetType() != biol::GetSSTypes().COIL
              )
              {
                BCL_MessageDbg
                (
                  "Insert into table: seq_dist=" + util::Format()( seq_euc_distance.First())
                  + "  euc_dist=" + util::Format()( seq_euc_distance.Second())
                  + "  norm_dist=" + util::Format()( score::Loop::NormalizeDistance( seq_euc_distance))
                );

                if( m_OutputOption == e_Table)
                {
                  const std::string sse_pair( "sse" + util::Format()( sse_a_number) + "_sse" + util::Format()( sse_b_number));
                  loop_dist_table.InsertRow
                  (
                    model_name + "_" + sse_pair,
                    storage::Vector< double>::Create
                    (
                      seq_euc_distance.First(),
                      seq_euc_distance.Second(),
                      score::Loop::NormalizeDistance( seq_euc_distance),
                      is_sse_a_b_consecutive
                    ),
                    true
                  );
                }

                BCL_MessageDbg( "table size=" + util::Format()( loop_dist_table.GetSize()));

                is_sse_a_b_consecutive = false;
                ++sse_b_number;
              }
            }

            if( ( **sse_a_itr).GetType() != biol::GetSSTypes().COIL)
            {
              ++sse_a_number;
            }
          }
        } // end of the current chain
      } // end of current model

      // write statistics
      std::ostringstream stream;
      if( m_OutputOption == e_Histogram || m_OutputOption == e_LogHistogram)
      {
        stream << loop_dist_histogram;

      }
      else if( m_OutputOption == e_Table)
      {
        loop_dist_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopDistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes loop distance statistics."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".loopdistance"
      );

      parameters.AddInitializer
      (
        "bin_size",
        "the bin size for the histogram",
        io::Serialization::GetAgent( &m_BinSize)
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

