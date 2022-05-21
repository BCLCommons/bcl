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
#include "scorestat/bcl_scorestat_loop_angle.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_aa.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_loop.h"
#include "score/bcl_score_loop_angle.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopAngle::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopAngle())
    );

    //! @brief OutputOption as string
    //! @param OutputOption the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopAngle::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_descriptors[] =
      {
        "NormalizedByTotalCount",
        "NormalizedByColumn",
        "Log",
        "Raw",
        "Table",
        GetStaticClassName< OutputOption>()
      };

      return s_descriptors[ OUTPUT_OPTION];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopAngle::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_file_extensions[] =
      {
        "loop_angle_count_norm.histograms",
        "loop_angle_column_norm.histograms",
        "loop_angle_log.histograms",
        "loop_angle_raw.histograms",
        "loop_angle.tbl",
        GetStaticClassName< OutputOption>()
      };

      return s_output_file_extensions[ OUTPUT_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopAngle::LoopAngle() :
      m_OutputOption( e_NormalizedByCount),
      m_Chains( "A"),
      m_VisualizationFlag( false)
    {
    }

    //! @brief virtual copy constructor
    LoopAngle *LoopAngle::Clone() const
    {
      return new LoopAngle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopAngle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopAngle::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopAngle::GetAlias() const
    {
      static const std::string s_Name( "LoopAngle");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string LoopAngle::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // Initialize variables
      size_t max_sequence_distance( score::LoopAngle::GetDefaultMaxmimumSequenceDistance());

      // min value for cosine
      const double min_cos( -1.0);

      // minimum loop size
      const double min_loop_len( 0.0);

      // loop size increment
      const double loop_size_increment( 1.0);

      // Number of bins as integer
      const int x_num_bins( s_DefaultNumberBinsX);
      const int y_num_bins( s_DefaultNumberBinsY);

      // compute bin size for y axis
      const double bin_size_y( 2.0 / double( y_num_bins));

      // loops with sequence distance <= m_MaxsSequenceDistance
      math::Histogram cos_angle_short_loops_histogram( min_cos, bin_size_y, s_DefaultNumberBinsY);

      // loops with sequence distance > m_MaxsSequenceDistance
      math::Histogram cos_angle_long_loops_histogram( min_cos, bin_size_y, s_DefaultNumberBinsY);

      math::Histogram2D euclidean_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      math::Histogram2D sequence_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      math::Histogram2D euclidean_over_sequence_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      storage::Table< double> loop_angle_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sequence_dist", "euclidean_dist", "de/log(ds)", "consecutive_sses", "cos_angle"
          )
        )
      );

      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get pdb filename
        const util::ShPtr< util::Wrapper< std::string> > &model_filename_ptr
        (
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );

        const std::string model_filename( model_filename_ptr.IsDefined() ? model_filename_ptr->GetData() : "");

        // create pdb factory for writing visualization pdbs
        pdb::Factory factory;

        //iterate over all chains
        const util::ShPtrVector< assemble::Chain> &chains( protein_model.GetChains());
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( chains.Begin()), chain_itr_end( chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // To simplify naming in the code
          const assemble::Chain &chain( **chain_itr);

          // skip chains that are not desired
          if( m_Chains.find( chain.GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( "Skip chain " + util::Format()( chain.GetChainID()) + std::string( ", not in chains to use: ") + m_Chains);
            continue;
          }

          // calculate radius of gyration
          const double radius_of_gyration( coord::RadiusOfGyration( chain.GetAtomCoordinates()));
          BCL_MessageDbg( "radius_of_gyration=" + util::Format()( radius_of_gyration));
          // calculate center of mass
          linal::Vector3D center_of_mass( coord::CenterOfMass( chain.GetAtomCoordinates(), true));
          BCL_MessageDbg( "center_of_mass=" + center_of_mass.ToString());

          // for each pair of (consecutive) sses: calculate euclidean distance and angles to center of gravity
          storage::Set< biol::SSType> ss_types;
          ss_types.InsertElement( biol::GetSSTypes().HELIX); // only collect helices and strands, no coil
          ss_types.InsertElement( biol::GetSSTypes().STRAND);
          util::SiPtrVector< const assemble::SSE> sses( chain.GetSSEs( ss_types));
          size_t sse_a_number( 0); // count the position of sse_a to clarify the output

          // Storage List appended atom lines
          util::ShPtrList< pdb::Line> vis_file_lines;

          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_a( sses.Begin()), sse_itr_end( sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a, ++sse_a_number
          )
          {
            if( !( *sse_itr_a)->IsDefined()) // skip sses with undefined body
            {
              continue;
            }

            // start itr for second sse one after the first
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a);
            size_t sse_b_number( sse_a_number); // count the position of sse_b to clarify the output
            bool sse_a_b_consecutive( false); // currently sse_itr_b points to the same sse as sse_itr_a
            if( sse_itr_b != sse_itr_end) // move itr to next sse if not at end
            {
              ++sse_itr_b;
              ++sse_b_number;
              sse_a_b_consecutive = true; // now it points to the next one after sse_itr_a
            }
            // iterate over all pairs of sses
            for( ; sse_itr_b != sse_itr_end; ++sse_itr_b, ++sse_b_number, sse_a_b_consecutive = false)
            {
              if( !( *sse_itr_b)->IsDefined()) // skip sses with undefined body
              {
                continue;
              }

              BCL_MessageDbg
              (
                "Calculate loop angle statistics for sse pair: sse" + util::Format()( sse_a_number) + "="
                + util::Format()( ( **sse_itr_a).GetIdentification())
                + "  sse" + util::Format()( sse_b_number) + "=" + util::Format()( ( **sse_itr_b).GetIdentification())
                + "  consecutive=" + util::Format()( sse_a_b_consecutive)
              );

              // calculate the sequence distance between the two SSEs
              const util::ShPtr< biol::AABase> last_aa_of_prev_sse( ( **sse_itr_a).GetLastAA());
              const util::ShPtr< biol::AABase> first_aa_of_next_sse( ( **sse_itr_b).GetFirstAA());
              BCL_MessageDbg( "seq_id_last_aa_of_prev_sse=" + util::Format()( last_aa_of_prev_sse->GetSeqID()));
              BCL_MessageDbg( "seq_id_first_aa_of_next_sse=" + util::Format()( first_aa_of_next_sse->GetSeqID()));

              // calculate the sequence distance
              const size_t sequence_distance( biol::SequenceSeparation( *last_aa_of_prev_sse, *first_aa_of_next_sse));
              BCL_MessageDbg( "sequence_distance=" + util::Format()( sequence_distance));

              // get positions of begin and end of z-axis of the two SSEs
              const linal::Vector3D prev_end_of_z( ( **sse_itr_a).EndOfZ());
              const linal::Vector3D next_begin_of_z( ( **sse_itr_b).BeginOfZ());
              BCL_MessageDbg( "prev_end_of_z=" + prev_end_of_z.ToString());
              BCL_MessageDbg( "next_begin_of_z=" + next_begin_of_z.ToString());

              // calculate the euclidean distance
              const double euclidean_distance( linal::Distance( next_begin_of_z, prev_end_of_z));
              BCL_MessageDbg( "euclidean_distance=" + util::Format()( euclidean_distance));

              // cosine of projection angle between center_of_mass->prev_end_of_z and center_of_mass->next_begin_of_z
              const double cosine_of_proj_angle( linal::ProjAngleCosinus( center_of_mass, prev_end_of_z, next_begin_of_z));
              BCL_MessageDbg( "cosine_of_proj_angle=" + util::Format()( cosine_of_proj_angle));

              if( m_OutputOption != e_Table)
              {
                // add data to histograms
                if( sequence_distance <= max_sequence_distance)
                {
                  cos_angle_short_loops_histogram.PushBack( cosine_of_proj_angle);
                }
                else
                {
                  cos_angle_long_loops_histogram.PushBack( cosine_of_proj_angle);
                }

                euclidean_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( euclidean_distance, cosine_of_proj_angle)
                );

                sequence_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( sequence_distance, cosine_of_proj_angle)
                );

                euclidean_over_sequence_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( euclidean_distance / std::log( sequence_distance), cosine_of_proj_angle)
                );
              }
              else
              {
                // add data to table
                const std::string sse_str( "sse" + util::Format()( sse_a_number) + "_sse" + util::Format()( sse_b_number));
                loop_angle_table.InsertRow
                (
                  model_filename + "_" + sse_str,
                  storage::Vector< double>::Create
                  (
                    sequence_distance,
                    euclidean_distance,
                    score::Loop::NormalizeDistance( storage::Pair< size_t, double>( sequence_distance, euclidean_distance)),
                    sse_a_b_consecutive,
                    cosine_of_proj_angle
                  ),
                  true
                );
              }

              // only perform if path was set from commandline
              if( m_VisualizationFlag)
              {
                // add atom lines for center_of_mass, prev_end_of_z, next_begin_of_z
                biol::Atom atom_center_of_mass( center_of_mass, biol::GetAtomTypes().CA);
                biol::Atom atom_prev_end_of_z( prev_end_of_z, biol::GetAtomTypes().CA);
                biol::Atom atom_next_begin_of_z( next_begin_of_z, biol::GetAtomTypes().CA);
                biol::AA amino_acid;

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_center_of_mass, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_prev_end_of_z, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_next_begin_of_z, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

              } // end if

            } // for sse_itr_b
          } // for sse_itr_a

          // Write the file for each protein in the ensemble
          if( m_VisualizationFlag)
          {
            std::string vis_filename( model_filename + GetOutFilePostfix() + ".pdb");

            // create handler and add lines
            pdb::Handler pdb_handler;
            pdb_handler.AppendLines( vis_file_lines);

            // write visualization pdb
            io::OFStream pdb_write_stream;
            io::File::MustOpenOFStream( pdb_write_stream, vis_filename);
            pdb_handler.WriteLines( pdb_write_stream);
            io::File::CloseClearFStream( pdb_write_stream);
          }

        } // for chain
      }// End Ensemble iteration

      std::stringstream ostream;
      if( m_OutputOption == e_Table)
      {
        loop_angle_table.WriteFormatted( ostream);
      }
      else
      {
        // normalization for histogram outputs
        switch( m_OutputOption)
        {
          case e_NormalizedByCount:
            euclidean_distance_cos_angle_histogram.Normalize();
            sequence_distance_cos_angle_histogram.Normalize();
            euclidean_over_sequence_distance_cos_angle_histogram.Normalize();
            break;
          case e_NormalizedByColumn:
            euclidean_distance_cos_angle_histogram.NormalizeY();
            sequence_distance_cos_angle_histogram.NormalizeY();
            euclidean_over_sequence_distance_cos_angle_histogram.NormalizeY();
            break;
          case e_Log:
            euclidean_distance_cos_angle_histogram.Log();
            sequence_distance_cos_angle_histogram.Log();
            euclidean_over_sequence_distance_cos_angle_histogram.Log();
            break;
          default:
            break;
        } // end normalization
        ostream << cos_angle_short_loops_histogram << '\n'
                << cos_angle_long_loops_histogram << '\n'
                << euclidean_distance_cos_angle_histogram << '\n'
                << sequence_distance_cos_angle_histogram << '\n'
                << euclidean_over_sequence_distance_cos_angle_histogram;
      }
      return ostream.str();
    } // end of operator()

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopAngle::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes loop distance statistics."
      );

      parameters.AddInitializer
      (
        "output option",
        "the form of output to create during analysis, Normalized by count, Normalized by Column, Log, or Raw",
        io::Serialization::GetAgent( &m_OutputOption),
        "NormalizedByTotalCount"
      );

      parameters.AddInitializer
      (
        "chains",
        "a string of chains to use for the analysis",
        io::Serialization::GetAgent( &m_Chains),
        "A"
      );

      parameters.AddInitializer
      (
        "visualize",
        "set path to visualize the model in pymol",
        io::Serialization::GetAgent( &m_VisualizationFlag),
        "False"
      );

      return parameters;
    } // end of GetSerializer function
  } // namespace scorestat
} // namespace bcl

