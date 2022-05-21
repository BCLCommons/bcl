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
#include "scorestat/bcl_scorestat_sse_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "score/bcl_score_loop.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPacking::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSEPacking())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSEPacking::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        "Histogram2D",
        GetStaticClassName< SSEPacking::OutputOption>()
      };

      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSEPacking::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_file_extentions[] =
      {
        "sse_packing.tbl",
        "sse_packing.histograms2D",
        GetStaticClassName< SSEPacking::OutputOption>()
      };

      return s_output_file_extentions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPacking::SSEPacking() :
      m_OutputOption( e_Table),
      m_SSEDistanceBinSize( 1.0),
      m_SSEAngleNumberBins( 24),
      m_StrandDistanceBinSize( 0.25),
      m_StrandAngleNumberBins( 24),
      m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    SSEPacking *SSEPacking::Clone() const
    {
      return new SSEPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSEPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSEPacking::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the sse distance bin size for the 2D histogram
    //! @return the sse distance bin size for the 2D histogram
    const double &SSEPacking::GetSSEDistanceBinSize() const
    {
      return m_SSEDistanceBinSize;
    }

    //! @brief returns the sse twist angle bin size for the 2D histogram
    //! @return the sse twist angle bin size for the 2D histogram
    const size_t &SSEPacking::GetSSEAngleNumberBins() const
    {
      return m_SSEAngleNumberBins;
    }

    //! @brief returns the strand distance bin size for the 2D histogram
    //! @return the strand distance bin size for the 2D histogram
    const double &SSEPacking::GetStrandDistanceBinSize() const
    {
      return m_StrandDistanceBinSize;
    }

    //! @brief returns the strand twist angle bin size for the 2D histogram
    //! @return the strand twist angle bin size for the 2D histogram
    const size_t &SSEPacking::GetStrandAngleNumberBins() const
    {
      return m_StrandAngleNumberBins;
    }

    //! @brief returns the fragment minimum interface length
    //! @return the fragment minimum interface length
    const double &SSEPacking::GetFragmentMinInterfaceLength() const
    {
      return m_FragmentMinInterfaceLength;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSEPacking::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSEPacking::GetAlias() const
    {
      static const std::string s_name( "SSEPacking");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSEPacking::operator ()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // initialize a table for storing sse packing angles and distances
      storage::Table< double> sse_packing_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sse_packing_type", "twist_angle", "shortest_distance"
          )
        )
      );

      // initialize a vector of 2D histograms for storing sse packing angles and distances
      storage::Vector< math::Histogram2D> sse_packing_angle_distance
      (
        contact::Types::s_NumberValidTypes,
        math::Histogram2D
        (
          storage::VectorND< 2, double>( -math::g_Pi, 0),
          storage::VectorND< 2, double>( 2 * math::g_Pi / m_SSEAngleNumberBins, m_SSEDistanceBinSize),
          storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_SSEDistanceBinSize))
        )
      );

      // initialize a 2D histogram for storing strand-strand packing angle and distance
      math::Histogram2D strand_strand_packing_angle_distance
      (
        storage::VectorND< 2, double>( -math::g_Pi, 0),
        storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_StrandDistanceBinSize),
        storage::VectorND< 2, size_t>( m_StrandAngleNumberBins, size_t( 20 / m_StrandDistanceBinSize))
      );

      // initialize a vector of 2D histograms for storing sse fragment packing angles
      storage::Vector< math::Histogram2D> sse_fragment_angle_distance
      (
        contact::Types::s_NumberValidTypes,
        math::Histogram2D
        (
          storage::VectorND< 2, double>( -math::g_Pi, 0),
          storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_SSEDistanceBinSize),
          storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_SSEDistanceBinSize))
        )
      );

      // initialize a 2D histogram for storing strand-strand fragment packing angle and distance
      math::Histogram2D strand_strand_fragment_angle_distance
      (
        storage::VectorND< 2, double>( -math::g_Pi, 0),
        storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_StrandDistanceBinSize),
        storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_StrandDistanceBinSize))
      );

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
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // get all sses in the current model
        const util::SiPtrVector< const assemble::SSE> all_sses( protein_model.GetSSEs());

        // iterate over all sses
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_a_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
          sse_a_itr != sse_itr_end;
          ++sse_a_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *sse_a_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // skip coils
          if( ( *sse_a_itr)->GetType() == biol::GetSSTypes().COIL)
          {
            continue;
          }

          // get the second sse
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_b_itr( sse_a_itr + 1);
             sse_b_itr != sse_itr_end;
            ++sse_b_itr
          )
          {
            // skip undesired chains
            if( !m_ChainIds.empty() && m_ChainIds.find( ( *sse_b_itr)->GetChainID()) == std::string::npos)
            {
              continue;
            }

            // skip coils
            if( ( *sse_b_itr)->GetType() == biol::GetSSTypes().COIL)
            {
              continue;
            }

            // skip sse pairs that are possibly from a single broken sse
            if( biol::CalculateSequenceDistance( **sse_a_itr, **sse_b_itr) <= 1)
            {
              continue;
            }

            // create an object of SSEGeometryPacking from sse_a and sse_b
            const assemble::SSEGeometryPacking sse_pack( **sse_a_itr, **sse_b_itr);

            // collect twist angle and distance and store them into corresponding histograms
            if( m_OutputOption == e_Histogram)
            {
              SSEPackingAngleDistance( sse_pack, sse_packing_angle_distance, strand_strand_packing_angle_distance);
            }
            else if( m_OutputOption == e_Table)
            {
              // collect twist angle and distance and store them into corresponding table
              sse_packing_table.InsertRow
              (
                model_name + ": " + ( *sse_a_itr)->GetChainID() + "_" + ( *sse_b_itr)->GetChainID(),
                storage::Vector< double>::Create
                (
                  sse_pack.GetContactType(),
                  sse_pack.GetTwistAngle(),
                  sse_pack.GetDistance()
                ),
                true
              );
            }

            // get the fragments for sse_a and sse_b
            storage::Vector< storage::List< assemble::SSEGeometryPacking> > fragment_packing_lists
            (
              assemble::SSEGeometryPackerAllFragmentPairs( m_FragmentMinInterfaceLength).operator ()
              (
                **sse_a_itr, **sse_b_itr
              )
            );

            // iterate over a list of fragment packing list
            for
            (
              storage::Vector< storage::List< assemble::SSEGeometryPacking> >::const_iterator
                fragment_packing_list_itr( fragment_packing_lists.Begin()), fragment_packing_list_itr_end( fragment_packing_lists.End());
              fragment_packing_list_itr != fragment_packing_list_itr_end;
              ++fragment_packing_list_itr
            )
            {
              // iterate over fragment packing list
              for
              (
                storage::List< assemble::SSEGeometryPacking>::const_iterator
                  fragment_packing_itr( fragment_packing_list_itr->Begin()), fragment_packing_itr_end( fragment_packing_list_itr->End());
                fragment_packing_itr != fragment_packing_itr_end;
                ++fragment_packing_itr
              )
              {
                if( m_OutputOption == e_Histogram)
                {
                  // collect twist angle and distance and store them into corresponding 2D histograms
                  SSEPackingAngleDistance( *fragment_packing_itr, sse_fragment_angle_distance, strand_strand_fragment_angle_distance);
                }
                else if( m_OutputOption == e_Table)
                {
                  // collect twist angle and distance and store them into corresponding table
                  sse_packing_table.InsertRow
                  (
                    model_name + ": " + ( *sse_a_itr)->GetChainID() + "_" + ( *sse_b_itr)->GetChainID(),
                    storage::Vector< double>::Create
                    (
                      ( *fragment_packing_itr).GetContactType(),
                      ( *fragment_packing_itr).GetTwistAngle(),
                      ( *fragment_packing_itr).GetDistance()
                    ),
                    true
                  );
                }
              }
            }
          } // end of iterating over the second sse
        } // end of iterating over all sses
      } // end of iterating over protein ensemble

      // write to files
      std::stringstream stream;
      io::OFStream write;

      if( m_OutputOption == e_Histogram)
      {
        // write sse_packing statistics
        io::File::MustOpenOFStream( write, "sse_angle_distance.histograms2D");
        write << assemble::SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength() << '\n';
        write << contact::GetTypes().HELIX_HELIX << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_HELIX);
        write << contact::GetTypes().HELIX_SHEET << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_SHEET);
        write << contact::GetTypes().HELIX_STRAND << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_STRAND);
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().STRAND_STRAND);
        write << contact::GetTypes().SHEET_SHEET << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().SHEET_SHEET);
        io::File::CloseClearFStream( write);

        // write sse_angle_distance_strand_strand
        io::File::MustOpenOFStream( write, "strand_angle_distance.histograms2D");
        write << assemble::SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength() << '\n';
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << strand_strand_packing_angle_distance << '\n';
        io::File::CloseClearFStream( write);

        // write sse_fragment statistics
        io::File::MustOpenOFStream( write, "sse_fragment_angle_distance.histograms2D");
        write << contact::GetTypes().HELIX_HELIX << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_HELIX) << "\n";
        write << contact::GetTypes().HELIX_SHEET << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_SHEET) << "\n";
        write << contact::GetTypes().HELIX_STRAND << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_STRAND) << "\n";
        write << contact::GetTypes().STRAND_STRAND << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().STRAND_STRAND) << "\n";
        write << contact::GetTypes().SHEET_SHEET << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().SHEET_SHEET) << "\n";
        io::File::CloseClearFStream( write);

        // write strand_fragment packing statistics
        io::File::MustOpenOFStream( write, "strand_fragment_angle_distance.histograms2D");
        write << m_SSEDistanceBinSize << '\n';
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << strand_strand_fragment_angle_distance;
        io::File::CloseClearFStream( write);

        // print instructions
        stream << "Refer to the following four files: " << '\n';
        stream << "sse_angle_distance.histograms2D" << '\n';
        stream << "strand_angle_distance.histograms2D" << '\n';
        stream << "sse_fragment_angle_distance.histograms2D" << '\n';
        stream << "strand_fragment_angle_distance.histograms2D";
      }
      else if( m_OutputOption == e_Table)
      {
        // write everything into one file
        sse_packing_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPacking::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes sse packing statistics."
      );

      parameters.AddInitializer
      (
        "sse_distance_bin_size",
        "the bin size of distance for the 2D histogram",
        io::Serialization::GetAgent( &m_SSEDistanceBinSize),
        "1.0"
      );

      parameters.AddInitializer
      (
        "sse_angle_number_bins",
        "number of bins for sse distance",
        io::Serialization::GetAgent( &m_SSEAngleNumberBins),
        "24"
      );

      parameters.AddInitializer
      (
        "strand_distance_bin_size",
        "the bin size of distance for the 2D histogram",
        io::Serialization::GetAgent( &m_StrandDistanceBinSize),
        "0.25"
      );

      parameters.AddInitializer
      (
        "strand_angle_number_bins",
        "the bin size of angle for the 2D histogram",
        io::Serialization::GetAgent( &m_StrandAngleNumberBins),
        "24"
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

    //! @brief collects twist angle and shortest distance for sse packing types
    //! @param SSE_PACK the sse packing type
    //!        SSEPACKING_ANGLE_DISTANCE a vector of 2D histograms that store twist angle and shortest distance
    //!        of sse packing types
    //!        STRAND_STRAND_ANGLE_DISTANCE a 2D histogram that stores the twist angle and shortest distance of
    //!        strand_strand packing
    //! @return a vector of 2D histograms that store twist angle and shortest distance of sse packing types
    storage::Vector< math::Histogram2D> &SSEPacking::SSEPackingAngleDistance
    (
      const assemble::SSEGeometryPacking &SSE_PACK,
      storage::Vector< math::Histogram2D> &SSEPACKING_ANGLE_DISTANCE,
      math::Histogram2D &STRAND_STRAND_ANGLE_DISTANCE
    ) const
    {
      // make pair of twist angle and shortest distance
      storage::VectorND< 2, double> angle_distance( SSE_PACK.GetTwistAngle(), SSE_PACK.GetDistance());

      // get packing type
      contact::Type contact_type( SSE_PACK.GetContactType());

      // consider all sorts of packing types, stores twist angle and shortest distance into 2D histogram for each sse packing type
      // strand-only packing
      if
      (
        contact_type == contact::GetTypes().STRAND_STRAND
        || contact_type == contact::GetTypes().SHEET_SHEET
        || contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().STRAND_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetStrandStrandPairingWeight()
        );
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().SHEET_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
        STRAND_STRAND_ANGLE_DISTANCE.PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetStrandStrandPairingWeight()
        );
      }
      // helix-only packing
      else if( contact_type == contact::GetTypes().HELIX_HELIX)
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_HELIX).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
      }
      // mixed packing
      else if
      (
        contact_type == contact::GetTypes().HELIX_STRAND
        || contact_type == contact::GetTypes().STRAND_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * ( 1 - SSE_PACK.GetRelativePositionWeight())
        );
      }
      else if
      (
        contact_type == contact::GetTypes().HELIX_SHEET
        || contact_type == contact::GetTypes().SHEET_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
      }
      // undefined packing
      else if
      (
        contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND
        || contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * ( 1 - SSE_PACK.GetRelativePositionWeight())
        );
      }

      // return the vector of 2D histograms
      return SSEPACKING_ANGLE_DISTANCE;
    }
  } // namespace scorestat
} // namespace bcl
