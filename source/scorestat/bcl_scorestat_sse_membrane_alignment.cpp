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
#include "scorestat/bcl_scorestat_sse_membrane_alignment.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "biol/bcl_biol_environment_types.h"
#include "biol/bcl_biol_membrane.h"
#include "biol/bcl_biol_ss_types.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_sse_membrane_alignment.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEMembraneAlignment::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSEMembraneAlignment())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSEMembraneAlignment::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< SSEMembraneAlignment::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSEMembraneAlignment::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "sse_membrane_alignment.histograms",
        GetStaticClassName< SSEMembraneAlignment::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEMembraneAlignment::SSEMembraneAlignment() :
      m_OutputOption( e_Histogram),
      m_NumberOfBins( 9),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSEMembraneAlignment *SSEMembraneAlignment::Clone() const
    {
      return new SSEMembraneAlignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSEMembraneAlignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSEMembraneAlignment::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the number of bins for the histogram
    //! @return the number of bins for the histogram
    const size_t &SSEMembraneAlignment::GetNumberOfBins() const
    {
      return m_NumberOfBins;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSEMembraneAlignment::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSEMembraneAlignment::GetAlias() const
    {
      static std::string s_name( "SSEMembraneAlignment");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSEMembraneAlignment::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics of sse membrane alignment
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::Histogram> > sse_membrane_alignment;
      storage::Map< biol::EnvironmentType, storage::Vector< math::Histogram> > strand_axis_alignment;

      // iterate over SSTypes of interest
      for
      (
        biol::SSTypes::const_iterator
          ss_type_itr( biol::GetSSTypes().Begin()),
          ss_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
        ss_type_itr != ss_type_itr_end;
        ++ss_type_itr
      )
      {
        // iterate over EnvironmentTypes
        for
        (
          biol::EnvironmentTypes::const_iterator
            env_type_itr( biol::GetEnvironmentTypes().Begin()), env_type_itr_end( biol::GetEnvironmentTypes().End());
          env_type_itr != env_type_itr_end;
          ++env_type_itr
        )
        {
          // initialize sse_membrane_alignment
          sse_membrane_alignment[ *ss_type_itr][ *env_type_itr] = math::Histogram
          (
            0,
            math::g_Pi / ( m_NumberOfBins * 2),
            m_NumberOfBins
          );

          // if the SSEType is Strand
          if( *ss_type_itr == biol::GetSSTypes().STRAND)
          {
            // initialize strand_axis_alignment
            strand_axis_alignment[ *env_type_itr] = storage::Vector< math::Histogram>
            (
              coord::GetAxes().GetEnumCount(),
              math::Histogram( 0, math::g_Pi / ( m_NumberOfBins * 2), m_NumberOfBins)
            );
          }
        } // end of iteration over EnvironmentTypes
      } // end of iteration over SSTypes

      // iterate over the protein ensemble
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

        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // get the membrane associated with the current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains in current protein model
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
            BCL_MessageStd( "Skip chain " + util::Format()( ( *chain_itr)->GetChainID()) + " in " + model_name);
            continue;
          }

          // get all sses in current chain
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
            // get all sse fragments of current sse
            const util::SiPtrVector< const assemble::SSEGeometryInterface> &all_fragments( ( *sse_itr)->GetSSEGeometries());

            // iterate over fragments of current sse
            for
            (
              util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
                fragment_itr( all_fragments.Begin()), fragment_itr_end( all_fragments.End());
              fragment_itr != fragment_itr_end;
              ++fragment_itr
            )
            {
              // get current sse
              const assemble::SSEGeometryInterface &current_sse( **fragment_itr);

              // determine the environment type in which the current sse can be found
              // get environment type of current sse
              const biol::EnvironmentType environment_type
              (
                sp_membrane.IsDefined() ?
                    sp_membrane->DetermineEnvironmentType( current_sse.GetCenter()) :
                      biol::GetEnvironmentTypes().e_Solution
              );

              if( !environment_type.IsDefined())
              {
                BCL_MessageCrt
                (
                  "environment for this sse is not defined " + current_sse.GetIdentification()
                );
                continue;
              }

              double align_weight( 1.0);
              score::SSEMembraneAlignment score_sse_membrane_alignment;

              // if the current sse is a strand
              if( current_sse.GetType() == biol::GetSSTypes().STRAND)
              {
                // change the align_weight for strand
                align_weight = score_sse_membrane_alignment.WeightXAxis( current_sse, *sp_membrane);

                // for strands the axes are important
                for
                (
                  coord::Axis::const_iterator
                    axis_itr( coord::GetAxes().Begin()), axis_itr_end( coord::GetAxes().End());
                  axis_itr != axis_itr_end;
                  ++axis_itr
                )
                {
                  strand_axis_alignment[ environment_type]( *axis_itr).PushBack
                  (
                    score_sse_membrane_alignment.AngleToMembranePlane( current_sse, *axis_itr, *sp_membrane)
                  );
                }
              }

              // insert the angle into the map with appropriate weight
              sse_membrane_alignment[ current_sse.GetType()][ environment_type].PushBack
              (
                score_sse_membrane_alignment.AngleToMembranePlane( current_sse, coord::GetAxes().e_Z, *sp_membrane),
                align_weight
              );
            } // end of itration over fragments of sse
          } // end of iteration over all sses in the current chain
        } // end of iteration over all chains in the current protein model
      } // end of iteration over the protein ensemble

      // write out the statistics
      std::stringstream stream;
      io::OFStream write;

      if( m_OutputOption == e_Histogram)
      {
        // write output files
        io::File::MustOpenOFStream( write, "sse_membrane_alignment.histograms");
        write << sse_membrane_alignment;
        io::File::CloseClearFStream( write);

        io::File::MustOpenOFStream( write, "strand_angle_membrane.histograms");
        write << strand_axis_alignment;
        io::File::CloseClearFStream( write);

        // print instructions
        stream << "Refer to the following two files: \n";
        stream << "sse_membrane_alignment.histograms" << '\n';
        stream << "strand_angle_membrane.histograms";
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEMembraneAlignment::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes secondary structure element (SSE) membrane alignment statistics."
      );

      parameters.AddInitializer
      (
        "number_bins",
        "number of bins for the histogram",
        io::Serialization::GetAgent( &m_NumberOfBins),
        "9"
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

