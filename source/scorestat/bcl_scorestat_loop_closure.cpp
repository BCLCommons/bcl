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
#include "scorestat/bcl_scorestat_loop_closure.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopClosure::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopClosure())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopClosure::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< LoopClosure::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopClosure::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "loop_closure.histogram",
        GetStaticClassName< LoopClosure::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopClosure::LoopClosure() :
      m_OutputOption( e_Histogram),
      m_DistanceBinSize( 0.1),
      m_MaxDistance( 50.001), //
      m_NumResidues( 30),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    LoopClosure *LoopClosure::Clone() const
    {
      return new LoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &LoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopClosure::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the bin size for the histogram
    //! @return the bin size for the histogram
    const double &LoopClosure::GetDistanceBinSize() const
    {
      return m_DistanceBinSize;
    }

    //! @brief returns the maximum distance
    //! @returns the maximum distance
    const double &LoopClosure::GetMaxDistance() const
    {
      return m_MaxDistance;
    }

    //! @brief returns the number of residues
    //! @returns the number of residues
    const size_t &LoopClosure::GetNumResidues() const
    {
      return m_NumResidues;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &LoopClosure::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopClosure::GetAlias() const
    {
      static std::string s_name( "LoopClosure");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string LoopClosure::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for loop closure
      storage::Vector< math::Histogram> loop_closure_histograms
      (
        m_NumResidues,
        math::Histogram( 0.0, m_DistanceBinSize, size_t( m_MaxDistance / m_DistanceBinSize))
      );

      // iterator over the protein ensemble
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
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

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
            continue;
          }

          // get all sses in current chain
          const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &all_sses( ( *chain_itr)->GetData());

          // iterate over all sses in current chain
          for
          (
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_itr_a( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a
          )
          {
            // get current sse type
            const biol::SSType &current_sse_type( ( *sse_itr_a)->GetType());

            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && current_sse_type == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !( *sse_itr_a)->IsDefined())
            {
              continue;
            }

            // iterate over the second sse
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator sse_itr_b( sse_itr_a);
            ++sse_itr_b;
            for
            (
              ;
              sse_itr_b != sse_itr_end;
              ++sse_itr_b
            )
            {
              // skip undefined sses
              if( !( *sse_itr_b)->IsDefined())
              {
                continue;
              }

              // get the sequence distance
              const size_t sequence_distance( biol::CalculateSequenceDistance( **sse_itr_a, **sse_itr_b));

              // if the distance is larger than needed
              if( sequence_distance > loop_closure_histograms.GetSize() - 1)
              {
                continue;
              }

              // get euclidean distance
              const double euclidean_distance
              (
                biol::Distance
                (
                  ( *sse_itr_a)->GetLastAA()->GetAtom( biol::GetAtomTypes().C),
                  ( *sse_itr_b)->GetFirstAA()->GetAtom( biol::GetAtomTypes().N)
                )
              );

              // skip undefined distance
              if( euclidean_distance == util::GetUndefinedDouble())
              {
                continue;
              }

              // push the distance into histogram
              if( m_OutputOption == e_Histogram)
              {
                loop_closure_histograms( sequence_distance).PushBack( euclidean_distance);
              }

            } // end of iteration over sses
          } // end of iteration over chains
        } // end of iteration over current protein model
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Histogram)
      {
        stream << loop_closure_histograms;
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopClosure::GetSerializer() const
    {
       io::Serializer parameters;
       parameters.SetClassDescription
       (
         "Computes loop closure statistics."
       );

       parameters.AddInitializer
       (
         "bin_size",
         "the bin size for the histogram",
         io::Serialization::GetAgent( &m_DistanceBinSize),
         "0.1"
       );

       parameters.AddInitializer
       (
         "maximum_distance",
         "maximum distance",
         io::Serialization::GetAgent( &m_MaxDistance),
         "50.001"
       );

       parameters.AddInitializer
       (
         "number_residues",
         "number of residues",
         io::Serialization::GetAgent( &m_NumResidues),
         "30"
       );

       parameters.AddInitializer
       (
         "chain_ids",
         "string of chain ids to be analyzed",
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
