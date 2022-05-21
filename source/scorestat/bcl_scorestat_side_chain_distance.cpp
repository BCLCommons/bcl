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
#include "scorestat/bcl_scorestat_side_chain_distance.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
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
    const util::SiPtr< const util::ObjectInterface> SideChainDistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SideChainDistance())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SideChainDistance::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< SideChainDistance::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SideChainDistance::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "side_chain_distance.histogram",
        GetStaticClassName< SideChainDistance::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SideChainDistance::SideChainDistance() :
      m_OutputOption( e_Histogram),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SideChainDistance *SideChainDistance::Clone() const
    {
      return new SideChainDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SideChainDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SideChainDistance::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SideChainDistance::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SideChainDistance::GetAlias() const
    {
      static std::string s_name( "SideChainDistance");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SideChainDistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for distance of side chain center of mass to sse main axis
      storage::Map< biol::SSType, storage::Map< biol::AAType, math::Histogram> > side_chain_distance;

      // initialize side_chain_distance
      for
      (
        biol::SSTypes::const_iterator
          sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
        sse_type_itr != sse_type_itr_end;
        ++sse_type_itr
      )
      {
        // iterate over aa types
        for
        (
          biol::AATypes::const_iterator
            aa_type_itr( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().End());
          aa_type_itr != aa_type_itr_end;
          ++aa_type_itr
        )
        {
          side_chain_distance[ *sse_type_itr][ *aa_type_itr] = math::Histogram( 0.0, 0.1, 100);
        }
      }

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

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

        // get all chains of current protein model
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
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get helices and strands the model
          const util::SiPtrVector< const assemble::SSE> &all_sses
          (
            current_protein_model.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
          );

          // the pdb file might be missing sse entries
          if( all_sses.IsEmpty())
          {
            BCL_MessageStd( "Warning: " + model_filename + " probably does not have sse entries in the pdb file, skipping");
            continue;
          }

          // iterate over all sses
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip sses that are shorter than minimum fragment length
            if( ( *sse_itr)->GetSize() <= ( *sse_itr)->GetType()->GetFragmentLength())
            {
              continue;
            }

            // create an iterator to the central amino acid residue
            biol::AASequence::const_iterator center_aa_itr
            (
              ( *sse_itr)->Begin() + ( ( ( *sse_itr)->GetType()->GetFragmentLength()) - 1) / 2
            );

            // create an iterator to the end amino acid residue
            biol::AASequence::const_iterator end_aa_itr( ( *sse_itr)->GetData().End());

            // iterate over fragments
            for
            (
              util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
                frag_itr( ( *sse_itr)->GetFragments().Begin()), frag_itr_end( ( *sse_itr)->GetFragments().End());
              frag_itr != frag_itr_end && center_aa_itr != end_aa_itr;
              ++frag_itr, ++center_aa_itr
            )
            {
              // get the center of mass of the center amino acid residue
              const linal::Vector3D com_center_aa( ( *center_aa_itr)->CalculateCenterOfMassOfSideChain());

              // get the line segment of current fragment
              const coord::LineSegment3D lineseg_from_fragment( ( *frag_itr)->GetMainAxis());

              // get distance from center aa to sse fragment
              const double distance
              (
                coord::CalculateDistancePointFromLineSegment( lineseg_from_fragment, com_center_aa).First()
              );

              // insert
              side_chain_distance[ ( *sse_itr)->GetType()][ ( *center_aa_itr)->GetType()].PushBack( distance);
            }
          }
        } // end of iteration over all chains
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Histogram)
      {
        // iterate over sse type
        for
        (
          biol::SSTypes::const_iterator
            sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
          sse_type_itr != sse_type_itr_end;
          ++sse_type_itr
        )
        {
          // write the name of sse type
          stream << ( *sse_type_itr)->GetName() << '\n';

          // iterate over aa types
          for
          (
            biol::AATypes::const_iterator
              aa_type_itr( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().End());
            aa_type_itr != aa_type_itr_end;
            ++aa_type_itr
          )
          {
            // write name of amino acid residue
            stream << ( *aa_type_itr)->GetName() << '\n';

            // write histogram
            stream << side_chain_distance[ *sse_type_itr][ *aa_type_itr];
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
    io::Serializer SideChainDistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics for side chain distance from main SSE axis."
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

