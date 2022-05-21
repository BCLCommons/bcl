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

// include header of this class
#include "fold/bcl_fold_locator_unconnected_segments.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> LocatorUnconnectedSegments::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorUnconnectedSegments())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorUnconnectedSegments::LocatorUnconnectedSegments()
    {
    }

    //! @brief clone function
    //! @return pointer to a new LocatorUnconnectedSegments
    LocatorUnconnectedSegments *LocatorUnconnectedSegments::Clone() const
    {
      return new LocatorUnconnectedSegments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LocatorUnconnectedSegments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates unconnected segments in a protein model and computes a parameterization of the segment
    //! @detail two segments are considered unconnected if there is at least one residue without coordinates in
    //! between
    //! @param MODEL protein model for which to find unconnected segments
    //! @return parameterizations of unconnected segments in the given protein model
    util::ShPtrVector< LoopParameters> LocatorUnconnectedSegments::Locate
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // get all loops in the model and filter out the terminal ones and fully defined ones
      const util::SiPtrVector< const assemble::SSE> loops( MODEL.GetSSEs( biol::GetSSTypes().COIL));
      util::SiPtrVector< const assemble::SSE> loops_filtered;
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        const assemble::SSE &sse( **loop_it);
        const int num_residues( MODEL.GetChain( sse.GetChainID())->GetNumberAAs());
        if( sse.GetFirstAA()->GetSeqID() != 1 && sse.GetLastAA()->GetSeqID() != num_residues)
        {
          // check if all coordinates in this loop are defined
          const util::ShPtrVector< biol::AABase> &residues( sse.GetData());
          bool partially_undefined( false);
          for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
          {
            const biol::AABase &current_res( **res_it);
            if( !current_res.HasDefinedCoordinates())
            {
              partially_undefined = true;
              break;
            }
          }
          if( partially_undefined)
          {
            loops_filtered.PushBack( *loop_it);
          }
        }
      }

      // compute the parameterizations for all remaining loops
      util::ShPtrVector< LoopParameters> segments;
      for
      (
        auto loop_it( loops_filtered.Begin()), loop_it_end( loops_filtered.End());
        loop_it != loop_it_end;
        ++loop_it
      )
      {
        // get the anchor residues for this loop
        const int anchor_1_seq_id( ( **loop_it).GetFirstAA()->GetSeqID() - 1);
        const int anchor_2_seq_id( ( **loop_it).GetLastAA()->GetSeqID() + 1);
        const char chain_id( ( **loop_it).GetChainID());
        const util::SiPtrVector< const biol::AABase> &residues( MODEL.GetChain( chain_id)->GetAminoAcids());
        const biol::AABase &anchor_1( *residues( anchor_1_seq_id - 1));
        const biol::AABase &anchor_2( *residues( anchor_2_seq_id - 1));

        // compute the parameterization of the loop
        if( anchor_2.GetSeqID() - anchor_1.GetSeqID() > 0)
        {
          segments.PushBack( LoopParameters::Create( anchor_1, anchor_2));
        }
      }

      return segments;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &LocatorUnconnectedSegments::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &LocatorUnconnectedSegments::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
