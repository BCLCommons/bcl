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
#include "fold/bcl_fold_loop_segment_sequence_order.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_compare.h"
#include "fold/bcl_fold_locator_loop_segment.h"
#include "fold/bcl_fold_loop_segment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //! @brief return true if SEGMENT_A comes before SEGMENT_B in sequence
    //! @param SEGMENT_A first LoopSegment
    //! @param SEGMENT_B second LoopSegment
    //! @return true if SEGMENT_A comes before SEGMENT_B in sequence
    bool LoopSegmentSequenceOrder::operator()
    (
      const LoopSegment &SEGMENT_A, const LoopSegment &SEGMENT_B
    ) const
    {
      // make sure the chain ids match
      BCL_Assert
      (
        SEGMENT_A.GetSSE()->GetChainID() == SEGMENT_B.GetSSE()->GetChainID(), "chain id of segments in are different"
      );

      // determine which segment comes before the other in sequence
      return assemble::SSELessThanNoOverlap().operator()( SEGMENT_A.GetSSE(), SEGMENT_B.GetSSE());
    }

    //! @brief return true if SEGMENT_LOCATOR_A comes before SEGMENT_LOCATOR_B in sequence
    //! @param SEGMENT_LOCATOR_A first LoopSegment
    //! @param SEGMENT_LOCATOR_B second LoopSegment
    //! @return true if SEGMENT_LOCATOR_A comes before SEGMENT_LOCATOR_B in sequence
    bool LoopSegmentSequenceOrder::operator()
    (
      const LocatorLoopSegment &SEGMENT_LOCATOR_A, const LocatorLoopSegment &SEGMENT_LOCATOR_B
    ) const
    {
      // make sure that the chain ids match
      BCL_Assert
      (
        SEGMENT_LOCATOR_A.GetLocatorSSE().GetChainID() == SEGMENT_LOCATOR_B.GetLocatorSSE().GetChainID(),
        "chain id of segment locators in are different"
      );

      // get the start and end seq ids of the first locator
      const storage::VectorND< 2, int> segment_a_start_end( SEGMENT_LOCATOR_A.GetLocatorSSE().GetSSEID());

      // get the start seq id of the second locator
      const int segment_b_start( SEGMENT_LOCATOR_B.GetLocatorSSE().GetSSEID().First());

      // true if the start and end seq id's of the firsts loop segment sse are less than "segment_b_start"
      if
      (
        ( segment_a_start_end.First() < segment_b_start)
        &&
        ( segment_a_start_end.Second() < segment_b_start)
      )
      {
        // return true since "SEGMENT_LOCATOR_A" is less than "SEGMENT_LOCATOR_B"
        return true;
      }

      // else return larger
      return false;
    }

  } // namespace fold
} // namespace bcl
