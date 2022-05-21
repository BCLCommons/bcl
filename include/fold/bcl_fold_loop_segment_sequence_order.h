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

#ifndef BCL_FOLD_LOOP_SEGMENT_SEQUENCE_ORDER_H_
#define BCL_FOLD_LOOP_SEGMENT_SEQUENCE_ORDER_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopSegmentSequenceOrder
    //! @brief This is a function class for comparing two loop segments so that the first one in terms of sequence
    //!        is determined to be less than the second one in terms of sequence
    //!
    //! @see @link example_fold_loop_segment_sequence_order.cpp @endlink
    //! @author alexanns
    //! @date August 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopSegmentSequenceOrder
    {

    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SEGMENT_A comes before SEGMENT_B in sequence
      //! @param SEGMENT_A first LoopSegment
      //! @param SEGMENT_B second LoopSegment
      //! @return true if SEGMENT_A comes before SEGMENT_B in sequence
      bool operator()( const LoopSegment &SEGMENT_A, const LoopSegment &SEGMENT_B) const;

      //! @brief return true if SEGMENT_LOCATOR_A comes before SEGMENT_LOCATOR_B in sequence
      //! @param SEGMENT_LOCATOR_A first LoopSegment
      //! @param SEGMENT_LOCATOR_B second LoopSegment
      //! @return true if SEGMENT_LOCATOR_A comes before SEGMENT_LOCATOR_B in sequence
      bool operator()
      (
        const LocatorLoopSegment &SEGMENT_LOCATOR_A, const LocatorLoopSegment &SEGMENT_LOCATOR_B
      ) const;

    }; // class LoopSegmentSequenceOrder

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_SEGMENT_SEQUENCE_ORDER_H_ 
