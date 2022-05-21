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

#ifndef BCL_ALIGN_MULTIPLE_ALIGNER_INTERFACE_H_
#define BCL_ALIGN_MULTIPLE_ALIGNER_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically
#include "score/bcl_score.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MultipleAlignerInterface
    //! @brief This is an interface class for defining a method of aligning multiple sequences.
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date Nov 12, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class MultipleAlignerInterface :
      virtual public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      virtual void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does aligns data together
      //! @param ALIGNMENTS is the list of Alignments to be aligned (no const ShPtrList, goes in new alignment as child)
      //! @return returns a pair of the Alignment and a double which is the score
      virtual storage::Pair< AlignmentNode< t_Member>, double> AlignMultiple
      (
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENTS
      ) const = 0;

    }; // template class MultipleAlignerInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_MULTIPLE_ALIGNER_INTERFACE_H_
