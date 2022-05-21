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

#ifndef BCL_ALIGN_PAIRWISE_ALIGNER_INTERFACE_H_
#define BCL_ALIGN_PAIRWISE_ALIGNER_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically
#include "score/bcl_score.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_node.h" // because Pair< T1, T2>::T1 type needs to be defined
#include "biol/bcl_biol_aa_base.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PairwiseAlignerInterface
    //! @brief This is an interface class for defining a method of aligning a pair of sequences.
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date Nov 12, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class PairwiseAlignerInterface :
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

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A the first alignment to be aligned (no const, it must go in AlignmentNode's m_ChildAlignments)
      //! @param ALIGNMENT_B the second alignment to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      virtual storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_A,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_B
      ) const = 0;

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT the alignment to be aligned
      //! @param ALIGNMENT_LIST the list of alignments ALIGNMENT is aligned with
      //! @return returns a pair of the Alignment and a double which is the score
      virtual storage::List< storage::Pair< AlignmentNode< t_Member>, double> > AlignPairwise
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT,
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENT_LIST
      ) const
      {
        storage::List< storage::Pair< AlignmentNode< t_Member>, double> > result_list;

        for
        (
          typename util::ShPtrList< AlignmentInterface< t_Member> >::iterator
            itr( ALIGNMENT_LIST.Begin()),
            itr_end( ALIGNMENT_LIST.End());
          itr != itr_end;
          ++itr
        )
        {
          result_list.Append( AlignPair( ALIGNMENT, *itr));
        }

        return result_list;
      }

    }; // template class PairwiseAlignerInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_PAIRWISE_ALIGNER_INTERFACE_H_ 
