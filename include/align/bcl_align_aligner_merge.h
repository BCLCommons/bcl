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

#ifndef BCL_ALIGN_ALIGNER_MERGE_H_
#define BCL_ALIGN_ALIGNER_MERGE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_alignment_node.h"
#include "bcl_align_multiple_aligner_interface.h"
#include "bcl_align_pairwise_aligner_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerMerge
    //! @brief This class creates an alignment by merging all given alignments
    //!
    //! @tparam t_Member type of object that the Assignment stores
    //!
    //! @see @link example_align_aligner_merge.cpp @endlink
    //! @author heinzes1
    //! @date Nov 1, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerMerge :
      public MultipleAlignerInterface< t_Member>,
      public PairwiseAlignerInterface< t_Member>
    {
    private:

    //////////////
    // typedefs //
    //////////////

      //! const iterator on list of assignments
      typedef typename util::ShPtrList< Assignment< t_Member> >::const_iterator sp_assignment_const_itr;

      //! pair of const iterators on lists of assignments
      typedef std::pair< sp_assignment_const_itr, sp_assignment_const_itr> pair_itrs;

    //////////
    // data //
    //////////

      score::AssignmentWithGap< t_Member> m_ScoreAssignment; //!< scoring function object

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AlignerMerge< t_Member>
      AlignerMerge< t_Member> *Clone() const
      {
        return new AlignerMerge< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE)
      {
        m_ScoreAssignment = SCORE;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A, ALIGNMENT_B are the Alignments to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_A,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_B
      ) const
      {
        // add both alignments to an alignment list and use AlignMultiple()
        util::ShPtrList< AlignmentInterface< t_Member> > alignment_list;
        alignment_list.Append( ALIGNMENT_A);
        alignment_list.Append( ALIGNMENT_B);
        return AlignMultiple( alignment_list);
      }

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENTS is the list of Alignments to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignMultiple
      (
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENTS
      ) const
      {
        if( ALIGNMENTS.GetSize() < 2) // no alignment possible
        {
          return storage::Pair< AlignmentNode< t_Member>, double>( AlignmentNode< t_Member>(), util::GetUndefinedDouble());
        }

        // create new alignment to hold the assignments
        AlignmentNode< t_Member> alignment( ALIGNMENTS);

        // preallocate and initialize pair of begin and end iterators and size for each alignment
        std::vector< pair_itrs> itrs;
        std::vector< size_t> depths;
        itrs.reserve( ALIGNMENTS.GetSize());
        depths.reserve( ALIGNMENTS.GetSize());
        for
        (
          typename util::ShPtrList< AlignmentInterface< t_Member> >::const_iterator
            itr( ALIGNMENTS.Begin()), itr_end( ALIGNMENTS.End());
          itr != itr_end;
          ++itr
        )
        {
          const util::ShPtr< AlignmentInterface< t_Member> > &current_alignment( *itr);
          pair_itrs pair_itr( current_alignment->GetAssignments().Begin(), current_alignment->GetAssignments().End());
          itrs.push_back( pair_itr);
          depths.push_back( current_alignment->GetDepth());
        }

        // insert members until the last alignment ends, fill shorter alignments with gaps
        bool done;
        do
        {
          // create new assignment
          util::ShPtr< Assignment< t_Member> > new_assignment( new Assignment< t_Member>());
          // add members to new assignment
          std::vector< size_t>::const_iterator depths_itr( depths.begin());
          for
          (
            typename std::vector< pair_itrs>::iterator itrs_itr( itrs.begin()), itrs_itr_end( itrs.end());
            itrs_itr != itrs_itr_end;
            ++itrs_itr, ++depths_itr
          )
          {
            // if current itr is at the end of the alignment, add gaps, otherwise members
            if( itrs_itr->first == itrs_itr->second)
            {
              new_assignment->Append( util::SiPtrList< const t_Member>( *depths_itr));
            }
            else
            {
              new_assignment->Append( ( *itrs_itr->first)->GetMembers());
            }
          }
          // add assignment to alignment
          alignment.Append( new_assignment);

          // move all itrs one assignment forward in all alignments unless they are the end of their alignment
          for
          (
            typename std::vector< pair_itrs>::iterator itrs_itr( itrs.begin()), itrs_itr_end( itrs.end());
            itrs_itr != itrs_itr_end;
            ++itrs_itr
          )
          {
            if( itrs_itr->first != itrs_itr->second)
            {
              ++( itrs_itr->first);
            }
          }

          // check if all itrs are at the end of their respective alignments
          done = true;
          for
          (
            typename std::vector< pair_itrs>::const_iterator itrs_itr( itrs.begin()), itrs_itr_end( itrs.end());
            itrs_itr != itrs_itr_end;
            ++itrs_itr
          )
          {
            done &= ( itrs_itr->first == itrs_itr->second);
          }
        }
        while( !done);

        return storage::Pair< AlignmentNode< t_Member>, double>( alignment, CalculateScore( alignment));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM; // return the stream
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM; // return the stream
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the score from individual assignments
      //! @param ALIGNMENT is the alignment whose score is calculated
      //! @return returns a double which is the score
      double CalculateScore( const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        double score( 0);

        //iterate over all
        typename AlignmentLeaf< t_Member>::const_iterator assignment_itr( ALIGNMENT.GetAssignments().Begin());
        typename AlignmentLeaf< t_Member>::const_iterator assignment_itr_end( ALIGNMENT.GetAssignments().End());
        for( ; assignment_itr != assignment_itr_end; ++assignment_itr)
        {
          score += m_ScoreAssignment( **assignment_itr);
        }

        return score;
      }

    }; // template class AlignerMerge

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignerMerge< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignerMerge< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_MERGE_H_
