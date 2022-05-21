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

#ifndef BCL_ALIGN_ALIGNMENT_INTERFACE_H_
#define BCL_ALIGN_ALIGNMENT_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_assignment.h"
#include "bcl_align_sequence_interface.h"
#include "iterate/bcl_iterate_generic.h"
#include "score/bcl_score_assignment_with_gap.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentInterface
    //! @brief interface for classes storing an alignment
    //! @details each class derived from the interface has to implement the private abstract methods; private methods
    //! can only be called from the interface, thus only from the public interface methods ensuring consistency:
    //! - adding or removing a sequence is not allowed as it makes all assignments invalid (construct a new alignment)
    //! - an assignment to be added has to have the same depth as the alignment
    //! - an assignment to be added must not be empty
    //! - TODO: each element in an assignment has to be from a different sequence; needs check in Assignment
    //!
    //! @tparam t_Member elements of the sequence used in the alignment
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date Sep 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentInterface :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! typedef for const_iterator; allow iterating over assignment, but not changes to them
      typedef typename util::ShPtrList< Assignment< t_Member> >::const_iterator const_iterator;
      //! typedef for const_reverse_iterator
      typedef typename util::ShPtrList< Assignment< t_Member> >::const_reverse_iterator const_reverse_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual empty constructor resetting only the assignment list
      //! @return pointer to new AlignmentInterface< t_Member>
      virtual AlignmentInterface< t_Member> *Empty() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the sequences the alignment is build from; this is a single sequence for leaf alignments;
      //! doesn't return reference as nodes assemble it which would make them return a reference to a temporary variable
      //! @return a list of ShPtr to SequenceInterface
      virtual util::ShPtrList< SequenceInterface< t_Member> > GetSequences() const = 0;

      //! @brief get the list of child alignments
      //! @return list of ShPtr to AlignmentInterface
      virtual const util::ShPtrList< AlignmentInterface< t_Member> > &GetChildAlignments() const = 0;

      //! @brief get an iterator to the list of child alignments; can be implemented even without shared ownership
      //! @return iterator to the child alignments
      virtual iterate::Generic< const AlignmentInterface< t_Member> > GetChildAlignmentsIterator() const
      {
        return
          iterate::Generic< const AlignmentInterface< t_Member> >
          (
            GetChildAlignments().Begin(),
            GetChildAlignments().End()
          );
      }

      //! @brief GetDepth gives the number of t_Members are in each Assignment i.e. the number of sequences
      //! @return size_t which is the number of t_Members that are in each of the Assignments
      virtual size_t GetDepth() const = 0;

      //! @brief returns the list of pointers to assignments
      //! @return the list of ShPtr to Assignment
      virtual const util::ShPtrList< Assignment< t_Member> > &GetAssignments() const = 0;

      //! @brief returns size of the assignment container
      //! @return size, i.e. number of assignments stored
      virtual size_t GetSize() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the identifiers of all sequences
      //! @return list of strings for the sequence identifiers
      storage::List< std::string> GetSequenceIds() const
      {
        const util::ShPtrList< SequenceInterface< t_Member> > sequences( GetSequences());
        storage::List< std::string> sequence_ids;

        // iterate over sequences and save sequence ids
        for
        (
          typename util::ShPtrList< SequenceInterface< t_Member> >::const_iterator
            itr( sequences.Begin()),
            itr_end( sequences.End());
          itr != itr_end;
          ++itr
        )
        {
          sequence_ids.Append( ( **itr).GetSequenceId());
        }

        return sequence_ids;
      }

      //! @brief checks whether any assignment are stored
      //! @return if the assignment container is empty or not
      virtual bool IsEmpty() const = 0;

      //! @brief reset the assignment storage container, removing all assignments
      virtual void ResetAssignments() = 0;

      //! @brief Prepend inserts SP_ASSIGNMENT at the front
      //! @param SP_ASSIGNMENT assignment to be inserted
      //! @return if operation was successful
      bool Prepend( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        if( IsValidAssignment( SP_ASSIGNMENT))
        {
          PrependImplementation( SP_ASSIGNMENT); // if assignment has same number sequences as alignment
          return true;
        }
        return false;
      }

      //! @brief Prepend inserts SP_LIST_ASSIGNMENT at the front; preserves ordering: prepending ab to x gives abx
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      //! @return if operation was successful and all assignment were inserted; none are inserted if unsuccessful
      bool Prepend( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        if( IsValidAssignmentList( SP_LIST_ASSIGNMENT))
        {
          PrependImplementation( SP_LIST_ASSIGNMENT);
          return true;
        }
        return false;
      }

      //! @brief Append inserts SP_ASSIGNMENT at the end
      //! @param SP_ASSIGNMENT assignment to be inserted
      //! @return if operation was successful
      bool Append( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        if( IsValidAssignment( SP_ASSIGNMENT))
        {
          AppendImplementation( SP_ASSIGNMENT); // if assignment has same number sequences as alignment
          return true;
        }
        return false;
      }

      //! @brief Append inserts SP_LIST_ASSIGNMENT at the end; preserves ordering: appending ab to x gives xab
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      //! @return if operation was successful and all assignment were inserted; none are inserted if unsuccessful
      bool Append( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        if( IsValidAssignmentList( SP_LIST_ASSIGNMENT))
        {
          AppendImplementation( SP_LIST_ASSIGNMENT);
          return true;
        }
        return false;
      }

      //! @brief scores the alignment with the given assignment scoring function
      //! @param SCORE_ASSIGNMENT the scoring function for scoring assignments
      //! @return the score
      double Score( const score::AssignmentWithGap< t_Member> SCORE_ASSIGNMENT) const
      {
        double score( 0);

        // iterate over all assignments
        typename util::ShPtrList< Assignment< t_Member> >::const_iterator itr( GetAssignments().Begin());
        typename util::ShPtrList< Assignment< t_Member> >::const_iterator itr_end( GetAssignments().End());
        for( ; itr != itr_end; ++itr)
        {
          score += SCORE_ASSIGNMENT( **itr);
        }

        return score;
      }

      //! @brief checks if the given alignments is a complete subalignment of *this alignment
      //! @param SUB_ALIGNMENT alignment that is supposed to be a subalignment of *this
      //! @return bool if SUB_ALIGNMENT is a subalignment
      bool IsSubAlignment( const AlignmentInterface< t_Member> &SUB_ALIGNMENT) const
      {
        // check that the sequences are the same and SUBALIGNMENT is not empty
        if( !( GetSequences() == SUB_ALIGNMENT.GetSequences()) || SUB_ALIGNMENT.IsEmpty())
        {
          return false;
        }

        // all assignments
        const util::ShPtrList< Assignment< t_Member> > &assignments( GetAssignments());
        const util::ShPtrList< Assignment< t_Member> > &sub_assignments( SUB_ALIGNMENT.GetAssignments());

        // iterator on first assignment in the ALIGNMENT and in the SUB_ALIGNMENT
        typename util::ShPtrList< Assignment< t_Member> >::const_iterator
          itr( assignments.Begin()), itr_end( assignments.End()),
          sub_itr( sub_assignments.Begin()), sub_itr_end( sub_assignments.End());

        // find the first agreement by looping over ALIGNMENT and find first assignment of SUB_ALIGNMENT
        while( true)
        {
          if( itr == itr_end) // return false, if reached the end of the alignment without finding
          {
            return false;
          }
          if( ( **itr).IsIdentical( **sub_itr))
          {
            break;
          }

          ++itr; // keep moving forward
        }

        // itr members and sub_itr members do agree at this point; now check all succeeding assignments agree
        while( sub_itr != sub_itr_end)
        {
          // return false, if itr reached the end of ALIGNMENT, but not SUB_ALIGNMENT; or if assignments disagree
          if( itr == itr_end || !( ( **itr).IsIdentical( **sub_itr)))
          {
            return false;
          }

          ++itr, ++sub_itr; // now move both itrs forward
        }

        return true; // if the while loop was exited b/c the end of the subalignment was reached
      }

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief checks if the given assignment is valid for this alignment
      //! @param SP_ASSIGNMENT assignment to check
      //! @return if valid
      bool IsValidAssignment( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT) const
      {
        // no empty assignments can are valid and no assignments with wrong size
        return !SP_ASSIGNMENT->GetMembers().IsEmpty() && GetDepth() == SP_ASSIGNMENT->GetMembers().GetSize();
      }

      //! @brief checks if the given assignment list is valid for this alignment
      //! @param SP_LIST_ASSIGNMENT list of assignments to check
      //! @return if valid
      bool IsValidAssignmentList( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT) const
      {
        // loop over all ShPtr<Assignment>
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator
            itr( SP_LIST_ASSIGNMENT.Begin()),
            itr_end( SP_LIST_ASSIGNMENT.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !IsValidAssignment( *itr))
          {
            return false;
          }
        }

        return true;
      }

      //! @brief private abstract method to insert SP_ASSIGNMENT at the front
      //! @param SP_ASSIGNMENT assignment to be inserted
      virtual void PrependImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT) = 0;

      //! @brief private abstract method to insert SP_LIST_ASSIGNMENT at the front
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      virtual void PrependImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT) = 0;

      //! @brief private abstract method to insert SP_ASSIGNMENT at the end
      //! @param SP_ASSIGNMENT assignment to be inserted
      virtual void AppendImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT) = 0;

      //! @brief private abstract method to insert SP_LIST_ASSIGNMENT at the end
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      virtual void AppendImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT) = 0;

    }; // template class AlignmentInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNMENT_INTERFACE_H_
