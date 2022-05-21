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

#ifndef BCL_ALIGN_ALIGNMENT_HIT_H_
#define BCL_ALIGN_ALIGNMENT_HIT_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_word.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentHit
    //! @brief class storing a hit consisting of a pair of alignment words
    //!
    //! @tparam t_Member elements of the sequence used in the alignment
    //!
    //! @see @link example_align_alignment_hit.cpp @endlink
    //! @author heinzes1
    //! @date Apr 12, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentHit :
      public AlignmentInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< AlignmentWord< t_Member> > m_AlignmentA; //!< first alignment word
      typename AlignmentInterface< t_Member>::const_iterator m_ItrAlignmentABegin; //!< itr to begin of first alignment
      typename AlignmentInterface< t_Member>::const_iterator m_ItrAlignmentAEnd; //!< itr to end of first alignment
      util::ShPtr< AlignmentWord< t_Member> > m_AlignmentB; //!< second alignment word
      typename AlignmentInterface< t_Member>::const_iterator m_ItrAlignmentBBegin; //!< itr to begin of second alignment
      typename AlignmentInterface< t_Member>::const_iterator m_ItrAlignmentBEnd; //!< itr to end of second alignment

      util::ShPtrList< AlignmentInterface< t_Member> > m_ChildAlignments; //!< list of pointers to child alignments
      util::ShPtrList< Assignment< t_Member> > m_Assignments; //!< list of pointers to assignments denoted by itrs

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AlignmentHit() :
        m_AlignmentA( new AlignmentWord< t_Member>()),
        m_ItrAlignmentABegin( m_AlignmentA->GetBeginItr()),
        m_ItrAlignmentAEnd( m_AlignmentA->GetEndItr()),
        m_AlignmentB( new AlignmentWord< t_Member>()),
        m_ItrAlignmentBBegin( m_AlignmentB->GetBeginItr()),
        m_ItrAlignmentBEnd( m_AlignmentB->GetEndItr()),
        m_ChildAlignments(),
        m_Assignments()
      {
        UpdateChildAlignments();
        UpdateAssignments();
      }

      //! @brief constructor taking two word alignments
      //! @param ALIGNMENT_A first word alignment
      //! @param ALIGNMENT_B second word alignment
      AlignmentHit
      (
        const util::ShPtr< AlignmentWord< t_Member> > &ALIGNMENT_A,
        const util::ShPtr< AlignmentWord< t_Member> > &ALIGNMENT_B
      ) :
        m_AlignmentA( ALIGNMENT_A),
        m_ItrAlignmentABegin( m_AlignmentA->GetBeginItr()),
        m_ItrAlignmentAEnd( m_AlignmentA->GetEndItr()),
        m_AlignmentB( ALIGNMENT_B),
        m_ItrAlignmentBBegin( m_AlignmentB->GetBeginItr()),
        m_ItrAlignmentBEnd( m_AlignmentB->GetEndItr()),
        m_ChildAlignments(),
        m_Assignments()
      {
        UpdateChildAlignments();
        UpdateAssignments();
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentHit< t_Member>
      AlignmentHit< t_Member> *Clone() const
      {
        return new AlignmentHit< t_Member>( *this);
      }

      //! @brief virtual empty constructor resetting only the assignment list
      //! @return pointer to new AlignmentInterface< t_Member>
      AlignmentInterface< t_Member> *Empty() const
      {
        return new AlignmentHit< t_Member>( m_AlignmentA, m_AlignmentB);
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

      //! @brief returns the sequences the alignment is build from; this is a single sequence for leaf alignments;
      //! doesn't return reference as nodes assemble it which would make them return a reference to a temporary variable
      //! @return a list of ShPtr to SequenceInterface
      util::ShPtrList< SequenceInterface< t_Member> > GetSequences() const
      {
        util::ShPtrList< SequenceInterface< t_Member> > sequence_list;
        sequence_list.Append( m_AlignmentA->GetSequences());
        sequence_list.Append( m_AlignmentB->GetSequences());
        return sequence_list;
      }

      //! @brief get the list of child alignments
      //! @return list of ShPtr to AlignmentInterface
      const util::ShPtrList< AlignmentInterface< t_Member> > &GetChildAlignments() const
      {
        return m_ChildAlignments;
      }

      //! @brief GetDepth gives the number of t_Members are in each Assignment i.e. the number of sequences
      //! @return size_t which is the number of t_Members that are in each of the Assignments
      size_t GetDepth() const
      {
        return m_AlignmentA->GetDepth() + m_AlignmentB->GetDepth();
      }

      //! @brief returns the list of pointers to assignments
      //! @return the list of ShPtr to Assignment
      const util::ShPtrList< Assignment< t_Member> > &GetAssignments() const
      {
        return m_Assignments;
      }

      //! @brief returns size of the assignment container
      //! @return size, i.e. number of assignments stored
      size_t GetSize() const
      {
        return m_Assignments.GetSize();
      }

      //! @brief returns the itr of the begin of the subalignment in the first child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetBeginItrA() const
      {
        return m_ItrAlignmentABegin;
      }

      //! @brief returns the itr of the begin of the subalignment in the second child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetBeginItrB() const
      {
        return m_ItrAlignmentBBegin;
      }

      //! @brief returns the itr of the end of the subalignment in the first child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetEndItrA() const
      {
        return m_ItrAlignmentAEnd;
      }

      //! @brief returns the itr of the end of the subalignment in the second child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetEndItrB() const
      {
        return m_ItrAlignmentBEnd;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether any assignment are stored
      //! @return if the assignment container is empty or not
      bool IsEmpty() const
      {
        return m_Assignments.IsEmpty();
      }

      //! @brief reset the assignment storage container to the assignments in the words
      void ResetAssignments()
      {
        m_ItrAlignmentABegin = m_AlignmentA->GetBeginItr();
        m_ItrAlignmentAEnd = m_AlignmentA->GetEndItr();
        m_ItrAlignmentBBegin = m_AlignmentB->GetBeginItr();
        m_ItrAlignmentBEnd = m_AlignmentB->GetEndItr();
        UpdateAssignments();
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
        // read members
        io::Serialize::Read( m_AlignmentA, ISTREAM);
        //io::Serialize::Read( m_ItrAlignmentABegin, ISTREAM); // cannot read itrs
        //io::Serialize::Read( m_ItrAlignmentAEnd, ISTREAM); // cannot read itrs
        io::Serialize::Read( m_AlignmentB, ISTREAM);
        //io::Serialize::Read( m_ItrAlignmentBBegin, ISTREAM); // cannot read itrs
        //io::Serialize::Read( m_ItrAlignmentBEnd, ISTREAM); // cannot read itrs

        io::Serialize::Read( m_ChildAlignments, ISTREAM);
        io::Serialize::Read( m_Assignments, ISTREAM);

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_AlignmentA, OSTREAM, INDENT);
        //io::Serialize::Write( m_ItrAlignmentABegin, OSTREAM, INDENT); // cannot write itrs
        //io::Serialize::Write( m_ItrAlignmentAEnd, OSTREAM, INDENT); // cannot write itrs
        io::Serialize::Write( m_AlignmentB, OSTREAM, INDENT);
        //io::Serialize::Write( m_ItrAlignmentBBegin, OSTREAM, INDENT); // cannot write itrs
        //io::Serialize::Write( m_ItrAlignmentBEnd, OSTREAM, INDENT); // cannot write itrs

        io::Serialize::Write( m_ChildAlignments, OSTREAM, INDENT);
        io::Serialize::Write( m_Assignments, OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief returns the assignment build from the two assignments pointed at by the itrs
      //! @param first itr
      //! @param second itr
      //! @return ShPtr to an assignment
      util::ShPtr< Assignment< t_Member> > GetAssignment
      (
        typename AlignmentInterface< t_Member>::const_iterator &ITR_A,
        typename AlignmentInterface< t_Member>::const_iterator &ITR_B
      ) const
      {
        util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>()); // create new assignment
        assignment->Append( ( **ITR_A).GetMembers()); // add members for each child alignment
        assignment->Append( ( **ITR_B).GetMembers());
        return assignment;
      }

    public:

      //! @brief get the assignment for the next prepend operation
      //! @return the assignment
      util::ShPtr< Assignment< t_Member> > GetPrependAssignment() const
      {
        // check if prepend is possible
        if
        (
          m_ItrAlignmentABegin == m_AlignmentA->GetCompleteAlignment()->GetAssignments().Begin()
            || m_ItrAlignmentBBegin == m_AlignmentB->GetCompleteAlignment()->GetAssignments().Begin()
        )
        {
          return util::ShPtr< Assignment< t_Member> >();
        }

        typename AlignmentInterface< t_Member>::const_iterator
          prepend_itr_a( m_ItrAlignmentABegin),
          prepend_itr_b( m_ItrAlignmentBBegin);
        --prepend_itr_a, --prepend_itr_b; // move itr to next assignment before the hit
        return GetAssignment( prepend_itr_a, prepend_itr_b);
      }

      //! @brief get the assignment for the next append operation
      //! @return the assignment
      util::ShPtr< Assignment< t_Member> > GetAppendAssignment() const
      {
        // check if append is possible
        if
        (
          m_ItrAlignmentAEnd == m_AlignmentA->GetCompleteAlignment()->GetAssignments().End()
            || m_ItrAlignmentBEnd == m_AlignmentB->GetCompleteAlignment()->GetAssignments().End()
        )
        {
          return util::ShPtr< Assignment< t_Member> >();
        }

        typename AlignmentInterface< t_Member>::const_iterator itr_a( m_ItrAlignmentAEnd), itr_b( m_ItrAlignmentBEnd);
        return GetAssignment( itr_a, itr_b);
      }

      //! @brief prepends the next assignment before the current hit if possible
      void PrependNextAssignment()
      {
        util::ShPtr< Assignment< t_Member> > assignment( GetPrependAssignment());
        if( assignment.IsDefined())
        {
          PrependNextAssignmentImplementation( assignment);
        }
      }

      //! @brief appends the next assignment after the current hit if possible
      void AppendNextAssignment()
      {
        util::ShPtr< Assignment< t_Member> > assignment( GetAppendAssignment());
        if( assignment.IsDefined())
        {
          AppendNextAssignmentImplementation( assignment);
        }
      }

    private:

      //! @brief prepends the given assignment to the hit and updates the itrs; no validity checks
      //! @param SP_ASSIGNMENT the assignment to add
      void PrependNextAssignmentImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        --m_ItrAlignmentABegin; // update itrs
        --m_ItrAlignmentBBegin;
        m_Assignments.PushFront( SP_ASSIGNMENT); // update assignment list
      }

      //! @brief appends the given assignment to the hit and updates the itrs; no validity checks
      //! @param SP_ASSIGNMENT the assignment to add
      void AppendNextAssignmentImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        ++m_ItrAlignmentAEnd; // update itrs
        ++m_ItrAlignmentBEnd;
        m_Assignments.PushBack( SP_ASSIGNMENT); // update assignment list
      }

      //! @brief private abstract method to insert SP_ASSIGNMENT at the front
      //! @param SP_ASSIGNMENT assignment to be inserted
      void PrependImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        // check if prepend is possible
        if
        (
          m_ItrAlignmentABegin == m_AlignmentA->GetCompleteAlignment()->GetAssignments().Begin()
            || m_ItrAlignmentBBegin == m_AlignmentB->GetCompleteAlignment()->GetAssignments().Begin()
        )
        {
          return;
        }

        if( GetPrependAssignment()->IsIdentical( *SP_ASSIGNMENT)) // compare assignments
        {
          PrependNextAssignmentImplementation( SP_ASSIGNMENT);
        }
      }

      //! @brief private abstract method to insert SP_LIST_ASSIGNMENT at the front
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void PrependImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_reverse_iterator // use reverse iterators
            itr( SP_LIST_ASSIGNMENT.ReverseBegin()),
            itr_end( SP_LIST_ASSIGNMENT.ReverseEnd());
          itr != itr_end;
          ++itr
        )
        {
          PrependImplementation( *itr); // append single assignment
        }
      }

      //! @brief private abstract method to insert SP_ASSIGNMENT at the end
      //! @param SP_ASSIGNMENT assignment to be inserted
      void AppendImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        // check if append is possible
        if
        (
          m_ItrAlignmentAEnd == m_AlignmentA->GetCompleteAlignment()->GetAssignments().End()
            || m_ItrAlignmentBEnd == m_AlignmentB->GetCompleteAlignment()->GetAssignments().End()
        )
        {
          return;
        }

        if( GetAppendAssignment()->IsIdentical( *SP_ASSIGNMENT)) // compare assignments
        {
          AppendNextAssignmentImplementation( SP_ASSIGNMENT);
        }
      }

      //! @brief private abstract method to insert SP_LIST_ASSIGNMENT at the end
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void AppendImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator
            itr( SP_LIST_ASSIGNMENT.Begin()),
            itr_end( SP_LIST_ASSIGNMENT.End());
          itr != itr_end;
          ++itr
        )
        {
          AppendImplementation( *itr);
        }
      }

      //! @brief update the child alignment list with the alignments m_AlignmentA and m_AlignmentB
      void UpdateChildAlignments()
      {
        m_ChildAlignments.Reset();
        m_ChildAlignments.PushBack( m_AlignmentA);
        m_ChildAlignments.PushBack( m_AlignmentB);
      }

      //! @brief update the assignments based on the itrs
      void UpdateAssignments()
      {
        m_Assignments.Reset();
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator
            itr_a( m_ItrAlignmentABegin),
            itr_b( m_ItrAlignmentBBegin);
          itr_a != m_ItrAlignmentAEnd && itr_b != m_ItrAlignmentBEnd;
          ++itr_a, ++itr_b
        )
        {
          util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>());
          assignment->Append( ( **itr_a).GetMembers());
          assignment->Append( ( **itr_b).GetMembers());
          m_Assignments.PushBack( assignment);
        }
      }

    }; // template class AlignmentHit

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignmentHit< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignmentHit< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNMENT_HIT_H_ 
