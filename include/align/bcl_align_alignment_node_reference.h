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

#ifndef BCL_ALIGN_ALIGNMENT_NODE_REFERENCE_H_
#define BCL_ALIGN_ALIGNMENT_NODE_REFERENCE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentNodeReference
    //! @brief a node alignment in the alignment tree; has no sequences and 2..n parent alignments, whose ownership
    //!        lies with any calling code
    //!
    //! @tparam t_Member type of elements of a sequence and an assignment
    //!
    //! @see @link example_align_alignment_node_reference.cpp @endlink
    //! @author mendenjl
    //! @date May 22, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentNodeReference :
      public AlignmentInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      util::SiPtrList< const AlignmentInterface< t_Member> > m_ChildAlignments; //!< list of pointers to child alignments
      util::ShPtrList< Assignment< t_Member> >         m_Assignments; //!< list of pointers to assignments

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, alignment node with an empty list of child alignments and assignments
      AlignmentNodeReference() :
        m_ChildAlignments(),
        m_Assignments()
      {
      }

      //! @brief constructor from a list of parent alignments
      //! @param SP_CHILD_ALIGNMENT_LIST list of ShPtr to parent alignments
      AlignmentNodeReference( const util::SiPtrList< const AlignmentInterface< t_Member> > &SP_CHILD_ALIGNMENT_LIST) :
        m_ChildAlignments( SP_CHILD_ALIGNMENT_LIST),
        m_Assignments()
      {
        // if less than two parent alignments, remove all; empty m_ParentAlignments -> no assignments can be inserted
        if( m_ChildAlignments.GetSize() == 1)
        {
          m_ChildAlignments.Reset();
        }
      }

      //! @brief constructor from two parent alignments
      //! @param SP_CHILD_ALIGNMENT_A ShPtr to first parent alignment
      //! @param SP_CHILD_ALIGNMENT_B ShPtr to second parent alignment
      AlignmentNodeReference
      (
        const util::SiPtr< const AlignmentInterface< t_Member> > &SP_CHILD_ALIGNMENT_A,
        const util::SiPtr< const AlignmentInterface< t_Member> > &SP_CHILD_ALIGNMENT_B
      ) :
        m_ChildAlignments(),
        m_Assignments()
      {
        m_ChildAlignments.PushBack( SP_CHILD_ALIGNMENT_A);
        m_ChildAlignments.PushBack( SP_CHILD_ALIGNMENT_B);
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentNodeReference< t_Member>
      AlignmentNodeReference< t_Member> *Clone() const
      {
        return new AlignmentNodeReference< t_Member>( *this);
      }

      //! @brief virtual empty constructor maintaining the child alignment list while resetting the assignment list
      //! @return pointer to new AlignmentNodeReference< t_Member>
      AlignmentNodeReference< t_Member> *Empty() const
      {
        return new AlignmentNodeReference< t_Member>( m_ChildAlignments);
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

      //! @brief returns the sequences the alignment is composed of, collected from all the parent alignments
      //! @return a list of ShPtr to SequenceInterface
      util::ShPtrList< SequenceInterface< t_Member> > GetSequences() const
      {
        util::ShPtrList< SequenceInterface< t_Member> > sequence_list;
        for
        (
          typename util::SiPtrList< const AlignmentInterface< t_Member> >::const_iterator
            itr( m_ChildAlignments.Begin()),
            itr_end( m_ChildAlignments.End());
          itr != itr_end;
          ++itr
        )
        {
          sequence_list.Append( ( *itr)->GetSequences()); // append sequences from each child alignment
        }
        return sequence_list;
      }

      //! @brief set the list of pointers to sequences
      //! @return list of ShPtr to AlignmentInterface
      const util::ShPtrList< AlignmentInterface< t_Member> > &GetChildAlignments() const
      {
        static const util::ShPtrList< AlignmentInterface< t_Member> > s_empty;
        BCL_Exit
        (
          this->GetClassIdentifier() + "::GetChildAlignments() cannot return a shared pointer list because it "
          "does not share ownership of child alignments",
          -1
        );
        return s_empty;
      }

      //! @brief get an iterator to the list of child alignments; can be implemented even without shared ownership
      //! @return iterator to the child alignments
      iterate::Generic< const AlignmentInterface< t_Member> > GetChildAlignmentsIterator() const
      {
        return
          iterate::Generic< const AlignmentInterface< t_Member> >
          (
            m_ChildAlignments.Begin(),
            m_ChildAlignments.End()
          );
      }

      //! @brief GetDepth gives the number of t_Members are in each Assignment i.e. the number of sequences
      //! @return size_t which is the number of t_Members that are in each of the Assignments
      size_t GetDepth() const
      {
        size_t depth( 0);
        for
        (
          typename util::SiPtrList< const AlignmentInterface< t_Member> >::const_iterator
            itr( m_ChildAlignments.Begin()),
            itr_end( m_ChildAlignments.End());
          itr != itr_end;
          ++itr
        )
        {
          depth += ( *itr)->GetDepth();
        }
        return depth;
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

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether any assignment are stored
      //! @return if the assignment container is empty or not
      bool IsEmpty() const
      {
        return m_Assignments.IsEmpty();
      }

      //! @brief reset the assignment storage container, removing all assignments
      void ResetAssignments()
      {
        m_Assignments.Reset();
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
        // Read member
        io::Serialize::Read( m_ChildAlignments, ISTREAM);
        io::Serialize::Read( m_Assignments    , ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_ChildAlignments, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Assignments    , OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief  method to insert SP_ASSIGNMENT at the front
      //! @param SP_ASSIGNMENT assignment to be inserted
      void PrependImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        m_Assignments.PushFront( SP_ASSIGNMENT);
      }

      //! @brief  method to insert SP_LIST_ASSIGNMENT at the front
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void PrependImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        m_Assignments.Prepend( SP_LIST_ASSIGNMENT);
      }

      //! @brief  method to insert SP_ASSIGNMENT at the end
      //! @param SP_ASSIGNMENT assignment to be inserted
      void AppendImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        m_Assignments.PushBack( SP_ASSIGNMENT);
      }

      //! @brief  method to insert SP_LIST_ASSIGNMENT at the end
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void AppendImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        m_Assignments.Append( SP_LIST_ASSIGNMENT);
      }

    }; // template class AlignmentNodeReference

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignmentNodeReference< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignmentNodeReference< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNMENT_NODE_REFERENCE_H_
