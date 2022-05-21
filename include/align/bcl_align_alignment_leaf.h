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

#ifndef BCL_ALIGN_ALIGNMENT_LEAF_H_
#define BCL_ALIGN_ALIGNMENT_LEAF_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_interface.h"
#include "bcl_align_sequence_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentLeaf
    //! @brief a leaf alignment in an alignment tree; has only a single sequence and no parent alignments
    //!
    //! @tparam t_Member type of elements of a sequence and an assignment
    //!
    //! @see @link example_align_alignment_leaf.cpp @endlink
    //! @author heinzes1
    //! @date Sep 28, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentLeaf :
      public AlignmentInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      //! list of pointers to sequences; multiple sequences possible e.g. for collapsed leafs or read in alignments
      //! m_Sequences has shared ownership of the t_Member, while assignments have only SiPtr to t_Member
      util::ShPtrList< SequenceInterface< t_Member> > m_Sequences;
      util::ShPtrList< Assignment< t_Member> > m_Assignments; //!< list of pointers to assignments

      static const util::ShPtrList< AlignmentInterface< t_Member> > s_EmptyChildAlignmentList;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a sequence; constructing without sequence not useful, so no default constructor
      //! @param SP_SEQUENCE ShPtr to sequence
      AlignmentLeaf( const util::ShPtr< SequenceInterface< t_Member> > &SP_SEQUENCE) :
        m_Sequences( 1, SP_SEQUENCE),
        m_Assignments()
      {
        // insert sequence members into alignment as assignments with a single element
        const util::SiPtrVector< const t_Member> members( SP_SEQUENCE->GetMembers());
        for
        (
          typename util::SiPtrVector< const t_Member>::const_iterator itr( members.Begin()), itr_end( members.End());
          itr != itr_end;
          ++itr
        )
        {
          // create ShPtr to new Assignment for adding to the end of alignment
          AlignmentInterface< t_Member>::Append( util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( *itr)));
        }
      }

      //! @brief constructor from a list of sequences
      //! @param SP_SEQUENCE ShPtr to sequence
      AlignmentLeaf( const util::ShPtrList< SequenceInterface< t_Member> > &SP_SEQUENCE_LIST) :
        m_Sequences( SP_SEQUENCE_LIST),
        m_Assignments()
      {
      }

      //! @brief constructor from AlignmentInterface; this collapses a AlignmentNode into an AlignmentLeaf
      //! @param ALIGNMENT alignment of type AlignmentInterface to be converted into an AlignmentLeaf type
      AlignmentLeaf( const AlignmentInterface< t_Member> &ALIGNMENT) :
        m_Sequences( ALIGNMENT.GetSequences()),
        m_Assignments( ALIGNMENT.GetAssignments())
      {
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentLeaf< t_Member>
      AlignmentLeaf< t_Member> *Clone() const
      {
        return new AlignmentLeaf< t_Member>( *this);
      }

      //! @brief virtual empty constructor maintaining the sequence list while resetting the assignment list
      //! @return pointer to new AlignmentLeaf< t_Member>
      AlignmentLeaf< t_Member> *Empty() const
      {
        return new AlignmentLeaf< t_Member>( this->m_Sequences);
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

      //! @brief returns the single sequence the leaf alignment has
      //! @return a list of ShPtr to SequenceInterface
      util::ShPtrList< SequenceInterface< t_Member> > GetSequences() const
      {
        return m_Sequences;
      }

      //! @brief set the list of pointers to sequences
      //! @return list of ShPtr to AlignmentInterface
      const util::ShPtrList< AlignmentInterface< t_Member> > &GetChildAlignments() const
      {
        return s_EmptyChildAlignmentList; // always returns an empty list, no children in a leaf
      }

      //! @brief GetDepth gives the number of t_Members are in each Assignment i.e. the number of sequences
      //! @return size_t which is the number of t_Members that are in each of the Assignments
      size_t GetDepth() const
      {
        return m_Sequences.GetSize();
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
        // read member
        io::Serialize::Read( m_Sequences  , ISTREAM);
        io::Serialize::Read( m_Assignments, ISTREAM);

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
        io::Serialize::Write( m_Sequences  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Assignments, OSTREAM, INDENT);

        // return
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

    }; // template class AlignmentLeaf

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignmentLeaf< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignmentLeaf< t_Member>())
    );

    // instantiate s_EmptyChildAlignmentList
    template< typename t_Member>
    const util::ShPtrList< AlignmentInterface< t_Member> > AlignmentLeaf< t_Member>::s_EmptyChildAlignmentList;

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNMENT_LEAF_H_
