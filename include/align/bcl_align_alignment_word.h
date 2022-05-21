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

#ifndef BCL_ALIGN_ALIGNMENT_WORD_H_
#define BCL_ALIGN_ALIGNMENT_WORD_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_alignment_node.h"
#include "biol/bcl_biol_aa_sequence.h"

// external includes - sorted alphabetically
#include <algorithm>

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentWord
    //! @brief class storing an word i.e. a piece of sequence of an alignment (a sub-alignment)
    //!
    //! @tparam t_Member elements of the sequence used in the alignment
    //!
    //! @see @link example_align_alignment_word.cpp @endlink
    //! @author heinzes1
    //! @date Mar 26, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentWord :
      public AlignmentInterface< t_Member>
    {
    private:

    //////////
    // data //
    //////////

      util::ShPtr< AlignmentInterface< t_Member> > m_Alignment; //!< alignment, on which the itrs are based upon
      typename AlignmentInterface< t_Member>::const_iterator m_ItrWordBegin; //!< itr to the begin of the word
      typename AlignmentInterface< t_Member>::const_iterator m_ItrWordEnd; //!< itr to the end of the word

      util::ShPtrList< AlignmentInterface< t_Member> > m_ChildAlignments; //!< list of pointers to child alignments
      util::ShPtrList< Assignment< t_Member> > m_Assignments; //!< list of pointers to assignments denoted by itrs

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor; creates an empty leaf and point itrs to assignments.end of that alignment
      AlignmentWord() :
        m_Alignment( util::ShPtr< AlignmentInterface< t_Member> >( new AlignmentNode< t_Member>())),
        m_ItrWordBegin( m_Alignment->GetAssignments().End()),
        m_ItrWordEnd( m_Alignment->GetAssignments().End()),
        m_Assignments()
      {
        UpdateChildAlignments();
        UpdateAssignments();
      }

      //! @brief constructor from ShPtr to AlignmentInterface; itrs will point to begin and end of that alignment
      //! @param ALIGNMENT alignment of type AlignmentInterface
      AlignmentWord( const util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT) :
        m_Alignment( ALIGNMENT),
        m_ItrWordBegin( m_Alignment->GetAssignments().End()),
        m_ItrWordEnd( m_Alignment->GetAssignments().End()),
        m_Assignments()
      {
        UpdateChildAlignments();
        UpdateAssignments();
      }

      //! @brief constructor from AlignmentInterface
      //! @param ALIGNMENT alignment of type ShPtr<AlignmentInterface>
      //! @param ITR_WORD_BEGIN itr to the begin of the word within the given alignment
      //! @param ITR_WORD_END itr to the end of the word within the given alignment
      AlignmentWord
      (
        const util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT,
        typename AlignmentInterface< t_Member>::const_iterator &ITR_WORD_BEGIN,
        typename AlignmentInterface< t_Member>::const_iterator &ITR_WORD_END
      ) :
        m_Alignment( ALIGNMENT),
        m_ItrWordBegin( ITR_WORD_BEGIN),
        m_ItrWordEnd( ITR_WORD_END),
        m_Assignments()
      {
        UpdateChildAlignments();
        UpdateAssignments();
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentWord< t_Member>
      AlignmentWord< t_Member> *Clone() const
      {
        return new AlignmentWord< t_Member>( *this);
      }

      //! @brief virtual empty constructor maintaining the sequence list while resetting the assignment itrs
      //! @return pointer to new AlignmentWord< t_Member>
      AlignmentWord< t_Member> *Empty() const
      {
        return new AlignmentWord< t_Member>( m_Alignment);
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
        return m_Alignment->GetSequences();
      }

      //! @brief set the list of pointers to sequences
      //! @return list of ShPtr to AlignmentInterface
      const util::ShPtrList< AlignmentInterface< t_Member> > &GetChildAlignments() const
      {
        return m_ChildAlignments;
      }

      //! @brief returns the complete alignment of which the word is the subalignment of
      //! @return the alignment
      const util::ShPtr< AlignmentInterface< t_Member> > &GetCompleteAlignment() const
      {
        return m_Alignment;
      }

      //! @brief GetDepth gives the number of t_Members are in each Assignment i.e. the number of sequences
      //! @return size_t which is the number of t_Members that are in each of the Assignments
      size_t GetDepth() const
      {
        return m_Alignment->GetDepth();
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

      //! @brief returns the itr of the begin of the subalignment in the child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetBeginItr() const
      {
        return m_ItrWordBegin;
      }

      //! @brief returns the itr of the end of the subalignment in the child alignment
      //! @return the itr
      const typename AlignmentInterface< t_Member>::const_iterator &GetEndItr() const
      {
        return m_ItrWordEnd;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether any assignment are stored
      //! @return if the alignment has any assignments or not
      bool IsEmpty() const
      {
        return m_Assignments.IsEmpty();
      }

      //! @brief reset the assignment storage container, removing all assignments
      void ResetAssignments()
      {
        m_ItrWordBegin = m_ItrWordEnd = m_Alignment->GetAssignments().End();
        UpdateAssignments();
      }

      //! @brief scores the alignment with itself using the given assignment scoring function
      //! @param ASSIGN_SCORE the scoring function for scoring assignments
      //! @return the score
      double ScoreSelf( const score::AssignmentWithGap< t_Member> &ASSIGN_SCORE) const
      {
        // create a second list for inserting,
        // on MacOS inserting into the list leads to an endless loop (memory allocation error)
        util::ShPtrList< SequenceInterface< t_Member> > sequence_list( GetSequences());
        util::ShPtrList< SequenceInterface< t_Member> > sequence_list_copy( GetSequences());

        sequence_list.InsertElements( sequence_list.End(), sequence_list_copy.Begin(), sequence_list_copy.End());
        AlignmentLeaf< t_Member> sequence_alignment( sequence_list); // create an alignment with the correct #sequences

        typename AlignmentLeaf< t_Member>::const_iterator itr_word( m_ItrWordBegin); // itrs to loop over query word
        for( ; itr_word != m_ItrWordEnd; ++itr_word)
        {
          util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>()); // create new assignment
          assignment->Append( ( **itr_word).GetMembers()); // add members twice
          assignment->Append( ( **itr_word).GetMembers());
          sequence_alignment.Append( assignment); // append assignment to sequence
        }

        return sequence_alignment.Score( ASSIGN_SCORE); // score alignment and return score
      }

      //! @brief scores a query word from an query alignment
      //! @param SEQUENCE a word sequence against which the query word is scored
      //! @param ASSIGN_SCORE the scoring function for scoring assignments
      //! @return the score
      double ScoreWith( const biol::AASequence &SEQUENCE, const score::AssignmentWithGap< t_Member> &ASSIGN_SCORE) const
      {
        util::ShPtrList< SequenceInterface< t_Member> > sequence_list( GetSequences());
        sequence_list.Append( util::ShPtr< biol::AASequence>( SEQUENCE.Clone())); // ensure correct #sequences
        AlignmentLeaf< t_Member> sequence_alignment( sequence_list); // create an alignment

        typename AlignmentLeaf< t_Member>::const_iterator itr_word( m_ItrWordBegin); // itrs to loop over query word
        typename util::ShPtrVector< t_Member>::const_iterator // itrs to loop over sequence
          itr_sequence( SEQUENCE.Begin()),
          itr_sequence_end( SEQUENCE.End());
        for( ; itr_word != m_ItrWordEnd && itr_sequence != itr_sequence_end; ++itr_word, ++itr_sequence)
        {
          util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>()); // create new assignment
          assignment->Append( ( **itr_word).GetMembers()); // add members
          assignment->Append( *itr_sequence);
          sequence_alignment.Append( assignment); // append assignment to sequence
        }

        return sequence_alignment.Score( ASSIGN_SCORE); // score alignment and return score
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
        io::Serialize::Read( m_Alignment, ISTREAM);
        //io::Serialize::Read( m_ItrWordBegin, ISTREAM); // cannot read itrs
        //io::Serialize::Read( m_ItrWordEnd, ISTREAM); // cannot read itrs
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
        io::Serialize::Write( m_Alignment, OSTREAM, INDENT);
        //io::Serialize::Write( m_ItrWordBegin, OSTREAM, INDENT); // cannot write itrs
        //io::Serialize::Write( m_ItrWordEnd, OSTREAM, INDENT); // cannot write itrs
        io::Serialize::Write( m_Assignments, OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief returns a simple string representation of the alignment word
      //! @return a string
      std::string ToString() const
      {
        std::string alignment_string;

        // write all assignments (as sequence of one letter codes for all members) separated by semicolons
        for( typename AlignmentInterface< t_Member>::const_iterator itr( m_ItrWordBegin); itr != m_ItrWordEnd; ++itr)
        {
          Assignment< t_Member> assignment( **itr);
          alignment_string.append( assignment.ToString().append( ";"));
        }

        return alignment_string;
      }

    private:

      //! @brief  method to insert SP_ASSIGNMENT at the front
      //! @param SP_ASSIGNMENT assignment to be inserted
      void PrependImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        if( m_ItrWordBegin != m_Alignment->GetAssignments().Begin()) // if something can be prepended
        {
          typename AlignmentInterface< t_Member>::const_iterator prepend_itr( m_ItrWordBegin);
          --prepend_itr; // get itr pointing to the assignment before the first assignment of the word
          if( *prepend_itr == SP_ASSIGNMENT) // if new assignment is this assignment, then add it
          {
            m_ItrWordBegin = prepend_itr; // update itr
            m_Assignments.PushFront( SP_ASSIGNMENT); // update assignment list
          }
        }
      }

      //! @brief  method to insert SP_LIST_ASSIGNMENT at the front
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void PrependImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_reverse_iterator // use reverse iterators
            itr( SP_LIST_ASSIGNMENT.ReverseBegin()),
            itr_end( SP_LIST_ASSIGNMENT.ReverseEnd());
          itr != itr_end && itr != m_Alignment->GetAssignments().ReverseEnd();
          ++itr
        )
        {
          PrependImplementation( *itr); // append single assignment
        }
      }

      //! @brief  method to insert SP_ASSIGNMENT at the end
      //! @param SP_ASSIGNMENT assignment to be inserted
      void AppendImplementation( const util::ShPtr< Assignment< t_Member> > &SP_ASSIGNMENT)
      {
        if( m_ItrWordEnd != m_Alignment->GetAssignments().End()) // if something can be appended
        {
          if( *m_ItrWordEnd == SP_ASSIGNMENT) // test if new assignment is pointed at by old word end iterator
          {
            ++m_ItrWordEnd; // update itr
            m_Assignments.PushBack( SP_ASSIGNMENT); // update assignment list
          }
        }
      }

      //! @brief  method to insert SP_LIST_ASSIGNMENT at the end
      //! @param SP_LIST_ASSIGNMENT assignment list to be inserted
      void AppendImplementation( const util::ShPtrList< Assignment< t_Member> > &SP_LIST_ASSIGNMENT)
      {
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator
            itr( SP_LIST_ASSIGNMENT.Begin()),
            itr_end( SP_LIST_ASSIGNMENT.End());
          itr != itr_end && itr != m_Alignment->GetAssignments().End();
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
        m_ChildAlignments.PushBack( m_Alignment);
      }

      //! @brief update the assignments based on the itrs
      void UpdateAssignments()
      {
        m_Assignments.Reset();
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator itr( m_ItrWordBegin);
          itr != m_ItrWordEnd;
          ++itr
        )
        {
          m_Assignments.PushBack( *itr);
        }
      }

    }; // template class AlignmentWord

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignmentWord< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignmentWord< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNMENT_WORD_H_ 
