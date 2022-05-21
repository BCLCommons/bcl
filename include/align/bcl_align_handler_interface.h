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

#ifndef BCL_ALIGN_HANDLER_INTERFACE_H_
#define BCL_ALIGN_HANDLER_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically
#include "function/bcl_function.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_align_aligner_merge.h"
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_alignment_node.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerInterface
    //! @brief interface class from which all methods to write and read an Alignment are derived
    //!
    //! @tparam t_Member type of elements of a sequence and an assignment
    //!
    //! @remarks example unnecessary
    //! @author alexanns, heinzes1
    //! @date 03/09/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerInterface :
      public util::ObjectInterface
    {

    private:

    /////////////////
    // data access //
    /////////////////

      //! @brief gap char for a specific alignment output format
      //! @return the gap char
      virtual char GetGapChar() const = 0;

    public:

      //! @brief set the scoring function for single assignment, if relevant for that handler
      //! @param SP_ASSIGNMENT_SCORE the assignment score
      virtual void SetAssignmentScore
      (
        const util::ShPtr< function::UnaryInterface< const Assignment< t_Member>, double> > &SP_ASSIGNMENT_SCORE
      )
      {
      }

      //! @brief get the default file extension for alignments of that handler
      //! @return string of default file extension
      virtual const std::string &GetFileExtension() const = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief ReadAlignment reads an alignment from a given stream
      //! @param ISTREAM the stream from which the alignment is read
      //! @param SEQUENCE_TYPE the implementation type of SequenceInterface
      //! @param CHAIN_ID chain id used to build the sequences
      //! @return the read alignment
      util::ShPtr< AlignmentNode< t_Member> > ReadAlignment
      (
        std::istream &ISTREAM,
        const SequenceInterface< t_Member> &SEQUENCE_TYPE,
        const char CHAIN_ID = 'A'
      ) const
      {
        // read istream into id-and-gapped-sequence format
        storage::Vector< storage::VectorND< 2, std::string> > gapped_strings( ReadAlignmentToStrings( ISTREAM));

        // create alignment list with one sequence each
        util::ShPtrList< AlignmentInterface< t_Member> > alignment_list;
        storage::Vector< storage::VectorND< 2, std::string> >::const_iterator
          itr( gapped_strings.Begin()),
          itr_end( gapped_strings.End());
        for( ; itr != itr_end; ++itr) // for all sequences read
        {
          // create an ShPtr< SequenceInterface> from an empty sequence string and create alignment from this sequence
          util::ShPtr< SequenceInterface< t_Member> > seq( SEQUENCE_TYPE.Construct( itr->First(), std::string(), CHAIN_ID));
          util::ShPtr< AlignmentInterface< t_Member> > alignment( new AlignmentLeaf< t_Member>( seq));

          // loop over all letters in gapped sequence string and add t_Members
          std::string current_gapped_string( itr->Second());
          for
          (
            std::string::const_iterator itr( current_gapped_string.begin()), itr_end( current_gapped_string.end());
            itr != itr_end;
            ++itr
          )
          {
            char letter( *itr);
            util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>());

            // if the letter encodes a gap
            if( letter == GetGapChar())
            {
              // dont add a gap to the sequence, just create an assignment with an empty SiPtr
              assignment->Append( util::SiPtr< const t_Member>());
            }
            else
            {
              // otherwise add member to sequence from character
              seq->AddMember( letter);
              // from the sequence get a SiPtr to the last member, and add it to the assignment
              util::SiPtr< const t_Member> member( seq->GetLastMember());
              assignment->Append( member);
            }

            alignment->Append( assignment);
          }

          // add alignment to list of alignments
          alignment_list.Append( alignment);
        }

        // merge
        AlignerMerge< t_Member> aligner;
        util::ShPtr< AlignmentNode< t_Member> > result_alignment( aligner.AlignMultiple( alignment_list).First().Clone());
        return result_alignment;
      }

      //! @brief ReadAlignment reads an alignment from a given stream using the given sequences
      //! @param ISTREAM the stream from which the alignment is read
      //! @param SEQUENCES
      //! @return the read alignment
      util::ShPtr< AlignmentNode< t_Member> > ReadAlignment
      (
        std::istream &ISTREAM,
        const util::ShPtrVector< SequenceInterface< t_Member> > &SEQUENCES
      ) const
      {
        // read istream into id-and-gapped-sequence format
        storage::Vector< storage::VectorND< 2, std::string> > gapped_strings( ReadAlignmentToStrings( ISTREAM));

        // create alignment list with one sequence each
        util::ShPtrList< AlignmentInterface< t_Member> > alignment_list;
        storage::Vector< storage::VectorND< 2, std::string> >::const_iterator
          itr( gapped_strings.Begin()),
          itr_end( gapped_strings.End());
        typename util::ShPtrVector< SequenceInterface< t_Member> >::const_iterator
          itr_sequence( SEQUENCES.Begin()),
          itr_sequence_end( SEQUENCES.End());
        for( ; itr != itr_end && itr_sequence != itr_sequence_end; ++itr, ++itr_sequence) // for all sequences given
        {
          // create alignment from sequence list
          // (do not use the constructor taking a sequence, it automatically creates assignments)
          const util::ShPtrList< SequenceInterface< t_Member> > sequence_list( 1, *itr_sequence);
          util::ShPtr< AlignmentInterface< t_Member> > alignment( new AlignmentLeaf< t_Member>( sequence_list));

          // loop over all letters in gapped sequence string and given sequence simultaneously and add t_Members
          std::string current_gapped_string( itr->Second());
          typename util::ShPtrVector< t_Member>::const_iterator itr_member( ( **itr_sequence).Begin());
          for
          (
            std::string::const_iterator itr( current_gapped_string.begin()), itr_end( current_gapped_string.end());
            itr != itr_end; // do not test for itr_member_end here, as itr might have gaps at the end that would be missed
            ++itr // do not increase itr_member here, it does not have gaps as itr does
          )
          {
            char letter( *itr);
            util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>());

            if( letter == GetGapChar()) // if the letter encodes a gap, create an assignment with an empty SiPtr
            {
              assignment->Append( util::SiPtr< const t_Member>());
            }
            else // otherwise add SiPtr to current itr_sequence to the assignment
            {
              util::SiPtr< const t_Member> member( *itr_member);
              assignment->Append( member);
              ++itr_member;
            }

            alignment->Append( assignment);
          }

          // add alignment to list of alignments
          alignment_list.Append( alignment);
        }

        // merge
        AlignerMerge< t_Member> aligner;
        util::ShPtr< AlignmentNode< t_Member> > result_alignment( aligner.AlignMultiple( alignment_list).First().Clone());
        return result_alignment;
      }

      //! @brief WriteAlignment prints an alignment to a given stream
      //! @param OSTREAM the stream which the alignment is written to
      //! @param ALIGNMENT the alignment which is output
      //! @return the passed std::ostream
      virtual std::ostream &WriteAlignment
      (
        std::ostream &OSTREAM,
        const AlignmentInterface< t_Member> &ALIGNMENT
      ) const = 0;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads the alignment from the istream into a vector of pairs of identifier, gapped sequence
      //! @param ISTREAM the stream to read from
      //! @return a vector of pairs of identifier, gapped sequence
      virtual storage::Vector< storage::VectorND< 2, std::string> > ReadAlignmentToStrings( std::istream &ISTREAM) const = 0;

    protected:

      //! @brief get all the aligned sequence from an alignment
      //! @param ALIGNMENT the alignment to extract the sequences from
      //! @param GAP_CHAR the gaps char to use when creating gapped strings
      //! @return a list of aligned sequence strings
      static storage::List< std::string> GetAlignedSequenceStrings
      (
        const AlignmentInterface< t_Member> &ALIGNMENT,
        const char GAP_CHAR
      )
      {
        // create as many empty strings as we have sequences in the alignment
        storage::List< std::string> aligned_sequences( ALIGNMENT.GetDepth());

        // get assignments and create iterators to loop over all assignments
        util::ShPtrList< Assignment< t_Member> > assignments( ALIGNMENT.GetAssignments());
        typename util::ShPtrList< Assignment< t_Member> >::const_iterator itr( assignments.Begin()), itr_end( assignments.End());
        for( ; itr != itr_end; ++itr)
        {
          // get members and create iterators to loop overall members and add them to their respective strings
          util::SiPtrList< const t_Member> members( ( *itr)->GetMembers());
          typename util::SiPtrList< const t_Member>::const_iterator itr( members.Begin()), itr_end( members.End());
          typename storage::List< std::string>::iterator string_itr( aligned_sequences.Begin());
          for( ; itr != itr_end; ++itr, ++string_itr)
          {
            const char char_id( itr->IsDefined() ? GetCharId< t_Member>( **itr) : GAP_CHAR);
            string_itr->append( 1, char_id);
          }
        }

        return aligned_sequences;
      }

    }; // template class HandlerInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_INTERFACE_H_
