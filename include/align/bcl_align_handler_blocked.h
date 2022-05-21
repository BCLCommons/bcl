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

#ifndef BCL_ALIGN_HANDLER_BLOCKED_H_
#define BCL_ALIGN_HANDLER_BLOCKED_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_handler_interface.h"
#include "score/bcl_score_assignment_with_gap.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerBlocked
    //! @brief class writes Alignments in blocked format
    //!
    //! @tparam t_Member is the type of object that the Assignment stores
    //!
    //! @see @link example_align_handler_blocked.cpp @endlink
    //! @author meilerj, alexanns
    //! @date 03/09/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerBlocked :
      public HandlerInterface< t_Member>
    {
    private:

    //////////
    // data //
    //////////

      bool   m_UnGapped;  //!< indicates whether the output should include gaps or not
      size_t m_BlockSize; //!< the number of characters on a line

      //! ScoreAssignment object
      score::AssignmentWithGap< t_Member> m_ScoreAssignment;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HandlerBlocked() :
        m_UnGapped(),
        m_BlockSize(),
        m_ScoreAssignment()
      {
      }

      //! @brief constructor with parameters of whether to output gaps or not and defining the width of the block
      //! @param UNGAPPED bool which indicates if gaps will be printed or not (true = no gaps; false = gaps)
      //! @param BLOCK_SIZE size_t which indicates the width of the block in which the Alignment will be printed
      //! @param SCORE_ASSIGNMENT ???
      HandlerBlocked
      (
        const bool UNGAPPED,
        const size_t BLOCK_SIZE,
        const score::AssignmentWithGap< t_Member> &SCORE_ASSIGNMENT
      ) :
        m_UnGapped( UNGAPPED),
        m_BlockSize( BLOCK_SIZE),
        m_ScoreAssignment( SCORE_ASSIGNMENT)
      {
      }

      //! @brief Clone is the virtual copy constructor
      HandlerBlocked< t_Member> *Clone() const
      {
        return new HandlerBlocked< t_Member>( *this);
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

      //! @brief get the default file extension for alignments of that handler
      //! @return string of default file extension
      const std::string &GetFileExtension() const
      {
        static const std::string s_file_extension( ".blocked");
        return s_file_extension;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief BlockSize gives copy of "m_BlockSize"
      //! @return returns copy of size_t "m_BlockSize"
      size_t BlockSize() const
      {
        return m_BlockSize;
      }

      //! @brief SetBlockSize "m_BlockSize" to a new size_t
      //! @@param BLOCK_SIZE size_t which "m_BlockSize" will be set to
      void SetBlockSize( const size_t BLOCK_SIZE)
      {
        m_BlockSize = BLOCK_SIZE;
      }

      //! @brief UnGapped gives copy of "m_UnGapped"
      //! @return returns copy of bool "m_UnGapped"
      bool UnGapped() const
      {
        return m_UnGapped;
      }

      //! @brief SetUnGapped changes "m_UnGapped" to new bool
      //! @param BOOL bool which "m_UnGapped" will be changed to
      void SetUnGapped( const bool BOOL)
      {
        m_UnGapped = BOOL;
      }

      //! @brief WriteAlignment prints an Alignment in the blocked format to a given stream
      //! @param OSTREAM is the stream which the Alignment is written to
      //! @param ALIGNMENT is the Alignment which is outputted in the blocked format
      //! @return returns a std::ostream
      //! TODO: make it output the correct score
      std::ostream &WriteAlignment( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
//        // write header
//        OSTREAM << ALIGNMENT.GetSize() << '\n' << '\n';
//
//        // create size_t "start_id", "cycle", "count", and "j"
//        size_t start_id( 0), //< initialize to 0
//               cycle( 0),    //< initialize to 0
//               count( 0),    //< initialize to 0, keeps track of how many characters have been printed on current line
//               j( start_id); //< initialize to start_id, keeps track of how many Assignments have been printed so far
//
//        // create Alignment const_iterator "itr_j" initialize to the beginning of the Alignment
//        typename AlignmentSimple< t_Member>::const_iterator itr_j( ALIGNMENT.GetAssignments().Begin());
//
//        // create Alignment const_iterator "itr_start_id"
//        typename AlignmentSimple< t_Member>::const_iterator itr_start_id( itr_j); //< initialize with "itr_j"
//
//        // write alignment
//        do
//        {
//          // output ruler
//          count = 0;
//          j = start_id;
//          itr_j = itr_start_id;
//
//          // start line at correct horizontal position
//          OSTREAM << "                        ";
//
//          // print 10 character segment of "." which ends in the number of characters printed so far
//          while( count < m_BlockSize && j < ALIGNMENT.GetSize())
//          {
//            // check if gaps can be printed or first SiPtr of the Assignment denoted by "itr_j" points to an amino acid
//            if( !m_UnGapped || ( *itr_j)->GetMembers().Begin()->IsDefined())
//            {
//              // check if "count" is evenly divisible by 10
//              if( !( ++count % 10))
//              {
//                // print 10 character segment of "." and the number of total spaces printed so far plus a blank space
//                OSTREAM << util::Format().W( 10).Fill( '.')( cycle * m_BlockSize + count) << ' ';
//              }
//            }
//
//            // increase "j" and "itr_j"
//            j++;
//            itr_j++;
//          }
//
//          // if "count" is not evenly divisible by 10 then print "." for the left out characters
//          if( count % 10)
//          {
//            OSTREAM << util::Format().W( count % 10).Fill( '.')( " ");
//
//          }
//
//          // print a new line
//          OSTREAM << '\n';
//
//          // reset "count", "j", and "itr_j" for outputing scores of the Assignments
//          count = 0;
//          j = start_id;
//          itr_j = itr_start_id;
//
//          // start line at correct horizontal position
//          OSTREAM << "                        ";
//
//          // output the scores of the Assignments
//          while( count < m_BlockSize && j < ALIGNMENT.GetSize())
//          {
//            // check if gaps can be printed or first SiPtr of the Assignment denoted by "itr_j" points to an amino acid
//            if( !m_UnGapped || ( *itr_j)->GetMembers().Begin()->IsDefined())
//            {
//              // create int "score" and initialize to the integer of the score of the Assignment denoted by "itr_j"
//              // TODO: make it output the correct score
//              int score( math::range( int( 2 * ( 1 + m_ScoreAssignment( **itr_j))), 0, 9)); //< ( 0 <= score <= 9)
//
//              // write score
//              OSTREAM << score;
//
//              // print space every 10 AAs
//              if( !( ++count % 10))
//              {
//                OSTREAM << ' ';
//              }
//            }
//
//            // increase "j" and "itr_j"
//            j++;
//            itr_j++;
//          }
//
//          // print end of line
//          OSTREAM << '\n';
//
//          // reset "itr_j" to "itr_start_id"
//          itr_j = itr_start_id;
//
//          // create Group const_iterator "itr_i" and initialize to the first t_GroupMember in the restraint::GroupCollection
//          // of the Assignment denoted by "itr_j"
//          typename util::SiPtrList< const t_Member>::const_iterator itr_i( ( *itr_j)->GetMembers().Begin());
//
//          // iterate over the depth of the Alignment to output the t_GroupMembers of the Assignment denoted by "itr_i"
//          for( size_t i( 0); i < ALIGNMENT.GetDepth(); ++i)
//          {
//            // reset "count" and "j" to zero and "start_id", respectively
//            count = 0;
//            j = start_id;
//
//            // reset "itr_j" to "itr_start_id" which is the beginning of the Alignment
//            itr_j = itr_start_id;
//
//            // set "itr_i" to the "i"th t_GroupMember of the restraint::GroupCollection
//            itr_i = (*itr_j)->GetNthMember( i);
//
//            // create size_ts "last" and "first" for holding the seq ids of the first and last aas in the Alignment
//            size_t last( util::GetUndefined< size_t>()),  //< initialize to undefined
//                   first( util::GetUndefined< size_t>()); //< initialize to undefined
//
//            // find "first" and "last" seq ids of the first and last amino acids in the Alignment
//            while( count < m_BlockSize && j < ALIGNMENT.GetSize())
//            {
//              // if the SiPtr denoted by "itr_i" points to an aa
//              if( itr_i->IsDefined())
//              {
//                // if "first" is not defined
//                if( !util::IsDefined( first))
//                {
//                  // "first" is the seq id of the aa denoted by "itr_i"
//                  first = ( *itr_i)->GetSeqID();
//                }
//
//                // set "last" to the seq id of the amino acid denoted by "itr_i"
//                last = ( *itr_i)->GetSeqID();
//
//                // increase "count"
//                count++;
//              }
//
//              // increase "j" and "itr_j"
//              j++;
//              itr_j++;
//
//              // if "itr_j" is not at the end of the Alignment
//              if( itr_j != ALIGNMENT.GetAssignments().End())
//              {
//                // set "itr_i" to the approprate Depth position within the Assignment denoted by "itr_j"
//                itr_i = ( *itr_j)->GetNthMember( i);
//              }
//            }
//
//            // reset "count" to zero, "j" to "start_id" and "itr_j" to "itr_start_id" for outputing the data
//            count = 0;
//            j = start_id;
//            itr_j = itr_start_id;
//
//            // set "itr_i" to the approprate Depth position within the Assignment denoted by "itr_j"
//            itr_i = ( *itr_j)->GetNthMember( i);
//
//            // print the current description; print empty description for now
//            OSTREAM << util::Format().W( 15).ForceW().L()( std::string()) << util::Format().W( 8).S().L()( first + 1);
//
//            // print the amino acids in the Alignment
//            while( count < m_BlockSize && j < ALIGNMENT.GetSize())
//            {
//              // check if gaps can be printed or 1st SiPtr of the Assignment denoted by "itr_j" points to an amino acid
//              if( !m_UnGapped || ( *itr_j)->GetMembers().Begin()->IsDefined())
//              {
//                // check if aa denoted by "itr_i" is defined
//                if( ( *itr_i).IsDefined())
//                {
//                  // print the one letter code for the amino acid currently denoted by "itr_i"
//                  OSTREAM << ( *itr_i)->GetType()->GetOneLetterCode();
//
//                  // set "last" to the seq id of the amino acid currently denoted by "itr_i"
//                  last = ( *itr_i)->GetSeqID();
//                }
//                else //< aa denoted by "itr_i" is not defined
//                {
//                  //  print dash
//                  OSTREAM << "-";
//                }
//
//                // print space every 10 AAs
//                if( !( ++count % 10))
//                {
//                  OSTREAM << ' ';
//                }
//              }
//
//              // increase "j" and "itr_j"
//              j++;
//              itr_j++;
//
//              // if "itr_j" is not at the end of the Alignment
//              if( itr_j != ALIGNMENT.GetAssignments().End())
//              {
//                // set "itr_i" to the appropriate Depth of the Assignment denoted by "itr_j"
//                itr_i = ( *itr_j)->GetNthMember( i);
//              }
//            }
//
//            // at the end of the block print the seq id of the last amino acid printed + 1
//            OSTREAM << util::Format().W( 8)( last + 1) << '\n';
//          }
//
//          // go to next block by printing end of line
//          OSTREAM << '\n';
//
//          // set "start_id" to "j", increase "cycle", and set "itr_start_id" to "itr_j"
//          start_id = j; cycle++;
//          itr_start_id = itr_j;
//        }
//        // do all this while "start_id" is not the size of the Alignment and "itr_start_id" is not at end of Alignment
//        while( start_id != ALIGNMENT.GetSize() && itr_start_id != ALIGNMENT.GetAssignments().End());

        // end
        return OSTREAM;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

    }; // template class HandlerBlocked

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_BLOCKED_H_
