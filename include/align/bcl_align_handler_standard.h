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

#ifndef BCL_ALIGN_HANDLER_STANDARD_H_
#define BCL_ALIGN_HANDLER_STANDARD_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_handler_interface.h"
#include "function/bcl_function_unary_interface.h"
#include "score/bcl_score_assignment_with_gap.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerStandard
    //! @brief class is used for printing Alignments in standard format;
    //! the printed secondary structure prediction information is ignored on reading in an alignment!
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @see @link example_align_handler_standard.cpp @endlink
    //! @author heinzes1
    //! @date Feb 19, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerStandard :
      public HandlerInterface< t_Member>
    {
    private:

    //////////
    // data //
    //////////

      //! ScoreAssignment object
      util::ShPtr< function::UnaryInterface< const Assignment< t_Member>, double> > m_ScoreAssignment;

      static const char        s_GapChar = '*'; //!< gap char
      static const std::string s_GapString; //!< gap string (non-integral types need to be defined outside of class
      static const size_t      s_ColumnWidth = 13; //!< width of column for input and output
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor taking a scoring function
      //! @param SCORE_ASSIGNMENT
      HandlerStandard( const util::ShPtr< function::UnaryInterface< const Assignment< t_Member>, double> > &SCORE_ASSIGNMENT) :
        m_ScoreAssignment()
      {
        SetAssignmentScore( SCORE_ASSIGNMENT);
      }

      //! @brief Clone is the virtual copy constructor
      //! @return pointer to new HandlerStandard< t_Member>
      HandlerStandard< t_Member> *Clone() const
      {
        return new HandlerStandard< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

    private:

      //! @brief gap char for a specific alignment output format
      //! @return the gap char
      char GetGapChar() const
      {
        return s_GapChar;
      }

    public:

      //! @brief get the default file extension for alignments of that handler
      //! @return string of default file extension
      const std::string &GetFileExtension() const
      {
        static const std::string s_file_extension( ".standard");
        return s_file_extension;
      }

      //! @brief set the scoring function for single assignment
      //! @param SP_ASSIGNMENT_SCORE the assignment score
      void SetAssignmentScore
      (
        const util::ShPtr< function::UnaryInterface< const Assignment< t_Member>, double> > &SP_ASSIGNMENT_SCORE
      )
      {
        if( SP_ASSIGNMENT_SCORE.IsDefined())
        {
          m_ScoreAssignment = SP_ASSIGNMENT_SCORE;
        }
        else
        {
          m_ScoreAssignment = util::CloneToShPtr( score::AssignmentWithGap< t_Member>());
        }
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief helper function printing an assignment
      //! @param OSTREAM the stream to write to
      //! @param ASSIGNMENT the assignment to write to the stream
      //! @return the stream
      std::ostream &WriteAssignment( std::ostream &OSTREAM, const Assignment< t_Member> &ASSIGNMENT) const
      {
        for
        (
          typename util::SiPtrList< const t_Member>::const_iterator
            itr( ASSIGNMENT.GetMembers().Begin()), itr_end( ASSIGNMENT.GetMembers().End());
          itr != itr_end;
          ++itr
        )
        {
          // write complete id or gap string
          OSTREAM << ' ' << ( itr->IsDefined() ? GetCompleteId< t_Member>( **itr) : s_GapString);
        }

        return OSTREAM;
      }

    public:

      //! @brief prints an Alignment in the standard format to a given stream
      //! @param OSTREAM is the stream which the Alignment is written to
      //! @param ALIGNMENT is the Alignment which is outputted in the standard format
      //! @return returns a std::ostream
      std::ostream &WriteAlignment( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        // write header and a number of space on the next line for the sequence_ids to be positioned correctly
        OSTREAM << ALIGNMENT.GetSize() << '\n' << util::Format().W( 12)( " ");

        // write sequence ids
        const util::Format seq_id_format( util::Format().L().W( s_ColumnWidth));
        const storage::List< std::string> sequence_ids( ALIGNMENT.GetSequenceIds());

        for
        (
          storage::List< std::string>::const_iterator itr( sequence_ids.Begin()), itr_end( sequence_ids.End());
          itr != itr_end;
          ++itr
        )
        {
          std::string short_sequence_id( *itr);
          // if too long, use substr and add '>' to indicate longer sequence_id
          if( short_sequence_id.size() > s_ColumnWidth)
          {
            short_sequence_id = short_sequence_id.substr( 0, s_ColumnWidth - 1) + '>';
          }
          OSTREAM << ' ' << seq_id_format( short_sequence_id);
        }

        // print new line
        OSTREAM << '\n';

        // create size_t "i" for keeping track of the number of Assignments that have been printed out so far
        size_t alignment_position( 0); // initialize to zero

        const util::Format align_pos_format( util::Format().W( 5));
        const util::Format assign_score_format( util::Format().W( 6).FFP( 3));

        // loop over ALIGNMENT, printing all the Assignments
        for
        (
          typename AlignmentLeaf< t_Member>::const_iterator
            assign_itr( ALIGNMENT.GetAssignments().Begin()), assign_itr_end( ALIGNMENT.GetAssignments().End());
          assign_itr != assign_itr_end;
          ++assign_itr, ++alignment_position
        )
        {
          // print the number of Assignments that have been printed so far
          OSTREAM << align_pos_format( alignment_position + 1) <<
            " " << assign_score_format( m_ScoreAssignment->operator()( **assign_itr));

          // print the Assignment currently denoted by "assign_itr"
          WriteAssignment( OSTREAM, **assign_itr);

          // print new line
          OSTREAM << '\n';
        }

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

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads the alignment from the istream into a vector of pairs of identifier, gapped sequence
      //! @param ISTREAM the stream to read from
      //! @return a vector of pairs of identifier, gapped sequence
      storage::Vector< storage::VectorND< 2, std::string> > ReadAlignmentToStrings( std::istream &ISTREAM) const
      {
        // create string to temporarily store content, and result vector
        std::string content;
        storage::Vector< storage::VectorND< 2, std::string> > gapped_sequences;

        // read in length of the alignment (use signed integer to detect reading more lines than specified in alignment)
        long long alignment_length;
        ISTREAM >> alignment_length;
        std::getline( ISTREAM, content); // read to the end of this line, throw it away

        // if nothing is left in the istream, return
        if( ISTREAM.eof())
        {
          return storage::Vector< storage::VectorND< 2, std::string> >();
        }

        // read first line containing identifiers of all sequences
        std::getline( ISTREAM, content);
        const size_t column_offset( s_ColumnWidth); // offset of the first column

        const storage::Vector< std::string> sequence_id_vec( util::SplitString( content));
        if( sequence_id_vec.IsEmpty() || ( s_ColumnWidth + 1) * ( sequence_id_vec.GetSize() - 1) + column_offset > content.size())
        {
          BCL_MessageCrt( "line with sequence identifier too short: " + content);
          return storage::Vector< storage::VectorND< 2, std::string> >();
        }

        for( size_t pos( 0), end_pos( sequence_id_vec.GetSize()); pos < end_pos; ++pos)
        {
          const size_t content_pos( ( s_ColumnWidth + 1) * pos + column_offset);
          std::string sequence_id( content.substr( content_pos, s_ColumnWidth));
          storage::VectorND< 2, std::string> gapped_sequence( sequence_id, ""); // identifier and empty gapped sequence
          gapped_sequences.PushBack( gapped_sequence);
        }

        const size_t col_olc( 4); // one letter code is every 4 colums
        const size_t col_offset_olc( 3); // first one letter code is in column #3
        // read all further lines containing the aligned sequences
        while( !ISTREAM.eof() && std::getline( ISTREAM, content).good())
        {
          --alignment_length; // to check reading all lines of the alignment
          if( alignment_length < 0) // stop processing lines if we already found the correct number
          {
            BCL_MessageCrt( "Ignore lines in alignment file over specified alignment length!");
            break;
          }

          storage::Vector< std::string> assignment_vec( util::SplitString( content));
          for( size_t seq_nr( 0), seq_nr_max( gapped_sequences.GetSize()); seq_nr < seq_nr_max; ++seq_nr)
          {
            // read one letter code from the column in assignment_vec for the respective sequence
            gapped_sequences( seq_nr).Second().append( assignment_vec( col_olc * seq_nr + col_offset_olc));
          }
        }

        if( alignment_length > 0) // print warning if not all lines were found
        {
          BCL_MessageCrt
          (
            "Incomplete alignment file: Found " + util::Format()( gapped_sequences( 0).Second().size())
              + ", expected " + util::Format()( gapped_sequences( 0).Second().size() + alignment_length) + " lines."
          );
        }

        return gapped_sequences;
      }

    }; // template class HandlerStandard

    // gap char for sequence pos, 1-letter code, 3-letter code, ss prediction
    template< typename t_Member>
    const std::string HandlerStandard< t_Member>::s_GapString =
      std::string( "    ") + s_GapChar + " " + s_GapChar + " " + s_GapChar + s_GapChar + s_GapChar + " " + s_GapChar;

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_STANDARD_H_
