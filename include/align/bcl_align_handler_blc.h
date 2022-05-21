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

#ifndef BCL_ALIGN_HANDLER_BLC_H_
#define BCL_ALIGN_HANDLER_BLC_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_handler_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerBLC
    //! @brief class writes Alignments in BLC format which can be read by the AMAS server located at
    //! http://www.compbio.dundee.ac.uk/manuals/amas/amas_manual.txt
    //!
    //! @tparam t_Member is the type of object that the Assignment stores
    //!
    //! @see @link example_align_handler_blc.cpp @endlink
    //! @author heinzes1
    //! @date 03/09/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerBLC :
      public HandlerInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      static const char s_GapChar = ' '; //!< gap char
      static const char s_IdentifierBeginChar = '>'; //!< char used the signal the begin of an identifier
      static const char s_IterationChar = '*'; //!< char to signal start and end of an iteration

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone is the virtual copy constructor
      //! @return pointer to new HandlerBLC< t_Member>
      HandlerBLC< t_Member> *Clone() const
      {
        return new HandlerBLC< t_Member>( *this);
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
        static const std::string s_file_extension( ".blc");
        return s_file_extension;
      }

    private:

      //! @brief gap char for a specific alignment output format
      //! @return the gap char
      char GetGapChar() const
      {
        return s_GapChar;
      }

    public:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief WriteAlignment prints an alignment to a given stream
      //! @param OSTREAM the stream which the alignment is written to
      //! @param ALIGNMENT the alignment which is outputted
      //! @return the passed std::ostream
      std::ostream &WriteAlignment( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        // write sequence identifiers
        WriteBLCHeader( OSTREAM, ALIGNMENT);

        // write iteration begin, alignment and iteration end
        OSTREAM << s_IterationChar << " iteration 1\n";
        WriteBLCAlignment( OSTREAM, ALIGNMENT);
        OSTREAM << s_IterationChar << '\n';

        return OSTREAM;
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief reads the alignment from the istream into a vector of pairs of identifier, gapped sequence
      //! @param ISTREAM the stream to read from
      //! @return a vector of pairs of identifier, gapped sequence
      storage::Vector< storage::VectorND< 2, std::string> > ReadAlignmentToStrings( std::istream &ISTREAM) const
      {
        // create string to temporarily store content, and result vector
        std::string content;
        storage::Vector< storage::VectorND< 2, std::string> > gapped_sequences;

        // read until the first fasta header '>'; everything before it is comment and ignored
        while( !ISTREAM.eof() && ISTREAM.peek() != s_IdentifierBeginChar)
        {
          std::getline( ISTREAM, content);
        }

        // read fasta header, one per line, starting with '>'
        while( !ISTREAM.eof() && ISTREAM.peek() == s_IdentifierBeginChar)
        {
          std::getline( ISTREAM, content);
          content.erase( 0, 1); // remove '>'
          storage::VectorND< 2, std::string> gapped_sequence( content, "");
          gapped_sequences.PushBack( gapped_sequence);
        }

        // read '* iteration N', return empty result if not found
        std::getline( ISTREAM, content);
        std::string search_string( 1, s_IterationChar);
        search_string.append( " iteration ");
        if( ISTREAM.eof() || content.find( search_string) == std::string::npos)
        {
          return storage::Vector< storage::VectorND< 2, std::string> >();
        }

        // read assignments until next '* iteration N+1' or end '*'
        while( !ISTREAM.eof() && ISTREAM.peek() != s_IterationChar)
        {
          std::getline( ISTREAM, content); // '\n' is not in content, so content.length() is #sequences aligned
          for( size_t seq_nr( 0), seq_nr_max( content.length()); seq_nr < seq_nr_max; ++seq_nr)
          {
            gapped_sequences( seq_nr).Second().append( 1, content[ seq_nr]);
          }
        }

        return gapped_sequences;
      }

      //! @brief writes the blc header to a ostream
      //! @param OSTREAM output stream to write to
      //! @param ALIGNMENT the alignment to extract the header from
      void WriteBLCHeader( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        const storage::List< std::string> sequence_ids( ALIGNMENT.GetSequenceIds());

        // print sequence identifiers
        for
        (
          storage::List< std::string>::const_iterator itr( sequence_ids.Begin()), itr_end( sequence_ids.End());
          itr != itr_end;
          ++itr
        )
        {
          OSTREAM << s_IdentifierBeginChar << *itr << '\n';
        }
      }

      //! @brief writes a core alignment in blc format
      //! @param OSTREAM output stream to write to
      //! @param ALIGNMENT the alignment to write out
      void WriteBLCAlignment( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        // write alignment, assignments are in a line, sequences are written vertically
        for
        (
          typename util::ShPtrList< Assignment< t_Member> >::const_iterator
            assignment_itr( ALIGNMENT.GetAssignments().Begin()),
            assignment_itr_end( ALIGNMENT.GetAssignments().End());
          assignment_itr != assignment_itr_end;
          ++assignment_itr
        )
        {
          // print the assignment as string and a new line afterwards
          OSTREAM << ( **assignment_itr).ToString( GetGapChar()) << '\n';
        }
      }

    }; // template class HandlerBLC

  } // namespace align
} // namespace bcl

#endif //BCL_ALIGN_HANDLER_BLC_H_
