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

#ifndef BCL_ALIGN_HANDLER_PIR_H_
#define BCL_ALIGN_HANDLER_PIR_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_handler_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerPIR
    //! @brief class is used for printing Alignments in the PIR format
    //! http://www.bioinformatics.nl/tools/crab_pir.html
    //! http://salilab.org/modeller/9v4/manual/node438.html
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @see @link example_align_handler_pir.cpp @endlink
    //! @author heinzes1, teixeipl
    //! @date 03/09/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerPIR :
      public HandlerInterface< t_Member>
    {
    private:

    //////////
    // data //
    //////////

      size_t m_BlockSize; //!< the number of characters on a line
      char m_Delim; //!< the delimiter used between lines during write/read
      static const char s_Delim = '\n'; //!< the default delimiter used, set to newline

      static const char   s_GapChar             = '-'; //!< gap char
      static const char   s_IdentifierBeginChar = '>'; //!< char used the signal the begin of an identifier
      static const char   s_IdentifierSplitChar = ' '; //!< char used to split the identifier in id and description
      static const char   s_SequenceEndChar     = '*'; //!< char to signal the end of a sequence
      static const size_t s_BlockSize           =  70; //!< default number of characters per line

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief (default) constructor with parameter for defining the width of the block
      //! @param BLOCK_SIZE the width of the block, default s_BlockSize characters
      HandlerPIR( const size_t BLOCK_SIZE = s_BlockSize) :
        m_BlockSize( BLOCK_SIZE > 0 ? BLOCK_SIZE : size_t( s_BlockSize)),
        m_Delim( s_Delim)
      {
      }

      //! @brief constructor with delimiter parameter
      //! @param DELIM the delimiter to be used when writing/reading, default is newline
      explicit HandlerPIR( const char DELIM) :
        m_Delim( DELIM)
      {
      }

      //! @brief Clone is the virtual copy constructor
      //! @return pointer to new HandlerPIR< t_Member>
      HandlerPIR< t_Member> *Clone() const
      {
        return new HandlerPIR< t_Member>( *this);
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
        static const std::string s_file_extension( ".pir");
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

      //! @brief returns the block size
      //! @return the block size
      size_t GetBlockSize() const
      {
        return m_BlockSize;
      }

      //! @brief sets the block size
      //! @param BLOCK_SIZE the new block size
      void SetBlockSize( const size_t BLOCK_SIZE)
      {
        m_BlockSize = BLOCK_SIZE > 0 ? BLOCK_SIZE : size_t( s_BlockSize);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief WriteAlignment prints an alignment to a given stream
      //! @param OSTREAM the stream which the alignment is written to
      //! @param ALIGNMENT the alignment which is outputted
      //! @return the passed std::ostream
      std::ostream &WriteAlignment( std::ostream &OSTREAM, const AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        const storage::List< std::string> sequence_ids( ALIGNMENT.GetSequenceIds());
        const storage::List< std::string> aligned_sequence
        (
          HandlerInterface< t_Member>::GetAlignedSequenceStrings( ALIGNMENT, GetGapChar())
        );

        // print every sequence
        storage::List< std::string>::const_iterator id_itr( sequence_ids.Begin()), id_itr_end( sequence_ids.End());
        storage::List< std::string>::const_iterator
          aligned_seq_itr( aligned_sequence.Begin()),
          aligned_seq_itr_end( aligned_sequence.End());
        for( ; id_itr != id_itr_end && aligned_seq_itr != aligned_seq_itr_end; ++id_itr, ++aligned_seq_itr)
        {
          WritePIRHeader( OSTREAM, *id_itr);
          WritePIRAlignedSequence( OSTREAM, *aligned_seq_itr);
        }

        return OSTREAM;
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
        io::Serialize::Read( m_BlockSize, ISTREAM);
        io::Serialize::Read( m_Delim, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_BlockSize, OSTREAM, INDENT);
        io::Serialize::Write( m_Delim, OSTREAM, INDENT);
        return OSTREAM;
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
        // create result vector
        storage::Vector< storage::VectorND< 2, std::string> > gapped_sequences;

        // read all sequences until end of file
        while( !ISTREAM.eof())
        {
          gapped_sequences.PushBack( ReadSingleAlignmentToString( ISTREAM));
        }

        return gapped_sequences;
      }

      //! @brief read a single aligned sequence, i.e. its identifier, description and sequence, from the given stream
      //! @param ISTREAM the stream to read from
      //! @return a vector containing id+description and sequence
      storage::VectorND< 2, std::string> ReadSingleAlignmentToString( std::istream &ISTREAM) const
      {
        // create string to temporarily store content
        std::string content;

        // read until start of identifier and description
        while
        (
          !ISTREAM.eof() && ISTREAM.peek() != s_IdentifierBeginChar && std::getline( ISTREAM, content, m_Delim).good()
        );

        // if nothing is left in the istream, return
        if( ISTREAM.eof())
        {
          return storage::VectorND< 2, std::string>();
        }

        // read identifier
        BCL_Assert
        (
          std::getline( ISTREAM, content).good(),
          "Could not read identifier of pir file!, got " + content + " instead of " + std::string( 1, m_Delim)
        );
        content.erase( 0, 1); // remove '>'
        storage::VectorND< 2, std::string> gapped_sequence( content, "");

        // read description
        BCL_Assert( std::getline( ISTREAM, content, m_Delim).good(), "Could not read description of pir file!");
        gapped_sequence.First().append( content);

        // read sequence until start of next sequence identifier: at least one line of sequence has to be in the stream
        while( !ISTREAM.eof() && ISTREAM.peek() != s_IdentifierBeginChar && std::getline( ISTREAM, content, m_Delim).good())
        {
          // delete '*' at the end if present
          if( !content.empty() && content[ content.length() - 1] == s_SequenceEndChar)
          {
            content.erase( content.length() - 1, 1);
          }
          // removed spaces and add to previous found sequence
          gapped_sequence.Second().append( util::RemoveSpacesFromString( content));
        }

        return gapped_sequence;
      }

      //! @brief write the pir header for one aligned sequence from the given sequence identifier
      //! @param OSTREAM the stream to write to
      //! @param ID the identifier to write
      void WritePIRHeader( std::ostream &OSTREAM, const std::string &ID) const
      {
        std::string sequence_id, sequence_description;

        // trim id and split into type+id and description
        const std::string trimmed_id( util::TrimString( ID));
        size_t pos( trimmed_id.find_first_of( s_IdentifierSplitChar));
        if( pos == std::string::npos)
        {
          sequence_id = trimmed_id;
        }
        else
        {
          sequence_id = trimmed_id.substr( 0, pos);
          sequence_description = trimmed_id.substr( pos + 1, trimmed_id.length() - pos - 1);
        }

        // if the type+id does not start with '>', write it
        if( sequence_id.length() > 0 && sequence_id[ 0] != s_IdentifierBeginChar)
        {
          OSTREAM << s_IdentifierBeginChar;
        }

        // write the complete identifier in two lines; hard code the sequence type as 'P1'
        OSTREAM //<< "P1;"
                << sequence_id << m_Delim << sequence_description << m_Delim;
      }

      //! @brief write the given aligned sequence in blocks
      //! @param OSTREAM the stream to write to
      //! @param ALIGNED_SEQUENCE the sequence to write out
      void WritePIRAlignedSequence( std::ostream &OSTREAM, const std::string &ALIGNED_SEQUENCE) const
      {
        // make a copy of the sequence to cut in blocks
        std::string sequence( ALIGNED_SEQUENCE);

        // write next line until the remaining sequence is shorter than blocksize - 1 (the end char has to fit too)
        while( sequence.length() >= m_BlockSize)
        {
          OSTREAM << sequence.substr( 0, m_BlockSize);
          sequence = sequence.substr( m_BlockSize);
          OSTREAM << m_Delim;
        }

        // write the rest and the end char
        OSTREAM << sequence << s_SequenceEndChar << m_Delim;
      }

    }; // template class HandlerPIR

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_PIR_H_
