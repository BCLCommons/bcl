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

#ifndef BCL_ALIGN_HANDLER_FASTA_H_
#define BCL_ALIGN_HANDLER_FASTA_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_handler_interface.h"
#include "biol/bcl_biol_aa_sequence_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerFasta
    //! @brief reads an multiple sequence alignment from a fasta file
    //! http://www.bioperl.org/wiki/FASTA_multiple_alignment_format
    //!
    //! @see @link example_align_handler_fasta.cpp @endlink
    //! @author heinzes1
    //! @date Apr 12, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerFasta :
      public HandlerInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      static const char   s_GapChar             = '-'; //!< gap char
      static const char   s_IdentifierBeginChar = '>'; //!< char used the signal the begin of an identifier
      static const size_t s_BlockSize           = 60;  //!< default number of characters per line

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new HandlerFasta
      HandlerFasta *Clone() const
      {
        return new HandlerFasta< t_Member>( *this);
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
        static const std::string s_file_extension( ".fasta");
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
          HandlerInterface< t_Member>::GetAlignedSequenceStrings( ALIGNMENT, s_GapChar)
        );

        // print every sequence
        storage::List< std::string>::const_iterator id_itr( sequence_ids.Begin()), id_itr_end( sequence_ids.End());
        storage::List< std::string>::const_iterator
          aligned_seq_itr( aligned_sequence.Begin()),
          aligned_seq_itr_end( aligned_sequence.End());
        for( ; id_itr != id_itr_end && aligned_seq_itr != aligned_seq_itr_end; ++id_itr, ++aligned_seq_itr)
        {
          WriteFastaHeader( OSTREAM, *id_itr);
          WriteFastaAlignedSequence( OSTREAM, *aligned_seq_itr);
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
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
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
        // create result vector
        storage::Vector< storage::VectorND< 2, std::string> > gapped_sequences;

        // read all sequences until end of file
        while( !ISTREAM.eof())
        {
          biol::AASequence sequence( biol::AASequenceFactory::BuildSequenceFromFASTA( ISTREAM));
          gapped_sequences.PushBack( storage::VectorND< 2, std::string>( sequence.GetSequenceId(), sequence.Sequence()));
        }

        return gapped_sequences;
      }

      //! @brief write the fasta header for one aligned sequence from the given sequence identifier
      //! @param OSTREAM the stream to write to
      //! @param ID the identifier to write
      void WriteFastaHeader( std::ostream &OSTREAM, const std::string &ID) const
      {
        // if the type+id does not start with '>', write it
        if( ID[ 0] != s_IdentifierBeginChar)
        {
          OSTREAM << s_IdentifierBeginChar;
        }

        // write the complete identifier
        OSTREAM // << "P1;"
                << ID << '\n';
      }

      //! @brief write the given aligned sequence in blocks
      //! @param OSTREAM the stream to write to
      //! @param ALIGNED_SEQUENCE the sequence to write out
      void WriteFastaAlignedSequence( std::ostream &OSTREAM, const std::string &ALIGNED_SEQUENCE) const
      {
        // make a copy of the sequence to cut in blocks
        std::string sequence( ALIGNED_SEQUENCE);

        // write next line until the remaining sequence is shorter than blocksize - 1 (the end char has to fit too)
        while( sequence.length() > s_BlockSize)
        {
          OSTREAM << sequence.substr( 0, s_BlockSize);
          sequence = sequence.substr( s_BlockSize);
          OSTREAM << '\n';
        }

        // write the rest and the end char
        OSTREAM << sequence << '\n';
      }

    }; // template class HandlerFasta

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_FASTA_H_
