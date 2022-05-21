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

#ifndef BCL_ALIGN_WORD_GENERATOR_SUBSEQUENCES_H_
#define BCL_ALIGN_WORD_GENERATOR_SUBSEQUENCES_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_word_generator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WordGeneratorSubsequences
    //! @brief generates a vector of words from an alignment by taking all possible subsequences
    //!
    //! @tparam t_Member the type of object that is used by the assignments in the alignment
    //!
    //! @see @link example_align_word_generator_subsequences.cpp @endlink
    //! @author heinzes1
    //! @date Mar 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class WordGeneratorSubsequences :
      public WordGeneratorInterface< t_Member>
    {

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WordGeneratorSubsequences()
      {
      }

      //! @brief Clone function
      //! @return pointer to new WordGeneratorSubsequences< t_Member>
      WordGeneratorSubsequences< t_Member> *Clone() const
      {
        return new WordGeneratorSubsequences< t_Member>( *this);
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

    ////////////////
    // operations //
    ////////////////

      //! @brief generates a vector of words from an alignment
      //! @param QUERY_ALIGNMENT the alignment
      //! @param WORD_LENGTH the length of the words to be generated
      //! @return a vector of words
      storage::Vector< AlignmentWord< t_Member> > Generate
      (
        const util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        const size_t WORD_LENGTH
      ) const
      {
        // check if ALIGNMENT is longer than m_WordLength, otherwise return empty word list
        if( QUERY_ALIGNMENT->GetSize() < WORD_LENGTH)
        {
          return storage::Vector< AlignmentWord< t_Member> >();
        }

        // create vector to store words
        storage::Vector< AlignmentWord< t_Member> > query_words;

        // create itrs marking begin and end of a word (the first word for now)
        typename AlignmentInterface< t_Member>::const_iterator
          itr_query_word_begin( QUERY_ALIGNMENT->GetAssignments().Begin()),
          itr_query_word_end( itr_query_word_begin);
        storage::AdvanceIterator( itr_query_word_end, QUERY_ALIGNMENT->GetAssignments().End(), WORD_LENGTH);

        // create first word before increasing itrs: create word, add assignments
        query_words.PushBack( AlignmentWord< t_Member>( QUERY_ALIGNMENT, itr_query_word_begin, itr_query_word_end));
        do
        {
          ++itr_query_word_begin; // move word by increasing itrs
          ++itr_query_word_end;

          query_words.PushBack( AlignmentWord< t_Member>( QUERY_ALIGNMENT, itr_query_word_begin, itr_query_word_end));
        }
        while( itr_query_word_end != QUERY_ALIGNMENT->GetAssignments().End());

        return query_words;
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
        // no members, just return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // no members, just return the stream
        return OSTREAM;
      }

    }; // template class WordGeneratorSubsequences

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> WordGeneratorSubsequences< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new WordGeneratorSubsequences< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_WORD_GENERATOR_SUBSEQUENCES_H_ 
