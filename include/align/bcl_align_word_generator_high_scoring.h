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

#ifndef BCL_ALIGN_WORD_GENERATOR_HIGH_SCORING_H_
#define BCL_ALIGN_WORD_GENERATOR_HIGH_SCORING_H_

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
    //! @class WordGeneratorHighScoring
    //! @brief generates words by comparing all possible sequences of word length to all subsequences of the alignment
    //!
    //! @tparam t_Member the type of object that is used by the assignments in the alignment
    //!
    //! @see @link example_align_word_generator_high_scoring.cpp @endlink
    //! @author heinzes1
    //! @date Mar 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class WordGeneratorHighScoring :
      public WordGeneratorInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      score::AssignmentWithGap< t_Member> m_ScoreAssignment; //!< ScoreAssignment object

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WordGeneratorHighScoring
      (
        const score::AssignmentWithGap< t_Member> &SCORE_ASSIGNMENT = score::AssignmentWithGap< t_Member>()
      ) :
        m_ScoreAssignment( SCORE_ASSIGNMENT)
      {
      }

      //! @brief Clone function
      //! @return pointer to new WordGeneratorHighScoring< t_Member>
      WordGeneratorHighScoring< t_Member> *Clone() const
      {
        return new WordGeneratorHighScoring< t_Member>( *this);
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
        // generate all possible words and their score, and sort according to score (negative/bad scoring words first)
        storage::Vector< storage::Pair< AlignmentWord< t_Member>, double> >
          word_score_pairs( GenerateWordAndScore( QUERY_ALIGNMENT, WORD_LENGTH));
        std::sort
        (
          word_score_pairs.Begin(), word_score_pairs.End(),
          storage::PairBinaryPredicateSecond< AlignmentWord< t_Member>, double>( **math::Comparisons< double>::GetEnums().e_Less)
        );

        // create vector to store words
        storage::Vector< AlignmentWord< t_Member> > query_words;
        size_t max( QUERY_ALIGNMENT->GetSize() * 10); // limit number of words returned
        for
        (
          typename storage::Vector< storage::Pair< AlignmentWord< t_Member>, double> >::const_reverse_iterator
            itr( word_score_pairs.ReverseBegin()), // use reverse itrs to get positive scores, good words
            itr_end( word_score_pairs.ReverseEnd());
          itr != itr_end && max > 0;
          ++itr, --max
        )
        {
          query_words.PushBack( itr->First());
        }

        return query_words;
      }

    public:

      //! @brief generates a vector of words from an alignment
      //! @param QUERY_ALIGNMENT the alignment
      //! @param WORD_LENGTH the length of the words to be generated
      //! @return a vector of words
      storage::Vector< storage::Pair< AlignmentWord< t_Member>, double> > GenerateWordAndScore
      (
        const util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        const size_t WORD_LENGTH
      ) const
      {
        // check if ALIGNMENT is longer than m_WordLength, otherwise return empty word list
        if( QUERY_ALIGNMENT->GetSize() < WORD_LENGTH)
        {
          return storage::Vector< storage::Pair< AlignmentWord< t_Member>, double> >();
        }

        // list of partial words and add an empty sequence to start the extension
        storage::List< biol::AASequence> sequence_list;
        sequence_list.Append( biol::AASequence());

        // find all sequences of m_WordLength
        for( size_t pos( 0); pos < WORD_LENGTH; ++pos)
        {
          storage::List< biol::AASequence> extended_sequence_list; // to store the extended sequences

          // for each partial word
          storage::List< biol::AASequence>::iterator itr( sequence_list.Begin()), itr_end( sequence_list.End());
          for( ; itr != itr_end; ++itr)
          {
            // add a new AA
            biol::AATypes::const_iterator
              aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( ( size_t)biol::AATypes::s_NumberStandardAATypes));
            for( ; aa_itr != aa_itr_end; ++aa_itr)
            {
              // extend partial word by one letter code denoted by itr
              biol::AASequence extended_sequence( *itr);
              extended_sequence.AddMember( ( **aa_itr).GetOneLetterCode());
              // add extended partial word to new list
              extended_sequence_list.Append( extended_sequence);
            }
          }

          sequence_list = extended_sequence_list; // save extended sequences
        }

        // create vector to store word and score pairs
        storage::Vector< storage::Pair< AlignmentWord< t_Member>, double> > word_score_pairs;
        word_score_pairs.AllocateMemory
        (
          math::Pow< size_t>( ( size_t)biol::AATypes::s_NumberStandardAATypes, WORD_LENGTH)
            * ( QUERY_ALIGNMENT->GetSize() - WORD_LENGTH + 1)
        );
        // for each sequence check with all the words from the query sequence
        for
        (
          storage::List< biol::AASequence>::const_iterator itr( sequence_list.Begin()), itr_end( sequence_list.End());
          itr != itr_end;
          ++itr
        )
        {
          // for each sequence and each word in the query alignment, create an alignment and score it
          // create itrs marking begin and end of a word (the first word for now)
          typename AlignmentInterface< t_Member>::const_iterator
            itr_query_word_begin( QUERY_ALIGNMENT->GetAssignments().Begin()),
            itr_query_word_end( itr_query_word_begin);
          storage::AdvanceIterator( itr_query_word_end, QUERY_ALIGNMENT->GetAssignments().End(), WORD_LENGTH);

          // create first word before increasing itrs: create word, add assignments
          AlignmentWord< t_Member> word( QUERY_ALIGNMENT, itr_query_word_begin, itr_query_word_end);
          double score( word.ScoreSelf( m_ScoreAssignment));
          word_score_pairs.PushBack( storage::Pair< AlignmentWord< t_Member>, double>( word, score));
          do
          {
            ++itr_query_word_begin; // move word by increasing itrs
            ++itr_query_word_end;

            AlignmentWord< t_Member> word( QUERY_ALIGNMENT, itr_query_word_begin, itr_query_word_end);
            double score( word.ScoreSelf( m_ScoreAssignment));
            word_score_pairs.PushBack( storage::Pair< AlignmentWord< t_Member>, double>( word, score));
          }
          while( itr_query_word_end != QUERY_ALIGNMENT->GetAssignments().End());
        }

        return word_score_pairs;
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
        io::Serialize::Read( m_ScoreAssignment, ISTREAM);

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
        io::Serialize::Write( m_ScoreAssignment, OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    }; // template class WordGeneratorHighScoring

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> WordGeneratorHighScoring< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new WordGeneratorHighScoring< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_WORD_GENERATOR_HIGH_SCORING_H_ 
