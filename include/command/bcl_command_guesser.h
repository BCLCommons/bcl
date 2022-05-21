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

#ifndef BCL_COMMAND_GUESSER_H_
#define BCL_COMMAND_GUESSER_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Guesser
    //! @brief a fuzzy string matcher with various heuristics to guess which option a human user intended
    //!
    //! @see @link example_command_guesser.cpp @endlink
    //! @author mendenjl, teixeipl
    //! @date Nov 17, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Guesser
    {

    public:

    //////////
    // enum //
    //////////

      //! Type of mismatch, ordered by increasing fuzziness (e.g. SomeStemsMatch will match a lot more than FirstLetters)
      enum MismatchType
      {
        e_CaseOrSpace,        //!< Case or space was wrong (e.g. ALPHA for option a_l_p_h_a)
        e_DefinedAlias,       //!< matches externally provided alias
        e_FirstLetters,       //!< Given argument is a prefix used in the options (e.g. Ch for Chemistry, ChemistryAtom)
        e_Suffix,             //!< Suffix differed (e.g. AwesomeFactories for option AwesomeFactory)
        e_SuffixCaseSpace,    //!< Suffix and case or space differed (AweSomeFactories for option Awesome Factory)
        e_Stems,              //!< all word stems in common, in the same order, and with the same case (e.g. AF or AFact for option AwesomeFactory)
        e_ReorderedWords,     //!< Reordered words (e.g. BetaAlpha for option AlphaBeta)
        e_ReorderedStems,     //!< Reordered stems (e.g. DogCatchers for option CatchDogs)
        e_Explicit,           //!< Inverse of abbreviation
        e_StrongAbbreviation, //!< Abbreviation, respecting case or other delimitation ( AlMgCa for option Aluminum Magnesium Carbon)
        e_WeakAbbreviation,   //!< Abbreviation, case / space insensitive (e.g. almgca for option Aluminum Magnesium Carbon)
        e_SomeWordsMatch,     //!< Some words were omitted, possibly also reordered (e.g. AbEf for option CdEfAb)
        e_SomeStemsMatch,     //!< Some stems match, and may also be reordered (e.g. ObscenelyLong for option LocatorForLongerObscenity)
        s_NumberTypes
      };

      //! @brief get the type name
      //! @param TYPE the type name
      //! @return string describing the type
      static const std::string &GetTypeName( const MismatchType &TYPE);

      //! @brief enum class wrapper for MismatchType
      typedef util::WrapperEnum< MismatchType, &GetTypeName, s_NumberTypes> TypeEnum;

    private:

    //////////
    // data //
    //////////

      //! map of defined aliases; aliases must match exactly (except for case and spacing) to be considered
      //! This should be used for cases when there is an alias for a string that do not share a common word stem
      storage::Map< std::string, std::string> m_Aliases;

      //! Vector of defined suffixes
      storage::Vector< storage::Map< std::string, std::string> > m_Suffixes;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief get the single instance of this class
      //! @return the single instance of this class
      static const Guesser &GetDefaultGuesser();

    ////////////////
    // operations //
    ////////////////

      //! @brief Guess what the user intended
      //! @param ARGUMENT the argument that was passed into the command line
      //! @param EXPECTED known, valid arguments
      //! @return a pair from best MismatchType to the strings that fell under that mismatch type
      storage::Pair< TypeEnum, storage::Vector< std::string> > Guess
      (
        const std::string &ARGUMENT,
        const storage::Vector< std::string> &EXPECTED
      ) const;

      //! @brief Guess what the user intended and write those guesses to STREAM
      //! @param ARGUMENT the argument that was passed into the command line
      //! @param EXPECTED known, valid arguments
      //! @param STREAM stream to write out the guesses to
      //! @param ARGUMENT_TYPE_NAME conceptual type of the argument
      void WriteGuesses
      (
        const std::string &ARGUMENT,
        const storage::Vector< std::string> &EXPECTED,
        std::ostream &STREAM,
        const std::string &ARGUMENT_TYPE_NAME
      ) const;

    ////////////////
    // operations //
    ////////////////

    public:

      //! @brief Registers the suffixes that have alternate forms that should be accepted
      //! @param SUFFIX to be normalized
      //! @param NORMALIZED suffix to create stem word
      //! @return a reference to the Guesser modified
      Guesser &RegisterSuffix( const std::string &SUFFIX, const std::string &NORMALIZED);

      //! @brief Protect a suffix from being changed if it is the longest matching suffix
      //! @param SUFFIX to be normalized
      //! @return a reference to the Guesser modified
      Guesser &ProtectSuffix( const std::string &SUFFIX);

      //! @brief Registers the aliases to be used
      //! @param ORIGINAL word
      //! @param ALIAS for original word
      //! @return a reference to the Guesser modified
      Guesser &RegisterAlias( const std::string &ORIGINAL, const std::string &ALIAS);

    private:

      //! @brief Applies suffix rules to get down to the word stem
      //! @param WORD to be stemmed
      //! @return word with the suffix standardized
      std::string GetStem( const std::string &WORD) const;

      //! @brief GetStem for an entire vector
      //! @param WORDS words to be stemmed
      //! @return the word stems
      storage::Vector< std::string> GetStems( const storage::Vector< std::string> &WORDS) const;

      //! @brief Splits input string into words based on white space characters, capitalization, and underscores
      //! @param STRING to be split into words
      //! @return Vector of split words
      storage::Vector< std::string> SplitIntoWords( const std::string &STRING) const;

      //! @brief Normalizes a string to remove capitalization, spacing, and underscores
      //! @param WORD to be normalized
      //! @return normalized word
      std::string NormalizeWord( const std::string &WORD) const;

      //! @brief Splits input string into words based on white space characters, capitalization, and underscores
      //! @param STRING to be split into words, letters of acronyms as separate entities
      //! @return Vector of split words, letters of acronyms as separate entities
      static storage::Vector< std::string> FullSplitIntoWords( const std::string &STRING);

      //! @brief Changes the capitalization of a word to lower
      //! @param WORD to be lower-cased
      static void ToLower( std::string &WORD);

      //! @brief Get the number of characters that match at the start of the string
      //! @param A, B the two strings of interest
      //! @return the number of letters at the start of the strings that match
      static size_t GetMatchingPrefixLength( const std::string &A, const std::string &B);

      //! @brief Test whether one string is the abbreviation of another
      //! @param A, B the two strings of interest
      //! @return true if A is an abbreviation prefix of B or B is an abbreviation of A
      static bool IsAbbreviation( const std::string &A, const std::string &B);

      //! @brief Determine whether a given string is obviously an acronym
      //! @param A the string of interest
      //! @return true if A is > 1 letter, and contains all vowels, or all consonants
      //! @note if y's are present, they count as a vowel only if all the rest of the word is consonants
      //! @note if the A is just a sequences of Y's, it counts as an acronym
      static bool IsObviousAcroynm( const std::string &A);

      //! @brief Test whether one string is an exact abbreviation / substring of another, with arbitrary number of gaps
      //! @param SUBSTR the posited substring
      //! @param STRING the string to test whether SUBSTR is a subsequence of
      //! @param START position to start at in STRING
      //! @return # of non-consecutive gaps, or undefined if there is no matching sub-sequence
      static size_t MatchingSequenceWithGaps
      (
        const std::string &WITH_GAPS,
        const std::string &STRING,
        const size_t &START = 0
      );

    }; //class Guesser

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_GUESSER_H_
