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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "command/bcl_command_guesser.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_csi_substructure.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace command
  {
    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &Guesser::GetTypeName( const MismatchType &TYPE)
    {
      static const std::string s_descriptors[ s_NumberTypes + 1] =
      {
        "CaseOrSpace",
        "DefinedAlias",
        "FirstLetters",
        "Suffix",
        "SuffixCaseOrSpace",
        "Stems",
        "ReorderedWords",
        "ReorderedStems",
        "Explict",
        "StrongAbbreviation",
        "WeakAbbreviation",
        "SomeWordsMatch",
        "SomeStemsMatch",
        GetStaticClassName< MismatchType>()
      };
      return s_descriptors[ size_t( TYPE)];
    }
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief get the single instance of this class
    //! @return the single instance of this class
    const Guesser &Guesser::GetDefaultGuesser()
    {
      // single instance of this class, initialized with the desired replacements
      static Guesser s_instance
      (
        Guesser()
        // suffixes that should never be changed
        .ProtectSuffix( "ss")     // ss is already a depluralized suffix
        .ProtectSuffix( "is")     // e.g. this
        .ProtectSuffix( "us")     // e.g. thus
        .ProtectSuffix( "se")     // e.g. these, portugese, geese
        .ProtectSuffix( "ite")    // e.g. favorite
        // suffixes that should be changed (order is irrelevant here)
        .RegisterSuffix( "s", "")        // Rods -> Rod
        .RegisterSuffix( "ed", "")       // Checked -> Check
        .RegisterSuffix( "ies", "y")     // Properties -> Property
        .RegisterSuffix( "ied", "y")     // Tried -> Try
        .RegisterSuffix( "oes", "o")     // potatoes -> potato
        .RegisterSuffix( "xes", "x")     // foxes -> fox
        .RegisterSuffix( "ves", "fe")    // lives -> life
        .RegisterSuffix( "lves", "lf")   // wolves -> wolf, calves -> calf
        .RegisterSuffix( "rices", "rix") // matrices -> matrix
        .RegisterSuffix( "sses", "ss")   // caresses  ->  caress
        .RegisterSuffix( "shes", "sh")   // bashes  ->  bash
        .RegisterSuffix( "ches", "ch")   // patches  ->  patch
        .RegisterSuffix( "xes", "x")     // boxes  ->  box
        .RegisterSuffix( "zzes", "zz")   // buzzes  ->  buzzes
        // These are usually suffixes, though there are many exceptions that a more sophisticated word stemming
        // algorithm would need to handle.  As the name of the class suggests, this is all heuristics anyway, so
        // if some app needs more precision, please create a separate word stemming class instead.
        .RegisterSuffix( "er", "")
        .RegisterSuffix( "ced", "ce")
        .RegisterSuffix( "sed", "se")
        .RegisterSuffix( "hed", "he")
        .RegisterSuffix( "ared", "are")
        .RegisterSuffix( "ered", "ere")
        .RegisterSuffix( "ored", "ore")
        .RegisterSuffix( "aned", "ane")
        .RegisterSuffix( "ened", "ene")
        .RegisterSuffix( "ined", "ine")
        .RegisterSuffix( "oned", "one")
        .RegisterSuffix( "ured", "ure")
        .RegisterSuffix( "ers", "")
        .RegisterSuffix( "or", "")
        .RegisterSuffix( "ing", "")
        .RegisterSuffix( "ors", "")
        .RegisterSuffix( "ion", "")
        .RegisterSuffix( "est", "")
        .RegisterSuffix( "able", "")
        .RegisterSuffix( "ible", "")
        .RegisterSuffix( "iest", "y")
        .RegisterSuffix( "ize", "")
        .RegisterSuffix( "ment", "")
        .RegisterSuffix( "ments", "")
        .RegisterSuffix( "able", "")
        .RegisterSuffix( "ful", "")
        .RegisterSuffix( "ions", "")
        .RegisterSuffix( "ness", "")
        .RegisterSuffix( "ly", "")
        .RegisterSuffix( "mic", "m")
        .RegisterSuffix( "nessly", "")
        .RegisterSuffix( "ity", "e")
        .RegisterSuffix( "ative", "ate")
        .RegisterSuffix( "ysis", "yze")
      );

      // ultimately word stemming is probably not optimal to determine user intent
      // An ideal approach would use alignments, with a score that decreases by distance from the start of the word

      return s_instance;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Guess what the user intended
    //! @param ARGUMENT the argument that was passed into the command line
    //! @param EXPECTED known, valid arguments
    //! @return a map from MismatchType to the strings that fell under that mismatch type
    storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> > Guesser::Guess
    (
      const std::string &ARGUMENT,
      const storage::Vector< std::string> &EXPECTED
    ) const
    {
      // create an object for the return value
      storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> > match_type_and_closest_matches;
      Guesser::TypeEnum &match_type( match_type_and_closest_matches.First());
      storage::Vector< std::string> &closest_matches( match_type_and_closest_matches.Second());

      // Normalize argument
      const std::string normalized_arg( NormalizeWord( ARGUMENT));

      // a vector for holding the normalized expected version of each word, to avoid recalculating
      storage::Vector< std::string> normalized_expected;
      const size_t expected_size( EXPECTED.GetSize());
      normalized_expected.AllocateMemory( expected_size);
      {
        // Check for CaseOrSpace case
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string normalized_option( NormalizeWord( *itr));
          if( normalized_option == normalized_arg)
          {
            // just differed in cases or spaces
            closest_matches.PushBack( *itr);
          }
          else if( closest_matches.IsEmpty())
          {
            // cache the normalized string for further use
            normalized_expected.PushBack( normalized_option);
          }
        }
        if( !closest_matches.IsEmpty())
        {
          // options differing in only case or spaces were found, so return them
          match_type = e_CaseOrSpace;
          return match_type_and_closest_matches;
        }
      }

      // Check for DefinedAlias case
      storage::Map< std::string, std::string>::const_iterator itr( m_Aliases.Find( normalized_arg));
      if( itr != m_Aliases.End() && EXPECTED.Find( itr->second) < expected_size)
      {
        // a defined alias was found and is also present in the expected strings, so return it
        closest_matches.PushBack( itr->second);
        match_type = e_DefinedAlias;
        return match_type_and_closest_matches;
      }

      // check for first letters case
      {
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          // check whether the given word is a prefix of any options
          if( GetMatchingPrefixLength( *itr_norm, normalized_arg) == normalized_arg.size())
          {
            closest_matches.PushBack( *itr);
          }
        }

        // check for first letters case
        if( !closest_matches.IsEmpty())
        {
          // users first letters all matched the beginning of something, so return that
          match_type = e_FirstLetters;
          return match_type_and_closest_matches;
        }
      }

      // Suffix case
      // check for matches that just differ in suffix by using the word stemming algorithm
      {
        const std::string arg_stem( GetStem( ARGUMENT));
        const std::string normalized_arg_stem( NormalizeWord( arg_stem));
        match_type = e_SuffixCaseSpace;
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          // if any case sensitive matches have been found, only check for case-sensitive matching strings
          if( match_type == e_Suffix)
          {
            if( arg_stem == GetStem( *itr))
            {
              closest_matches.PushBack( *itr);
            }
          }
          else
          {
            // get the expected stem; possibly with case removed
            const std::string expected_stem( GetStem( *itr));

            // standardize case
            const std::string normalized_expected_stem( NormalizeWord( expected_stem));

            if( normalized_arg_stem == normalized_expected_stem)
            {
              // case insensitive match, check for case sensitive match
              if( arg_stem == expected_stem)
              {
                // remove inferior case insensitive matches
                match_type = e_Suffix;
                closest_matches.Reset();
              }
              closest_matches.PushBack( *itr);
            }
          }
        }

        // return matches that differed in suffix and possibly case/space
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      const storage::Vector< std::string> arg_split_into_words( SplitIntoWords( ARGUMENT));
      const storage::Vector< std::string> arg_split_into_stems( GetStems( arg_split_into_words));

      // cache each expected option split into words and stems for later use
      storage::Vector< storage::Vector< std::string> > expected_split_into_words;
      storage::Vector< storage::Vector< std::string> > expected_split_into_stems;
      expected_split_into_words.AllocateMemory( expected_size);
      expected_split_into_stems.AllocateMemory( expected_size);

      // stems case
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          expected_split_into_words.PushBack( SplitIntoWords( *itr));
          expected_split_into_stems.PushBack( GetStems( expected_split_into_words.LastElement()));
          if( arg_split_into_stems == expected_split_into_stems.LastElement())
          {
            closest_matches.PushBack( *itr);
          }
        }

        // check for same stems case
        if( !closest_matches.IsEmpty())
        {
          // words had the same stem set, so return the given words
          match_type = e_Stems;
          return match_type_and_closest_matches;
        }
      }

      // check for reordered words or stems
      {
        match_type = e_ReorderedStems;
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr_split_words( expected_split_into_words.Begin()), itr_split_stems( expected_split_into_stems.Begin());
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_split_words, ++itr_split_stems
        )
        {
          if( arg_split_into_words.GetSize() != itr_split_words->GetSize())
          {
            // skip items with wrong #s of words
            continue;
          }
          if( graph::CSISubstructure::IsContainedIn( arg_split_into_words, *itr_split_words))
          {
            if( match_type == e_ReorderedStems)
            {
              // remove the inferior reordered stems matches
              closest_matches.Reset();
              match_type = e_ReorderedWords;
            }
            // found this set of words as a subset of one of the options
            closest_matches.PushBack( *itr);
          }
          else if( match_type == e_ReorderedStems)
          {
            // no reordered words have been found, look for reordered stems instead
            if( graph::CSISubstructure::IsContainedIn( arg_split_into_stems, *itr_split_stems))
            {
              closest_matches.PushBack( *itr);
            }
          }
        }

        // check whether any reordered words or stems were found
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for explicit cases, where the ARGUMENT is a fully written out version of the option
      {
        size_t best_abbreviation_gaps( FullSplitIntoWords( ARGUMENT).GetSize());
        match_type = e_Explicit;
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          // determine the number of gaps, considering the words to be sequence alignments
          const size_t gaps_size( MatchingSequenceWithGaps( *itr, ARGUMENT));
          if( gaps_size <= best_abbreviation_gaps)
          {
            if( gaps_size < best_abbreviation_gaps)
            {
              // discard inferior matches
              best_abbreviation_gaps = gaps_size;
              closest_matches.Reset();
            }
            closest_matches.PushBack( *itr);
          }
        }
        // if case/space sensitive abbreviations were found, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for case/space-sensitive abbreviations; prefer case sensitive abbreviations
      {
        // threshold; if more than this # of non-consecutive gaps are found, then do not consider it a valid abbreviation
        size_t best_abbreviation_gaps( FullSplitIntoWords( ARGUMENT).GetSize());
        match_type = e_WeakAbbreviation;
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          // determine the number of gaps, considering the words to be sequence alignments
          const size_t gaps_size_strong( MatchingSequenceWithGaps( ARGUMENT, *itr));
          if( util::IsDefined( gaps_size_strong) && match_type == e_WeakAbbreviation)
          {
            // case sensitive abbreviation matched, prior matches were case insensitive, which is inferior
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size_strong;

            match_type = e_StrongAbbreviation;
            closest_matches.PushBack( *itr);
          }
          else if( gaps_size_strong <= best_abbreviation_gaps)
          {
            if( gaps_size_strong < best_abbreviation_gaps)
            {
              // fewer gaps than seen before, remove existing inferior matches
              closest_matches.Reset();

              // update best abbreviation gaps and add *itr as the first match
              best_abbreviation_gaps = gaps_size_strong;
            }

            // add another match
            closest_matches.PushBack( *itr);
          }

          // skip weak abbreviations, if strong abbreviations have already been found
          if( match_type == e_StrongAbbreviation)
          {
            continue;
          }

          // test gaps size, ignoring case
          const size_t gaps_size_weak( MatchingSequenceWithGaps( normalized_arg, *itr_norm));

          // ignore if gaps size is unreasonable
          if( gaps_size_weak < best_abbreviation_gaps)
          {
            // fewer gaps than seen before, remove existing inferior matches
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size_weak;
            closest_matches.PushBack( *itr);
          }
          else if( gaps_size_weak == best_abbreviation_gaps)
          {
            // add another match
            closest_matches.PushBack( *itr);
          }
        }
        // if case/space sensitive abbreviations were found, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for matching abbreviations, regardless of case
      {
        size_t best_abbreviation_gaps( util::GetUndefined< size_t>());
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          const size_t gaps_size( MatchingSequenceWithGaps( normalized_arg, *itr_norm));
          if( gaps_size < best_abbreviation_gaps)
          {
            // fewer gaps than seen before, remove existing inferior matches
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size;
            closest_matches.PushBack( *itr);
          }
          else if( util::IsDefined( gaps_size) && gaps_size == best_abbreviation_gaps)
          {
            // add another match
            closest_matches.PushBack( *itr);
          }
        }
        if( !closest_matches.IsEmpty())
        {
          match_type = e_WeakAbbreviation;
          return match_type_and_closest_matches;
        }
      }

      // check for maximum overlapping words and stems
      // If there are more overlapping stems than words, return overlapping stems, otherwise,
      // return the overlapping words
      {
        // track the max overlap of words or stems; start at 1 so to avoid finding all cases where nothing
        // matches initially
        size_t max_overlap( 1);

        // track the type (words or stems) with the greatest overlap; in case of a tie, prefer words since
        // they consitute a better matching criteria
        match_type = e_SomeStemsMatch;

        // iterate over words and split words
        storage::Vector< storage::Vector< std::string> >::iterator
          itr_option_words( expected_split_into_words.Begin()), itr_option_stems( expected_split_into_stems.Begin());
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_option_words, ++itr_option_stems
        )
        {
          // compute # of equal words
          const size_t word_overlap( graph::CSISubstructure::GetOverlap( arg_split_into_words, *itr_option_words));

          // if there were at least as many overlapping words as were seen before
          if( word_overlap >= max_overlap)
          {
            // # of overlapping words > than previously seen # of overlapping words -> Remove existing (inferior) matches
            if( word_overlap > max_overlap || match_type == e_SomeStemsMatch)
            {
              // type == e_SomeStemsMatch catches when there are as many overlapping words as the best # of overlapping stems,
              // Because overlapping words are better than overlapping stems, remove the existing inferior matches
              closest_matches.Reset();
              max_overlap = word_overlap;
              match_type = e_SomeWordsMatch;
            }

            // add the new match
            closest_matches.PushBack( *itr);
          }

          // now check for overlap of stems
          const size_t stem_overlap( graph::CSISubstructure::GetOverlap( arg_split_into_stems, *itr_option_stems));
          if( stem_overlap >= max_overlap)
          {
            if( stem_overlap > max_overlap)
            {
              // more matching stems than ever seen, remove all existing inferior matches
              closest_matches.Reset();
              max_overlap = stem_overlap;
              match_type = e_SomeStemsMatch;
            }
            if( match_type == e_SomeStemsMatch)
            {
              // add additional stem-based match, unless equivalent overlap has already been found with word-based
              // overlap
              closest_matches.PushBack( *itr);
            }
          }
        }
        // if there were any matches based on common stems or words, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // no good matching mechanism, just return all the words
      match_type = s_NumberTypes;
      closest_matches = EXPECTED;
      return match_type_and_closest_matches;
    }

    //! @brief Guess what the user intended and write those guesses to STREAM
    //! @param ARGUMENT the argument that was passed into the command line
    //! @param EXPECTED known, valid arguments
    //! @param STREAM stream to write out the guesses to
    //! @param ARGUMENT_TYPE_NAME conceptual type of the argument
    void Guesser::WriteGuesses
    (
      const std::string &ARGUMENT,
      const storage::Vector< std::string> &EXPECTED,
      std::ostream &STREAM,
      const std::string &ARGUMENT_TYPE_NAME
    ) const
    {
      const storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> >
        guess_type_guesses( Guess( ARGUMENT, EXPECTED));

      // write out the argument type name, if one was given, or parameter otherwise
      if( ARGUMENT_TYPE_NAME.empty())
      {
        STREAM << "Given parameter \"";
      }
      else
      {
        STREAM << "Given " << ARGUMENT_TYPE_NAME << " \"";
      }
      STREAM << ARGUMENT << "\" is unknown";

      // if no remotely good matches were found, print the full list
      if( guess_type_guesses.First() == Guesser::s_NumberTypes)
      {
        STREAM << ", here are the allowed values {";
        if( guess_type_guesses.Second().GetSize())
        {
          std::copy( EXPECTED.Begin(), EXPECTED.End() - 1, std::ostream_iterator< std::string>( STREAM, ", "));
          STREAM << EXPECTED.LastElement();
        }
        STREAM << "}";
      }
      else if( guess_type_guesses.Second().GetSize() == size_t( 1))
      {
        // write out the matches
        STREAM << ", did you mean \"" << guess_type_guesses.Second()( 0) << "\"";
      }
      else
      {
        STREAM << ", perhaps you meant one of these: { ";
        std::copy
        (
          guess_type_guesses.Second().Begin(),
          guess_type_guesses.Second().End() - 1,
          std::ostream_iterator< std::string>( STREAM, ", ")
        );
        STREAM << guess_type_guesses.Second().LastElement();
        STREAM << "}";
      }
      STREAM << '\n';
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Registers the affixes for the prefix and suffix list
    //! @param AFFIX to be normalized
    //! @param NORMALIZED affix to create stem word
    //! @return a reference to the Guesser modified
    Guesser &Guesser::RegisterSuffix( const std::string &AFFIX, const std::string &NORMALIZED)
    {
      if( AFFIX.size() < NORMALIZED.size())
      {
        return RegisterSuffix( NORMALIZED, AFFIX);
      }

      const std::string normalized_affix( NormalizeWord( AFFIX));
      const std::string normalized_word( NormalizeWord( NORMALIZED));
      if( m_Suffixes.GetSize() <= AFFIX.size())
      {
        m_Suffixes.Resize( AFFIX.size() + 1);
      }

      // add the normalized word in at the normal place
      m_Suffixes( normalized_affix.size())[ normalized_affix] = normalized_word;

      // make the suffix map to itself to prevent non-specific matching in non-trivial cases
      if( !normalized_word.empty() && normalized_affix != normalized_word)
      {
        m_Suffixes( normalized_word.size())[ normalized_word] = normalized_word;
      }

      return *this;
    }

    //! @brief Protect a suffix from being changed if it is the longest matching suffix
    //! @param SUFFIX to be normalized
    //! @return a reference to the Guesser modified
    Guesser &Guesser::ProtectSuffix( const std::string &SUFFIX)
    {
      return RegisterSuffix( SUFFIX, SUFFIX);
    }

    //! @brief Registers the aliases to be used
    //! @param ORIGINAL word
    //! @param ALIAS for original word
    //! @return a reference to the Guesser modified
    Guesser &Guesser::RegisterAlias( const std::string &ORIGINAL, const std::string &ALIAS)
    {
      m_Aliases[ NormalizeWord( ORIGINAL)] = ALIAS;
      return *this;
    }

    //! @brief Applies suffix rules to get down to the word stem
    //! @param WORD to be stemmed
    //! @return word with the suffix standardized
    std::string Guesser::GetStem( const std::string &WORD) const
    {
      // determine the word size; this may change so do not make it constant
      size_t word_size( WORD.size());

      // skip words with size <= 2, they cannot have a prefix or affix
      if( word_size <= 2)
      {
        return WORD;
      }

      // do not try to match affixes that would leave only a single letter on the word, e.g., if the affix is our,
      // do not match four.
      const size_t max_affix_size( word_size - 2);

      if( !m_Suffixes.IsEmpty())
      {
        // start by fixing the suffix, if possible
        // start with the largest suffixes
        for( size_t suffix_size( std::min( max_affix_size, m_Suffixes.GetSize() - 1)); suffix_size > 0; --suffix_size)
        {
          // get the suffix of the given size from the given word
          const std::string actual_suffix( WORD.substr( word_size - suffix_size));

          // check for it in the appropriate map
          storage::Map< std::string, std::string>::const_iterator itr( m_Suffixes( suffix_size).Find( actual_suffix));
          if( itr != m_Suffixes( suffix_size).End())
          {
            // replace the given suffix
            if( itr->first != itr->second)
            {
              return WORD.substr( 0, word_size - suffix_size) + itr->second;
            }
            else
            {
              return WORD;
            }
          }
        }
      }
      return WORD;
    }

    //! @brief GetStem for an entire vector
    //! @param WORDS words to be stemmed
    //! @return the word stems
    storage::Vector< std::string> Guesser::GetStems( const storage::Vector< std::string> &WORDS) const
    {
      storage::Vector< std::string> stems;
      stems.AllocateMemory( WORDS.GetSize());
      // call standardize suffix on each element in the vector
      for
      (
        storage::Vector< std::string>::const_iterator itr( WORDS.Begin()), itr_end( WORDS.End());
        itr != itr_end;
        ++itr
      )
      {
        stems.PushBack( Guesser::GetStem( *itr));
      }
      return stems;
    }

    // @brief Splits input string into words based on white space characters, capitalization, and underscores
    // @param STRING to be split into words
    // @return Vector of split words
    storage::Vector< std::string> Guesser::SplitIntoWords( const std::string &STRING) const
    {
      // test whether the trimmed string appears to be an object data label, if so, just take the first word
      size_t label_delimiter_pos( STRING.find_first_of( "(="));
      if( label_delimiter_pos && label_delimiter_pos != std::string::npos)
      {
        // if an object data label delimiter was found, just split the first word into words
        return SplitIntoWords( STRING.substr( 0, label_delimiter_pos));
      }

      // get the trimmed string
      const std::string trimmed( util::TrimString( STRING));

      storage::Vector< std::string> split;

      // check whether spaces or underscores are present
      if( trimmed.find_first_of( "_ \n\t\r") != std::string::npos)
      {
        // just split the string as normal
        split = util::SplitString( STRING, "_ \n\t\r");
      }
      else
      {
        // split based on capitalization
        // a word is defined by the regex: [A-Z0-9]+[a-z0-9]*
        // also, a word may not be composed entirely of numbers
        size_t pos( 0), size( trimmed.size());
        while( pos < size)
        {
          // skip non-alphanumeric characters
          if( !isalnum( trimmed[ pos]))
          {
            ++pos;
            continue;
          }

          // add caps and digits
          const size_t original_pos( pos);
          ++pos;
          while( pos < size && ( isupper( trimmed[ pos]) || isdigit( trimmed[ pos])))
          {
            ++pos;
          }

          size_t number_letters( pos - original_pos);
          // if there was only one initial letter, then we have a normal word
          // 2+ initial letters means an abbreviation sequence, followed by the end, or a normal word
          // If it was an abbreviation sequence followed by a word, then we need to backup one
          if( number_letters > size_t( 1) && pos < size)
          {
            --pos;
          }
          else
          {
            // typical case, the initial letter sequence is complete
            while( pos < size && ( islower( trimmed[ pos]) || isdigit( trimmed[ pos])))
            {
              ++pos;
            }

            // if the last letter was a digit, and we are not at the end, then that digit belongs with the next word
            if( pos < size && isdigit( trimmed[ pos - 1]))
            {
              --pos;
            }
          }

          // safety catch, if there are no letters, increment pos
          if( pos == original_pos)
          {
            ++pos;
          }
          split.PushBack( trimmed.substr( original_pos, pos - original_pos));
        }
      }
      std::for_each( split.Begin(), split.End(), ToLower);
      return split;
    }

    // @brief Normalizes a string to remove capitalization, spacing, and underscores
    // @param WORD to be normalized
    // @return normalized word
    std::string Guesser::NormalizeWord( const std::string &WORD) const
    {
      // form a new string without the spaces, _, and with all lower case letters
      std::string normalized_word;
      normalized_word.reserve( WORD.size());
      for( size_t i( 0), size( WORD.size()); i < size; ++i)
      {
        if( !isspace( WORD[ i]) && WORD[ i] != '_')
        {
          normalized_word += tolower( WORD[ i]);
        }
      }
      return normalized_word;
    }

    //! @brief Splits input string into words based on white space characters, capitalization, and underscores
    //! @param STRING to be split into words, letters of acronyms as separate entities
    //! @return Vector of split words, letters of acronyms as separate entities
    storage::Vector< std::string> Guesser::FullSplitIntoWords( const std::string &STRING)
    {
      // test whether the trimmed string appears to be an object data label, if so, just take the first word
      size_t label_delimiter_pos( STRING.find_first_of( "(="));
      if( label_delimiter_pos != std::string::npos)
      {
        // if an object data label delimiter was found, just split the first word into words
        return FullSplitIntoWords( STRING.substr( 0, label_delimiter_pos));
      }

      // first, get the trimmed string
      const std::string trimmed( util::TrimString( STRING));

      storage::Vector< std::string> split;

      // check whether spaces or underscores are present
      if( trimmed.find_first_of( "_ \n\t\r") != std::string::npos)
      {
        // just split the string as normal
        const storage::Vector< std::string> partially_split( util::SplitString( STRING, "_ \n\t\r"));

        for( size_t j( 0), size( partially_split.GetSize()); j < size; ++j)
        {
          if( !partially_split( j).empty())
          {
            // split this word into multiple words, if necessary
            split.Append( FullSplitIntoWords( partially_split( j)));
          }
        }
      }
      else
      {
        // split based on capitalization
        // a word is defined by the regex: [A-Z0-9]+[a-z0-9]*
        // also, a word may not be composed entirely of numbers
        size_t pos( 0), size( trimmed.size());
        while( pos < size)
        {
          // add caps and digits
          const size_t original_pos( pos);

          // the first letter is automatically a new word
          ++pos;

          // typical case, the initial letter sequence is complete
          while( pos < size && islower( trimmed[ pos]))
          {
            ++pos;
          }

          std::string new_string( trimmed.substr( original_pos, pos - original_pos));
          ToLower( new_string);

          // no further splitting is possible on single characters
          if( !IsObviousAcroynm( new_string))
          {
            split.PushBack( new_string);
          }
          else
          {
            // fully split the word
            std::string letter( size_t( 1), ' ');
            // words where all letters are either consonants or vowels are generally acronyms (aa, ss, pdb, sdf)
            for( size_t i( 0), new_size( new_string.size()); i < new_size; ++i)
            {
              letter[ 0] = new_string[ i];
              split.PushBack( letter);
            }
          }
        }
      }
      return split;
    }

    // @brief Just changes the capitalization of a word to lower
    // @param WORD to be lower-cased
    void Guesser::ToLower( std::string &WORD)
    {
      // form a new string without the spaces, _, and with all lower case letters
      for( size_t i( 0), size( WORD.size()); i < size; ++i)
      {
        WORD[ i] = tolower( WORD[ i]);
      }
    }

    //! @brief Get the number of characters that match at the start of the string
    //! @param A, B the two strings of interest
    //! @return the number of letters at the start of the strings that match
    size_t Guesser::GetMatchingPrefixLength( const std::string &A, const std::string &B)
    {
      if( A.size() <= B.size())
      {
        return std::mismatch( A.begin(), A.end(), B.begin()).first - A.begin();
      }
      return std::mismatch( B.begin(), B.end(), A.begin()).first - B.begin();
    }

    //! @brief Test whether one string is the abbreviation of another
    //! @param A, B the two strings of interest
    //! @return true if A is an abbreviation prefix of B or B is an abbreviation of A
    bool Guesser::IsAbbreviation( const std::string &A, const std::string &B)
    {
      if( A.size() <= B.size())
      {
        return std::mismatch( A.begin(), A.end(), B.begin()).first == A.end();
      }
      return std::mismatch( B.begin(), B.end(), A.begin()).first == B.end();
    }

    //! @brief Determine whether a given string is obviously an acronym
    //! @param A the string of interest
    //! @return true if A is > 1 letter, and contains all vowels, or all consonants
    //! @note if y's are present, they count as a vowel only if all the rest of the word is consonants
    //! @note if the A is just a sequences of Y's, it counts as an acronym
    bool Guesser::IsObviousAcroynm( const std::string &A)
    {
      static const std::string s_vowels( "aeiou");

      // test for 1 letter words
      if( A.size() <= size_t( 1))
      {
        return false;
      }

      size_t vowels_count( 0);
      // match the first letter independent of case
      if( s_vowels.find( tolower( A[ 0])) != std::string::npos)
      {
        ++vowels_count;
      }

      // count vowels for the rest of the string
      for( std::string::const_iterator itr( A.begin() + 1), itr_end( A.end()); itr != itr_end; ++itr)
      {
        if( s_vowels.find( *itr) != std::string::npos)
        {
          ++vowels_count;
        }
      }

      // if there appeared to be no vowels, look for a y
      // ignore the first character, since Y is always a consonant when it is the first letter
      if( !vowels_count && A.find( 'y', 1) != std::string::npos)
      {
        // make sure that the word was not all y's
        if( tolower( A[ 0]) != 'y' || A.find_first_not_of( 'y', 1) != std::string::npos)
        {
          ++vowels_count;
        }
      }

      // return true if there were no vowels or no consonants
      return !vowels_count || vowels_count == A.size();
    }

    //! @brief Test whether one string is an exact abbreviation / substring of another, with arbitrary number of gaps
    //! @param SUBSTR the posited substring
    //! @param STRING the string to test whether SUBSTR is a subsequence of
    //! @param START position to start at
    //! @return # of non-consecutive gaps, or undefined if there is no matching sub-sequence
    size_t Guesser::MatchingSequenceWithGaps
    (
      const std::string &WITH_GAPS,
      const std::string &STRING,
      const size_t &START
    )
    {
      // check that the alignment could exist
      if( STRING.size() < START + WITH_GAPS.size() || STRING.empty())
      {
        return util::GetUndefined< size_t>();
      }
      else if( WITH_GAPS.empty())
      {
        return 0;
      }

      // find the first matching letter
      const size_t first_match( STRING.find( WITH_GAPS[ 0], START));

      // no match: return
      if( first_match == std::string::npos || STRING.size() < first_match + WITH_GAPS.size())
      {
        return util::GetUndefined< size_t>();
      }

      bool last_was_gap( false);
      size_t number_non_consecutive_gaps( 0);

      // track # of available gaps left in string
      size_t remaining_gaps( STRING.size() - first_match - WITH_GAPS.size());

      // get iterators onto both strings
      std::string::const_iterator
        itr_gap( WITH_GAPS.begin()), itr_gap_end( WITH_GAPS.end()), itr_full( STRING.begin() + first_match);

      while( itr_gap != itr_gap_end)
      {
        if( *itr_gap != *itr_full)
        {
          // test whether any gaps remain
          if( !remaining_gaps)
          {
            return util::GetUndefined< size_t>();
          }
          // do not allow gaps to traverse upper letters except at the beginning and end of the string
          else if( itr_gap != WITH_GAPS.begin() && isupper( *itr_full))
          {
            // look for a match further on in the string
            return MatchingSequenceWithGaps( WITH_GAPS, STRING, first_match + 1);
          }
          --remaining_gaps;
          if( !last_was_gap)
          {
            ++number_non_consecutive_gaps;
          }
          last_was_gap = true;
        }
        else
        {
          ++itr_gap;
          last_was_gap = false;
        }
        ++itr_full;
      }

      if( !number_non_consecutive_gaps)
      {
        return 0;
      }
      return std::min( number_non_consecutive_gaps, MatchingSequenceWithGaps( WITH_GAPS, STRING, first_match + 1));
    }

  } // namespace command
} // namespace bcl
