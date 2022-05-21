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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "util/bcl_util_string_replacement.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_string_replacement.cpp
  //!
  //! @author mendenjl
  //! @date Apr 09, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilStringReplacement :
    public ExampleInterface
  {

  public:

    ExampleUtilStringReplacement *Clone() const
    {
      return new ExampleUtilStringReplacement( *this);
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

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // StringReplacement()
      util::StringReplacement default_string_replacement;

      // StringReplacement( const MatchContext &CONTEXT, const std::string & MATCH, const std::string & REPLACE)
      util::StringReplacement replace_fast_with_faster_any( util::StringReplacement::e_Any, "fast", "faster");
      util::StringReplacement replace_fast_with_faster_exact( util::StringReplacement::e_Exact, "fast", "faster");
      util::StringReplacement replace_fast_with_faster_prefix( util::StringReplacement::e_Prefix, "fast", "faster");
      util::StringReplacement replace_fast_with_faster_suffix( util::StringReplacement::e_Suffix, "fast", "faster");
      util::StringReplacement replace_fast_with_faster_word( util::StringReplacement::e_Word, "fast", "faster");
      util::StringReplacement replace_faster_with_fastest_any( util::StringReplacement::e_Any, "faster", "fastest");
      util::StringReplacement remove_faster( util::StringReplacement::e_Any, "faster");

    /////////////////
    // data access //
    /////////////////

      // IsDefined() const
      BCL_ExampleCheck( default_string_replacement.IsDefined(), false);
      BCL_ExampleCheck( replace_fast_with_faster_suffix.IsDefined(), true);
      BCL_ExampleCheck( remove_faster.IsDefined(), true);

      // GetMatch() const
      BCL_ExampleCheck( replace_fast_with_faster_suffix.GetMatch(), "fast");
      BCL_ExampleCheck( remove_faster.GetMatch(), "faster");

      // GetMatchContext() const
      BCL_ExampleCheck( replace_fast_with_faster_suffix.GetMatchContext(), util::StringReplacement::e_Suffix);
      BCL_ExampleCheck( replace_fast_with_faster_exact.GetMatchContext(), util::StringReplacement::e_Exact);

      // GetMatchContextDescription( const MatchContext &ENUM)
      BCL_ExampleCheck
      (
        util::StringReplacement::GetMatchContextDescription( util::StringReplacement::e_Suffix),
        "Suffix"
      );

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck // 5 matches of fast anywhere
      (
        replace_fast_with_faster_any.FindAllMatches( "fast faster fastest breakfast breakfaster").size(),
        size_t( 5)
      );
      BCL_ExampleCheck // 3 words begin with fast
      (
        replace_fast_with_faster_prefix.FindAllMatches( "fast faster fastest breakfast breakfaster").size(),
        size_t( 3)
      );
      BCL_ExampleCheck // 2 words end with fast
      (
        replace_fast_with_faster_suffix.FindAllMatches( "fast faster fastest breakfast breakfaster").size(),
        size_t( 2)
      );
      BCL_ExampleCheck // the string does not exactly match fast
      (
        replace_fast_with_faster_exact.FindAllMatches( "fast faster fastest breakfast breakfaster").size(),
        size_t( 0)
      );
      BCL_ExampleCheck // there is only word "fast" in the string
      (
        replace_fast_with_faster_word.FindAllMatches( "fast faster fastest breakfast breakfaster").size(),
        size_t( 1)
      );

      // FindNextMatch( const std::string &STRING_TO_SEARCH, const size_t &START_POSITION ) const
      BCL_ExampleCheck // matches at first
      (
        replace_fast_with_faster_any.FindNextMatch( "fast faster fastest breakfast breakfaster", 0),
        size_t( 0)
      );
      BCL_ExampleCheck // matches at first
      (
        replace_fast_with_faster_any.FindNextMatch( "fast faster fastest breakfast breakfaster", 1),
        size_t( 5)
      );

      const std::string test_string( "fast faster fastest breakfast breakfaster");

      // ReplaceEachIn( std::string &STRING) const
      {
        std::string test_string_copy( test_string);
        replace_fast_with_faster_any.ReplaceEachIn( test_string_copy);
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          "faster fasterer fasterest breakfaster breakfasterer",
          "ReplaceEachIn " + test_string + " using " + util::Format()( replace_fast_with_faster_any)
        );
      }
      {
        std::string test_string_copy( test_string);
        replace_fast_with_faster_exact.ReplaceEachIn( test_string_copy);
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          test_string,
          "ReplaceEachIn " + test_string + " using " + util::Format()( replace_fast_with_faster_exact)
        );
      }
      {
        std::string test_string_copy( test_string);
        replace_fast_with_faster_prefix.ReplaceEachIn( test_string_copy);
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          "faster fasterer fasterest breakfast breakfaster",
          "ReplaceEachIn " + test_string + " using " + util::Format()( replace_fast_with_faster_prefix)
        );
      }
      {
        std::string test_string_copy( test_string);
        replace_fast_with_faster_suffix.ReplaceEachIn( test_string_copy);
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          "faster faster fastest breakfaster breakfaster",
          "ReplaceEachIn " + test_string + " using " + util::Format()( replace_fast_with_faster_suffix)
        );
      }
      {
        std::string test_string_copy( test_string);
        replace_fast_with_faster_word.ReplaceEachIn( test_string_copy);
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          "faster faster fastest breakfast breakfaster",
          "ReplaceEachIn " + test_string + " using " + util::Format()( replace_fast_with_faster_word)
        );
      }

      // ReplaceAllIn( std::string &STRING) const should call ReplaceEachIn until no more replacements are found.
      // create some more interesting replacements for this case
      util::StringReplacement remove_duplicate_fast( util::StringReplacement::e_Word, "fast fast", "fast");
      const std::string fast5_and_faster_original( "fast fast fast fast fast and faster");
      std::string fast5_and_faster_replaced( fast5_and_faster_original);
      remove_duplicate_fast.ReplaceAllIn( fast5_and_faster_replaced);
      BCL_ExampleIndirectCheck
      (
        fast5_and_faster_replaced,
        "fast and faster",
        "ReplaceAllIn " + fast5_and_faster_original + " using " + util::Format()( remove_duplicate_fast)
      );

      // ReplaceEachWithExclusions( std::string &STRING, const storage::List< StringReplacement> &EXCLUSIONS ) const
      {
        std::string test_string_copy( test_string);
        storage::List< util::StringReplacement> keep_faster( 1, remove_faster);
        replace_fast_with_faster_prefix.ReplaceEachWithExclusions
        (
          test_string_copy,
          keep_faster
        );
        BCL_ExampleIndirectCheck
        (
          test_string_copy,
          "faster faster fasterest breakfast breakfaster",
          "ReplaceEachWithExclusions, exclusion is not removing faster, from "
          + test_string + " using " + util::Format()( replace_fast_with_faster_prefix)
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      // IsNonVariableCharacter( const char &CHAR)
      BCL_ExampleCheck( util::StringReplacement::IsNonVariableCharacter( 'c'), false);
      BCL_ExampleCheck( util::StringReplacement::IsNonVariableCharacter( 'C'), false);
      BCL_ExampleCheck( util::StringReplacement::IsNonVariableCharacter( '+'), true);
      BCL_ExampleCheck( util::StringReplacement::IsNonVariableCharacter( ' '), true);
      BCL_ExampleCheck( util::StringReplacement::IsNonVariableCharacter( '_'), false);

      // SafeSubstr( const std::string &STRING, const size_t &POSITION, const size_t &SIZE)
      BCL_ExampleCheck( util::StringReplacement::SafeSubstr( "", 5, 1), "");
      BCL_ExampleCheck( util::StringReplacement::SafeSubstr( "hello", 5, 3), "");
      BCL_ExampleCheck( util::StringReplacement::SafeSubstr( "hello", 0, 8), "hello");
      BCL_ExampleCheck( util::StringReplacement::SafeSubstr( "hello", 3, 1), "l");

      // SafeReplaceAt( const std::string &STRING, const size_t &POSITION, const size_t &FOUND_STRING_SIZE, const std::string &REPLACE_STRING)
      BCL_ExampleCheck( util::StringReplacement::SafeReplaceAt( "", 5, 1, "hello "), "hello ");
      BCL_ExampleCheck( util::StringReplacement::SafeReplaceAt( "hello", 5, 3, " tokyo"), "hello tokyo");
      BCL_ExampleCheck( util::StringReplacement::SafeReplaceAt( "hello", 0, 8, " tokyo"), " tokyo");
      BCL_ExampleCheck( util::StringReplacement::SafeReplaceAt( "hello", 3, 1, " tokyo"), "hel tokyoo");

    ////////////////
    // comparison //
    ////////////////

      BCL_ExampleCheck( replace_fast_with_faster_exact < replace_fast_with_faster_any, true);
      BCL_ExampleCheck( replace_fast_with_faster_any < replace_fast_with_faster_exact, false);
      BCL_ExampleCheck( replace_fast_with_faster_exact < remove_faster, true);
      BCL_ExampleCheck( replace_fast_with_faster_suffix < replace_fast_with_faster_prefix, true);
      BCL_ExampleCheck( replace_fast_with_faster_prefix < replace_fast_with_faster_suffix, false);

      return 0;
    } //end ExampleUtilStringReplacement

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleUtilStringReplacement

  const ExampleClass::EnumType ExampleUtilStringReplacement::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilStringReplacement())
  );
} // namespace bcl

