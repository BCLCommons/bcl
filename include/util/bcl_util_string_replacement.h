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

#ifndef BCL_UTIL_STRING_REPLACEMENT_H_
#define BCL_UTIL_STRING_REPLACEMENT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"
#include "bcl_util_wrapper_enum.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StringReplacement
    //! @brief a class used to replace arbitrary string patterns
    //!
    //! @see @link example_util_string_replacement.cpp @endlink
    //! @author mendenjl
    //! @date   04/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StringReplacement :
      public ObjectInterface
    {

    public:

      //! Context in which to look for a match
      enum MatchContext
      {
        e_Any,    //!< match regardless of what is around the string
        e_Prefix, //!< match only if the string is preceded by neither letters, numbers, or _
        e_Suffix, //!< match only if the string is succeeded by neither letters, numbers, or _
        e_Word,   //!< combination of e_Prefix and e_Suffix
        e_Exact,  //!< match only if the strings match exactly,
        s_NumberContexts //!< last element of enum class
      };

      //! @brief GetMatchContextDescription provides the name of ENUM
      //! @param ENUM - the context for which a name is desired
      static const std::string &GetMatchContextDescription( const MatchContext &ENUM);

      //! @brief enum class wrapper for match context
      typedef WrapperEnum< MatchContext, &GetMatchContextDescription, s_NumberContexts> MatchContextEnum;

    private:

    //////////
    // data //
    //////////

      MatchContextEnum m_Context;     //!< indicates what can surround a match in a test string
      std::string      m_Match;       //!< string that this replacer will match, if in m_Context
      std::string      m_Replacement; //!< string to replace m_Match with if found in m_Context of a given string

    public:

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! an empty default constructor
      StringReplacement()
      {
      }

      //! @brief construct from data members
      //! @param CONTEXT
      //! @param MATCH
      //! @param REPLACE
      StringReplacement
      (
        const MatchContext &CONTEXT,
        const std::string  &MATCH,
        const std::string  &REPLACE = std::string()
      );

      //! @brief Clone function
      //! @return pointer to new StringReplacer
      StringReplacement *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief determine if this string replacement has a defined match
      //! @return true if this string replacement is defined
      bool IsDefined() const
      {
        return !m_Match.empty();
      }

      //! @brief return the match context
      //! @return the match context
      MatchContextEnum GetMatchContext() const
      {
        return m_Context;
      }

      //! @brief return the string that this replacement will look for
      const std::string &GetMatch() const
      {
        return m_Match;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief finds all matching positions for a given search string and context
      //! @param STRING_TO_SEARCH the string we are searching inside
      //! @return position for all occurrences
      std::vector< size_t> FindAllMatches
      (
        const std::string &STRING_TO_SEARCH
      ) const;

      //! @brief finds all matching positions for a given search string and context
      //! @param STRING_TO_SEARCH the string we are searching inside
      //! @param START_POSITION the position to test for the match
      //! @return position of the next match
      size_t FindNextMatch
      (
        const std::string &STRING_TO_SEARCH,
        const size_t &START_POSITION
      ) const;

      //! @brief ReplaceEachIn makes the replacement in a given string
      //! @param STRING the string to perform the replacements in
      //! @return the number of replacements made
      //! Non-recursive, so
      //! StringReplacement x( e_Any, "ab", "abab");
      //! string y( "abcd");
      //! x.ReplaceEachIn( y); // returns 1; now y == "ababcd"
      size_t ReplaceEachIn( std::string &STRING) const;

      //! @brief ReplaceAllIn recursively makes the replacement
      //! @param STRING the string to perform the replacements in
      //! Recursive, so
      //! StringReplacement x( e_Any, "abc", "ab");
      //! string y( "abcc");
      //! ReplaceAllIn( y); // y == "ab"
      void ReplaceAllIn( std::string &STRING) const;

      //! @brief ReplaceEachWithExclusions makes the replacement in a given string, excluding things that are matched
      //!        by other StringReplacements
      //! @param STRING the string to perform the replacements in
      //! @param EXCLUSIONS a list of string replacements.  If the string at a given position matches any exclusion, the
      //!        replacement will not happen
      //! @return the number of replacements made
      size_t ReplaceEachWithExclusions
      (
        std::string &STRING,
        const storage::List< StringReplacement> &EXCLUSIONS
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief SafeSubstr Replacement function for std::string.substr()
      //! std::string.substr(position, size) fails (throws an exception) if the position parameter is at or
      //! past the end of the string, even if called with 0 size.  This function returns an empty string instead
      //! @param STRING the string to get a substring of
      //! @param POSITION the index to start the substring with
      //! @param SIZE how big to make the substring
      //! @return STRING.substr( POSITION, SIZE) if position and size were both valid.  An empty string otherwise
      static std::string SafeSubstr
      (
        const std::string &STRING,
        const size_t &POSITION,
        const size_t &SIZE = std::string::npos
      )
      {
        return POSITION >= STRING.size()
               ? std::string()
               : STRING.substr( POSITION, SIZE);
      }

      //! @brief SafeReplaceAt uses SafeSubstr to replace one string with position and size with another
      //! @param STRING the string to get a substring of
      //! @param POSITION the position where we want to begin inserting REPLACE_STRING at
      //! @param FOUND_STRING_SIZE the size of string to replace
      //! @param REPLACE_STRING the string to replace FOUND_STRING_SIZE with
      //! @return the new string, after replacement
      static std::string SafeReplaceAt
      (
        const std::string &STRING,
        const size_t &POSITION,
        const size_t &FOUND_STRING_SIZE,
        const std::string &REPLACE_STRING
      )
      {
        return SafeSubstr( STRING, 0, POSITION) + REPLACE_STRING + SafeSubstr( STRING, POSITION + FOUND_STRING_SIZE);
      }

      //! @brief IsNonVariableCharacter Helper function for determining whether we are in a given context
      //! @param CHAR the character to check
      //! @return true iff CHAR is an alphanumeric character or an _
      static bool IsNonVariableCharacter( const char &CHAR)
      {
        return !isalnum( CHAR) && CHAR != '_';
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Replace makes the replacement in a given string at a given position
      //! @param STRING the string we are interested in
      //! @param POSITION the position to test for the match.  This is moved forward to the end of the replaced string
      //!        if a replacement was performed
      void Replace( std::string &STRING, size_t &POSITION) const;

      //! @brief Matches check whether a given string matches this string replacement at a given position
      //! @param STRING the string we are interested in
      //! @param POSITION the position to test for the match
      //! @return true if STRING matches this string replacement at POSITION
      bool Matches( const std::string &STRING, const size_t &POSITION) const;

    ////////////////
    // comparison //
    ////////////////

    public:

      //! @brief operator < test StringReplacements for specificity
      //! @param A a string replacement
      //! @param B another string replacement
      //! @return true iff A should be performed before B
      //! A string replacement has priority over another if it either has a more specific context, or an
      //! equivalent context but matches a larger string
      friend bool operator <( const StringReplacement &A, const StringReplacement &B);

    }; // class StringReplacement

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_STRING_REPLACEMENT_H_

