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
#include "util/bcl_util_string_replacement.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <cstring>

namespace bcl
{
  namespace util
  {

    // instantiate s_Instance
    const SiPtr< const ObjectInterface> StringReplacement::s_Instance
    (
      GetObjectInstances().AddInstance( new StringReplacement())
    );

    //! @brief GetMatchContextDescription provides the name of ENUM
    //! @param ENUM - the context for which a name is desired
    const std::string &StringReplacement::GetMatchContextDescription( const StringReplacement::MatchContext &ENUM)
    {
      static const std::string s_descriptors[] =
      {
        "Any",
        "Prefix",
        "Suffix",
        "Word",
        "Exact",
        GetStaticClassName< MatchContext>()
      };
      return s_descriptors[ size_t( ENUM)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! Construct from data members
    StringReplacement::StringReplacement
    (
      const MatchContext &CONTEXT,
      const std::string  &MATCH,
      const std::string  &REPLACE
    ) :
      m_Context( CONTEXT),
      m_Match( MATCH),
      m_Replacement( REPLACE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StringReplacer
    StringReplacement *StringReplacement::Clone() const
    {
      return new StringReplacement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &StringReplacement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief ReplaceEachIn makes the replacement in a given string
    //! @param STRING the string to perform the replacements in
    //! @return the number of replacements made
    //! Non-recursive, so
    //! StringReplacement x( e_Any, "ab", "abab");
    //! string y( "abcd");
    //! x.ReplaceEachIn( y); // returns 1; now y == "ababcd"
    size_t StringReplacement::ReplaceEachIn( std::string &STRING) const
    {
      size_t count( 0);

      // for each letter in the (non-constant) STRING
      // string size is non-constant, so don't move the .size() out of the for loop
      for( size_t i( 0); i < STRING.size(); i++)
      {
        if( Matches( STRING, i))
        {
          ++count;
          Replace( STRING, i);
        }
      }

      return count;
    }

    //! @brief ReplaceAllIn recursively makes the replacement
    //! @param STRING the string to perform the replacements in
    //! Recursive, so
    //! StringReplacement x( e_Any, "abc", "ab");
    //! string y( "abcc");
    //! ReplaceAllIn( y); // y == "ab"
    void StringReplacement::ReplaceAllIn( std::string &STRING) const
    {
      // must ensure that replacement string does not match the match string or else this would go on forever
      BCL_Assert
      (
        FindNextMatch( m_Replacement, 0) == std::string::npos,
        "Call to ReplaceAllIn( " + STRING + ") is recursive " + GetClassIdentifier()
        + " = " + util::Format()( *this)
      );

      // keep calling the non-recursive version until no more replacements are found
      while( ReplaceEachIn( STRING));
    }

    //! @brief ReplaceEachWithExclusions makes the replacement in a given string, excluding things that are matched
    //!        by other StringReplacements
    //! @param STRING the string to perform the replacements in
    //! @param EXCLUSIONS a list of string replacements.  If the string at a given position matches any exclusion, the
    //!        replacement will not happen
    //! @return the number of replacements made
    //! Non-recursive, so
    //! StringReplacement x( e_Any, "ab", "abab");
    //! string y( "abcd");
    //! x.ReplaceEachIn( y); // returns 1; now y == "ababcd"
    size_t StringReplacement::ReplaceEachWithExclusions
    (
      std::string &STRING,
      const storage::List< StringReplacement> &EXCLUSIONS
    ) const
    {
      size_t count( 0);

      storage::List< StringReplacement>::const_iterator itr_end( EXCLUSIONS.End());

      // for each letter in the (non-constant) STRING
      // string size is non-constant, so don't move the .size() out of the for loop
      for( size_t i( 0); i < STRING.size(); i++)
      {
        if( Matches( STRING, i)) // move to the next position because this position was excluded
        {
          bool matched_exclusion( false);
          for( storage::List< StringReplacement>::const_iterator itr( EXCLUSIONS.Begin()); itr != itr_end; ++itr)
          {
            if( itr->Matches( STRING, i))
            {
              i += itr->GetMatch().size() - 1;
              matched_exclusion = true;
              break;
            }
          }
          if( !matched_exclusion)
          {
            ++count;
            Replace( STRING, i);
          }
        }
      }

      return count;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief finds all matching positions for a given search string and context
    //! @param STRING_TO_SEARCH the string we are searching inside
    //! @return position for all occurrences
    std::vector< size_t> StringReplacement::FindAllMatches
    (
      const std::string &STRING_TO_SEARCH
    ) const
    {
      std::vector< size_t> matching_positions;
      for( size_t pos( 0); pos < STRING_TO_SEARCH.size(); ++pos)
      {
        if( Matches( STRING_TO_SEARCH, pos))
        {
          matching_positions.push_back( pos);
        }
      }

      return matching_positions;
    }

    //! @brief finds all matching positions for a given search string and context
    //! @param STRING_TO_SEARCH the string we are searching inside
    //! @param START_POSITION the position to test for the match
    //! @return position of the next match
    size_t StringReplacement::FindNextMatch
    (
      const std::string &STRING_TO_SEARCH,
      const size_t &START_POSITION
    ) const
    {
      for
      (
        size_t pos( START_POSITION), last_pos( STRING_TO_SEARCH.size());
        pos < last_pos;
        pos++
      )
      {
        if( Matches( STRING_TO_SEARCH, pos))
        {
          return pos;
        }
      }

      return std::string::npos;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StringReplacement::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Context, ISTREAM);
      io::Serialize::Read( m_Match, ISTREAM);
      io::Serialize::Read( m_Replacement, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StringReplacement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Context, OSTREAM, INDENT);
      io::Serialize::Write( m_Match, OSTREAM, INDENT);
      io::Serialize::Write( m_Replacement, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Matches check whether a given string matches this string replacement at a given position
    //! @param STRING the string we are interested in
    //! @param POSITION the position to test for the match
    //! @return true if STRING matches this string replacement at POSITION
    bool StringReplacement::Matches( const std::string &STRING, const size_t &POSITION) const
    {
      if( m_Context == e_Exact)
      {
        // only match an exact context at the first position if the strings match exactly
        return ( POSITION == 0 && STRING == m_Match);
      }
      else if( POSITION + m_Match.size() <= STRING.size())
      {
        // STRING could contain the string to be matched
        if( m_Context == e_Word || m_Context == e_Prefix) // check the start context
        {
          // Check the letter at POSITION-1, it if exists, make sure it is a valid variable character
          if( POSITION > 0 && !IsNonVariableCharacter( STRING[ POSITION - 1]))
          {
            return false; // start context did not match
          }
        }

        if( m_Context == e_Word || m_Context == e_Suffix)
        {
          // Check the letter just after where the match would end, it if exists
          if
          (
            POSITION + m_Match.size() < STRING.size()
            && !IsNonVariableCharacter( STRING[ POSITION + m_Match.size()])
          )
          {
            return false; // end context did not match
          }
        }

        return !strncmp( STRING.c_str() + POSITION, m_Match.c_str(), m_Match.size());
      }

      return false; // string was not big enough
    }

    //! @brief Replace makes the replacement in a given string at a given position
    //! @param STRING the string we are interested in
    //! @param POSITION the position to test for the match.  This is moved forward to the end of the replaced string
    //!        if a replacement was performed
    void StringReplacement::Replace( std::string &STRING, size_t &POSITION) const
    {
      if( m_Replacement.size() == m_Match.size())
      {
        // same size, no need to take substrings, just overwrite what is already there
        strncpy( &STRING[ POSITION], m_Replacement.c_str(), m_Replacement.size());
      }
      else
      {
        if( POSITION == 0)
        {
          STRING = m_Replacement + SafeSubstr( STRING, POSITION + m_Match.size());
        }
        else
        {
          // different sizes, use substrings
          STRING = STRING.substr( 0, POSITION) + m_Replacement + SafeSubstr( STRING, POSITION + m_Match.size());
        }
      }

      POSITION += m_Replacement.size() - 1;
    }

  ////////////////
  // comparison //
  ////////////////

    //! @brief operator < test StringReplacements for specificity
    //! @param A a string replacement
    //! @param B another string replacement
    //! @return true iff A should be performed before B
    //! A string replacement has priority over another if it either has a more specific context, or an
    //! equivalent context but matches a larger string
    bool operator <( const StringReplacement &A, const StringReplacement &B)
    {
      return size_t( A.m_Context) > size_t( B.m_Context)
             ||
             (
               A.m_Match.size() > B.m_Match.size()
               &&
               (
                 A.m_Context == B.m_Context
                 ||
                 ( A.m_Context == StringReplacement::e_Prefix && B.m_Context == StringReplacement::e_Suffix)
               )
             );
    }

  } // namespace util
} // namespace bcl
