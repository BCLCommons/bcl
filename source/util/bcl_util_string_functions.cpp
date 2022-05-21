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
#include "util/bcl_util_string_functions.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief test whether a given string begins with another string
    //! @param STRING the full string
    //! @param TEST_PREFIX the prefix
    //! @return true if STRING begins with TEST_PREFIX
    bool StartsWith( const std::string &STRING, const std::string &TEST_PREFIX)
    {
      return STRING.size() >= TEST_PREFIX.size()
             && std::equal( STRING.begin(), STRING.begin() + TEST_PREFIX.size(), TEST_PREFIX.begin());
    }

    //! @brief test whether a given string ends with another string
    //! @param STRING the full string
    //! @param TEST_SUFFIX the prefix
    //! @return true if STRING ends with TEST_SUFFIX
    bool EndsWith( const std::string &STRING, const std::string &TEST_SUFFIX)
    {
      return STRING.size() >= TEST_SUFFIX.size()
             && std::equal( STRING.end() - TEST_SUFFIX.size(), STRING.end(), TEST_SUFFIX.begin());
    }

    //! @brief Convert characters ]A-Z] in a string converted to lower case
    //! @param STRING the full string
    //! @return the string with all upper-case characters converted to lower case
    std::string ToLower( const std::string &STRING)
    {
      std::string lowered( STRING);
      for( std::string::iterator itr( lowered.begin()), itr_end( lowered.end()); itr != itr_end; ++itr)
      {
        *itr = tolower( *itr);
      }
      return lowered;
    }

    //! @brief Convert characters [a-z] in a string converted to upper case
    //! @param STRING the full string
    //! @return the string with all lower-case characters converted to upper case
    std::string ToUpper( const std::string &STRING)
    {
      std::string raised( STRING);
      for( std::string::iterator itr( raised.begin()), itr_end( raised.end()); itr != itr_end; ++itr)
      {
        *itr = toupper( *itr);
      }
      return raised;
    }

    //! @brief Remove a specified set of starting and trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning and end of the string
    //! @return the string stripped of starting and ending characters found in CHARS_TO_STRIP
    std::string Strip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for first character not to be stripped
      const std::string::size_type pos1( STRING.find_first_not_of( CHARS_TO_STRIP));

      if( pos1 == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }

      //searches for the last character that will not be stripped
      const std::string::size_type pos2( STRING.find_last_not_of( CHARS_TO_STRIP));

      //returns substring from pos1 of length pos2 - pos1 + 1
      return STRING.substr( pos1, pos2 - pos1 + 1);
    }

    //! @brief Remove a specified set of starting characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning of the string
    //! @return the string stripped of starting characters found in CHARS_TO_STRIP
    std::string LStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for first character not to be stripped
      const std::string::size_type pos( STRING.find_first_not_of( CHARS_TO_STRIP));

      if( pos == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }
      return STRING.substr( pos);
    }

    //! @brief Remove a specified set of trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the end of the string
    //! @return the string stripped of ending characters found in CHARS_TO_STRIP
    std::string RStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for the last character that will not be stripped
      const std::string::size_type pos( STRING.find_last_not_of( CHARS_TO_STRIP));

      if( pos == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }

      //returns substring from first character to the last non-stripped character
      return STRING.substr( 0, pos + 1);
    }

    //! @brief repeat a string a particular number of times
    //! @param STRING the string to repeat
    //! @param REPETITIONS the number of times to repeat the string
    //! @return the string, repeated REPETITIONS times
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    std::string Repeat( const std::string &STRING, const size_t &REPETITIONS)
    {
      size_t repetitions_binary( REPETITIONS);
      std::string repeated( REPETITIONS & 1 ? STRING : std::string());
      repeated.reserve( STRING.size() * REPETITIONS);

      // at each iteration in the loop below, string_doubler doubles in size
      std::string string_doubler( STRING);
      string_doubler.reserve( STRING.size() * REPETITIONS);
      while( repetitions_binary >>= 1)
      {
        string_doubler += string_doubler;

        // if the remaining # of repeats is odd, add the doubled string to the repeated string
        if( repetitions_binary & 1)
        {
          repeated += string_doubler;
        }
      }
      return repeated;
    }

    //! @brief join a vector of strings with a different internal string
    //! @param JOINER the string to insert between consecutive elements of the vector of strings
    //! @param STRINGS array of strings to join with the JOINER
    //! @return STRINGS(0){JOINER}{STRINGS(1)}{JOINER}...{STRINGS(N)}
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    std::string Join( const std::string &JOINER, const storage::Vector< std::string> &STRINGS)
    {
      std::string joined;
      if( STRINGS.IsEmpty())
      {
        return joined;
      }

      // initialize with the 1st element
      storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
      joined = *itr;

      // join the remaining elements
      for( ++itr; itr != itr_end; ++itr)
      {
        joined += JOINER;
        joined += *itr;
      }
      return joined;
    }

    //! returns a vector of std::string from an array of char *
    storage::Vector< std::string> StringListFromCharacterArray( const int NUMBER_ARGUMENTS, const char *ARGUMENTS[])
    {
      //list of strings
      storage::Vector< std::string> string_list;
      string_list.AllocateMemory( NUMBER_ARGUMENTS);

      // copy arguments
      for( int i( 0); i < NUMBER_ARGUMENTS; ++i)
      {
        string_list.PushBack( std::string( ARGUMENTS[ i]));
      }

      //end
      return string_list;
    }

    //! reads all strings from an istream and stores them consecutively in a storage vector of strings
    storage::Vector< std::string> StringListFromIStream( std::istream &ISTREAM)
    {
      //initialize string_list
      storage::Vector< std::string> string_list;

      //read strings from the ISTREAM till the end
      while( !ISTREAM.eof())
      {
        // get the next character
        char c( ISTREAM.get());
        if( !isspace( c) && !ISTREAM.eof())
        {
          // add a new string at the end of the vector and get a reference on it
          string_list.PushBack();
          std::string &current_argument( string_list.LastElement());

          // keep track of whether we are in quotes, escapes, etc. quotes and escapes
          bool in_quotes( false);
          bool in_escape( false);
          char quotes_char( '"'); // if in_quotes is true, this is the type of quotes that we are presently in

          do
          {
            if( in_escape)
            {
              in_escape = false;
              current_argument += c;
            }
            else if( c == '\\')
            {
              in_escape = true;
            }
            else if( in_quotes)
            {
              if( c == quotes_char)
              {
                in_quotes = false;
              }
              else
              {
                current_argument += c;
              }
            }
            else if( c == '\'' || c == '"')
            {
              in_quotes = true;
              quotes_char = c;
            }
            else
            {
              current_argument += c;
            }

            c = char( ISTREAM.get());
          } while( !ISTREAM.eof() && ( !isspace( c) || in_quotes));
        }
      }

      //end
      return string_list;
    }

    //! @brief SplittedStringLineListFromIStream reads all lines from an istream, splits them, and then stores them
    //! @param ISTREAM the stream from which the lines will be read
    //! @param SPLITTER the character to be used to split the line strings
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    storage::Vector< storage::Vector< std::string> > SplittedStringLineListFromIStream
    (
      std::istream &ISTREAM, const std::string &SPLITTER
    )
    {
      // create const storage::Vector "lines" and initialize with all the unsplitted lines in "ISTREAM"
      const storage::Vector< std::string> lines( StringLineListFromIStream( ISTREAM));

      // create storage::Vector of storage::Vector of strings "splitted_lines" which will hold all the splitted lines
      // of "lines"
      storage::Vector< storage::Vector< std::string> > splitted_lines;

      // iterate through "lines" in order to split the lines into the individual strings that make up each line
      for
      (
        storage::Vector< std::string>::const_iterator iter( lines.Begin()), iter_end( lines.End());
        iter != iter_end;
        ++iter
      )
      {
        // add the strings of the line currently denoted by "itr" to "splitted_lines"
        if( !iter->empty())
        {
          splitted_lines.PushBack( SplitString( *iter, SPLITTER));
        }
      }

      // return the lines of "ISTREAM" which have been split into the strings that make them up
      return splitted_lines;
    }

    //! returns a vector of substrings of STRING, that don't contain SPLITTER
    storage::Vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER)
    {
      // initialize results
      storage::Vector< std::string> result;
      // Skip delimiters at beginning.
      std::string::size_type last_pos( STRING.find_first_not_of( SPLITTER, 0));
      // Find first "non-delimiter".
      std::string::size_type pos( STRING.find_first_of( SPLITTER, last_pos));

      // as long as something was found
      while( std::string::npos != pos || std::string::npos != last_pos)
      {
        // Found a token, add it to the vector.
        result.PushBack( STRING.substr( last_pos, pos - last_pos));
        // Skip delimiters.  Note the "not_of"
        last_pos = STRING.find_first_not_of( SPLITTER, pos);
        // Find next "non-delimiter"
        pos = STRING.find_first_of( SPLITTER, last_pos);
      }

      // end
      return result;
    }

    //! returns a string which has all REPLACE replaced with REPLACEMENT in ORIGINAL
    std::string ReplaceString( const std::string &ORIGINAL, const std::string &REPLACE, const std::string &REPLACEMENT)
    {
      // split class identifier at REPLACE
      const storage::Vector< std::string> split_original( SplitString( ORIGINAL, REPLACE));

      // is there anything to replace
      if( split_original.GetSize() <= 1)
      {
        return ORIGINAL;
      }

      // initialize filename with the first part without REPLACEMENT in front
      std::string result( split_original.FirstElement());
      // loop over remaining parts and add them to filename after REPLACEMENT
      for
      (
        storage::Vector< std::string>::const_iterator name_itr( split_original.Begin() + 1),
          name_itr_end( split_original.End());
        name_itr != name_itr_end;
        ++name_itr
      )
      {
        result += REPLACEMENT + *name_itr;
      }

      return result;
    }

    //! trims spaces from beginning and end of a copied string
    std::string TrimString( const std::string &STRING)
    {
      return Strip( STRING, " \n\t\r\0");
    }

    //! test whether string is numerical (double, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    bool IsNumerical( const std::string &STRING)
    {
      return LengthOfFloatingPointType( STRING) == STRING.size();
    }

    //! @brief Find the number of characters that belong to a character string that can be converted to a numerical type
    //! @brief (e.g. float, double, int, size_t, etc.), including whitespace
    //! @param STRING the string to search
    //! @param START is where number (or whitespace) should begin in the string
    //! @return the number of characters in the number starting at STRING[ START]
    //! @note returns 0 if STRING[ START]...STRING[START+X]
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators,
    size_t LengthOfFloatingPointType( const std::string &STRING, const size_t &START)
    {
      size_t i( START);

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a double,
      // none of the loops would continue if they reached it.
      const char *my_string = STRING.c_str();

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces
      {
        ++i;
      }

      if( my_string[ i] == '-' || my_string[ i] == '+') // Optional sign after spaces
      {
        ++i;
      }

      if( isdigit( my_string[ i]))
      {
        ++i; // MYSTRING[ i] is a digit; no need to check it again

        while( isdigit( my_string[ i])) // go through the digits
        {
          ++i;
        }

        if( my_string[ i] == '.') // allow a decimal
        {
          ++i; // step over the decimal
          while( isdigit( my_string[ i])) // and then any digits that follow
          {
            ++i;
          }
        }
      }
      else if( my_string[ i] == '.' && isdigit( my_string[ i + 1])) // number that starts with decimal, e.g. .9
      {
        ++i; // step over the decimal
        while( isdigit( my_string[ i])) // go through the numbers after the decimal point
        {
          ++i;
        }
      }
      else // something other than a number
      {
        return size_t( 0);
      }

      // lastly, account for the optional exponent
      // 'i' will move forward iff there is a sign or number after the e/E
      if( my_string[ i] == 'e' || my_string[ i] == 'E') //
      {
        if( isdigit( my_string[ i + 1])) // implied positive sign after e or E, e.g. 1e9 == 1e+9
        {
          i += 2;
        }
        else if( ( my_string[ i + 1] == '-' || my_string[ i + 1] == '+') && isdigit( my_string[ i + 2])) // [eE][+-][:digit] after the number
        { // [eE]?[+-][[:digits:]+] after the mantissa
          i += 3;
        }

        while( isdigit( my_string[ i])) // walk past the remainder of the exponent, if any
        {
          ++i;
        }
      }

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
      {
        ++i;
      }

      // i is now 1 past the last character that could be considered part of the number, so
      // returning i-START gives the number of characters in the array.
      return ( i - START);
    }

    //! Find the number of consecutive characters in a string that can be converted to a single integer type
    //! (e.g. long, int, short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long/integer/short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    size_t LengthOfIntegerType( const std::string &STRING, const size_t &START)
    {
      size_t i = START;

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a numeric type,
      // none of the loops would continue if they reached it.
      const char *const &my_string( STRING.c_str());

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Optional spaces
      {
        i++;
      }

      if( my_string[ i] == '+' || my_string[ i] == '-') // Optional sign after spaces
      {
        i++;
      }

      if( isdigit( my_string[ i]))
      {
        ++i;
        while( isdigit( my_string[ i]))
        {
          ++i;
        }
        while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
        {
          i++;
        }
        return ( i - START);
      }
      else // no digits found; go back to start
      {
        return size_t( 0);
      }
    }

    //! Find the number of consecutive characters in a string that can be converted to a single unsigned integer type
    //! (e.g. size_t or unsigned long/int/short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long/integer/short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    size_t LengthOfUnsignedIntegerType( const std::string &STRING, const size_t &START)
    {
      size_t i = START;

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a numeric type,
      // none of the loops would continue if they reached it.
      const char *const &my_string = STRING.c_str();

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Optional spaces
      {
        ++i;
      }

      if( my_string[ i] == '+') // Optional sign after spaces
      {
        ++i;
      }

      if( isdigit( my_string[ i])) // now at the first digit.
      {
        ++i;
        while( isdigit( my_string[ i]))
        {
          ++i;
        }
        while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
        {
          i++;
        }
        return ( i - START);
      }
      else // no digits found
      {
        return size_t( 0);
      }
    }

    //! searches for spaces and removes them from the string
    std::string RemoveSpacesFromString( const std::string &STRING)
    {
      std::string cleaned_string( STRING);

      // erase every space
      cleaned_string.erase
      (
        std::remove( cleaned_string.begin(), cleaned_string.end(), ' '),
        cleaned_string.end()
      );

      return cleaned_string;
    }

    //! @brief ConvertStringToBoolean takes a string and converts it to a boolean value
    //! @param STRING the string that will be converted to a boolean, should be either "true" or "false"
    //! @return boolean either true of false depending on STRING
    bool
    ConvertStringToBoolean( const std::string &STRING)
    {
      // remove tabs, spaces, new lines, etc.
      const std::string string( TrimString( STRING));

      // make sure the string is either true or false
      BCL_Assert
      (
        string == "true" || string == "false",
        "string to convert to boolean is neither \"true\" nor \"false\" but is \"" + string + "\""
      );

      // return whether the string was true
      return string == "true";
    }

    //! @brief wrap a string for a given line length into multiple lines
    //! @details like split string, only that it will also split, if the resulting string is longer that LINE_LENGTH
    //! @param STRING the string to be wrapped
    //! @param LINE_LENGTH length of lines
    //! @param WRAPPER_CHARS all characters to wrap at
    //! @return vector of lines
    storage::Vector< std::string> WrapString( const std::string &STRING, const size_t LINE_LENGTH, const std::string &WRAPPER_CHARS)
    {
      // stores the wrapped lines
      storage::Vector< std::string> wrapped_lines;

      const std::string::size_type string_length( STRING.length());
      std::string::size_type current_start( 0);
      std::string::size_type current_end( 0);
      std::string::size_type skip( 0); // if wrap is at space, remove that space

      // as long as end was not hit
      while( current_start != string_length)
      {
        // potentially go till end of line
        current_end = current_start + LINE_LENGTH;

        // current end beyond the string
        if( current_end > string_length)
        {
          current_end = string_length;
        }
        // not at the end yet, try to wrap at space
        else
        {
          // find the last space before the max line length
          const std::string::size_type current_space_pos( STRING.rfind( WRAPPER_CHARS, current_end));

          // try to wrap the string at the space position
          if( current_space_pos != std::string::npos && current_space_pos > current_start)
          {
            current_end = current_space_pos;
            skip = 1;
          }
        }

        // print the substring from start with length
        wrapped_lines.PushBack( STRING.substr( current_start, current_end - current_start));

        // update the start to the current end
        current_start = current_end + skip;
        skip = 0;
      }

      // end
      return wrapped_lines;
    }

    //! ignore lines that start with '#' and empty lines
    void ChopHeader( std::istream &ISTREAM)
    {
      std::string line;
      while( ISTREAM.peek() == '#' || ISTREAM.peek() == '\n')
      {
        std::getline( ISTREAM, line);
      }
    }

    //! returns the maximum matching string in terms of pair< size_t, size_t>( index, length) ~ substr( index, length) within first and second string
    storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
    MaximumMatchingSubstring( const std::string &STRING_A, const std::string &STRING_B, const size_t &MINIMUM_SIZE)
    {
      // search second within first string
      // start with maximum string and reduce size iteratively
      size_t size( std::min( STRING_B.size(), STRING_A.size()));
      while( size >= MINIMUM_SIZE)
      {
        // shift window through first sequence
        for( size_t i = 0; i < STRING_B.size() - size + 1; ++i)
        {
          size_t index( STRING_A.find( STRING_B.substr( i, size)));
          if( index < std::string::npos)
          {
            return storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
            (
              STRING_B.substr( i, size),
              storage::Pair< size_t, size_t>( index, size),
              storage::Pair< size_t, size_t>( i, size)
            );
          }
        }
        --size;
      }
      return storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
      (
        std::string( ""),
        storage::Pair< size_t, size_t>( GetUndefined< size_t>(), GetUndefined< size_t>()),
        storage::Pair< size_t, size_t>( GetUndefined< size_t>(), GetUndefined< size_t>())
      );
    }

    //! @brief find the matching brackets revers
    //! @param STRING the string of interest
    //! @param POSITION the position of the closing bracket
    //! @param BRACKET_OPEN the character for opening bracket
    //! @param BRACKET_CLOSE the character for closing bracket
    //! @return the Position of the matching opening bracket
    size_t
    RFindMatchingBracket
    (
      const std::string &STRING,
      const size_t POSITION,
      const char BRACKET_OPEN,
      const char BRACKET_CLOSE
    )
    {
      if( STRING[ POSITION] != BRACKET_CLOSE)
      {
        return POSITION;
      }

      //number of brackets that need to be closed and the current pos in the string
      size_t number_brackets_to_close( 1), current_pos( POSITION);
      //iterate as long as there is a a bracket to close or the rend is reached
      while( number_brackets_to_close > 0 && current_pos != std::string::npos)
      {
        //step back
        --current_pos;
        //increase or decrease number_brackets_to_close
        if( STRING[ current_pos] == BRACKET_CLOSE)
        {
          ++number_brackets_to_close;
        }
        else if( STRING[ current_pos] == BRACKET_OPEN)
        {
          --number_brackets_to_close;
        }
      }

      //check that every pair of brackets are closed
      BCL_Assert
      (
        number_brackets_to_close == 0,
        std::string( "could not find every matching pair of \'")
          + BRACKET_OPEN + "\' and \'" + BRACKET_CLOSE + "\' in \""
          + STRING + "\" starting at " + util::Format()( POSITION)
      );

      //return position pointing to the matching openeing bracket
      return current_pos;
    }

    //! @brief StringLineListFromIStream reads all lines from an istream and stores them consecutively as strings
    //! @param ISTREAM the stream from which the lines will be read
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    storage::Vector< std::string> StringLineListFromIStream( std::istream &ISTREAM)
    {
      // create storage::Vector "string_lines" to hold the consecutive lines being read in from "ISTREAM"
      storage::Vector< std::string> string_lines;

      // create std::string "current_line" to hold the current line of the accessibility file
      std::string current_line;

      // read in all the line of "ISTREAM"
      while( std::getline( ISTREAM, current_line))
      {
        // add "current_line" to "string_lines"
        string_lines.PushBack( TrimString( current_line));
      }

      // return "string_lines" which has all the lines which were in "ISTREAM"
      return string_lines;
    }

  } // namespace util
} // namespace bcl
