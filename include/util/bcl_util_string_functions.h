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

#ifndef BCL_UTIL_STRING_FUNCTIONS_H_
#define BCL_UTIL_STRING_FUNCTIONS_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_class_descriptor.h"
#include "bcl_util_logger_interface.h"
#include "bcl_util_message.h"
#include "bcl_util_string_numeric_conversion.h"
#include "bcl_util_undefined.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically
#include <string>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_util_string_functions.h
  //! @brief functions for splitting/parsing strings and converting strings to numeric values
  //!
  //! @see @link example_util_string_functions.cpp @endlink
  //! @author woetzen, karakam, mendenjl
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace util
  {

    //! @brief test whether a given string begins with another string
    //! @param STRING the full string
    //! @param TEST_PREFIX the prefix
    //! @return true if STRING begins with TEST_PREFIX
    BCL_API bool StartsWith( const std::string &STRING, const std::string &TEST_PREFIX);

    //! @brief test whether a given string ends with another string
    //! @param STRING the full string
    //! @param TEST_SUFFIX the prefix
    //! @return true if STRING ends with TEST_SUFFIX
    BCL_API bool EndsWith( const std::string &STRING, const std::string &TEST_SUFFIX);

    //! @brief Convert characters [A-Z] in a string converted to lower case
    //! @param STRING the full string
    //! @return the string with all upper-case characters converted to lower case
    BCL_API std::string ToLower( const std::string &STRING);

    //! @brief Convert characters [a-z] in a string converted to upper case
    //! @param STRING the full string
    //! @return the string with all lower-case characters converted to upper case
    BCL_API std::string ToUpper( const std::string &STRING);

    //! @brief Remove a specified set of starting and trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning and end of the string
    //! @return the string stripped of starting and ending characters found in CHARS_TO_STRIP
    BCL_API std::string Strip( const std::string &STRING, const std::string &CHARS_TO_STRIP);

    //! @brief Remove a specified set of starting characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning of the string
    //! @return the string stripped of starting characters found in CHARS_TO_STRIP
    BCL_API std::string LStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP);

    //! @brief Remove a specified set of trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the end of the string
    //! @return the string stripped of ending characters found in CHARS_TO_STRIP
    BCL_API std::string RStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP);

    //! @brief repeat a string a particular number of times
    //! @param STRING the string to repeat
    //! @param REPETITIONS the number of times to repeat the string
    //! @return the string, repeated REPETITIONS times
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    BCL_API std::string Repeat( const std::string &STRING, const size_t &REPETITIONS);

    //! @brief join a vector of strings with a different internal string
    //! @param JOINER the string to insert between consecutive elements of the vector of strings
    //! @param STRINGS array of strings to join with the JOINER
    //! @return STRINGS(0){JOINER}{STRINGS(1)}{JOINER}...{STRINGS(N)}
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    BCL_API std::string Join( const std::string &JOINER, const storage::Vector< std::string> &STRINGS);

    //! returns a vector of std::string from an array of char *
    BCL_API
    storage::Vector< std::string>
    StringListFromCharacterArray( const int NUMBER_ARGUMENTS, const char *ARGUMENTS[]);

    //! reads all strings from an istream and stores them consecutively in a stroagevector of strings
    BCL_API
    storage::Vector< std::string>
    StringListFromIStream( std::istream &ISTREAM);

    //! @brief StringLineListFromIStream reads all lines from an istream and stores them consecutively as strings
    //! @param ISTREAM the stream from which the lines will be read
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    BCL_API
    storage::Vector< std::string>
    StringLineListFromIStream( std::istream &ISTREAM);

    //! @brief SplittedStringLineListFromIStream reads all lines from an istream, splits them, and then stores them
    //! @param ISTREAM the stream from which the lines will be read
    //! @param SPLITTER the character to be used to split the line strings
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    BCL_API
    storage::Vector< storage::Vector< std::string> >
    SplittedStringLineListFromIStream
    (
      std::istream &ISTREAM, const std::string &SPLITTER = " \n\t\r\0"
    );

    //! returns a vector of substrings of STRING, that don't contain SPLITTER
    BCL_API
    storage::Vector< std::string>
    SplitString( const std::string &STRING, const std::string &SPLITTER = " \n\t\r\0");

    //! returns a string which has all REPLACE replaced with REPLACEMENT in ORIGINAL
    BCL_API
    std::string
    ReplaceString( const std::string &ORIGINAL, const std::string &REPLACE, const std::string &REPLACEMENT);

    //! trims spaces from beginning and end of a copied string
    BCL_API
    std::string
    TrimString( const std::string &STRING);

    //! test whether string is numerical (double, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    BCL_API
    bool
    IsNumerical( const std::string &STRING);

    //! Find the number of consecutive characters in a string that can be converted to a single floating-point type
    //! (e.g. float, double), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a double or float
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    BCL_API
    size_t
    LengthOfFloatingPointType( const std::string &STRING, const size_t &START = 0);

    //! Find the number of consecutive characters in a string that can be converted to a single integer type
    //! (e.g. long, int, short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long, integer, or short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    BCL_API
    size_t
    LengthOfIntegerType( const std::string &STRING, const size_t &START = 0);

    //! Find the number of consecutive characters in a string that can be converted to a single unsigned integer type
    //! (e.g. size_t or unsigned long/int/short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long/integer/short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    BCL_API
    size_t
    LengthOfUnsignedIntegerType( const std::string &STRING, const size_t &START = 0);

    //! searches for spaces and removes them from the string
    BCL_API
    std::string
    RemoveSpacesFromString( const std::string &STRING);

    //! @brief ConvertStringToBoolean takes a string and converts it to a boolean value
    //! @param STRING the string that will be converted to a boolean, should be either "true" or "false"
    //! @return boolean either true of false depending on STRING
    BCL_API
    bool
    ConvertStringToBoolean( const std::string &STRING);

    //! converts std::string to object of template typename t_DataType
    template< typename t_DataType>
    inline
    t_DataType
    ConvertStringToNumericalValue( const std::string &STRING)
    {
      // remove tabs, spaces, new lines, etc.
      const std::string string( TrimString( STRING));

      // if string is empty an undefined is returned
      if( string.empty())
      {
        return GetUndefined< t_DataType>();
      }
      t_DataType new_t;
      BCL_Assert( TryConvertFromString( new_t, string, GetLogger()), "Numeric conversion of " + string + " failed");
      return new_t;
    }

    //! reads the next t_DataType from the stream
    template< typename t_DataType>
    t_DataType ReadNumberFromStream( std::istream &ISTREAM)
    {
      // make sure the stream is still good
      BCL_Assert( ISTREAM.good(), "Missing number before end of file / stream");

      // skip over any whitespace
      ISTREAM >> std::ws;

      // make sure the stream is still good
      BCL_Assert( ISTREAM.good(), "Missing number before end of file / stream");

      // if the next character is an 'n', then we're looking at a nan
      if( ISTREAM.peek() == 'n')
      {
        // read the next 3 characters to make sure they spell "nan"
        static char buffer[ 4];
        ISTREAM.read( buffer, 3);

        // ensure that the string really was nan
        BCL_Assert( std::string( buffer) == "nan", "Expected nan, but found: " + std::string( buffer));

        // string was nan, so return the undefined t_DataType
        return GetUndefined< t_DataType>();
      }

      // make sure we have not reached the end of file
      BCL_Assert( ISTREAM.good(), "Missing number before end of file / stream");

      // the value is not an nan, so read it in
      t_DataType value;
      ISTREAM >> value;

      return value;
    }

    //! returns a vector of t_DataType of STRING, that don't contain SPLITTER
    template< typename t_DataType>
    inline
    storage::Vector< t_DataType>
    SplitStringToNumerical( const std::string &STRING, const std::string &SPLITTER = " \n\t\r\0")
    {
      // remove tabs, spaces, new lines, etc.
      const storage::Vector< std::string> result( SplitString( TrimString( STRING), SPLITTER));
      storage::Vector< t_DataType> converted;

      // iterate over all strings and convert them to t_DataType
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_str( result.Begin()),
          itr_str_end( result.End());
        itr_str != itr_str_end; ++itr_str
      )
      {
        // skip empty
        if( itr_str->empty())
        {
          continue;
        }

        converted.PushBack( ConvertStringToNumericalValue< t_DataType>( *itr_str));
      }

      // end
      return converted;
    }

    //! @brief wrap a string for a given line length into multiple lines
    //! @details like split string, only that it will also split, if the resulting string is longer that LINE_LENGTH
    //! @param STRING the string to be wrapped
    //! @param LINE_LENGTH length of lines
    //! @param WRAPPER_CHARS all characters to wrap at
    //! @return vector of lines
    BCL_API
    storage::Vector< std::string>
    WrapString( const std::string &STRING, const size_t LINE_LENGTH, const std::string &WRAPPER_CHARS = " \n\t\r\0");

    //! ignore lines that start with '#' and empty lines
    BCL_API
    void
    ChopHeader( std::istream &ISTREAM);

    //! returns the maximum matching string in terms of pair< size_t, size_t>( index, length) ~ substr( index, length) within first and second string
    BCL_API
    storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
    MaximumMatchingSubstring( const std::string &STRING_A, const std::string &STRING_B, const size_t &MINIMUM_SIZE = 5);

    //! @brief find the matching brackets reverse
    //! @param STRING the string of interest
    //! @param POSITION the position of the closing bracket
    //! @param BRACKET_OPEN the character for opening bracket
    //! @param BRACKET_CLOSE the character for closing bracket
    //! @return the Position of the matching opening bracket
    BCL_API
    size_t
    RFindMatchingBracket
    (
      const std::string &STRING,
      const size_t POSITION,
      const char BRACKET_OPEN,
      const char BRACKET_CLOSE
    );

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_STRING_FUNCTIONS_H_
