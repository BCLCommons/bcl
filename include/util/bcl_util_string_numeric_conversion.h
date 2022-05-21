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

#ifndef BCL_UTIL_STRING_NUMERIC_CONVERSION_H_
#define BCL_UTIL_STRING_NUMERIC_CONVERSION_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_util_string_numeric_conversion.h
  //! @brief functions for splitting/parsing strings and converting strings to numeric values
  //! @detail These functions are higher performance and more error tolerant than similar functions in util_string_functions,
  //!         and may eventually replace all numeric conversion functions in that header
  //! @see @link example_util_string_numeric_conversion.cpp @endlink
  //! @author mendenjl
  //! @date Oct 16, 2012
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace util
  {
    // floating point types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( double &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( float &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( unsigned char &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( char &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( signed char &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( short &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( int &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( long &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( unsigned short &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( unsigned int &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( unsigned long &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( unsigned long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    // bool

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString( bool &TYPE, const std::string &STR, std::ostream &ERR_STREAM);

    // Overloads including desired range
    // floating point types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      double &TYPE,
      const double &MIN,
      const double &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      float &TYPE,
      const float &MIN,
      const float &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      unsigned char &TYPE,
      const unsigned char &MIN,
      const unsigned char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      char &TYPE,
      const char &MIN,
      const char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      signed char &TYPE,
      const signed char &MIN,
      const signed char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      short &TYPE,
      const short &MIN,
      const short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      int &TYPE,
      const int &MIN,
      const int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      long &TYPE,
      const long &MIN,
      const long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      long long &TYPE,
      const long long &MIN,
      const long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      unsigned short &TYPE,
      const unsigned short &MIN,
      const unsigned short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      unsigned int &TYPE,
      const unsigned int &MIN,
      const unsigned int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      unsigned long &TYPE,
      const unsigned long &MIN,
      const unsigned long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    BCL_API bool TryConvertFromString
    (
      unsigned long long &TYPE,
      const unsigned long long &MIN,
      const unsigned long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    );

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_STRING_NUMERIC_CONVERSION_H_
