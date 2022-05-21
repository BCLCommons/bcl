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
#include "util/bcl_util_string_numeric_conversion.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <iostream>
#include <map>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

namespace bcl
{
  namespace util
  {
    // floating point types

    //! @brief a helper function to test for a string that is nan, inf, +inf, or -inf (all case-insensitive)
    //! @param TYPE variable to try to convert STRING into
    //! @param STR string of interest
    //! @return true on success
    template< typename t_DataType>
    bool TryConvertNanInfString( t_DataType &TYPE, const std::string &STR)
    {
      if( STR.size() == size_t( 3))
      {
        if( !strcasecmp( STR.c_str(), "nan"))
        {
          TYPE = GetUndefined< t_DataType>();
          return true;
        }
        else if( !strcasecmp( STR.c_str(), "inf"))
        {
          TYPE = std::numeric_limits< t_DataType>::infinity();
          return true;
        }
      }
      else if( STR.size() == size_t( 4) && !strcasecmp( STR.c_str() + 1, "inf"))
      {
        // check for +/- inf as well
        if( STR[ 0] == '+')
        {
          TYPE = std::numeric_limits< t_DataType>::infinity();
          return true;
        }
        else if( STR[ 0] == '-')
        {
          TYPE = -std::numeric_limits< t_DataType>::infinity();
          return true;
        }
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( double &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // convert to value
      TYPE = strtod( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length == size_t( 0) && TryConvertNanInfString( TYPE, STR))
        {
          return true;
        }
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( float &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // declare a variable to hold the end position for string, which strtof will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // convert to value
      TYPE = strtof( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        // handle the underflow condition; not really a problem, just set to 0
        if( STR.find( "e-") != std::string::npos)
        {
          ERR_STREAM << "Warning; interpreting " << STR << " as 0 because it is too small to fit in a float\n";
          errno = current_err;
          return true;
        }
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nan's specifically, since the current mingw cross-compiler does not handle them
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length == size_t( 0) && TryConvertNanInfString( TYPE, STR))
        {
          return true;
        }
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }
      return true;
    }

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      TYPE = STR[ 0];
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      char &TYPE,
      const char &MIN,
      const char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      TYPE = STR[ 0];
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( signed char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      return true;
    }

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( short &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          long( std::numeric_limits< short>::min()),
          long( std::numeric_limits< short>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = short( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( int &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          long( std::numeric_limits< int>::min()),
          long( std::numeric_limits< int>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = int( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtol( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }
      else if( TYPE < 0)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoll( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }
      else if( TYPE < 0)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned short &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          (unsigned long)( std::numeric_limits< unsigned short>::min()),
          (unsigned long)( std::numeric_limits< unsigned short>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = ( unsigned short)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned int &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          (unsigned long)( std::numeric_limits< unsigned int>::min()),
          (unsigned long)( std::numeric_limits< unsigned int>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = ( unsigned int)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoul( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoull( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    // bool

    //! @brief create a static map containing all valid boolean strings
    std::map< std::string, bool> GetBoolStrMap()
    {
      // these are from http://yaml.org/type/bool.html, except 0/1, which are bcl extensions
      std::map< std::string, bool> bool_map;
      bool_map[ "0"]     = bool_map[ "n"]     = bool_map[ "N"]     = false;
      bool_map[ "no"]    = bool_map[ "No"]    = bool_map[ "NO"]    = false;
      bool_map[ "false"] = bool_map[ "False"] = bool_map[ "FALSE"] = false;
      bool_map[ "off"]   = bool_map[ "Off"]   = bool_map[ "OFF"]   = false;
      bool_map[ "1"]     = bool_map[ "y"]     = bool_map[ "Y"]     = true;
      bool_map[ "yes"]   = bool_map[ "Yes"]   = bool_map[ "YES"]   = true;
      bool_map[ "true"]  = bool_map[ "True"]  = bool_map[ "TRUE"]  = true;
      bool_map[ "on"]    = bool_map[ "On"]    = bool_map[ "ON"]    = true;
      return bool_map;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( bool &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      static const std::map< std::string, bool> s_bool_map( GetBoolStrMap());
      std::map< std::string, bool>::const_iterator itr( s_bool_map.find( STR));
      if( itr == s_bool_map.end())
      {
        ERR_STREAM << "Expected a bool value (0/n/N/no/No/No/false/False/FALSE/off/Off/OFF/1/y/Y/yes/Yes/YES/true/True/TRUE/on/On/ON) but received: " << STR;
        return false;
      }
      TYPE = itr->second;
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      double &TYPE,
      const double &MIN,
      const double &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtod( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length != size_t( 0) || !TryConvertNanInfString( TYPE, STR))
        {
          ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }

      if( IsNaN( TYPE) || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      float &TYPE,
      const float &MIN,
      const float &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtof( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length != size_t( 0) || !TryConvertNanInfString( TYPE, STR))
        {
          ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }

      if( IsNaN( TYPE) || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }
      return true;
    }

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned char &TYPE,
      const unsigned char &MIN,
      const unsigned char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      if( TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << "Expected a character in the range " << MIN << '-' << MAX << ", but got \"" << STR << "\"";
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      signed char &TYPE,
      const signed char &MIN,
      const signed char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      if( TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << "Expected a character in the range " << MIN << '-' << MAX << ", but got \"" << STR << "\"";
        return false;
      }
      return true;
    }

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      short &TYPE,
      const short &MIN,
      const short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, long( MIN), long( MAX), STR, ERR_STREAM))
      {
        TYPE = short( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      int &TYPE,
      const int &MIN,
      const int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, long( MIN), long( MAX), STR, ERR_STREAM))
      {
        TYPE = int( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      long &TYPE,
      const long &MIN,
      const long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtol( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
          return false;
        }
      }
      else if( TYPE < 0 || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      long long &TYPE,
      const long long &MIN,
      const long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoll( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
          return false;
        }
      }
      else if( TYPE < 0 || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned short &TYPE,
      const unsigned short &MIN,
      const unsigned short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, ( unsigned long)( MIN), ( unsigned long)( MAX), STR, ERR_STREAM))
      {
        TYPE = ( unsigned short)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned int &TYPE,
      const unsigned int &MIN,
      const unsigned int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, ( unsigned long)( MIN), ( unsigned long)( MAX), STR, ERR_STREAM))
      {
        TYPE = ( unsigned int)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned long &TYPE,
      const unsigned long &MIN,
      const unsigned long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoul( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-' || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned long long &TYPE,
      const unsigned long long &MIN,
      const unsigned long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoull( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-' || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

  } // namespace util
} // namespace bcl
