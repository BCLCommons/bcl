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
#include "util/bcl_util_undefined.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#if defined (__APPLE__)
  #include <math.h>
#elif defined (__GNUC__)
  #include <math.h>     // _isnan function in Linux located in math.h
#elif defined (_MSC_VER) //this block exist to avoid unneccesarry warnings with Visual C++
  #include <float.h>    // isnan function in Windows located in float.h
  #define isnan _isnan  // in float.h isnan is called _isnan
  #define isinf !_finite
#endif
//this block is necessary since for mingw and cygwin, isnan and isinf is defined in namespace std
#if defined (__MINGW32__) || defined(__CYGWIN__)
  #include <cmath> // either include <math.h> and use isnan or include <cmath> and use std::isnan
  using std::isnan;
  using std::isinf;
#endif

namespace bcl
{
  namespace util
  {

    //! @brief specialization of template GetUndefined for float
    template<>
    const float &GetUndefined< float>()
    {
      static const float s_undefined( std::numeric_limits< float>::quiet_NaN());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for double
    template<>
    const double &GetUndefined< double>()
    {
      static const double s_undefined( std::numeric_limits< double>::quiet_NaN());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for short
    template<>
    const short &GetUndefined< short>()
    {
      static const short s_undefined( std::numeric_limits< short>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned short
    template<>
    const unsigned short &GetUndefined< unsigned short>()
    {
      static const unsigned short s_undefined( std::numeric_limits< unsigned short>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for int
    template<>
    const int &GetUndefined< int>()
    {
      static const int s_undefined( std::numeric_limits< int>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned int
    template<>
    const unsigned int &GetUndefined< unsigned int>()
    {
      static const unsigned int s_undefined( std::numeric_limits< unsigned int>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for long
    template<>
    const long &GetUndefined< long>()
    {
      static const long s_undefined( std::numeric_limits< long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned long
    template<>
    const unsigned long &GetUndefined< unsigned long>()
    {
      static const unsigned long s_undefined( std::numeric_limits< unsigned long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for long long
    template<>
    const long long &GetUndefined< long long>()
    {
      static const long long s_undefined( std::numeric_limits< long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned long
    template<>
    const unsigned long long &GetUndefined< unsigned long long>()
    {
      static const unsigned long long s_undefined( std::numeric_limits< unsigned long long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for char
    template<>
    const char &GetUndefined< char>()
    {
      static const char s_undefined( std::numeric_limits< char>::max());
      return s_undefined;
    }

    const double &GetUndefinedDouble()
    {
      static const double s_undefined( std::numeric_limits< double>::quiet_NaN());
      return s_undefined;
    }

    const size_t &GetUndefinedSize_t()
    {
      static const size_t s_undefined( std::numeric_limits< size_t>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for bool
    template<>
    BCL_API
    const bool &GetUndefined< bool>()
    {
      // it is not really possible to define an undefined bool; but in most contexts, false will do
      static const bool s_undefined( false);
      return s_undefined;
    }

    //! function to return whether the supplied UndefinedObject is defined or not
    bool IsDefined( const UndefinedObject &OBJECT)
    {
      return false;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for double
    bool IsDefined( const double DOUBLE)
    {
      return isnan( DOUBLE) == 0 && isinf( DOUBLE) == 0;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for float
    bool IsDefined( const float FLOAT)
    {
      return isnan( FLOAT) == 0 && isinf( FLOAT) == 0;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned int
    bool IsDefined( const unsigned int INT)
    {
      return INT != GetUndefined< unsigned int>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for int
    bool IsDefined( const int INT)
    {
      return INT != GetUndefined< int>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long
    bool IsDefined( const unsigned long LONG)
    {
      return LONG != GetUndefined< unsigned long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for long
    bool IsDefined( const long LONG)
    {
      return LONG != GetUndefined< long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned short
    bool IsDefined( const unsigned short SHORT)
    {
      return SHORT != GetUndefined< unsigned short>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for short
    bool IsDefined( const short SHORT)
    {
      return SHORT != GetUndefined< short>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long long
    bool IsDefined( const unsigned long long LONG)
    {
      return LONG != GetUndefined< unsigned long long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for long long
    bool IsDefined( const long long LONG)
    {
      return LONG != GetUndefined< long long>();
    }

    //! function to return whether the supplied value is not a number
    bool IsNaN( const double &DOUBLE)
    {
      return DOUBLE != DOUBLE;
    }

    //! function to return whether the supplied value is not a number
    bool IsNaN( const float &FLOAT)
    {
      return FLOAT != FLOAT;
    }

  } // namespace util
} // namespace bcl
