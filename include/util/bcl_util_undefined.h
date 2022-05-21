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

#ifndef BCL_UTIL_UNDEFINED_H_
#define BCL_UTIL_UNDEFINED_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "type/bcl_type_enable_if.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UndefinedObject
    //! @brief Trivial class used in constructors to indicate that the object should be initialized in an undefined state
    //!
    //! @see @link example_util_undefined.cpp @endlink
    //! @author woetzen, karakam
    //! @date Nov 20, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API UndefinedObject
    {
    }; // class UndefinedObject

    //! @brief template function, that returns the undefined object on access
    //! a function that has a static variable, that is only constructued once, and only if the function is called.
    //! By default it assumes, that there is a constructor, that takes an UndefinedObject::Definitionstate and class
    //! that one with e_Undefined;
    //! @return an undefined object of t_DataType, that should give the Answer to IsDefined( OBJECT) as false

    template< typename t_DataType>
    inline
    const t_DataType &GetUndefined()
    {
      static const t_DataType s_undefined = t_DataType( UndefinedObject());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for double
    template<>
    BCL_API
    const double &GetUndefined< double>();

    //! @brief specialization of template GetUndefined for float
    template<>
    BCL_API
    const float &GetUndefined< float>();

    //! @brief specialization of template GetUndefined for short
    template<>
    BCL_API
    const short &GetUndefined< short>();

    //! @brief specialization of template GetUndefined for unsigned short
    template<>
    BCL_API
    const unsigned short &GetUndefined< unsigned short>();

    //! @brief specialization of template GetUndefined for int
    template<>
    BCL_API
    const int &GetUndefined< int>();

    //! @brief specialization of template GetUndefined for unsigned int
    template<>
    BCL_API
    const unsigned int &GetUndefined< unsigned int>();

    //! @brief specialization of template GetUndefined for long
    template<>
    BCL_API
    const long &GetUndefined< long>();

    //! @brief specialization of template GetUndefined for unsigned long
    template<>
    BCL_API
    const unsigned long &GetUndefined< unsigned long>();

    //! @brief specialization of template GetUndefined for long long
    template<>
    BCL_API
    const long long &GetUndefined< long long>();

    //! @brief specialization of template GetUndefined for unsigned long long
    template<>
    BCL_API
    const unsigned long long &GetUndefined< unsigned long long>();

    //! @brief specialization of template GetUndefined for char
    template<>
    BCL_API
    const char &GetUndefined< char>();

    //! @brief specialization of template GetUndefined for bool
    template<>
    BCL_API
    const bool &GetUndefined< bool>();

    BCL_API
    const double &GetUndefinedDouble();

    BCL_API
    const size_t &GetUndefinedSize_t();

    //! function to return whether the supplied UndefinedObject is defined or not
    BCL_API
    bool IsDefined( const UndefinedObject &OBJECT);

    //! function to return whether the supplied UndefinedObject is defined or not for double
    BCL_API
    bool IsDefined( const double DOUBLE);

    //! function to return whether the supplied UndefinedObject is defined or not for float
    BCL_API
    bool IsDefined( const float FLOAT);

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned int
    BCL_API
    bool IsDefined( const unsigned int SIZE_T);

    //! function to return whether the supplied UndefinedObject is defined or not for int
    BCL_API
    bool IsDefined( const int INT);

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned short
    BCL_API
    bool IsDefined( const unsigned short SIZE_T);

    //! function to return whether the supplied UndefinedObject is defined or not for short
    BCL_API
    bool IsDefined( const short SHORT);

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long
    BCL_API
    bool IsDefined( const unsigned long SIZE_T);

    //! function to return whether the supplied UndefinedObject is defined or not for long
    BCL_API
    bool IsDefined( const long LONG);

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long long
    BCL_API
    bool IsDefined( const unsigned long long SIZE_T);

    //! function to return whether the supplied UndefinedObject is defined or not for long long
    BCL_API
    bool IsDefined( const long long LONG);

    //! function to return whether the supplied value is not a number
    BCL_API
    bool IsNaN( const double &DOUBLE);

    //! function to return whether the supplied value is not a number
    BCL_API
    bool IsNaN( const float &FLOAT);

    //! function to return whether the supplied value is not a number
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, bool>::Type
      IsNaN( const t_DataType &INTEGER)
    {
      return false;
    }

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_UNDEFINED_H_
