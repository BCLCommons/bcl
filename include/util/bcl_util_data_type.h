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

#ifndef BCL_UTIL_DATA_TYPE_H_
#define BCL_UTIL_DATA_TYPE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataType
    //! @brief language independent data types used in the bcl based loosely off @link http://yaml.org/type/ @endlink
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Feb 24, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API DataType
    {
    public:
      //! enum for what the type of an object is, modeled after YAML's set of language-independent types
      enum Type
      {
        e_Bool,              //!< True or False
        e_Char,              //!< Single character
        e_Int,               //!< [+-]?[0-9]+
        e_UnsignedInt,       //!< [0-9]+
        e_Float,             //!< Floating point number
        e_String,            //!< Any string
        s_NumberScalarTypes, //!< Generic: Any scalar type; note, all types after this should be non-scalar
        e_Sequence,          //!< A sequence of items, which themselves may be any type
        e_Map,               //!< Map of name / value pairs, which themselves may be any type
        e_Set,               //!< An ordered sequence; behaves similarly to sequence, except in comparison
        e_StaticObject,      //!< A name / alias, followed by multiple parameters
        e_DynamicObject,     //!< An object whose underlying type cannot be determined statically
        e_Pointer,           //!< Null, or a reference or any of the value types
        s_NumberTypes
      };

      //! @brief get the type name
      //! @param TYPE the type name
      //! @return string describing the type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief detect whether a given type must be scalar, e.g. not have any sub-arguments
      //! @param TYPE the type
      //! @return bool - true if the type must be scalar
      static bool TestMustBeScalar( const Type &TYPE);

      //! @brief test whether the underlying cpp type for a given object is known
      //! @param TYPE the type
      //! @return bool - true if the type is known
      static bool IsUnderlyingTypeKnown( const Type &TYPE);

      //! @brief convenience function for containers to write out their size requirements
      //! @param OSTREAM stream to write out size requirements to
      //! @param MIN_SIZE the minimum size of the container
      //! @param MAX_SIZE the maximum size of the container
      static void WriteSizeRequirements( std::ostream &OSTREAM, const size_t &MIN_SIZE, const size_t &MAX_SIZE);
    };
  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_DATA_TYPE_H_

