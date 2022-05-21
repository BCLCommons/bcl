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
#include "util/bcl_util_data_type.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace util
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &DataType::GetTypeName( const DataType::Type &TYPE)
    {
      static const std::string s_names[ s_NumberTypes + 1] =
      {
        "Bool",
        "Char",
        "Int",
        "UnsignedInt",
        "Float",
        "String",
        "Scalar",
        "Sequence",
        "Map",
        "Set",
        "StaticObject",
        "DynamicObject",
        "Pointer",
        GetStaticClassName< DataType::Type>()
      };
      return s_names[ TYPE];
    }

    //! @brief detect whether a given type must be scalar, e.g. not have any sub-arguments
    //! @param TYPE the type
    //! @return bool - true if the type must be scalar
    bool DataType::TestMustBeScalar( const DataType::Type &TYPE)
    {
      return TYPE <= s_NumberScalarTypes;
    }

    //! @brief test whether the underlying cpp type for a given object is known
    //! @param TYPE the type
    //! @return bool - true if the type is known
    bool DataType::IsUnderlyingTypeKnown( const Type &TYPE)
    {
      return TYPE < e_DynamicObject;
    }

    //! @brief convenience function for containers to write out their size requirements
    //! @param OSTREAM stream to write out size requirements to
    //! @param MIN_SIZE the minimum size of the container
    //! @param MAX_SIZE the maximum size of the container
    void DataType::WriteSizeRequirements
    (
      std::ostream &OSTREAM,
      const size_t &MIN_SIZE,
      const size_t &MAX_SIZE
    )
    {
      if( MIN_SIZE > size_t( 0))
      {
        OSTREAM << " with ";
        if( MAX_SIZE == std::numeric_limits< size_t>::max())
        {
          OSTREAM << "at least " << MIN_SIZE;
        }
        else if( MIN_SIZE == MAX_SIZE)
        {
          OSTREAM << MIN_SIZE;
        }
        else
        {
          OSTREAM << "between " << MIN_SIZE << " and " << MAX_SIZE;
        }
      }
      else if( MAX_SIZE < std::numeric_limits< size_t>::max())
      {
        OSTREAM << " with at most " << MAX_SIZE;
      }
      OSTREAM << ' ';
    }

  } // namespace util
} // namespace bcl
