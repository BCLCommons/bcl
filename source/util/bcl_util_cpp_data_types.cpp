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
#include "util/bcl_util_cpp_data_types.h"

// includes from bcl - sorted alphabetically
#include "type/bcl_type_compare.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief conversion to a string from a Types
    //! @param DATA_TYPE the data type to get a string for
    //! @return a string representing that data type
    const std::string &CPPDataTypes::GetCPPDatatypeName( const CPPDataTypes::Types &DATA_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "DT_UNKNOWN",
        "DT_CHAR",
        "DT_DOUBLE",
        "DT_FLOAT",
        "DT_INT",
        "DT_UNSIGNED_INT",
        "DT_SIZE_T",
        "DT_STRING",
        "DT_BOOL",
        GetStaticClassName< CPPDataTypes::Types>()
      };
      return s_descriptors[ size_t( DATA_TYPE)];
    }

    //! @brief c++ string for datatype
    //! @param DATA_TYPE the data type to get a string for
    //! @return the string representing that data type as it would be written in the source code
    const std::string &CPPDataTypes::GetCPPString( const CPPDataTypes::Types &DATA_TYPE)
    {
      static const std::string s_strings[] =
      {
        "",
        "char",
        "double",
        "float",
        "int",
        "unsigned int",
        "size_t",
        "std::string",
        "bool",
        ""
      };
      return s_strings[ size_t( DATA_TYPE)];
    }

    //! @brief Types from template parameter for char
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< char>()
    {
      return e_Char;
    }

    //! @brief Types from template parameter for double
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< double>()
    {
      return e_Double;
    }

    //! @brief Types from template parameter for float
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< float>()
    {
      return e_Float;
    }

    //! @brief Types from template parameter for int
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< int>()
    {
      return e_Int;
    }

    //! @brief Types from template parameter for unsigned ints
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< unsigned int>()
    {
      // Prefer the real type (unsigned int) to the typedef-ed size_t,
      //   unsigned int is often the same as size_t on 32-bit architectures
      return e_UnsignedInt;
    }

    //! @brief Types from template parameter for unsigned longs
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< unsigned long>()
    {
      return type::Compare< unsigned long, size_t>::e_Same ? e_SizeT : e_UnsignedInt;
    }

    //! @brief Types from template parameter for string
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< std::string>()
    {
      return e_String;
    }

    //! @brief Types from template parameter for bool
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< bool>()
    {
      return e_Bool;
    }

  } // namespace util

} // namespace bcl
