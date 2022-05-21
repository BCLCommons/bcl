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

#ifndef BCL_UTIL_CPP_DATA_TYPES_H_
#define BCL_UTIL_CPP_DATA_TYPES_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_util_cpp_data_types.h
  //! @brief standard c++ datatypes
  //!
  //! @see @link example_util_cpp_data_types.cpp @endlink
  //! @author woetzen
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace util
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CPPDataTypes
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_util_cpp_data_types.cpp @endlink
    //! @author woetzen
    //! @date Apr 30, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CPPDataTypes
    {
    public:

    ///////////
    // types //
    ///////////

      //! enum standard c++ datatypes
      enum Types
      {
        e_Unknown,        //!< unknown data type
        e_Char,           //!< char data type
        e_Double,         //!< double data type
        e_Float,          //!< float data type
        e_Int,            //!< int data type
        e_UnsignedInt,    //!< unsigned int data type
        e_SizeT,          //!< size_t data type
        e_String,         //!< string data type
        e_Bool,           //!< bool data type
        s_NumberDataTypes //!< # of data types
      };

      //! @brief conversion to a string from a Types
      //! @param DATA_TYPE the data type to get a string for
      //! @return a string representing that data type
      static const std::string &GetCPPDatatypeName( const Types &DATA_TYPE);

      //! @brief c++ string for datatype
      //! @param DATA_TYPE the data type to get a string for
      //! @return the string representing that data type as it would be written in the source code
      static const std::string &GetCPPString( const Types &DATA_TYPE);

      //! @brief enum class wrapper for Types
      typedef WrapperEnum< Types, &GetCPPDatatypeName, s_NumberDataTypes> TypeEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CPPDataTypes();

    public:

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief Types from template parameter
      //! @return Types for given t_DataType
      template< typename t_DataType>
      static TypeEnum DataTypeFromTemplate();

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class CPPDataTypes

  } // namespace util

} // namespace bcl

#endif //BCL_UTIL_CPP_DATA_TYPES_H_
