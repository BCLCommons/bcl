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

#ifndef BCL_UTIL_BINARY_FUNCTION_INTERFACE_SERIALIZABLE_H_
#define BCL_UTIL_BINARY_FUNCTION_INTERFACE_SERIALIZABLE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_binary_function_interface.h"
#include "bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionInterfaceSerializable
    //! @brief This class is a Binary Function Interface class
    //! @details It is supposed to be used a a base class for every mathematical purposes, if you want to minimize, or
    //! sum up functions
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Sep 16, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterfaceSerializable :
      public BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>,
      public SerializableInterface
    {
    public:

      //! @brief virtual copy constructor
      virtual BinaryFunctionInterfaceSerializable *Clone() const = 0;

    }; //end template class BinaryFunctionInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_BINARY_FUNCTION_INTERFACE_SERIALIZABLE_H_
