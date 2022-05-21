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

#ifndef BCL_UTIL_BINARY_FUNCTION_INTERFACE_H_
#define BCL_UTIL_BINARY_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_class_descriptor.h"
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionInterface
    //! @brief This class is a Binary Function Interface class
    //! @details It is supposed to be used a a base class for every mathematical purposes, if you want to minimize, or
    //! sum up functions
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterface :
      virtual public ObjectInterface
    {
    public:

      //! Typedefs with standard STL names
      //! These allow other functions to know the argument and result types at compile time
      //! These types can/should be used as template parameters to other templates when it will
      //! enhance clarity or generality (typically by reducing the number of template parameters)
      //! For example of how they should be used, see util::BinaryFunctionSTLWrapper
      //! @note This is equivalent to deriving from std::binary_function<>, but keeps us from inheriting
      //! @note from non-virtual bases and unnecessary multiple inheritance
      typedef const t_ArgumentType1 first_argument_type;  //!< type of the first argument
      typedef const t_ArgumentType2 second_argument_type; //!< type of the first argument
      typedef t_ResultType          result_type;          //!< type of the return type

    ///////////////
    // operators //
    ///////////////

      //! @brief operator the implements the binary operation on the two arguments returning a result
      //! @param ARGUMENT1 argument 1
      //! @param ARGUMENT2 argument 2
      //! @return the Result of the binary operation
      virtual result_type operator()( first_argument_type &ARGUMENT1, second_argument_type &ARGUMENT2) const = 0;

    }; //end template class BinaryFunctionInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_BINARY_FUNCTION_INTERFACE_H_
