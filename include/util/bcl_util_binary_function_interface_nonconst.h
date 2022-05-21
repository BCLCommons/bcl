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

#ifndef BCL_UTIL_BINARY_FUNCTION_INTERFACE_NONCONST_H_
#define BCL_UTIL_BINARY_FUNCTION_INTERFACE_NONCONST_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionInterfaceNonConst
    //! @brief This class is a Binary Function Interface class without constants
    //! @details This should be able to replace bcl_util_binary_function_interface since
    //! constness of arguments/function can be delegated to the design of those functions.
    //!
    //! @remarks example unnecessary
    //! @author riddeljs
    //! @date 02.11.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterfaceNonConst :
      public ObjectInterface
    {
    public:

      //! Typedefs with standard STL names to allow functors to use their argument and result types at compile time
      //! These types can/should be used as template parameters to other templates when it will
      //! enhance clarity or generality (typically by reducing the number of template parameters)
      //! For example of how they should be used, see util::BinaryFunctionSTLWrapper
      //! @note This is equivalent to deriving from std::binary_function<>, but keeps us from inheriting
      //! @note from non-virtual bases and unnecessary multiple inheritance
      typedef t_ArgumentType1 first_argument_type;  //!< type of the first argument
      typedef t_ArgumentType2 second_argument_type; //!< type of the first argument
      typedef t_ResultType    result_type;          //!< type of the return type

    ///////////////
    // operators //
    ///////////////

      //! @brief operator the implements the binary operation on the two arguments returning a result
      //! @param ARGUMENT1 argument 1
      //! @param ARGUMENT2 argument 2
      //! @return the Result of the binary operation
      virtual t_ResultType operator()( t_ArgumentType1 &ARGUMENT1, t_ArgumentType2 &ARGUMENT2) const = 0;

    }; //end template class BinaryFunctionInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_BINARY_FUNCTION_INTERFACE_NONCONST_H_
