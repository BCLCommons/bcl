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

#ifndef BCL_UTIL_FOUR_INPUT_FUNCTION_INTERFACE_H_
#define BCL_UTIL_FOUR_INPUT_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <functional>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FourInputFunctionInterface
    //! @brief This class is a Function Interface class for functions with Four inputs.
    //! @details It is supposed to be used a a base class for every mathematical purposes, if you want to minimize, or sum up functions
    //!
    //! @remarks example unnecessary
    //! @author riddeljs
    //! @date 12.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ArgumentType4, typename t_ResultType>
    class FourInputFunctionInterface :
      public ObjectInterface
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      virtual FourInputFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ArgumentType4, t_ResultType> *Clone() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator the implements the operation on the Four arguments returning a result
      //! @param ARGUMENT1 argument 1
      //! @param ARGUMENT2 argument 2
      //! @param ARGUMENT3 argument 3
      //! @param ARGUMENT4 argument 4
      //! @return the Result of the operation
      virtual t_ResultType operator()
      (
        t_ArgumentType1 &ARGUMENT1,
        t_ArgumentType2 &ARGUMENT2,
        t_ArgumentType3 &ARGUMENT3,
        t_ArgumentType4 &ARGUMENT4
      ) = 0;

    }; //end template class FourInputFunctionInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_FOUR_INPUT_FUNCTION_INTERFACE_H_
