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

#ifndef BCL_UTIL_FUNCTION_INTERFACE_H_
#define BCL_UTIL_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_function_interface_nonconst.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionInterface
    //! @brief This class is a Function Interface class
    //! @details It is supposed to be used as a base class for every mathematical purposes,
    //! if you want to minimize, or sum up functions
    //!
    //! @remarks example unnecessary
    //!
    //! @author woetzen
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterface :
      public FunctionInterfaceNonConst< const t_ArgumentType, t_ResultType>
    {
      // FunctionInterface is an alias for a unary function with a const argument

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual FunctionInterface *Clone() const = 0;

    }; // template class FunctionInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_FUNCTION_INTERFACE_H_
