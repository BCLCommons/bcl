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

#ifndef BCL_UTIL_FUNCTION_INTERFACE_NONCONST_H_
#define BCL_UTIL_FUNCTION_INTERFACE_NONCONST_H_

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
    //! @class FunctionInterfaceNonConst
    //! @brief This class is a Function Interface class, but the const types are removed
    //! @details It is used for functions which may mutate the argument
    //!
    //! @remarks example unnecessary
    //! @author woetzen, riddeljs, mendenjl
    //! @date 12.31.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceNonConst :
      public virtual ObjectInterface
    {

    public:

      //! Typedefs with standard STL names to allow functors to use their argument and result types at compile time
      //! These types can/should be used as template parameters to other templates when it will
      //! enhance clarity or generality (typically by reducing the number of template parameters)
      //! For example of how they should be used, see util::BinaryFunctionSTLWrapper
      //! @note This is equivalent to deriving from std::unary_function<>, but keeps us from inheriting
      //! @note from non-virtual bases and unnecessary multiple inheritance
      //! Derived classes will inherit these typedefs, so there is no need to redeclare them elsewhere
      typedef t_ArgumentType argument_type; //!< argument type of operator()
      typedef t_ResultType   return_type;   //!< result type of operator()

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual FunctionInterfaceNonConst *Clone() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an argument (possibly void) and returning a t_ResultType object
      //! @param ARGUMENT argument of interest
      //! @return the resultant object
      virtual return_type operator()( t_ArgumentType &ARGUMENT) const = 0;

    }; //end template class util::FunctionInterfaceNonCnst

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_FUNCTION_INTERFACE_NONCONST_H_
