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

#ifndef BCL_MATH_BINARY_FUNCTION_INTERFACE_H_
#define BCL_MATH_BINARY_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_binary_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionInterface
    //! @brief is the base class for all binary mathematical functions.
    //! @details This util::BinaryFunctionInterface derived class provides the interface for functions that work
    //! on two arguments to derive from.
    //!
    //! @tparam t_ArgumentType1 Type of the first argument to the function
    //! @tparam t_ArgumentType2 Type of the second argument to the function
    //! @tparam t_ResulType     Type of the result of the function
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Jun 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterface :
      public util::BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    {
    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      virtual const std::string &GetScheme() const
      {
        return this->GetClassIdentifier();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking two ARGUMENTs and returning a t_ResultType object
      //! @param ARGUMENT1 First argument to be used to evaluate the function
      //! @param ARGUMENT2 Second argument to be used to evaluate the function
      //! @return function value of the given argument
      virtual t_ResultType operator()( const t_ArgumentType1 &ARGUMENT1, const t_ArgumentType2 &ARGUMENT2) const = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT1 First argument to be used to evaluate the function
      //! @param ARGUMENT2 Second argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      virtual std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType1 &ARGUMENT1,
        const t_ArgumentType2 &ARGUMENT2,
        std::ostream &OSTREAM
      ) const
      {
        // end
        return OSTREAM;
      }

    }; // template class BinaryFunctionInterface

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_BINARY_FUNCTION_INTERFACE_H_
