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

#ifndef BCL_MATH_FUNCTION_INTERFACE_H_
#define BCL_MATH_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionInterface
    //! @brief This class is a Function Interface class
    //! @details It is supposed to be used as a base class for every mathematical function.
    //! It requires to overwrite an operator(), which takes t_ArgumentType and returns t_ResultType.
    //! So one can basically define a function Square like:
    //! class Square :
    //!   public FunctionInterfaceSerializable< double, double>
    //! {
    //!   ...
    //!   double operator()( const double &ARGUMENT) const
    //!   {
    //!     return ARGUMENT * ARGUMENT;
    //!   }
    //! }
    //!
    //! it also has three useful functions, that can be overwritten, which are:
    //! Scheme - one can output some information about the function like "square"
    //! WriteDetailedSchemeAndValues - one can write even more details like "square of number 3 = 3 * 3 = 9"
    //!
    //! @remarks Interface may be removed/replaced by function::Interface in the future
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterface :
      public util::FunctionInterface< t_ArgumentType, t_ResultType>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual FunctionInterface *Clone() const = 0;

      //! @brief destructor
      virtual ~FunctionInterface()
      {
      }

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

      //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      virtual t_ResultType operator()( const t_ArgumentType &ARGUMENT) const = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      virtual std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType &ARGUMENT,
        std::ostream &OSTREAM
      ) const
      {
        // end
        return OSTREAM;
      }

    }; // template class FunctionInterface

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_FUNCTION_INTERFACE_H_
