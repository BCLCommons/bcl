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

#ifndef BCL_FUNCTION_UNARY_INTERFACE_H_
#define BCL_FUNCTION_UNARY_INTERFACE_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace function
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UnaryInterface
    //! @brief this class is a interface for functions, that takes one argument and returns an argument
    //! it has the for result F( x)
    //!
    //! @tparam t_ArgumentType is the type of x
    //! @tparam t_ResultType is the type of Result
    //!
    //! It is supposed to be used a a base class for every function class.
    //! It requires to overwrite an operator(), which takes t_ArgumentType and returns t_ResultType.
    //! So one can basically define a function Square like:
    //! class Square :
    //!   public UnaryInterface< const double, double>
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
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 08/01/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class UnaryInterface :
      public util::ObjectInterface
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
      virtual UnaryInterface< t_ArgumentType, t_ResultType> *Clone() const = 0;

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
      virtual t_ResultType operator()( t_ArgumentType &ARGUMENT) const = 0;

    }; // template class UnaryInterface

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_UNARY_INTERFACE_H_ 
