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

#ifndef BCL_MATH_SUM_FUNCTION_H_
#define BCL_MATH_SUM_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_sum_function_mixin.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SumFunction
    //! @brief This class is a template for t_ResultType = Function( t_ArgumentType&) allowing addition
    //! of functions with coefficients.
    //! @details Assuming y being element of t_ResultType, x being element of t_ArgumentType and s being scalar
    //! this class can represent the following functions:
    //! y1 = s1 * g1(x)
    //! y2 = s2 * g2(x)
    //!
    //! f(x) = y0 + y1 + y2 = y0 + s1 * g1(x) + s2 * g2(x)
    //!
    //! operators provide the functionality to divide the sumfunction, add an absolute to y0, add coefficients
    //! s3 * g3(x), add another sumfunction and to create sumfunctions from function interfaces and coefficients.
    //!
    //! @see @link example_math_sum_function.cpp @endlink
    //! @author meilerj, woetzen
    //! @date 27.08.2004
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class SumFunction :
      public SumFunctionMixin< FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> >
    {
    public:

      //! @typedefs
      typedef FunctionInterfaceSerializable< t_ArgumentType, t_ResultType>  t_Interface;
      typedef SumFunctionMixin< t_Interface>                                Base;
      typedef typename Base::Term Term;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SumFunction( const std::string &SCHEME = "") :
        Base( SCHEME)
      {
      }

      //! @brief constructor from base class
      SumFunction( const Base &BASE) :
        Base( BASE)
      {
      }

      //! @brief construct from ABS_CONST
      //! @param ABS_CONST y0 for the sumfunction
      explicit SumFunction( const t_ResultType &ABS_CONST) :
        Base( ABS_CONST)
      {
      }

      //! @brief construct from ShPtr on Function, optional coefficient, and absolute value
      //! y = y0 + s1*g1(x)
      //! @param SP_FUNCTION g1(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunction
      (
        const util::Implementation< t_Interface> &SP_FUNCTION,
        const double &COEFFICIENT,
        const t_ResultType &ABS_CONST
      ) :
        Base( SP_FUNCTION, COEFFICIENT, ABS_CONST)
      {
      }

      //! @brief construct from Function, optional coefficient, and optional absolute value
      //! y = y0 + s1*g1(x)
      //! @param FUNCTION    g1(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunction
      (
        const t_Interface &FUNCTION,
        const double &COEFFICIENT,
        const t_ResultType &ABS_CONST
      ) :
        Base( FUNCTION, COEFFICIENT, ABS_CONST)
      {
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x) + s2*g2(x)
      //! @param SP_FUNCTION_A g1(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param SP_FUNCTION_B g2(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param COEFFICIENT_A s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param COEFFICIENT_B s2(x) as coefficient for g2(x) -> s2 * g2(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunction
      (
        const util::Implementation< t_Interface> &SP_FUNCTION_A,
        const util::Implementation< t_Interface> &SP_FUNCTION_B,
        const double &COEFFICIENT_A,
        const double &COEFFICIENT_B,
        const t_ResultType &ABS_CONST
      ) :
        Base( SP_FUNCTION_A, SP_FUNCTION_B, COEFFICIENT_A, COEFFICIENT_B, ABS_CONST)
      {
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x) + s2*g2(x)
      //! @param FUNCTION_A    g1(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param FUNCTION_B    g2(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param COEFFICIENT_A s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param COEFFICIENT_B s2(x) as coefficient for g2(x) -> s2 * g2(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunction
      (
        const t_Interface &FUNCTION_A,
        const t_Interface &FUNCTION_B,
        const double &COEFFICIENT_A,
        const double &COEFFICIENT_B,
        const t_ResultType &ABS_CONST
      ) :
        Base( FUNCTION_A, FUNCTION_B, COEFFICIENT_A, COEFFICIENT_B, ABS_CONST)
      {
      }

      //! @brief copy constructor
      //! @return pointer to new SumFunction which is a copy of this
      SumFunction *Clone() const
      {
        return new SumFunction( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assign Function
      //! @param FUNCTION g1( x) with s1 = 1
      //! @return reference to this SumFunction
      SumFunction &operator =( const t_Interface &FUNCTION)
      {
        Base::operator =( FUNCTION);

        // return this
        return *this;
      }

      //! @brief assignment operator
      //! @param SUM_FUNCTION sum function this should be assigned to
      //! @return reference to this SumFunction
      SumFunction &operator =( const Base &SUM_FUNCTION)
      {
        // assign data member
        Base::operator =( SUM_FUNCTION);

        // return this
        return *this;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator += Function
      //! @param FUNCTION function g(x) to be added to the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      SumFunction &operator +=( const t_Interface &FUNCTION)
      {
        Base::operator +=( FUNCTION);
        return *this;
      }

      //! @brief operator -= Function
      //! @param FUNCTION function g(x) to be subtracted from the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      SumFunction &operator -=( const t_Interface &FUNCTION)
      {
        Base::operator -=( FUNCTION);
        return *this;
      }

      //! @brief operator += SumFunction
      //! @param SUM_FUNCTION sum function to be added to the sum function as reference
      //! @return reference to this altered sum function
      SumFunction &
      operator +=
      (
        const Base &SUM_FUNCTION
      )
      {
        Base::operator +=( SUM_FUNCTION);

        // return
        return *this;
      }

      //! @brief operator -= SumFunction
      //! @param SUM_FUNCTION sum function to be subtracted from the sum function as reference
      //! @return reference to this altered sum function
      SumFunction &
      operator -=
      (
        const Base &SUM_FUNCTION
      )
      {
        Base::operator -=( SUM_FUNCTION);
        return *this;
      }

      //! @brief operator += VALUE adds y to y0
      //! @param VALUE y to be added to y0
      SumFunction &
      operator +=
      (
        const t_ResultType &VALUE
      )
      {
        // add to absolute
        Base::operator +=( VALUE);

        // return
        return *this;
      }

      //! @brief operator -= VALUE subtract y from y0
      //! @param VALUE y to be subtracted from y0
      SumFunction &
      operator -=
      (
        const t_ResultType &VALUE
      )
      {
        // subtract from absolute
        Base::operator -=( VALUE);

        // return
        return *this;
      }

      //! @brief operator *= SCALAR
      //! mutlipy absolute and each coefficient with the given scalar: s*f(x) = s*y0 + s*s1*g1(x) + s*s2*g2(x)...
      //! @param SCALAR scalar to be used to multiply every term
      //! @return reference to altered sum function
      SumFunction &
      operator *=
      (
        const double &SCALAR
      )
      {
        Base::operator *=( SCALAR);

        // return this
        return *this;
      }

      //! @brief operator /= SCALAR
      //! devide absolute and each coefficient with the given scalar: f(x)/s = y0/s + s1/s*g1(x) + s2/s*g2(x)...
      //! @param SCALAR scalar to be used to devide every term
      //! @return reference to altered sum function
      SumFunction &
      operator /=
      (
        const double &SCALAR
      )
      {
        Base::operator /=( SCALAR);

        return *this;
      }

    }; // template class SumFunction

  //////////
  // data //
  //////////

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> SumFunction< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new SumFunction< t_ArgumentType, t_ResultType>())
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_SUM_FUNCTION_H_
