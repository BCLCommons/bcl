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

#ifndef BCL_MATH_H_
#define BCL_MATH_H_

// include the namespace forward header
#include "bcl_math.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "type/bcl_type_chooser.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_sequence.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_math.h
  //! @brief namespace for math classes and functions in the biochemistry library
  //! @details This class contains the constants: Pi, Euler, Permeability of Free Space ( magnetic constant), Boltzmann,
  //! @details and the plank constant.  It also defines basic math operations on data
  //! @see @link example_math.cpp @endlink
  //! @author mueller
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace math
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    const double g_Pi(        3.14159265359); //!< The constant Pi [1]
    const double g_Euler(     2.71828182846); //!< The constant  e [1]
    const double g_My(        1.256637e-6  ); //!< Permeability of Free Space my [T1^2 * m^3 * J^-1]   = 4 * PI * 10e-7
    const double g_Boltzmann( 1.380652e-23 ); //!< Boltzmann constant b [J * K^-1]
    const double g_Plank(     6.626176e-34 ); //!< Planks constant h [J * s]
    const double g_Gas(       8.31446261815324); //!< Gas constant R [J * K^-1 * mol^-1]

    //! defines "absolute" for int
    inline
    int Absolute( const int A)
    {
#ifdef ANDROID
      return std::fabs( A);
#else
      return std::abs( A);
#endif
    }

    //! defines "absolute" for size_t
    inline
    size_t
    Absolute( const size_t A)
    {
      return( A);
    }

    //! defines "absolute" for double
    inline
    double
    Absolute( const double A)
    {
      return std::abs( A);
    }

    //! defines "absolute" for complex< double>
    inline
    std::complex< double>
    Absolute( const std::complex< double> &A)
    {
      return std::abs( A);
    }

    //! set all elements to their absolute numbers
    template< typename t_Container>
    inline typename type::EnableIf< type::IsSequence< t_Container>::value, void>::Type
    Absolute( t_Container &CONTAINER)
    {
      //absolute each element
      for( typename t_Container::iterator itr( CONTAINER.Begin()), itr_end( CONTAINER.End()); itr != itr_end; ++itr)
      {
        *itr = std::abs( *itr);
      }
    }

    //! defines "absolute" for std::pair< t_DataType1, t_DataType2>
    template< typename t_DataType1, typename t_DataType2>
    inline std::pair< t_DataType1, t_DataType2> Absolute( const std::pair< t_DataType1, t_DataType2> &A)
    {
      return std::pair< t_DataType1, t_DataType2>( Absolute( A.first), Absolute( A.second));
    }

    //! shortcut to A*A analogous to sqrt
    template< typename t_DataType>
    inline t_DataType Sqr( const t_DataType &A)
    {
      return A * A;
    }

    //! returns sqrt of t_DataType
    //! This version is for any datatype other than any non floating point data type.  This casts the data type as a
    //! double prior to taking the square root
    template< typename t_DataType>
    inline
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, double>::Type Sqrt( const t_DataType &X)
    {
      // convert argument to double, so that only one of the alternative for std::sqrt (double, long double and float)
      // is called
      return std::sqrt( double( X));
    }

    //! This version is for floating data point type.
    template< typename t_DataType>
    inline
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type Sqrt( const t_DataType &X)
    {
      return std::sqrt( X);
    }

    template< typename t_DataType>
    inline t_DataType Pow( const t_DataType &VALUE, const t_DataType &EXPONENT)
    {
      return std::pow( VALUE, EXPONENT);
    }

    //! returns X to the power of EXPONENT analogous to pow but for size_t
    template<>
    BCL_API
    size_t Pow< size_t>( const size_t &VALUE, const size_t &EXPONENT);

    //! copies sign of B to A for signed case
    template< typename t_DataType, typename t_OtherDataType>
    inline
    typename type::EnableIf< std::numeric_limits< t_OtherDataType>::is_signed, t_DataType>::Type
      Sign( const t_DataType &A, const t_OtherDataType &B)
    {
      return ( B) >= t_OtherDataType( 0) ? Absolute( A) : -Absolute( A);
    }

    //! copies sign of B to A for unsigned case
    template< typename t_DataType, typename t_OtherDataType>
    inline
    typename type::EnableIf< !std::numeric_limits< t_OtherDataType>::is_signed, t_DataType>::Type
      Sign( const t_DataType &A, const t_OtherDataType &B)
    {
      return Absolute( A);
    }

    //! computes sqrt(A*A+B*B) in a numerically stable manner
    template< typename t_DataType>
    inline typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type
    Pythag( const t_DataType &A, const t_DataType &B)
    {
      // numerically stable version
      const t_DataType abs_a( Absolute( A)), abs_b( Absolute( B));
      if( abs_a > abs_b)
      {
        // A > B, compute A * Sqrt( 1 + (B/A)^2)
        return abs_a * Sqrt( t_DataType( 1) + Sqr( abs_b / abs_a));
      }
      else if( abs_a < abs_b)
      {
        // A < B, compute B * Sqrt( 1 + (A/B)^2)
        return abs_b * Sqrt( t_DataType( 1) + Sqr( abs_a / abs_b));
      }
      // A == B, compute A * Sqrt(2)
      static const double sqrt_two( std::sqrt( 2.0));
      return abs_a * sqrt_two;
    }

    //! computes sqrt(A*A+B*B) in a numerically stable manner
    template< typename t_DataType>
    inline typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type
    Pythag( const t_DataType &A, const t_DataType &B)
    {
      // numerically stable version
      const double abs_a( Absolute( A)), abs_b( Absolute( B));
      if( abs_a > abs_b)
      {
        // A > B, compute A * Sqrt( 1 + (B/A)^2)
        return abs_a * Sqrt( 1.0 + Sqr( abs_b / abs_a));
      }
      else if( abs_a < abs_b)
      {
        // A < B, compute B * Sqrt( 1 + (A/B)^2)
        return abs_b * Sqrt( 1.0 + Sqr( abs_a / abs_b));
      }
      // A == B, compute A * Sqrt(2)
      static const double sqrt_two( std::sqrt( 2.0));
      return abs_a * sqrt_two;
    }

    //! defines "<" for complex< t_DataType>
    template< typename t_DataType>
    inline
    bool
    operator <( const std::complex< t_DataType> &A, const std::complex< t_DataType> &B)
    {
      return std::abs( A) < std::abs( B);
    }

    //! defines "<=" for complex< t_DataType>
    template< typename t_DataType>
    inline
    bool
    operator <=( const std::complex< t_DataType> &A, const std::complex< t_DataType> &B)
    {
      return std::abs( A) <= std::abs( B);
    }

    //! defines ">" for complex< t_DataType>
    template< typename t_DataType>
    inline
    bool
    operator >( const std::complex< t_DataType> &A, const std::complex< t_DataType> &B)
    {
      return std::abs( A) > std::abs( B);
    }

    //! defines ">=" for complex< t_DataType>
    template< typename t_DataType>
    inline
    bool
    operator >=( const std::complex< t_DataType> &A, const std::complex< t_DataType> &B)
    {
      return std::abs( A) >= std::abs( B);
    }

    //! test if argument is power of 2
    inline
    bool
    IsPowerOfTwo( const size_t &ARGUMENT)
    {
      // http://en.wikipedia.org/wiki/Power_of_two#Fast_algorithm_to_check_if_a_positive_number_is_a_power_of_two
      return ( ARGUMENT != 0) && ( ( ARGUMENT & ( ARGUMENT - 1)) == 0);
    }

    //! test if ARGUMENT_A is equal within absolute Tolerance
    inline
    bool
    EqualWithinAbsoluteTolerance
    (
      const double &TARGET_VALUE, const double &TEST_VALUE, const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon()
    )
    {
      //check that difference between Arguments is smaller than tolerance
      return Absolute( TARGET_VALUE - TEST_VALUE) <= ( ABSOLUTE_TOLERANCE);
    }

    //! test if ARGUMENT_A is almost equal ARGUMENT_B (default TOLERANCE is one per mill of larger ARGUMENT)
    template< typename t_DataType, typename t_DataTypeB>
    inline typename type::EnableIf< !type::IsSequence< t_DataType>::value, bool>::Type
    EqualWithinTolerance
    (
      const t_DataType &TARGET_VALUE, const t_DataTypeB &TEST_VALUE, const double &RELATIVE_TOLERANCE = 0.001,
      const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon()
    )
    {
      // check whether absolute difference is smaller TRAGET_VALUE * RELATIVE_TOLERANCE + ABSOLUTE_TOLERANCE
      return EqualWithinAbsoluteTolerance( TARGET_VALUE, TEST_VALUE, Absolute( TARGET_VALUE * RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE);
    }

    template< typename t_ContainerType, typename t_ContainerTypeB>
    inline typename type::EnableIf< type::IsSequence< t_ContainerType>::value, bool>::Type
    EqualWithinTolerance
    (
      const t_ContainerType &TARGET_VECTOR,
      const t_ContainerTypeB &TEST_VECTOR,
      const double &RELATIVE_TOLERANCE = 0.001,
      const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon()
    )
    {
      // GetSize is not defined for certain types, like matrices, so use std::distance to compute the size instead
      // On vector types, std::distance(x,y) is O(1) (e.g. x-y)
      typename t_ContainerTypeB::const_iterator itr_b( TEST_VECTOR.Begin()), itr_b_end( TEST_VECTOR.End());
      typename t_ContainerType::const_iterator itr_a( TARGET_VECTOR.Begin()), itr_a_end( TARGET_VECTOR.End());
      if( std::distance( itr_b, itr_b_end) != std::distance( itr_a, itr_a_end))
      {
        return false;
      }
      // loop over each element in the containers to test that they are equal within tolerance
      for( ; itr_a != itr_a_end; ++itr_a, ++itr_b)
      {
        // check if those values are within tolerance
        if( !EqualWithinTolerance( *itr_a, *itr_b, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE))
        {
          // if not, return false
          return false;
        }
      }

      // end
      return true;
    }

    //! test if ARGUMENT_A is almost equal ARGUMENT_B within machine tolerance
    template< typename t_DataType, typename t_DataTypeB>
    inline typename type::EnableIf< !type::IsSequence< t_DataType>::value && !std::numeric_limits< t_DataType>::is_integer, bool>::Type
    EqualWithinMachineTolerance
    (
      const t_DataType &TARGET_VALUE, const t_DataTypeB &TEST_VALUE
    )
    {
      // determine which of the two types is smaller and use it's epsilon; unless type B is an integer, then use Type A's
      // epsilon
      typedef typename type::Chooser
      <
        sizeof( t_DataType) <= sizeof( t_DataTypeB) || std::numeric_limits< t_DataTypeB>::is_integer,
        t_DataType,
        t_DataTypeB
      >::Type t_Smaller;
      // check whether absolute difference is smaller TRAGET_VALUE * RELATIVE_TOLERANCE + ABSOLUTE_TOLERANCE
      // always check first for exact equality, since this common case is much faster to test than computing differences
      return
        TARGET_VALUE == TEST_VALUE
        || EqualWithinAbsoluteTolerance( TARGET_VALUE, TEST_VALUE, Absolute( TARGET_VALUE * std::numeric_limits< t_Smaller>::epsilon()) + std::numeric_limits< t_Smaller>::min());
    }

    //! test if ARGUMENT_A is almost equal ARGUMENT_B (default TOLERANCE is one per mill of larger ARGUMENT)
    template< typename t_DataType, typename t_DataTypeB>
    inline typename type::EnableIf< !type::IsSequence< t_DataType>::value && std::numeric_limits< t_DataType>::is_integer, bool>::Type
    EqualWithinMachineTolerance
    (
      const t_DataType &TARGET_VALUE, const t_DataTypeB &TEST_VALUE
    )
    {
      // integral types have 0 epsilon, so expect exact values
      return TARGET_VALUE == TEST_VALUE;
    }

    template< typename t_ContainerType, typename t_ContainerTypeB>
    inline typename type::EnableIf< type::IsSequence< t_ContainerType>::value, bool>::Type
    EqualWithinMachineTolerance
    (
      const t_ContainerType &TARGET_VECTOR,
      const t_ContainerTypeB &TEST_VECTOR
    )
    {
      // GetSize is not defined for certain types, like matrices, so use std::distance to compute the size instead
      // On vector types, std::distance(x,y) is O(1) (e.g. x-y)
      typename t_ContainerTypeB::const_iterator itr_b( TEST_VECTOR.Begin()), itr_b_end( TEST_VECTOR.End());
      typename t_ContainerType::const_iterator itr_a( TARGET_VECTOR.Begin()), itr_a_end( TARGET_VECTOR.End());
      if( std::distance( itr_b, itr_b_end) != std::distance( itr_a, itr_a_end))
      {
        return false;
      }
      // loop over each element in the containers to test that they are equal within tolerance
      for( ; itr_a != itr_a_end; ++itr_a, ++itr_b)
      {
        // check if those values are within tolerance
        if( !EqualWithinMachineTolerance( *itr_a, *itr_b))
        {
          // if not, return false
          return false;
        }
      }

      // end
      return true;
    }

    //! this function returns a weight according to an angle between between 0 .. pi :
    //! 0 .. pi/3 = 1; pi / 3 .. 2 * pi / 3: cosine function from 1 to 0; 2 * pi / 3 .. 1 = 0
    BCL_API double WeightBetweenZeroAndPi_ThreeSections( const double &ANGLE_RAD);

    //! this function returns a weight according to an angle between between 0 .. pi : cosfunction from 1 to 0;
    BCL_API double WeightBetweenZeroAndPi( const double &ANGLE_RAD);

    //! returns VALUE if larger than LIMIT, otherwise 0
    template< typename t_DataType>
    t_DataType FilterValuesSmallerEqualLimit( const t_DataType &VALUE, const t_DataType &LIMIT)
    {
      return ( ( VALUE > LIMIT) ? VALUE : t_DataType( 0));
    }

    //! @brief converts the given boolean to a sign ( false -> negative, true -> positive)
    //! @param BOOLEAN boolean value to be converted to a sign
    //! @return sign converted from boolean
    inline
    int ConvertBooleanToSign( const bool BOOLEAN)
    {
      return BOOLEAN ? 1 : -1;
    }

    //! @brief compute the binomial coefficient n | k
    //! @param N size of set from which to choose
    //! @param K the number to choose
    //! @return the number of ways to choose K unique elements from a set of size N
    BCL_API size_t BinomialCoefficient( const size_t &N, const size_t &K);

    //! @brief compute the factorial
    //! @param N size of set from which to choose
    //! @return the number of permutations (orderings) of a set of size N = product( 1...N)
    BCL_API uint64_t Factorial( const size_t &N);

    //! @brief relative error in error function - calculate 12 significant figures
    //! don't ask for > 15 figures (assuming usual 52 bit mantissa in a double)
    const double g_ErfRelativeError( 1e-12);

    //! @brief error function
    //! http://en.wikipedia.org/wiki/Error_function
    //! erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
    //!        = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
    //!        = 1-erfc(x)
    //! @param x argument to function
    //! @return the integral
    BCL_API double Erf( const double x);

    //! @brief complimanetary error function 1-erf(x)
    //! erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
    //!         = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
    //!         = 1-erf(x)
    //! expression inside [] is a continued fraction so '+' means add to denominator only
    BCL_API double Erfc( const double x);

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_H_
