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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! returns X to the power of EXPONENT analogous to pow but for size_t
    // fastest possible x ^ n Russian peasant algorithm, Integer needs to be an integral type
    template<>
    size_t Pow< size_t>( const size_t &VALUE, const size_t &EXPONENT)
    {
      // initialize p to one
      size_t p( 1);

      // Create local variables to operate on the data
      size_t value( VALUE);
      size_t exponent( EXPONENT);

      if( p != value && exponent)
      {
        // if exponent is an odd number set p equal to value
        if( exponent & 1)
        {
          p = value;
        }
        // This loop computes the value
        while( exponent)
        {
          value *= value;
          exponent >>= 1;
          if( exponent & 1)
          {
            p *= value;
          }
        } // end while
      } // end if
      return p;
    }

    //! this function returns a weight according to an angle between between 0 .. pi :
    //! 0 .. pi/3 = 1; pi / 3 .. 2 * pi / 3: cosfunction from 1 to 0; 2 * pi / 3 .. 1 = 0
    double WeightBetweenZeroAndPi_ThreeSections( const double &ANGLE_RAD)
    {
      // lower boundary - return
      static const double lower_angle(     g_Pi / 3);
      static const double upper_angle( 2 * g_Pi / 3);

      //if smaller than lower angle return 1
      if( ANGLE_RAD <= lower_angle)
      {
        return double( 1);
      }

      //if larger than upper angle return 1
      if( ANGLE_RAD >= upper_angle)
      {
        return double( 0);
      }

      //otherwise apply sine transition
      return ( cos( ( ANGLE_RAD - lower_angle) * 3) + 1) / 2;
    }

    //! this function returns a weight according to an angle between between 0 .. pi : cosine function from 1 to 0;
    double WeightBetweenZeroAndPi( const double &ANGLE_RAD)
    {
      //otherwise apply sine transition
      return ( cos( ANGLE_RAD) + 1) / 2;
    }

    //! @brief compute the binomial coefficient n | k
    //! @param N size of set from which to choose
    //! @param K the number to choose
    //! @return the number of ways to choose K unique elements from a set of size N
    size_t BinomialCoefficient( const size_t &N, const size_t &K)
    {
      // handle the case where K > N-K
      if( K > N - K)
      {
        return BinomialCoefficient( N, N - K);
      }
      size_t combinations( 1);
      for( size_t numerator( N), denominator( 1); denominator <= K; ++denominator, --numerator)
      {
        combinations = combinations * numerator / denominator;
      }
      return combinations;
    }

    //! @brief compute the factorial
    //! @param N size of set from which to choose
    //! @return the number of permutations (orderings) of a set of size N = product( 1...N)
    uint64_t Factorial( const size_t &N)
    {
      // 21! is larger than 1 << 64, so only store 0-20!
      static const size_t max_factorial( 20);
      static const uint64_t s_factorials[ max_factorial + 1] =
      {
        UINT64_C( 1),
        UINT64_C( 1),
        UINT64_C( 2),
        UINT64_C( 6),
        UINT64_C( 24),
        UINT64_C( 120),
        UINT64_C( 720),
        UINT64_C( 5040),
        UINT64_C( 40320),
        UINT64_C( 362880),
        UINT64_C( 3628800),
        UINT64_C( 39916800),
        UINT64_C( 479001600),
        UINT64_C( 6227020800),
        UINT64_C( 87178291200),
        UINT64_C( 1307674368000),
        UINT64_C( 20922789888000),
        UINT64_C( 355687428096000),
        UINT64_C( 6402373705728000),
        UINT64_C( 121645100408832000),
        UINT64_C( 2432902008176640000)
      };
      if( N > max_factorial)
      {
        // return the max 64 bit integer
        return std::numeric_limits< uint64_t>::max();
      }
      return s_factorials[ N];
    }

    //! @brief error function
    //! http://en.wikipedia.org/wiki/Error_function
    //! http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf
    //! only accurate to better than 4 * 10**(-4) in relative error
    //! @param x argument to function
    //! @return the integral
    double Erf( const double x)
    {
      const double a( 0.140012288686666606); // approximation for the factor used in the calculation, further information see link above
      const double xsqr( Sqr( x)); // used frequently in the equation
      const double integral( Sqrt( 1.0 - exp( -xsqr * ( 4.0 / g_Pi + a * xsqr) / ( 1.0 + a * xsqr))));
      if( x >= 0.0)
      {
        return integral;
      }
      else//  if ( x < 0)
      {
        return -integral;
      }
    }

    //! @brief complimentary error function 1-erf(x)
    double Erfc( const double x)
    {
      return 1.0 - Erf( x);
    }

  } // namespace math
} // namespace bcl
