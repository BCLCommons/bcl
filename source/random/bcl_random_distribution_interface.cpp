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
#include "random/bcl_random_distribution_interface.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"

// external includes - sorted alphabetically
#include <time.h>

namespace bcl
{
  namespace random
  {

  //////////
  // data //
  //////////

    //! @brief enum, that adds RandomSeed flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_RandomSeedFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag
      (
        DistributionInterface::GetFlagRandomSeed(),
        command::e_Random
      )
    );

    //! @brief the default seed
    //! @return the default seed
    uint64_t DistributionInterface::GetDefaultSeed()
    {
      return uint64_t( 5489);
    }

    //! @brief the default double range
    //! @return default double range [0,1)
    const math::Range< double> &DistributionInterface::GetDefaultDoubleRange()
    {
      static const math::Range< double> s_default_range( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
      return s_default_range;
    }

    //! @brief commandline flag to be used to set MessageLevel over the commandline
    util::ShPtr< command::FlagInterface> &DistributionInterface::GetFlagRandomSeed()
    {
      static util::ShPtr< command::FlagInterface> s_flag_random_seed
      (
        new command::FlagStatic
        (
          "random_seed",
          "adjust the random seed; if flag is used, system time is used as seed, if additional parameter is passed, the given number will be used, otherwise default will be used",
          command::Parameter
          (
            "seed",
            "random seed for the random number generator",
            command::ParameterCheckRanged< uint64_t>( 0, std::numeric_limits< uint64_t>::max()),
            util::Format()( DistributionInterface::GetDefaultSeed())
          ),
          &DistributionInterface::SetGlobalSeedFromCommandlineFlag
        )
      );

      return s_flag_random_seed;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set seed randomly
    //! @return seed that it was set too
    uint64_t DistributionInterface::RandomizeSeed()
    {
      // create time objec and retrieve time from system
      time_t t;
      time( &t);

      return SetSeed( uint64_t( t));
    }

    //! @brief set seed to the one given in the command line flag
    //! id the flag is given, without an actual numerival value, the system time is used as seed
    uint64_t DistributionInterface::SetSeedFromCommandlineFlag()
    {
      // if flag was set but no seed was given, use a random seed
      if( GetFlagRandomSeed()->GetFlag() && !GetFlagRandomSeed()->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // randomize the std seed and the bcl seed.  The std seed is used in certain std functions
        std::srand( time( 0));
        return RandomizeSeed();
      }
      // if flag was not use the default seed or when falg was set and seed was given, use that
      else
      {
        const size_t seed( GetFlagRandomSeed()->GetFirstParameter()->GetNumericalValue< uint64_t>());
        // randomize the std seed and the normal seed
        std::srand( size_t( seed));
        return SetSeed( seed);
      }
    }

    //! @brief set seed to the one given in the command line flag
    //! id the flag is given, without an actual numerical value, the system time is used as seed
    void DistributionInterface::SetGlobalSeedFromCommandlineFlag()
    {
      uint64_t seed( GetGlobalRandom().SetSeedFromCommandlineFlag());
      BCL_MessageCrt( "Seed was set to " + util::Format()( seed) + "!");
    }

  /////////////////////////////////////////////////
  // default ranges for generated random numbers //
  /////////////////////////////////////////////////

    //! @brief default range for random float
    //! @return range, in which Float() generates a random number
    const math::Range< float> &DistributionInterface::GetFloatRange()
    {
      static const math::Range< float> s_default_range( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
      return s_default_range;
    }

    //! @brief default range for random int
    //! @return range, in which Integer() generates a random number
    const math::Range< int> &DistributionInterface::GetIntegerRange()
    {
      static const math::Range< int> s_default_range( -std::numeric_limits< int>::max(), std::numeric_limits< int>::max());
      return s_default_range;
    }

    //! @brief default range for random unsigned int
    //! @return range, in which unsigned int generates a random number
    const math::Range< unsigned int> &DistributionInterface::GetUnsignedIntRange()
    {
      static const math::Range< unsigned int> s_default_range( 0, std::numeric_limits< unsigned int>::max());
      return s_default_range;
    }

    //! @brief default range for random unsigned long
    //! @return range, in which unsigned long generates a random number
    const math::Range< unsigned long> &DistributionInterface::GetUnsignedLongRange()
    {
      static const math::Range< unsigned long> s_default_range( 0, std::numeric_limits< unsigned long>::max());
      return s_default_range;
    }

    //! @brief default range for random unsigned long long
    //! @return range, in which unsigned long long generates a random number
    const math::Range< unsigned long long> &DistributionInterface::GetUnsignedLongLongRange()
    {
      static const math::Range< unsigned long long> s_default_range( 0, std::numeric_limits< unsigned long long>::max());
      return s_default_range;
    }

    //! @brief default range for random bool
    //! @return range, in which Double() generates a random number
    const math::Range< bool> &DistributionInterface::GetBooleanRange()
    {
      static const math::Range< bool> s_default_range( false, true);
      return s_default_range;
    }

    //! @brief default range for random char
    //! @return range, in which char generates a random number
    const math::Range< char> &DistributionInterface::GetCharRange()
    {
      static const math::Range< char> s_default_range
      (
        std::numeric_limits< char>::min(),
        std::numeric_limits< char>::max()
      );
      return s_default_range;
    }

    //! @brief templated function to default range double specialization
    //! @return default range for double
    template<>
    const math::Range< double> &DistributionInterface::GetRange< double>() const
    {
      return GetDoubleRange();
    }

    //! @brief templated function to default range float specialization
    //! @return default range for float
    template<>
    const math::Range< float> &DistributionInterface::GetRange< float>() const
    {
      return GetFloatRange();
    }

    //! @brief templated function to default range int specialization
    //! @return default range for int
    template<>
    const math::Range< int> &DistributionInterface::GetRange< int>() const
    {
      return GetIntegerRange();
    }

    //! @brief templated function to default range unsigned int specialization
    //! @return default range for unsigned int specialization
    template<>
    const math::Range< unsigned int> &DistributionInterface::GetRange< unsigned int>() const
    {
      return GetUnsignedIntRange();
    }

    //! @brief templated function to default range unsigned long specialization
    //! @return default range for unsigned long specialization
    template<>
    const math::Range< unsigned long> &DistributionInterface::GetRange< unsigned long>() const
    {
      return GetUnsignedLongRange();
    }

    //! @brief templated function to default range unsigned long long specialization
    //! @return default range for unsigned long long specialization
    template<>
    const math::Range< unsigned long long> &DistributionInterface::GetRange< unsigned long long>() const
    {
      return GetUnsignedLongLongRange();
    }

    //! @brief templated function to default range bool specialization
    //! @return default range for bool
    template<>
    const math::Range< bool> &DistributionInterface::GetRange< bool>() const
    {
      return GetBooleanRange();
    }

    //! @brief templated function to default range char specialization
    //! @return default range for char
    template<>
    const math::Range< char> &DistributionInterface::GetRange< char>() const
    {
      return GetCharRange();
    }

  /////////////////////////////////////
  // random numbers in default range //
  /////////////////////////////////////

    //! @brief random double
    //! @return random number in double range
    double DistributionInterface::Double() const
    {
      // get the range for doubles
      const math::Range< double> &range( GetDoubleRange());

      // offset to the left
      const double left_offset( double( range.GetLeftCondition()));

      // calculate the fraction of one double on the entire uint64_t interval
      const double fraction( 1.0 / ( double( std::numeric_limits< uint64_t>::max()) + left_offset + double( range.GetRightCondition())));

      // generate the random 64 bit int on the interval [0,64 bit int max]
      const double random_64_bit_uint( Unsigned64BitInt());

      // calculate the double
      const double random_double( ( random_64_bit_uint + left_offset) * fraction);

      return random_double;
    }

    //! @brief random float
    //! @return random number in range [0,1)
    float DistributionInterface::Float() const
    {
      return float( Double( GetFloatRange().Convert< double>()));
    }

  ///////////////////////////////////
  // random numbers in given range //
  ///////////////////////////////////

    //! @brief random double in given range
    //! @return random number in range RANGE
    double DistributionInterface::Double( const math::Range< double> &RANGE) const
    {
      return RANGE.Rescale( Double(), GetDoubleRange());
    }

    //! @brief random float in given range
    //! @return random number in range RANGE
    float DistributionInterface::Float( const math::Range< float> &RANGE) const
    {
      return float( RANGE.Convert< double>().Rescale( Double(), GetDoubleRange()));
    }

    //! @brief random integer in given range
    //! @return random integer in range RANGE
    int DistributionInterface::Integer( const math::Range< int> &RANGE) const
    {
      return int( RANGE.StandardizeRange().Convert< double>().Rescale( Double(), GetDoubleRange()));
    }

    //! @brief random unsigned integer in given range
    //! @return random unsigned bit integer in range RANGE
    unsigned int DistributionInterface::SizeT( const math::Range< unsigned int> &RANGE) const
    {
      return ( unsigned int)( RANGE.StandardizeRange().Convert< double>().Rescale( Double(), GetDoubleRange()));
    }

    //! @brief random unsigned long in given range
    //! @return random unsigned long in range RANGE
    unsigned long DistributionInterface::SizeT( const math::Range< unsigned long> &RANGE) const
    {
      return ( unsigned long)( RANGE.StandardizeRange().Convert< double>().Rescale( Double(), GetDoubleRange()));
    }

    //! @brief random unsigned long long in given range
    //! @return random unsigned long long in range RANGE
    unsigned long long DistributionInterface::SizeT( const math::Range< unsigned long long> &RANGE) const
    {
      return ( unsigned long long)( RANGE.StandardizeRange().Convert< double>().Rescale( Double(), GetDoubleRange()));
    }

    //! @brief generate a random boolean
    //! @return random true or false
    bool DistributionInterface::Boolean() const
    {
      // interpret highest-order bit as a bool
      // technically we could choose any bit, but the lowest order bits of numbers of several popular algorithms are
      // have a short period (e.g. 2 or 64), so choose the highest order bit instead
      return ( Unsigned64BitInt() & UINT64_C( 0x8000000000000000)) == 0;
    }

    //! @brief templated function to generate random number in given range specialized for double
    //! @return random double in RANGE
    template<>
    double DistributionInterface::Random< double>( const math::Range< double> &RANGE) const
    {
      return Double( RANGE);
    }

    //! @brief templated function to generate random number in given range specialized for float
    //! @return random float in RANGE
    template<>
    float DistributionInterface::Random< float>( const math::Range< float> &RANGE) const
    {
      return Float( RANGE);
    }

    //! @brief templated function to generate random number in given range specialized for int
    //! @return random int in RANGE
    template<>
    int DistributionInterface::Random< int>( const math::Range< int> &RANGE) const
    {
      return Integer( RANGE);
    }

    //! @brief templated function to generate random number in given range specialized for int
    //! @return random char in RANGE
    template<>
    char DistributionInterface::Random< char>( const math::Range< char> &RANGE) const
    {
      // note that we cannot use SizeT here because some implementations make normal char type signed, while others make
      // it unsigned
      return Integer( math::Range< int>( int( RANGE.GetMin()), int( RANGE.GetMax())));
    }

    //! @brief templated function to generate random number in default range specialized for unsigned ints
    //! @return random unsigned int in default range
    template<>
    unsigned int DistributionInterface::Random< unsigned int>( const math::Range< unsigned int> &RANGE) const
    {
      return SizeT( RANGE);
    }

    //! @brief templated function to generate random number in default range specialized for unsigned longs
    //! @return random unsigned long in default range
    template<>
    unsigned long DistributionInterface::Random< unsigned long>( const math::Range< unsigned long> &RANGE) const
    {
      return SizeT( RANGE);
    }

    //! @brief templated function to generate random number in default range specialized for unsigned long longs
    //! @return random unsigned long long in default range
    template<>
    unsigned long long DistributionInterface::Random< unsigned long long>
    (
      const math::Range< unsigned long long> &RANGE
    ) const
    {
      return SizeT( RANGE);
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for double
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    double DistributionInterface::Random( const double &MIN, const double &MAX) const
    {
      return Random< double>( math::Range< double>( MIN, MAX));
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for float
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    float DistributionInterface::Random( const float &MIN, const float &MAX) const
    {
      return Random< float>( math::Range< float>( MIN, MAX));
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for int
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    int DistributionInterface::Random( const int &MIN, const int &MAX) const
    {
      return Random< int>( math::Range< int>( MIN, MAX));
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for unsigned ints
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    unsigned int DistributionInterface::Random( const unsigned int &MIN, const unsigned int &MAX) const
    {
      return Random< unsigned int>( math::Range< unsigned int>( MIN, MAX));
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for unsigned longs
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    unsigned long DistributionInterface::Random( const unsigned long &MIN, const unsigned long &MAX) const
    {
      return Random< unsigned long>( math::Range< unsigned long>( MIN, MAX));
    }

    //! @brief templated function to generate random number if range [MIN, MAX] specialized for unsigned long longs
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right open range
    //! @return random number in range [MIN, MAX]
    template<>
    unsigned long long DistributionInterface::Random( const unsigned long long &MIN, const unsigned long long &MAX) const
    {
      return Random< unsigned long long>( math::Range< unsigned long long>( MIN, MAX));
    }

    //! @brief templated function to generate random char if range [MIN, MAX] specialized for int
    //! @param MIN minimum on left closed range
    //! @param MAX maximum on right closed range
    //! @return random char in range [MIN, MAX]
    template<>
    char DistributionInterface::Random( const char &MIN, const char &MAX) const
    {
      return Random< char>( math::Range< char>( MIN, MAX));
    }

    //! @brief generate a sign represented by an integer
    //! @return -1 or 1
    int DistributionInterface::Sign() const
    {
      return Boolean() ? 1 : -1;
    }

    //! @brief random number from a gaussian distribution
    //! uses the polar form of the Box-Mueller distribution
    //! http://www.taygeta.com/random/gaussian.html
    //! @param MEAN mean of the gaussian distribution
    //! @param STANDARD_DEVIATION standard deviation of the gaussian distribution
    //! @return a random number, which distribution is a gaussian distribution with mean MEAN and sd STANDARD_DEVIATON
    double DistributionInterface::RandomGaussian
    (
      const double MEAN,
      const double STANDARD_DEVIATION
    ) const
    {
      double x1, x2, w, y1;

      // use value from previous call
      if( m_OtherGaussianValueAvailable)
      {
        y1 = m_OtherGaussianValue;
        m_OtherGaussianValueAvailable = false;
      }
      else
      {
        do
        {
          x1 = double( 2.0) * Double() - double( 1.0);
          x2 = double( 2.0) * Double() - double( 1.0);
          w = x1 * x1 + x2 * x2;
        }
        while( w >= double( 1.0));

        w = sqrt( ( double( -2.0) * log( w)) / w);
        y1 = x1 * w;
        m_OtherGaussianValue = x2 * w;
        m_OtherGaussianValueAvailable = true;
      }

      return ( MEAN + y1 * STANDARD_DEVIATION);
    }

    //! @brief random discrete number from a poisson distribution
    //! This approach is from The transformed rejection method for generating Poisson random variables W. Hormann
    //! This is what is used in the boost library.  For lambda less than 10, a simple inversion techniqe is used
    //! @param LAMBDA expected Value of the poisson distribution
    //! @return a random number, which distribution is a poisson distribution with mean LAMBDA
    double DistributionInterface::RandomPoisson( const int LAMBDA) const
    {
      // This function is adapted from a version in Boost
      //   Distributed under the Boost Software License, Version 1.0.
      // Boost Software License - Version 1.0 - August 17th, 2003
      //
      // Permission is hereby granted, free of charge, to any person or organization
      // obtaining a copy of the software and accompanying documentation covered by
      // this license (the "Software") to use, reproduce, display, distribute,
      // execute, and transmit the Software, and to prepare derivative works of the
      // Software, and to permit third-parties to whom the Software is furnished to
      // do so, all subject to the following:
      //
      // The copyright notices in the Software and this entire statement, including
      // the above license grant, this restriction and the following disclaimer,
      // must be included in all copies of the Software, in whole or in part, and
      // all derivative works of the Software, unless such copies or derivative
      // works are solely in the form of machine-executable object code generated by
      // a source language processor.
      //
      // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
      // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
      // FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
      // SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
      // FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
      // ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
      // DEALINGS IN THE SOFTWARE.
      BCL_Assert( LAMBDA >= 0 && LAMBDA <= 100000000, "Lambda is outside of the usable range [ 0, 100,000,000]");

      // For Lambda less than 10 use this inversion algorithm from the boost library
      if( LAMBDA < 10)
      {
        double exp_mean( exp( -LAMBDA));
        double p( exp_mean);
        double x( 0);

        // generate a random number from a uniform distribution between 0 and 1
        double u( Double());

        // inverse transform method for low value of LAMBDA
        while( u > p)
        {
            u = u - p;
            ++x;
            p = LAMBDA * p / x;
        }
        return x;
      }
      else
      {
        // Prepare Table of log(k!), k =0,...,9
        static storage::Vector< double> poisson_table
        (
          storage::Vector< double>::Create
          (
            log( math::Factorial( 0)),
            log( math::Factorial( 1)),
            log( math::Factorial( 2)),
            log( math::Factorial( 3)),
            log( math::Factorial( 4)),
            log( math::Factorial( 5)),
            log( math::Factorial( 6)),
            log( math::Factorial( 7)),
            log( math::Factorial( 8)),
            log( math::Factorial( 9))
          )
        );

        // Initialize parameters
        double smu( math::Sqrt( LAMBDA));
        double   b( 0.931 + 2.53 * smu);
        double   a( -0.059 + 0.02483 * b);
        double inv_alpha( 1.1239 + 1.1328 / ( b - 3.4));
        double v_r( 0.9277 - 3.6224 / ( b - 2));

        // This will be executed until a return condition is met
        while( true)
        {
          // Step 1 of algorithm
          double u;
          double v( Double());

          if( v <= 0.86 * v_r)
          {
            u = v / v_r - 0.43;

            return std::floor( ( 2 * a / ( 0.5 - math::Absolute( u)) + b) * u + LAMBDA + 0.445);
          }
          // Step 2 of algorithm
          if( v >= v_r)
          {
            u = Double() - 0.5;
          }
          else
          {
            u = v / v_r - 0.93;
            u = ( ( u < 0) ? -0.5 : 0.5) - u;
            v = Double() * v_r;
          }
          // Step 3 of algorithm
          double us( 0.5 - math::Absolute( u));
          if( us < 0.013 && v > us)
          {
            continue;
          }
          double k( std::floor( ( 2 * a / us + b) * u + LAMBDA + 0.445));
          v = v * inv_alpha / ( a / ( us * us) + b);

          static double log_sqrt_2pi( log( math::Sqrt( 2 * math::g_Pi)));

          if( k >= 10)
          {
            if
            (
              log( v * smu) <=
              ( k + 0.5) * log( LAMBDA / k) -
              LAMBDA -
              log_sqrt_2pi +
              k -
              ( 1.0 / 12.0 - ( 1.0 / 360.0 - 1.0 / ( 1260.0 * k * k)) / ( k * k)) / k
            )
            {
              return k;
            }
          }
          else if( k >= 0)
          {
            if( log( v) <= k * log( LAMBDA) - LAMBDA - poisson_table( k))
            {
              return k;
            }
          }
        } // end while Loop
      } // end outer if Statement
      return 0.0;
    } // end poisson function
  } // namespace random
} // namespace bcl
