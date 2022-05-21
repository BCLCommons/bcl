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
#include "random/bcl_random.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace random
} // namespace bcl
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
#include "random/bcl_random_histogram_1d_distribution.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Histogram1DDistribution::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram1DDistribution())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Histogram1DDistribution::Histogram1DDistribution() :
      m_Distribution1D(),
      m_DistributionSum( 0.0),
      m_RandomNumberGenerator( GetGlobalRandom())
    {
    }

    //! @brief constructor taking histogram and random number generator
    //! @param HISTOGRAM the histogram providing the probability distribution
    //! @param RNG the random number generator used to help select entry from the distribution
    Histogram1DDistribution::Histogram1DDistribution
    (
      const math::Histogram &HISTOGRAM, const DistributionInterface &RNG
    ) :
      m_Distribution1D( HISTOGRAM.GetHistogram()),
      m_DistributionSum( m_Distribution1D.Sum()),
      m_RandomNumberGenerator( RNG)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Histogram1DDistribution
    Histogram1DDistribution *Histogram1DDistribution::Clone() const
    {
      return new Histogram1DDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram1DDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetDistributionVector gives the vector containing the probability distribution
    //! @return m_Distribution1D which is the vector containing the probability distribution
    const linal::Vector< double> &Histogram1DDistribution::GetDistributionVector() const
    {
      return m_Distribution1D;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Return a random object according to the supplied probabilities
    //! @param LOWER_BOUNDARY the lower boundary to be used
    //! @param BIN_SIZE the size of the bins
    //! @return double which is the actual bin label that has been selected
    const double Histogram1DDistribution::DetermineRandomCase
    (
      const double &LOWER_BOUNDARY,
      const double &BIN_SIZE
    ) const
    {
      // get a random index of m_Distribution1D
      const size_t random_case( DetermineRandomCase());

      // get the name of the bin that "random_case" denotes
      const double bin( LOWER_BOUNDARY + ( random_case + 0.5) * BIN_SIZE);

      // return the bin that was determined
      return bin;
    }

    //! @brief DetermineRandomCase gives random position in probability distribution according to the probabilities
    //! @return size_t which is a random position in m_Distribution1D biased by the distribution probabilities
    const size_t Histogram1DDistribution::DetermineRandomCase() const
    {
      // initiate a random number from 0 to 1.0
      double temp( m_RandomNumberGenerator.Random< double>( m_DistributionSum));

      // initialize the position holder
      size_t pos( 0);

      for
      (
        // initialize the pointers inside the distribution
        const double *ptr( m_Distribution1D.Begin()), *ptr_end( m_Distribution1D.End());
        ptr != ptr_end;
        ++ptr, ++pos
      )
      {
        // subtract the random number generated from the value pointed to in m_Distribution1D
        temp -= *ptr;

        // if the difference is negative or zero, then return the current position in the vector
        // because the larger the probability the greater the chance of causing this difference to be negative
        if( temp <= 0.0)
        {
          return pos;
        }
      }

      // should not reach this point
      BCL_Exit( "The probabilities are not set correctly or distribution is empty!!! " + util::Format()( temp) + " " + util::Format()( m_Distribution1D), -1);
      return pos;
    }

    //! @brief DetermineMostLikelyCase gives random most likely case
    //! @return size_t which is a random position in m_Distribution1D biased by the distribution probabilities
    const size_t Histogram1DDistribution::DetermineMostLikelyCase() const
    {
      storage::Vector< size_t> max_indices;
      double max_val( 0.0);
      // initialize the position holder
      size_t pos( 0);

      for
      (
        // initialize the pointers inside the distribution
        const double *ptr( m_Distribution1D.Begin()), *ptr_end( m_Distribution1D.End());
        ptr != ptr_end;
        ++ptr, ++pos
      )
      {
        if( *ptr > max_val)
        {
          max_val = *ptr;
          max_indices.Reset();
          max_indices.PushBack( pos);
        }
        else if( *ptr >= max_val * 0.95)
        {
          max_indices.PushBack( pos);
        }
      }
      return max_indices.GetSize() > size_t( 1)
             ? max_indices( GetRandomSizeT( max_indices.GetSize()))
             : max_indices( 0);
    }

    //! @brief Set the probability of a particular entry
    //! @param ENTRY entry to set the probability for
    //! @param PROB new probability
    //! @return true unless the distribution is now undefined - e.g. all events have a probability of 0
    bool Histogram1DDistribution::SetProbabilityOfEvent( const size_t &EVENT, const double &PROB)
    {
      BCL_Assert( PROB >= 0.0, "Negative probabilities are not allowed");
      m_DistributionSum += ( PROB - m_Distribution1D( EVENT));
      m_Distribution1D( EVENT) = PROB;
      return m_DistributionSum > 0.0;
    }

  ///////////////
  // operators //
  ///////////////

    const Histogram1DDistribution Histogram1DDistribution::operator =( const Histogram1DDistribution &HISTOGRAM)
    {
      m_Distribution1D = HISTOGRAM.m_Distribution1D;
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Histogram1DDistribution::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Distribution1D, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Histogram1DDistribution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Distribution1D, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace random

} // namespace bcl
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
#include "random/bcl_random_histogram_2d_distribution.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram2DDistribution::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram2DDistribution())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor which takes a pre-generated histogram2D
    //! @param HISTOGRAM a Histogram2D object which contains the statistics
    Histogram2DDistribution::Histogram2DDistribution( const math::Histogram2D &HISTOGRAM) :
      m_ProbabilityMatrix( HISTOGRAM.GetHistogram())
    {
      m_ProbabilityMatrix.AsVector().SetToSum( 1.0);
    }

    //! @brief virtual copy constructor
    Histogram2DDistribution *Histogram2DDistribution::Clone() const
    {
      return new Histogram2DDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @brief the class name as const &std::string
    const std::string &Histogram2DDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief  determines a biased random case from a probability distribution
    //! @return a randomly drawn 1D position from the probability distribution
    const size_t Histogram2DDistribution::DetermineRandomCase() const
    {
      // initiate a random number from 0 to 1.0
      double temp( GetGlobalRandom().Random< double>( 1.0));

      // initialize the position holder
      size_t pos( 0);

      for
      (
        // initialize the pointers inside the distribution
        const double *ptr( m_ProbabilityMatrix.Begin()), *ptr_end( m_ProbabilityMatrix.End());
        ptr != ptr_end;
        ++ptr, ++pos
      )
      {
        // subtract the random number generated from the value pointed to in m_ProbabilityMatrix
        temp -= *ptr;

        // if the difference is negative or zero, then return the current position in the matrix
        // because the larger the probability the greater the chance of causing this difference to be negative
        if( temp <= 0)
        {
          return pos;
        }
      }
      BCL_Exit( "The probabilities are not set correctly or distribution is empty!!!", -1);
      return pos;
    }

    //! @brief Determines a random 2D object according to the supplied probabilities
    //! @return VectorND< 2, size_t> of the determined random case
    const storage::VectorND< 2, size_t> Histogram2DDistribution::DetermineRandomCase2D() const
    {
      // initialize the 1D case to return a position
      const size_t position( DetermineRandomCase());

      // calculate the 2D coordinates given the 1D position
      const size_t x( position % m_ProbabilityMatrix.GetNumberCols());
      const size_t y( position / m_ProbabilityMatrix.GetNumberRows());

      return storage::VectorND< 2, size_t>( x, y);
    }

    //! @brief Return a random object according to the supplied probabilities
    //! @param X_BOUNDARY the first boundary in the x-direction
    //! @param Y_BOUNDARY the first boundary in the y-direction
    //! @param X_BINSIZE the size of the bin in the x-direction
    //! @param Y_BINSIZE the size of the bin in the y-direction
    //! @return VectorND<2, double> containing the coordinates of the determined random case transformed into the
    //! @return proper units
    const storage::VectorND< 2, double> Histogram2DDistribution::DetermineRandomCase2D
    (
      const double &X_BOUNDARY,
      const double &Y_BOUNDARY,
      const double &X_BINSIZE,
      const double &Y_BINSIZE
    ) const
    {
      // invoke the DetermineRandomCase2D to be used for the conversion
      const storage::VectorND< 2, size_t> position_2d( DetermineRandomCase2D());

      // given the 2D row/column position, calculate the "numerical" position
      const double x( X_BOUNDARY + ( position_2d.First() + 0.5) * X_BINSIZE);
      const double y( Y_BOUNDARY + ( position_2d.Second() + 0.5) * Y_BINSIZE);

      return storage::VectorND< 2, double>( x, y);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read the ProbabilityMatrix from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Histogram2DDistribution::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ProbabilityMatrix, ISTREAM);

      //end
      return ISTREAM;
      }

    //! @brief write ProbabilityMatrix to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Histogram2DDistribution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ProbabilityMatrix, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  } // namespace random
} // namespace bcl
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
#include "random/bcl_random_uniform_distribution.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

// disable vs warning c4351 new behavior of array elements being default initialized
#if defined (_MSC_VER)
  #pragma warning( disable : 4351)
#endif

namespace bcl
{
  namespace random
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor initialized with s_DefaultSeed
    UniformDistribution::UniformDistribution() :
      m_Index( 0)
    {
      SetSeedFromCommandlineFlag();
    }

    //! construct from seed
    UniformDistribution::UniformDistribution( const uint64_t SEED) :
      m_Index()
    {
      SetSeed( SEED);
    }

    //! copy constructor
    UniformDistribution::UniformDistribution( const UniformDistribution &UNIFORM_DISTRIBUTION) :
      DistributionInterface( UNIFORM_DISTRIBUTION),
      m_Index()
    {
      *this = UNIFORM_DISTRIBUTION;
    }

    //! clone
    UniformDistribution *UniformDistribution::Clone() const
    {
      return new UniformDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name if this function is overwritten
    const std::string &UniformDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return current seed
    uint64_t UniformDistribution::GetSeed() const
    {
      return m_RandomSeed;
    }

    //! @brief set the seed and initialize the first few values
    uint64_t UniformDistribution::SetSeed( const uint64_t SEED)
    {
      DistributionInterface::SetSeed( SEED);
      m_RandomSeed = SEED;
      m_StateVector[ 0] = SEED;
      // generate an initial set of random numbers using a linear congruential generator
      for( m_Index = 1; m_Index < s_N; ++m_Index)
      {
        m_StateVector[ m_Index] = UINT64_C( 6364136223846793005) * ( m_StateVector[ m_Index - 1] ^ ( m_StateVector[ m_Index - 1] >> 62)) + m_Index;
      }
      return SEED;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate random unsigned 64-bit integer
    //! @return integer uniformly distributed in range [0,2^64-1]
    uint64_t UniformDistribution::Unsigned64BitInt() const
    {
      uint64_t y;

      if( m_Index >= s_N) // generate s_N words at one time
      {
        GenerateNew64BitWords();
      }

      y = m_StateVector[ m_Index++];

      // Tempering
      y ^= ( y >> 29) & UINT64_C( 0x5555555555555555);
      y ^= ( y << 17) & UINT64_C( 0x71D67FFFEDA60000);
      y ^= ( y << 37) & UINT64_C( 0xFFF7EEE000000000);
      y ^= y >> 43;

      return y;
    }

    //! generate new words if end of pregenerated word list has been hit
    void UniformDistribution::GenerateNew64BitWords() const
    {
      // this is a constant array of two numbers used by the mersenne twister algorithm
      static const uint64_t mag01[ 2] = { UINT64_C( 0), UINT64_C( 0xB5026F5AA96619E9)};
      // Most significant 33 bits
      static const uint64_t s_upper_mask( UINT64_C( 0xFFFFFFFF80000000));
      // Least significant 31 bits
      static const uint64_t s_lower_mask( UINT64_C( 0x7FFFFFFF));

      uint64_t i, x;

      for( i = 0; i < s_N - s_M; ++i)
      {
        x = ( m_StateVector[ i] & s_upper_mask) | ( m_StateVector[ i + 1] & s_lower_mask);
        m_StateVector[ i] = m_StateVector[ i + s_M] ^ ( x >> 1) ^ mag01[ x & UINT64_C( 1)];
      }

      for( ; i < s_N - 1; ++i)
      {
        x = ( m_StateVector[ i] & s_upper_mask) | ( m_StateVector[ i + 1] & s_lower_mask);
        m_StateVector[ i] = m_StateVector[ i + ( s_M - s_N)] ^ ( x >> 1) ^ mag01[ x & UINT64_C( 1)];
      }
      x = ( m_StateVector[ s_N - 1] & s_upper_mask) | ( m_StateVector[ 0] & s_lower_mask);
      m_StateVector[ s_N - 1] = m_StateVector[ s_M - 1] ^ ( x >> 1) ^ mag01[ x & UINT64_C( 1)];

      m_Index = 0;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief equal operator
    UniformDistribution &UniformDistribution::operator =( const UniformDistribution &UNIFORM_DISTRIBUTION)
    {
      // assign members
      m_RandomSeed = UNIFORM_DISTRIBUTION.m_RandomSeed;
      std::copy( UNIFORM_DISTRIBUTION.m_StateVector, UNIFORM_DISTRIBUTION.m_StateVector + s_N, m_StateVector);
      m_Index      = UNIFORM_DISTRIBUTION.m_Index;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &UniformDistribution::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_RandomSeed, ISTREAM);
      for( uint64_t *ptr( m_StateVector), *ptr_end( m_StateVector + s_N); ptr != ptr_end; ++ptr)
      {
        io::Serialize::Read( *ptr, ISTREAM);
      }
      io::Serialize::Read( m_Index, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &UniformDistribution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_RandomSeed, OSTREAM, INDENT) << '\n';
      // first element
      io::Serialize::Write( *m_StateVector, OSTREAM, INDENT);
      for( const uint64_t *ptr( m_StateVector + 1), *ptr_end( m_StateVector + s_N); ptr != ptr_end; ++ptr)
      {
        OSTREAM << '\t';
        io::Serialize::Write( *ptr, OSTREAM, INDENT);
      }
      OSTREAM << '\n';

      io::Serialize::Write( m_Index, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! return a global RandomNumberGenerator
    DistributionInterface &GetGlobalRandom()
    {
      // static random object
      static DistributionInterface *s_random( new UniformDistribution());

      // return random object
      return *s_random;
    }

    //! @brief helper function to avoid creating a bound member function pointer every time during Shuffle
    //! @param MAX the maximum value to return'
    std::iterator_traits< storage::Vector< size_t>::const_iterator>::difference_type GetRandomSizeT
    (
      const std::iterator_traits< storage::Vector< size_t>::const_iterator>::difference_type &MAX
    )
    {
      return GetGlobalRandom().Random< size_t>( size_t( MAX - 1));
    }

  } // namespace random
} // namespace bcl

