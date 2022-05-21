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

