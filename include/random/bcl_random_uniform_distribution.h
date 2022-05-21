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

#ifndef BCL_RANDOM_UNIFORM_DISTRIBUTION_H_
#define BCL_RANDOM_UNIFORM_DISTRIBUTION_H_

// include the namespace header
#include "bcl_random.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.h"

// includes from bcl - sorted alphabetically
#include "bcl_random_distribution_interface.h"

// external includes - sorted alphabetically
#include <vector>

namespace bcl
{
  namespace random
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UniformDistribution
    //! @brief This class generates random long, double numbers
    //! @details the MT19937 variant of the mersenne twister algorithm is used
    //!
    //! @see @link example_random_uniform_distribution.cpp @endlink
    //! @author woetzen, mendenjl
    //! @date Mar 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API UniformDistribution :
      public DistributionInterface
    {

    private:

    //////////
    // data //
    //////////

      static const uint64_t s_N = 312;       //!< Period parameter N for 64-bit random number generation
      static const uint64_t s_M = 156;       //!< Period parameter M for 64-bit random number generation

      uint64_t          m_RandomSeed;        //!< global seed for random number generator
      mutable uint64_t  m_StateVector[ s_N]; //!< state vector for 64-bit Mersenne twister
      mutable uint64_t  m_Index;             //!< current position in m_StateVector64

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor initialized with DefaultSeed
      UniformDistribution();

      //! construct from seed
      UniformDistribution( const uint64_t SEED);

      //! copy constructor
      UniformDistribution( const UniformDistribution &UNIFORM_DISTRIBUTION);

      //! clone
      UniformDistribution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name if this function is overwritten
      const std::string &GetClassIdentifier() const;

      //! seeds to a defined state
      uint64_t SetSeed( const uint64_t SEED);

      //! return current seed
      uint64_t GetSeed() const;

      //! @brief default range for random double
      //! @return range, in which Double() generates a random number
      const math::Range< double> &GetDoubleRange() const
      {
        return GetDefaultDoubleRange();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief generate random unsigned 64-bit integer
      //! @return integer uniformly distributed in range [0,2^64-1]
      uint64_t Unsigned64BitInt() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      UniformDistribution &operator =( const UniformDistribution &UNIFORM_DISTRIBUTION);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! generate new words if end of pregenerated word list has been hit
      void GenerateNew64BitWords() const;

    }; // class UniformDistribution

    //! return a global RandomNumberGenerator
    BCL_API
    DistributionInterface &GetGlobalRandom();

    //! @brief helper function to avoid creating a bound member function pointer every time during Shuffle
    //! @param MAX the maximum value to return'
    BCL_API std::iterator_traits< std::vector< size_t>::const_iterator>::difference_type GetRandomSizeT
    (
      const std::iterator_traits< std::vector< size_t>::const_iterator>::difference_type &MAX
    );

    //! @brief helper function to randomly permute elements that is standard across compilers (unlike std::random_shuffle)
    template< typename t_DataType>
    void RandomShuffle( t_DataType FIRST, t_DataType LAST)
    {
      typename std::iterator_traits< t_DataType>::difference_type i, n;
      n = LAST - FIRST;
      for( i = n - 1; i > 0; --i)
      {
        std::swap( FIRST[ i], FIRST[ GetGlobalRandom().Random< size_t>( i)]);
      }
    }

  } // namespace random
} // namespace bcl

#endif // BCL_RANDOM_UNIFORM_DISTRIBUTION_H_

