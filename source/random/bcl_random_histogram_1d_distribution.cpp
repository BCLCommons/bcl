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
