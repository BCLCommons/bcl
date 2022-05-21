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

#ifndef BCL_RANDOM_HISTOGRAM_1D_DISTRIBUTION_H_
#define BCL_RANDOM_HISTOGRAM_1D_DISTRIBUTION_H_

// include the namespace header
#include "bcl_random.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_random_uniform_distribution.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Histogram1DDistribution
    //! @brief Class is for getting a random value from a distribution, where the random value is biased by the
    //!        probabilities observed in the distribution.
    //!
    //! @see @link example_random_histogram_1d_distribution.cpp @endlink
    //! @author alexanns
    //! @date Oct 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Histogram1DDistribution :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the distribution from which a random entry will be selected biased by the frequencies
      linal::Vector< double> m_Distribution1D;

      //! sum of m_Distribution1D
      double m_DistributionSum;

      //! the random number generator used to help select entry from the distribution
      const DistributionInterface &m_RandomNumberGenerator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Histogram1DDistribution();

      //! @brief constructor taking histogram and random number generator
      //! @param HISTOGRAM the histogram providing the probability distribution
      //! @param RNG the random number generator used to help select entry from the distribution
      Histogram1DDistribution( const math::Histogram &HISTOGRAM, const DistributionInterface &RNG = GetGlobalRandom());

      //! @brief Clone function
      //! @return pointer to new Histogram1DDistribution
      Histogram1DDistribution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetDistributionVector gives the vector containing the probability distribution
      //! @return m_Distribution1D which is the vector containing the probability distribution
      const linal::Vector< double> &GetDistributionVector() const;

      //! @brief Set the probability of a particular entry
      //! @param ENTRY entry to set the probability for
      //! @param PROB new probability
      //! @return true unless the distribution is now undefined - e.g. all events have a probability of 0
      bool SetProbabilityOfEvent( const size_t &EVENT, const double &PROB);

      //! @brief get the current sum of probability
      double GetDistributionSum() const
      {
        return m_DistributionSum;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Return a random object according to the supplied probabilities
      //! @param LOWER_BOUNDARY the lower boundary to be used
      //! @param BIN_SIZE the size of the bins
      //! @return double which is the actual bin label that has been selected
      const double DetermineRandomCase
      (
        const double &LOWER_BOUNDARY,
        const double &BIN_SIZE
      ) const;

      //! @brief DetermineRandomCase gives random position in probability distribution according to the probabilities
      //! @return size_t which is a random position in m_Distribution1D biased by the distribution probabilities
      const size_t DetermineRandomCase() const;

      //! @brief DetermineMostLikelyCase gives random most likely case
      //! @return size_t which is a random position in m_Distribution1D biased by the distribution probabilities
      const size_t DetermineMostLikelyCase() const;

    ///////////////
    // operators //
    ///////////////

      const Histogram1DDistribution operator =( const Histogram1DDistribution &HISTOGRAM);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class Histogram1DDistribution

  } // namespace random
} // namespace bcl

#endif // BCL_RANDOM_HISTOGRAM_1D_DISTRIBUTION_H_
