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

#ifndef BCL_RANDOM_HISTOGRAM_2D_DISTRIBUTION_H_
#define BCL_RANDOM_HISTOGRAM_2D_DISTRIBUTION_H_

// include the namespace header
#include "bcl_random.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_histogram_2d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Histogram2DDistribution
    //! @brief class for determining random cases from a 2D probability distribution
    //! @details This class is used when statistics in form of a histogram2D (i.e amino acid phi/psi angles) needs to be used in
    //! biased sampling.
    //!
    //! @see @link example_random_histogram_2d_distribution.cpp @endlink
    //! @author rouvelgh, woetzen
    //! @date May 18, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Histogram2DDistribution :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! matrix which contains the probabilities
      linal::Matrix< double> m_ProbabilityMatrix;

    public:

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Histogram2DDistribution()
      {
      }

      //! @brief constructor which takes a pre-generated histogram2D
      //! @param HISTOGRAM is a Histogram2D object which contains the statistical information used for determining the random cases
      Histogram2DDistribution( const math::Histogram2D &HISTOGRAM);

      //! @brief virtual copy constructor
      Histogram2DDistribution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the matrix associated with the object
      const linal::Matrix< double> &GetMatrix() const
      {
        return m_ProbabilityMatrix;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Determines a random 2D object according to the supplied probabilities
      //! @return VectorND< 2, size_t> of the determined random case
      const storage::VectorND< 2, size_t> DetermineRandomCase2D() const;

      //! @brief Return a random object according to the supplied probabilities
      //! @param X_BOUNDARY the boundary to be used in the x-direction
      //! @param Y_BOUNDARY the boundary to be used in the y-direction
      //! @param X_BINSIZE the size of the bin in the x-direction
      //! @param Y_BINSIZE the size of the bin in the y-direction
      //! @return VectorND<2, double> containing the coordinates of the determined random case transformed into the
      //! @return proper units
      const storage::VectorND< 2, double> DetermineRandomCase2D
      (
        const double &X_BOUNDARY,
        const double &Y_BOUNDARY,
        const double &X_BINSIZE,
        const double &Y_BINSIZE
      ) const;

    private:

      //! @brief Return a random object according to the supplied probabilities
      const size_t DetermineRandomCase() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ProbabilityMatrix from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProbabilityMatrix to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Histogram2DDistribution

  } // namespace random
} // namespace bcl

#endif // BCL_RANDOM_HISTOGRAM_2D_DISTRIBUTION_H_
