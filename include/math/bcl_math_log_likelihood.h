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

#ifndef BCL_MATH_LOG_LIKELIHOOD_H_
#define BCL_MATH_LOG_LIKELIHOOD_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_z_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LogLikelihood
    //! @brief Scores how well a given value is, compared to the mean and SD of a set of these values
    //! @details This class converts a score into its log likelihood based on a transformation of the provided score to its
    //! normalized value
    //!
    //! @see @link example_math_log_likelihood.cpp @endlink
    //! @author bitterd
    //! @date Jun 29, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LogLikelihood :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      //! mean and standard deviation of the distribution we score in, also used for standardization
      ZScore m_ZScore;

      //! an arbitary offset to the score, that shifts the confidence level, above which scores are returned negative
      //! in Tyka , DiMaio et.al. 2009 it was set to 0.5, which gives negative values for z-score > -1.0
      static const double s_ArbitaryScoreConfidenceOffset;

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
      LogLikelihood();

      //! @brief construct from mean and sd
      //! @param MEAN mean value of the considered distribution
      //! @param SIGMA standard deviation of the considered distribution
      LogLikelihood
      (
        const double MEAN,
        const double SIGMA
      );

      //! @brief Clone function
      //! @return pointer to new LogLikelihood
      LogLikelihood *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to Mean
      //! @return the mean used in the z-score
      double GetMean() const;

      //! @brief access to Standard Deviation
      //! @return the standard deviation used in the z-score
      double GetSigma() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates the score of the value in the normal distribution with the given mean and sigma
      //! @param ARGUMENT value which z-score shall be calculated from
      //! @return score of the given value
      double operator()( const double &ARGUMENT) const;

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

    }; // class LogLikelihood

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_LOG_LIKELIHOOD_H_
