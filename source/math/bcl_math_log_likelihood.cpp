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
#include "math/bcl_math_log_likelihood.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LogLikelihood::s_Instance
    (
      GetObjectInstances().AddInstance( new LogLikelihood())
    );

    //! an arbitary offset to the score, that shifts the confidence level, above which scores are returned negative
    //! in Tyka , DiMaio et.al. 2009 it was set to 0.5, which gives negative values for z-score > -1.0
    const double LogLikelihood::s_ArbitaryScoreConfidenceOffset( 0.5);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LogLikelihood::LogLikelihood() :
      m_ZScore()
    {
    }

    //! @brief construct with Mean and Standard Deviation
    //! @param MEAN mean value of the considered distribution
    //! @param SIGMA standard deviation of the considered distribution
    LogLikelihood::LogLikelihood
    (
      const double MEAN,
      const double SIGMA
    ) :
      m_ZScore( MEAN, SIGMA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LogLikelihood
    LogLikelihood *LogLikelihood::Clone() const
    {
      return new LogLikelihood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LogLikelihood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &LogLikelihood::GetScheme() const
    {
      static const std::string s_scheme( "log(0.5*(1-erf((X-MEAN)/SD)))");
      return s_scheme;
    }

    //! @brief access to Mean
    //! @return the mean used
    double LogLikelihood::GetMean() const
    {
      return m_ZScore.GetMean();
    };

    //! @brief access to Standard Deviation
    //! @return the sigma used
    double LogLikelihood::GetSigma() const
    {
      return m_ZScore.GetSigma();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the score of the value in the normal distribution with the given mean and sigma
    //! @param ARGUMENT value which z-score shall be calculated from
    //! @return score of the given value
    double LogLikelihood::operator()( const double &ARGUMENT) const
    {
      // make sure, x never reaches 1 to avoid -inf as a result
      const double x( erf( m_ZScore( ARGUMENT)) - 0.000001);
      return log( s_ArbitaryScoreConfidenceOffset * ( 1.0 - x));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LogLikelihood::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ZScore, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LogLikelihood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ZScore, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
