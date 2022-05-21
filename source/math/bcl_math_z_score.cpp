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
#include "math/bcl_math_z_score.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ZScore::s_Instance
    (
      GetObjectInstances().AddInstance( new ZScore())
    );

  //////////
  // data //
  //////////

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ZScore::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "z_score");

      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ZScore::ZScore() :
      m_Mean( 0),
      m_Sigma( 0)
    {
    }

    //! @brief construct with Mean and Standard Deviation
    //! @param MEAN mean value of the considered distribution
    //! @param SIGMA standard deviation of the considered distribution
    //! @param SCHEME
    ZScore::ZScore
    (
      const double MEAN,
      const double SIGMA,
      const std::string &SCHEME
    ) :
      m_Mean( MEAN),
      m_Sigma( SIGMA),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ZScore
    ZScore *ZScore::Clone() const
    {
      return new ZScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ZScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to Mean
    //! @return the mean used in the z-score
    double ZScore::GetMean() const
    {
      return m_Mean;
    };

    //! @brief access to Standard Deviation
    //! @return the sigma used in the z-score
    double ZScore::GetSigma() const
    {
      return m_Sigma;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ZScore::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the z-score of the value in the normal distribution with the given mean and sigma
    //! @param ARGUMENT value which z-score shall be calculated from
    //! @return z-score of the given value
    double ZScore::operator()( const double &ARGUMENT) const
    {
      return m_Sigma ? ( ARGUMENT - m_Mean) / m_Sigma : ARGUMENT - m_Mean;
    }

    //! @brief multiply z-score by given weight
    //! @param WEIGHT weight to multiply by
    //! @return reference to this
    ZScore &ZScore::operator *=( const double WEIGHT)
    {
      m_Sigma /= WEIGHT;
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ZScore::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Mean  , ISTREAM);
      io::Serialize::Read( m_Sigma , ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ZScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Mean  , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Sigma , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_Scheme, OSTREAM,      0);

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
