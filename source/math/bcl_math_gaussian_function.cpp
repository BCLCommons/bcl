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
#include "math/bcl_math_gaussian_function.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> GaussianFunction::s_Instance
    (
      util::Enumerated< FunctionInterfaceSerializable< double, double> >::AddInstance( new GaussianFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GaussianFunction::GaussianFunction() :
      m_Mean( 1),
      m_Sigma( 1),
      m_OneOverTwoSquareSigma( 1 / ( 2 * Sqr( 1)))
    {
    }

    //! @brief construct from mean and sigma
    GaussianFunction::GaussianFunction
    (
      const double MEAN,
      const double SIGMA
    ) :
      m_Mean( MEAN),
      m_Sigma( SIGMA),
      m_OneOverTwoSquareSigma( 1.0 / ( 2.0 * Sqr( SIGMA)))
    {
    }

    //! @brief Clone function
    //! @return pointer to new GaussianFunction
    GaussianFunction *GaussianFunction::Clone() const
    {
      return new GaussianFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GaussianFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to mean
    //! @return the mean of the gaussian
    double GaussianFunction::GetMean() const
    {
      return m_Mean;
    }

    //! @brief access to sigma
    //! @return the sigma of the guassian
    double GaussianFunction::GetSigma() const
    {
      return m_Sigma;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &GaussianFunction::GetAlias() const
    {
      static const std::string s_Name( "GaussianFunction");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get symmetric range relative to the mean for a given multiple of sigma
    //! @param MULTIPLE_SIGMA n multiple of sigma
    //! @return Range[ mean - n * sigma, mean + n * sigma]
    Range< double> GaussianFunction::GetRange( const size_t MULTIPLE_SIGMA) const
    {
      const double half_width( MULTIPLE_SIGMA * m_Sigma);
      return Range< double>( m_Mean - half_width, m_Mean + half_width);
    }

    //! @brief get normalization factor
    //! @return a factor so that integral f(x) -inf to +inf = 1
    double GaussianFunction::GetNormalizationParameter() const
    {
      return double( 1) / ( m_Sigma * Sqrt( 2 * g_Pi));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator f( x) = y = e^( -(x-mean)^2 / ( 2*sigma^2))
    //! @param ARGUMENT the x
    //! @return f( x)
    double GaussianFunction::operator()( const double &ARGUMENT) const
    {
      return exp( -( Sqr( ARGUMENT - m_Mean) * m_OneOverTwoSquareSigma));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GaussianFunction::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Mean, ISTREAM);
      io::Serialize::Read( m_Sigma, ISTREAM);

      // calculate members
      m_OneOverTwoSquareSigma = 1 / ( 2 * Sqr( m_Sigma));

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &GaussianFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Mean, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Sigma, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GaussianFunction::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Represents a Gaussian function. A function of the form f(x) = exp( -( x - mean)^2 / ( 2*sigma^2))"
      );

      parameters.AddInitializer
      (
        "mean",
        "mean of the gaussian distribution",
        io::Serialization::GetAgent( &m_Mean),
        "1"
      );

      parameters.AddInitializer
      (
        "sigma",
        "standard deviation or sigma of the distribution",
        io::Serialization::GetAgent( &m_Sigma),
        "1"
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool GaussianFunction::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_OneOverTwoSquareSigma = double( 1.0) / double( 2.0 * m_Sigma * m_Sigma);
      return true;
    }

  } // namespace math
} // namespace bcl
