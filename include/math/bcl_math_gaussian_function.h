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

#ifndef BCL_MATH_GAUSSIAN_FUNCTION_H_
#define BCL_MATH_GAUSSIAN_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GaussianFunction
    //! @brief represents a Gaussian function
    //! @details a function of the form f(x) = exp( -( x - mean)^2 / ( 2*sigma^2))
    //!
    //! @see @link example_math_gaussian_function.cpp @endlink
    //! @author woetzen
    //! @date Apr 26, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GaussianFunction :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_Mean;                  //!< mean of the gaussian distribution
      double m_Sigma;                 //!< standard deviation or sigma of the distribution
      double m_OneOverTwoSquareSigma; //! 1 / ( 2 * sigma^2) - used to save computation time

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
      GaussianFunction();

      //! @brief construct from mean and sigma
      GaussianFunction
      (
        const double MEAN,
        const double SIGMA
      );

      //! @brief Clone function
      //! @return pointer to new GaussianFunction
      GaussianFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief access to mean
      //! @return the mean of the gaussian
      double GetMean() const;

      //! @brief access to sigma
      //! @return the sigma of the guassian
      double GetSigma() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get symmetric range relative to the mean for a given multiple of sigma
      //! @param MULTIPLE_SIGMA n multiple of sigma
      //! @return Range[ mean - n * sigma, mean + n * sigma]
      Range< double> GetRange( const size_t MULTIPLE_SIGMA) const;

      //! @brief get normalization factor
      //! 1 / ( sigma * sqrt( 2*Pi))
      //! @return a factor so that integral f(x) -inf to +inf = 1
      double GetNormalizationParameter() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator f( x) = y = e^( (x-mean)^2 / ( 2*sigma^2))
      //! @param ARGUMENT the x
      //! @return f( x)
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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class GaussianFunction

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_GAUSSIAN_FUNCTION_H_
