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

#ifndef BCL_MATH_Z_SCORE_H_
#define BCL_MATH_Z_SCORE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "function/bcl_function_unary_interface.h"

//external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ZScore
    //! @brief is a class to calculate the z-score of a value in a given normal distribution
    //! @details The z-score indicates how many standard deviations an value is above or below the mean.
    //! For further information go to http://en.wikipedia.org/wiki/Standard_score
    //!
    //! @see @link example_math_z_score.cpp @endlink
    //! @author bitterd, woetzen
    //! @date Jun 9, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ZScore :
      public function::UnaryInterface< const double, double>
    {

    private:

    //////////
    // data //
    //////////

      double      m_Mean;   //!< mean of the distribution
      double      m_Sigma;  //!< sigma of the distribution
      std::string m_Scheme; //! scheme to be used in outputting schemes

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ZScore();

      //! @brief construct from mean and sd
      //! @param MEAN mean value of the considered distribution
      //! @param SIGMA standard deviation of the considered distribution
      ZScore
      (
        const double MEAN,
        const double SIGMA,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new ZScore
      ZScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to Mean
      //! @return the mean used in the z-score
      double GetMean() const;

      //! @brief access to Standard Deviation
      //! @return the standard deviation used in the z-score
      double GetSigma() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates the z-score of the value in the normal distribution with the given mean and sigma
      //! @param ARGUMENT value which z-score shall be calculated from
      //! @return z-score of the given value
      double operator()( const double &ARGUMENT) const;

      //! @brief multiply z-score by given weight
      //! @param WEIGHT weight to multiply by
      //! @return reference to this
      ZScore &operator *=( const double WEIGHT);

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

    }; // class ZScore

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_Z_SCORE_H_
