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

#ifndef BCL_MATH_POLYNOMIAL_H_
#define BCL_MATH_POLYNOMIAL_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Polynomial
    //! @brief represents a polynomial of the form f(x) = a0+a1*x+...+an*x^n
    //!
    //! @see @link example_math_polynomial.cpp @endlink
    //! @author mendenjl
    //! @date July 20, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Polynomial :
      public FunctionInterfaceSerializable< double, double>
    {
    private:

    //////////
    // data //
    //////////

      linal::Vector< double> m_Coefficients; //!< the coefficients of the vector, stored in order of increasing degree
      std::string     m_Scheme;       //!< the scheme of the function as string like: "f(x) = 5.01 + 3.51*x + -0.42*x^2"

    public:

      typedef linal::Vector< double>::const_iterator const_iterator;
      typedef linal::Vector< double>::iterator       iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Polynomial();

      //! @brief create a polynomial from an input iterator range over an arbitrary container
      template< typename t_InputIterator>
      Polynomial( const t_InputIterator &BEGIN, const t_InputIterator &END) :
        m_Coefficients( BEGIN, END)
      {
        UpdateScheme();
      }

      //! @brief Clone function
      //! @return pointer to new Polynomial
      Polynomial *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetCoefficients returns the coefficients of the polynomial
      //! @return returns vector of coefficients
      const linal::Vector< double> &GetCoefficients() const;

      //! @brief SetCoefficients set the coefficients of the polynomial
      //! @param COEFFICIENTS the new coefficients to use
      void SetCoefficients( const linal::Vector< double> &COEFFICIENTS);

      //! @brief GetDegree get the highest power of x that this polynomial calculates
      //! @return highest power of x that this polynomial calculates
      size_t GetDegree() const;

      //! @brief get the polynomial's formula
      //! @return a string like "math::Polynomial(x) = 5.01 + 3.51*x + -0.42*x^2"
      const std::string &GetScheme() const;

      //! @brief Get an iterator to the first coefficient
      //! @return an iterator to the first coefficient
      const_iterator Begin() const;

      //! @brief Get an iterator to the end of the coefficients
      //! @return an iterator to the end of the coefficients
      const_iterator End() const;

      //! @brief Get an iterator to the first coefficient
      //! @return an iterator to the first coefficient
      iterator Begin();

      //! @brief Get an iterator to the end of the coefficients
      //! @return an iterator to the end of the coefficients
      iterator End();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the derivative function of the polynomial
      //! @return the derivative of this polynomial function
      Polynomial Derivative() const;

      //! @brief create a polynomial up to x ^ DEGREE that is the linear least squares fit to the given coordinates
      //! @param X_COORDINATES the x coordinates to fit the polynomial to
      //! @param Y_COORDINATES the y coordinates to fit the polynomial to
      //! @param DEGREE the maximum degree that the polynomial should be of
      //! @return the resulting polynomial
      static Polynomial MakeFromCoordinates
      (
        const linal::Vector< double> &X_COORDINATES,
        const linal::Vector< double> &Y_COORDINATES,
        const size_t &DEGREE
      );

      //! @brief create a polynomial with COEFFICIENTS
      //! @param COEFFICIENTS the coefficients of the new polynomial
      //! @return the resulting polynomial
      static Polynomial MakeFromCoefficients( const linal::Vector< double> &COEFFICIENTS);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() compute the polynomial using horner's method
      //! @param ARGUMENT x for the polynomial
      //! @return returns the result of evaluating the polynomial at x=ARGUMENT
      double operator ()( const double &ARGUMENT) const;

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

      //! @brief update the scheme
      void UpdateScheme();

    }; // class Polynomial

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_POLYNOMIAL_H_

