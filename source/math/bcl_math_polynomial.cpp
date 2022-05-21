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
#include "math/bcl_math_polynomial.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialize.h"
#include "math/bcl_math_linear_least_squares.h"
#include "storage/bcl_storage_pair.h"
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
    const util::SiPtr< const util::ObjectInterface> Polynomial::s_Instance
    (
      util::Enumerated< Polynomial>::AddInstance( new Polynomial())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Polynomial::Polynomial()
    {
      UpdateScheme();
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Polynomial
    Polynomial *Polynomial::Clone() const
    {
      return new Polynomial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Polynomial::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetCoefficients returns the coefficients of the polynomial
    //! @return returns vector of coefficients
    const linal::Vector< double> &Polynomial::GetCoefficients() const
    {
      return m_Coefficients;
    }

    //! @brief SetCoefficients set the coefficients of the polynomial
    //! @param COEFFICIENTS the new coefficients to use
    void Polynomial::SetCoefficients( const linal::Vector< double> &COEFFICIENTS)
    {
      m_Coefficients = COEFFICIENTS;
      UpdateScheme();
    }

    //! @brief GetDegree get the highest power of x that this polynomial calculates
    //! @return highest power of x that this polynomial calculates
    size_t Polynomial::GetDegree() const
    {
      return m_Coefficients.IsEmpty() ? 0 : m_Coefficients.GetSize() - 1;
    }

    //! @brief get the polynomial's formula
    //! @return a string like "math::Polynomial(x) = 5.01 + 3.51*x + -0.42*x^2"
    const std::string &Polynomial::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief Get an iterator to the first coefficient
    //! @return an iterator to the first coefficient
    Polynomial::const_iterator Polynomial::Begin() const
    {
      return m_Coefficients.Begin();
    }

    //! @brief Get an iterator to the end of the coefficients
    //! @return an iterator to the end of the coefficients
    Polynomial::const_iterator Polynomial::End() const
    {
      return m_Coefficients.End();
    }

    //! @brief Get an iterator to the first coefficient
    //! @return an iterator to the first coefficient
    Polynomial::iterator Polynomial::Begin()
    {
      return m_Coefficients.Begin();
    }

    //! @brief Get an iterator to the end of the coefficients
    //! @return an iterator to the end of the coefficients
    Polynomial::iterator Polynomial::End()
    {
      return m_Coefficients.End();
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Polynomial::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Represents a polynomial of the form f(x) = a0+a1*x+...+an*x^n");
      serializer.AddInitializer
      (
        "coefficients",
        "coefficients of the terms",
        io::Serialization::GetAgent( &m_Coefficients)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the derivative of the polynomial
    //! @return the derivative of the polynomial function
    Polynomial Polynomial::Derivative() const
    {
      if( GetDegree() == 0) // if polynomial is constant, derivative is zero, so just return the empty polynomial
      {
        return Polynomial();
      }

      // return a polynomial with 1
      linal::Vector< double> derivative_coefficients( GetDegree());

      // power law: d(an*x^n+an-1*x^(n-1)+...+a0) = n*an*x^(n-1)+(n-1)*an-1*x^(n-2)+...+a1
      for
      (
        size_t degree( 1), max_degree( GetDegree());
        degree <= max_degree;
        ++degree
      )
      {
        // d(m_Coefficients(degree)*x^degree)=degree*m_Coefficients(degree)*x^(degree-1)
        derivative_coefficients( degree - 1) = double( degree) * m_Coefficients( degree);
      }

      return MakeFromCoefficients( derivative_coefficients);
    }

    //! @brief create a polynomial with x ^ DEGREE that is the linear least squares fit to given coordinates
    //! @param X_COORDINATES the x coordinates to fit the polynomial to
    //! @param Y_COORDINATES the y coordinates to fit the polynomial to
    //! @param DEGREE the maximum degree that the polynomial should be of
    //! @return the resulting polynomial
    Polynomial Polynomial::MakeFromCoordinates
    (
      const linal::Vector< double> &X_COORDINATES,
      const linal::Vector< double> &Y_COORDINATES,
      const size_t &DEGREE
    )
    {
      BCL_Assert
      (
        X_COORDINATES.GetSize() == Y_COORDINATES.GetSize(),
        "given vectors of non-equal size"
      );

      const size_t number_coordinates( X_COORDINATES.GetSize());

      if( X_COORDINATES.IsEmpty()) // no coordinates to satisfy
      {
        return Polynomial(); // return a polynomial that always returns 0
      }
      else if( X_COORDINATES.GetSize() <= DEGREE) // under-determined system
      {
        // return a smaller polynomial that will pass through the points
        return MakeFromCoordinates( X_COORDINATES, Y_COORDINATES, number_coordinates - 1);
      }

      // given coordinates X(0...m), compute matrix
      // M =
      // {
      //   1,X(0),X(0)^2,X(0)^3,...
      //   1,X(1),X(1)^2,X(1)^3,...
      //   ........................
      //   1,X(m),X(m)^2,X(m)^3,...
      // }
      linal::Matrix< double> coordinate_powers( number_coordinates, DEGREE + 1);

      for( size_t coordinate_index( 0); coordinate_index < number_coordinates; ++coordinate_index)
      {
        double *const matrix_row( coordinate_powers[ coordinate_index]);

        double coordinate( X_COORDINATES( coordinate_index));
        double coordinate_to_degree_power( 1.0); // coordinate ^ degree
        for( size_t degree( 0); degree <= DEGREE; ++degree, coordinate_to_degree_power *= coordinate)
        {
          matrix_row[ degree] = coordinate_to_degree_power;
        }
      }

      // compute the least-squares approximation to the polynomial running through those points
      return MakeFromCoefficients
             (
               // just take the solutions, ignore the chi squared
               LinearLeastSquares::SolutionAndChiSquared( coordinate_powers, Y_COORDINATES).First()
             );
    }

    //! @brief create a polynomial with COEFFICIENTS
    //! @param COEFFICIENTS the coefficients of the new polynomial
    //! @return the resulting polynomial
    Polynomial Polynomial::MakeFromCoefficients( const linal::Vector< double> &COEFFICIENTS)
    {
      Polynomial poly;
      poly.SetCoefficients( COEFFICIENTS);

      return poly;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() compute the polynomial using horner's method
    //! @param ARGUMENT x for the polynomial
    //! @return returns the result of evaluating the polynomial at x=ARGUMENT
    double Polynomial::operator ()( const double &ARGUMENT) const
    {
      if( m_Coefficients.IsEmpty()) // handle the trivial case of no coefficients
      {
        return 0.0;
      }

      double result( m_Coefficients.Last()); // start off with the highest-degree coefficient

      // evaluate polynomial using horner's method, which is faster and more numerically stable than computing powers
      // ((an*x+an-1)x+an-2)+...
      for
      (
        // set itr just before the last element, which we've already used
        const double *itr( m_Coefficients.End() - 2),
                     *itr_end( m_Coefficients.Begin() - 1); // reverse end
        itr != itr_end;
        --itr // iterate from highest degree to lowest degree
      )
      {
        result = result * ARGUMENT + *itr;
      }

      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Polynomial::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Coefficients, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Polynomial::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Coefficients, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief update the scheme
    void Polynomial::UpdateScheme()
    {
      std::ostringstream formula;

      //write scheme
      formula << "F(x) = ";

      //write the formula
      if( !m_Coefficients.IsEmpty())
      {
        formula << m_Coefficients.First(); // constant term

        for( size_t power( 1), degree( GetDegree()); power <= degree; ++power)
        {
          formula << " + " << m_Coefficients( power) << "x"; // successive powers of x

          if( power > 1) // write x^power if the power > 1 (writing x^1 looks sloppy)
          {
            formula << "^" << power;
          }
        }
      }
      else
      {
        formula << "0.0";
      }

      m_Scheme = formula.str();
    }

  } // namespace math
} // namespace bcl

