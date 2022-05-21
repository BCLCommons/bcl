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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "math/bcl_math_polynomial.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_polynomial.cpp
  //!
  //! @author mendenjl
  //! @date August 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathPolynomial :
    public ExampleInterface
  {
  public:

    ExampleMathPolynomial *Clone() const
    {
      return new ExampleMathPolynomial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_ExampleCheck( math::Polynomial().GetCoefficients().IsEmpty(), true);

      // test clone
      BCL_ExampleCheck( util::ShPtr< math::Polynomial>( math::Polynomial().Clone()).IsDefined(), true);

      // now make a polynomial f(x)=4+3x+2x^2
      linal::Vector< double> coefficients( 3);
      coefficients( 0) = 4.0;
      coefficients( 1) = 3.0;
      coefficients( 2) = 2.0;

      // can either use via MakeFromCoefficients
      math::Polynomial made_polynomial( math::Polynomial::MakeFromCoefficients( coefficients));

      // or call SetCoefficients
      math::Polynomial set_polynomial;
      set_polynomial.SetCoefficients( coefficients);

      // make sure the coefficients are equal to one another
      BCL_ExampleCheck( made_polynomial.GetCoefficients(), set_polynomial.GetCoefficients());

      // assert that there are 3 coefficients
      BCL_ExampleAssert( made_polynomial.GetCoefficients().GetSize(), size_t( 3));

      // make sure the polynomial has the right coefficients.
      // Its okay to directly check for equality of the doubles in this case:
      //   IEEE guarantees that setting doubles that were set to the same value must compare equal
      //   (except nan)
      BCL_ExampleIndirectCheck
      (
        made_polynomial.GetCoefficients()( 0),
        coefficients( 0),
        "math::Polynomial::SetCoefficients()"
      );
      BCL_ExampleIndirectCheck
      (
        made_polynomial.GetCoefficients()( 1),
        coefficients( 1),
        "math::Polynomial::SetCoefficients()"
      );
      BCL_ExampleIndirectCheck
      (
        made_polynomial.GetCoefficients()( 2),
        coefficients( 2),
        "math::Polynomial::SetCoefficients()"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( math::Polynomial().GetClassIdentifier(), GetStaticClassName< math::Polynomial>());

      // test GetDegree function
      BCL_ExampleCheck( math::Polynomial().GetDegree(), 0);
      BCL_ExampleCheck( made_polynomial.GetDegree(), 2);

    ///////////////
    // operators //
    ///////////////

      // check operator()
      const double x( 6.0);
      BCL_ExampleCheck( made_polynomial( x), ( 2.0 * x + 3.0) * x + 4.0);

    ////////////////
    // operations //
    ////////////////

      // test GetDerivative
      math::Polynomial made_polynomial_derivative( made_polynomial.Derivative());
      BCL_ExampleAssert( made_polynomial_derivative.GetDegree(), size_t( 1));
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance( made_polynomial_derivative.GetCoefficients()( 0), 3.0, std::numeric_limits< double>::epsilon()),
        true
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance( made_polynomial_derivative.GetCoefficients()( 1), 4.0, std::numeric_limits< double>::epsilon()),
        true
      );

      // test MakeFromCoordinates
      // passing in no coordinates should give a polynomial of 0 degree
      BCL_ExampleCheck
      (
        math::Polynomial::MakeFromCoordinates( linal::Vector< double>(), linal::Vector< double>(), 10).GetCoefficients().GetSize(),
        0
      );

      linal::Vector< double> x_coordinates( 2);
      x_coordinates( 0) = 0.0;
      x_coordinates( 1) = 2.0;

      math::Polynomial y_equals_x( math::Polynomial::MakeFromCoordinates( x_coordinates, x_coordinates, 1));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( y_equals_x( 5.0), 5.0, 1e-8),
        true,
        "math::Polynomial::MakeFromCoordinates"
      );

      BCL_MessageStd
      (
        "y_equals_x( 5.0) == " + util::Format()( y_equals_x( 5.0))
        + " diff from 5.0= " + util::Format()( y_equals_x( 5.0) - 5.0)
      );

      // now add one to the y-coordinates to make sure we still get the right line
      math::Polynomial y_equals_x_plus_one
      (
        math::Polynomial::MakeFromCoordinates( x_coordinates, x_coordinates + 1.0, 1)
      );

      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( y_equals_x_plus_one( 5.0), 6.0, 1.0e-8),
        true,
        "math::Polynomial::MakeFromCoordinates"
      );

      // now sample some points from made_polynomial (above), to see if we recover roughly the same coefficients
      linal::Vector< double> x_coordinates_from_quadratic( 5), y_coordinates_from_quadratic( 5);

      x_coordinates_from_quadratic( 0) = 1.0;
      x_coordinates_from_quadratic( 1) = 3.0;
      x_coordinates_from_quadratic( 2) = 2.0;
      x_coordinates_from_quadratic( 3) = 7.0;
      x_coordinates_from_quadratic( 4) = -1.0; // the points need not be sorted

      double *y_itr( y_coordinates_from_quadratic.Begin());
      for
      (
        const double *x_itr( x_coordinates_from_quadratic.Begin()), *x_itr_end( x_coordinates_from_quadratic.End());
        x_itr != x_itr_end;
        ++x_itr, ++y_itr
      )
      {
        *y_itr = made_polynomial( *x_itr);
      }

      // make a new polynomial that should pass through the same points as made_polynomial
      math::Polynomial reverse_engineered_made_polynomial
      (
        math::Polynomial::MakeFromCoordinates( x_coordinates_from_quadratic, y_coordinates_from_quadratic, 2)
      );

      // check that the coefficients from reverse_engineered_made_polynomial are essentially the same as made_polynomial
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          made_polynomial.GetCoefficients(),
          reverse_engineered_made_polynomial.GetCoefficients(),
          0.0,
          0.01
        ),
        true,
        "math::Polynomial::MakeFromCoordinates with excess coordinates"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test that everything written out is read back in
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( made_polynomial, math::Polynomial()), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleMathPolynomial::s_Instance
  (
    GetExamples().AddEnum( ExampleMathPolynomial())
  );

} // namespace bcl
