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
#include "linal/bcl_linal_vector_2d_operations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_2d_operations.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVector2DOperations :
    public ExampleInterface
  {
  public:

    ExampleLinalVector2DOperations *Clone() const
    {
      return new ExampleLinalVector2DOperations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    /////////////////
    // preparation //
    /////////////////

      // construct a few vectors to use in operators
      const linal::Vector2D vector_a( 1.0, 1.0);
      const linal::Vector2D vector_b( -2.0, 0.0);
      const linal::Vector2D vector_c( 0.5, 0.5);
      const linal::Vector2D vector_d( 10.0, 1.0);
      const linal::Vector2D vector_e( 10.0, 10.0);

    ///////////////
    // operators //
    ///////////////

      // test == operator first
      BCL_ExampleCheck( vector_c, linal::Vector2D( 0.5, 0.5));

      // test != operator
      BCL_ExampleCheck( vector_c != linal::Vector2D( 0.55, -0.5), true);

      // test operator +
      BCL_ExampleCheck( +vector_a, vector_a);
      BCL_ExampleCheck( vector_a + vector_b, linal::Vector2D( -1.0, 1.0));
      BCL_ExampleCheck( vector_a + 5.5, linal::Vector2D( 6.5, 6.5));
      BCL_ExampleCheckWithinAbsTolerance( -1.2 + vector_d, linal::Vector2D( 8.8, -0.2), 0.001);

      // test all the operator -
      BCL_ExampleCheck( -vector_a, linal::Vector2D( -1.0, -1.0));
      BCL_ExampleCheck( vector_a - vector_b, linal::Vector2D( 3.0, 1.0));
      BCL_ExampleCheck( vector_a - 5.5, linal::Vector2D( -4.5, -4.5));
      BCL_ExampleCheckWithinAbsTolerance( -1.2 - vector_b, linal::Vector2D( 0.8, -1.2), 0.001);

      // test all the operator *
      BCL_ExampleCheck( vector_a * vector_d, 11.0);
      BCL_ExampleCheck( vector_a * 5.5, linal::Vector2D( 5.5, 5.5));
      BCL_ExampleCheck( -1.2 * vector_d, linal::Vector2D( -12.0, -1.2));

      // test the operator /
      BCL_ExampleCheck( vector_d / -2.5, linal::Vector2D( -4.0, -0.4));

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // Vector functions //
    //////////////////////

      // test EqualWithinTolerance function
      BCL_ExampleCheckWithinAbsTolerance( vector_a, linal::Vector2D( 1.001, 1.0), 0.01);

      // test SquareDistance function
      BCL_ExampleCheckWithinAbsTolerance( linal::SquareDistance( vector_c, vector_d), 90.5, 0.001);

      // test Distance function
      BCL_ExampleCheckWithinAbsTolerance( linal::Distance( vector_c, vector_d), math::Sqrt( 90.5), 0.001);

      // test ProjAngle functions
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::ProjAngle( vector_a, vector_b, vector_c, vector_d),
        acos( -29.0 / math::Sqrt( 90.5) / math::Sqrt( 10.0)),
        0.001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::ProjAngle( vector_a, vector_b, vector_c),
        acos( 2.0 / math::Sqrt( 0.5) / math::Sqrt( 10.0)),
        0.001
      );

      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::ProjAngle( vector_a, vector_b),
        acos( -2.0 / math::Sqrt( 2.0) / 2.0),
        0.001
      );

      // test ScalarProduct function
      BCL_ExampleCheckWithinAbsTolerance( linal::ScalarProduct( vector_a, vector_c), 1.0, 0.0001);

      // test CrossProduct functions
      BCL_ExampleCheckWithinAbsTolerance( linal::CrossProduct( vector_a, vector_b, vector_c, vector_d), 8.0, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( linal::CrossProduct( vector_a, vector_b, vector_c), 1.0, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( linal::CrossProduct( vector_a, vector_b), 2.0, 0.001);

      // test CalculateFootpoint and CalculateDistancePointFromLine functions
      // initialize point that is going to be projected onto line
      const linal::Vector2D point_to_be_projected( 1.0, 1.0);

      // define line by point origin and vector direction
      const linal::Vector2D origin( 0.0, 0.0);
      const linal::Vector2D direction( 1.0, 0.25);

      // calculate distance of point_to_be_projected from line (defined by origin and direction)
      const double distance( linal::CalculateDistancePointFromLine( point_to_be_projected, origin, direction));

      // check whether distance was calculated correctly
      BCL_ExampleCheckWithinAbsTolerance( distance, 0.727670, 0.001);

      // check whether footpoint was calculated correctly
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::CalculateFootpoint( point_to_be_projected, origin, direction),
        linal::Vector2D( 1.17647, 0.29418),
        0.0001
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVector2DOperations

  const ExampleClass::EnumType ExampleLinalVector2DOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVector2DOperations())
  );

} // namespace bcl

