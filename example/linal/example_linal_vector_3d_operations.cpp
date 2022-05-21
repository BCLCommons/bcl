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
#include "linal/bcl_linal_vector_3d_operations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_3d_operations.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVector3DOperations :
    public ExampleInterface
  {
  public:

    ExampleLinalVector3DOperations *Clone() const
    {
      return new ExampleLinalVector3DOperations( *this);
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
      const linal::Vector3D vector_a( 1.0, 1.0, 1.0);
      const linal::Vector3D vector_b( -2.0, 0.0, 1.0);
      const linal::Vector3D vector_c( 0.5, 0.5, -1.0);
      const linal::Vector3D vector_d( 10.0, 10.0, -5.0);
      const linal::Vector3D vector_e( 10.0, 10.0, 10.0);

    ///////////////
    // operators //
    ///////////////

      // test == operator first
      BCL_Example_Check
      (
        vector_c == linal::Vector3D( 0.5, 0.5, -1.0),
        "operator == isn't working correctly"
      );

      // test != operator first
      BCL_Example_Check
      (
        vector_c != linal::Vector3D( 0.55, -0.5, -1.0),
        "operator != isn't working correctly"
      );

      // test all the operator +
      const linal::Vector3D operator_plus_one_argument( +vector_a);
      const linal::Vector3D operator_plus_two_arguments( vector_a + vector_b);
      const linal::Vector3D operator_vector_plus_double( vector_a + 5.5);
      const linal::Vector3D operator_double_plus_vector( -1.2 + vector_d);

      BCL_Example_Check
      (
        operator_plus_one_argument == vector_a,
        "operator plus with one argument isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_plus_two_arguments == linal::Vector3D( -1.0, 1.0, 2.0),
        "operator plus with two arguments isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_vector_plus_double == linal::Vector3D( 6.5, 6.5, 6.5),
        "operator vector plus double isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_double_plus_vector == linal::Vector3D( 8.8, 8.8, -6.2),
        "operator double plus vector isn't working correctly"
      );

      // test all the operator -
      const linal::Vector3D operator_minus_one_argument( -vector_a);
      const linal::Vector3D operator_minus_two_arguments( vector_a - vector_b);
      const linal::Vector3D operator_vector_minus_double( vector_a - 5.5);
      const linal::Vector3D operator_double_minus_vector( -1.2 - vector_b);

      BCL_Example_Check
      (
        operator_minus_one_argument == linal::Vector3D( -1.0, -1.0, -1.0),
        "operator minus with one argument isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_minus_two_arguments == linal::Vector3D( 3.0, 1.0, 0.0),
        "operator minus with two arguments isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_vector_minus_double == linal::Vector3D( -4.5, -4.5, -4.5),
        "operator vector minus double isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_double_minus_vector == linal::Vector3D( 0.8, -1.2, -2.2),
        "operator double minus vector isn't working correctly"
      );

      // test all the operator *
      const linal::Vector3D operator_times_two_arguments( vector_a * vector_d);
      const linal::Vector3D operator_vector_times_double( vector_a * 5.5);
      const linal::Vector3D operator_double_times_vector( -1.2 * vector_d);

      BCL_Example_Check
      (
        operator_times_two_arguments == linal::Vector3D( 15.0),
        "operator times with two arguments isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_vector_times_double == linal::Vector3D( 5.5, 5.5, 5.5),
        "operator vector times double isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_double_times_vector == linal::Vector3D( -12.0, -12.0, 6.0),
        "operator double times vector isn't working correctly"
      );

      // test the operator /
      const linal::Vector3D operator_vector_divided_by_double( vector_d / -2.5);
      BCL_Example_Check
      (
        operator_vector_divided_by_double == linal::Vector3D( -4.0, -4.0, 2.0),
        "operator vector divided by double isn't working correctly"
      );

    ///////////////
    // operators //
    ///////////////

      BCL_Example_Check
      (
        operator_minus_one_argument == linal::Vector3D( -1.0, -1.0, -1.0),
        "operator minus with one argument isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_minus_two_arguments == linal::Vector3D( 3.0, 1.0, 0.0),
        "operator minus with two arguments isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_vector_minus_double == linal::Vector3D( -4.5, -4.5, -4.5),
        "operator vector minus double isn't working correctly"
      );
      BCL_Example_Check
      (
        operator_double_minus_vector == linal::Vector3D( 0.8, -1.2, -2.2),
        "operator double minus vector isn't working correctly"
      );

    //////////////////////
    // Vector functions //
    //////////////////////

      // test EqualWithinTolerance function
      BCL_Example_Check
      (
        math::EqualWithinTolerance( vector_a, linal::Vector3D( 1.001, 1.0, 1.0), 0.01, 0.01),
        "function EqualWithinTolerance isn't working correctly"
      );
      BCL_Example_Check
      (
        !math::EqualWithinTolerance( vector_a, linal::Vector3D( 1.025, 1.0, 1.0), 0.01, 0.01),
        "function EqualWithinTolerance isn't working correctly"
      );

      // test SquareDistance function
      const double square_distance_vector_c_vector_d( linal::SquareDistance( vector_c, vector_d));
      BCL_Example_Check
      (
        square_distance_vector_c_vector_d == 196.5,
        "function SquareDistance isn't working correctly"
      );

      // test Distance function
      const double distance_vector_c_vector_d( linal::Distance( vector_c, vector_d));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( distance_vector_c_vector_d, math::Sqrt( 196.5)),
        "function Distance isn't working correctly"
      );

      // test ProjAngle functions
      const double angle_vector_a_vector_b_vector_c_vector_d( linal::ProjAngle( vector_a, vector_b, vector_c, vector_d));
      const double angle_vector_a_vector_b_vector_c( linal::ProjAngle( vector_a, vector_b, vector_c));
      const double angle_vector_a_vector_b( linal::ProjAngle( vector_a, vector_b));

      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          angle_vector_a_vector_b_vector_c_vector_d, acos( -38.0 / math::Sqrt( 10.0) / math::Sqrt( 196.5))
        ),
        "function ProjAngle with 4 vectors isn't working correctly"
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          angle_vector_a_vector_b_vector_c, acos( 2.0 / math::Sqrt( 10.0) / math::Sqrt( 4.5))
        ),
        "function ProjAngle with 3 vectors isn't working correctly"
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          angle_vector_a_vector_b, acos( -1.0 / math::Sqrt( 3.0) / math::Sqrt( 5.0))
        ),
        "function ProjAngle with 2 vectors isn't working correctly"
      );

      // test ScalarProduct function
      const double scalar_product_vector_a_vector_c( linal::ScalarProduct( vector_a, vector_c));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( scalar_product_vector_a_vector_c, 0.0),
        "function ScalarProduct isn't working correctly"
      );

      // test CrossProduct functions
      const linal::Vector3D cross_product_vector_a_vector_b_vector_c_vector_d
      (
        linal::CrossProduct( vector_a, vector_b, vector_c, vector_d)
      );
      const linal::Vector3D cross_product_vector_a_vector_b_vector_c( linal::CrossProduct( vector_a, vector_b, vector_c));
      const linal::Vector3D cross_product_vector_a_vector_b( linal::CrossProduct( vector_a, vector_b));
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          cross_product_vector_a_vector_b_vector_c_vector_d, linal::Vector3D( 4.0, -12.0, -19.0)
        ),
        "function CrossProduct with 4 vectors isn't working correctly"
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          cross_product_vector_a_vector_b_vector_c, linal::Vector3D( 2.0, -6.0, 1.0)
        ),
        "function CrossProduct with 3 vectors isn't working correctly"
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          cross_product_vector_a_vector_b, linal::Vector3D( 1.0, -3.0, 2.0)
        ),
        "function CrossProduct with 2 vectors isn't working correctly"
      );

      // test Dihedral function
      const double dihedral_vector_a_vector_b_vector_c_vector_d
      (
        linal::Dihedral( vector_a, vector_b, vector_c, vector_d)
      );
      const linal::Vector3D cross_1( linal::CrossProduct( vector_a, vector_b, vector_b, vector_c));
      const linal::Vector3D cross_2( linal::CrossProduct( vector_b, vector_c, vector_c, vector_d));
      const double right_dihedral
      (
        std::atan2( ( vector_c - vector_b).Norm() * ( vector_b - vector_a) * cross_2, cross_1 * cross_2)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          dihedral_vector_a_vector_b_vector_c_vector_d,
          right_dihedral
        ),
        "function Dihedral isn't working correctly"
      );
      const double dihedral_vector_d_vector_c_vector_b_vector_a
      (
        linal::Dihedral( vector_d, vector_c, vector_b, vector_a)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        dihedral_vector_a_vector_b_vector_c_vector_d,
        dihedral_vector_d_vector_c_vector_b_vector_a,
        0.01
      );

      // test CoordinatesDihedral function
      const double correct_dihedral_angle( math::g_Pi);
      const linal::Vector3D dihedral_coordinates_vector_a_vector_b_vector_c
      (
        linal::CoordinatesDihedral( vector_a, vector_b, vector_c, 1.5, math::g_Pi / 2, correct_dihedral_angle)
      );
      const linal::Vector3D correct_dihedral_coordinates( linal::Vector3D( 0.92592, 1.22224, 2.48159));
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          linal::Dihedral( vector_c, vector_b, vector_a, dihedral_coordinates_vector_a_vector_b_vector_c),
          -correct_dihedral_angle
        ),
        "dihedral angle should be " + util::Format()( correct_dihedral_angle) + " but is " +
        util::Format()( linal::Dihedral( vector_c, vector_b, vector_a, dihedral_coordinates_vector_a_vector_b_vector_c))
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          dihedral_coordinates_vector_a_vector_b_vector_c,
          correct_dihedral_coordinates
        ),
        "function CoordinatesDihedral isn't working correctly "
      );

      // test CoordinatesAngle function
      const linal::Vector3D coordinates_angle_vector_a_vector_b_vector_c
      (
        linal::CoordinatesAngle( vector_a, vector_b, vector_c, 1.5, math::g_Pi / 2, math::g_Pi / 2)
      );
      const linal::Vector3D correct_angle_coordinates( linal::Vector3D( 1.46852, -0.405564, 1.23426));
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          coordinates_angle_vector_a_vector_b_vector_c,
          correct_angle_coordinates
        ),
        "function CoordinatesAngle isn't working correctly "
      );

      // test CoordinatesLinear function
      const linal::Vector3D coordinates_x_linear( linal::CoordinatesLinear( vector_a, vector_e, math::Sqrt( 3.0)));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( coordinates_x_linear, linal::Vector3D( 0.0, 0.0, 0.0)),
        "function CoordinatesLinear isn't working correctly " + util::Format()( coordinates_x_linear)
      );

      // test CoordinatesTrigonal function
      const linal::Vector3D coordinates_x_trigonal( linal::CoordinatesTrigonal( vector_a, vector_b, vector_c, 5.0));
      const linal::Vector3D correct_coordinates_x_trigonal( linal::Vector3D( 4.67525, 2.71269, 3.92562));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( coordinates_x_trigonal, correct_coordinates_x_trigonal),
        "function CoordinatesTrigonal isn't working correctly "
      );

      // test CoordinatesTetrahedral function
      BCL_ExampleCheckWithinTolerance
      (
        linal::CoordinatesTetrahedral( vector_a, vector_b, vector_c, vector_d, 5.0),
        linal::Vector3D( 4.14384, -0.664385, 4.5137),
        0.001
      );

      // test CalculateFootpoint and CalculateDistancePointFromLine functions
      // initialize point that is going to be projected onto line
      const linal::Vector3D point_to_be_projected( 1.0, 1.0, 1.0);

      // define line by point origin and vector direction
      const linal::Vector3D origin( 0.0, 0.0, 0.0);
      const linal::Vector3D direction( 1.0, 0.0, 0.0);

      // calculate distance of point_to_be_projected from line (defined by origin and direction)
      const double distance( linal::CalculateDistancePointFromLine( point_to_be_projected, origin, direction));

      // check whether distance was calculated correctly
      BCL_ExampleCheck( distance, math::Sqrt( double( 2.0)));

      // calculate footpoint of point_to_be_projected onto line (defined by origin and direction)
      const linal::Vector3D footpoint( linal::CalculateFootpoint( point_to_be_projected, origin, direction));

      // check whether distance was calculated correctly
      BCL_ExampleCheck( footpoint, linal::Vector3D( 1.0, 0.0, 0.0));

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVector3DOperations

  const ExampleClass::EnumType ExampleLinalVector3DOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVector3DOperations())
  );

} // namespace bcl

