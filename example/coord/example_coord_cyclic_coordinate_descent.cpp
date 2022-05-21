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
#include "coord/bcl_coord_cyclic_coordinate_descent.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_cyclic_coordinate_descent.cpp
  //!
  //! @author alexanns
  //! @date Oct 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordCyclicCoordinateDescent :
    public ExampleInterface
  {
  public:

    ExampleCoordCyclicCoordinateDescent *Clone() const
    {
      return new ExampleCoordCyclicCoordinateDescent( *this);
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

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      {
        const linal::Vector3D rotation_axis_origin( -2.378, 0.493, -5.464);
        const linal::Vector3D rotation_axis_end( -1.400, 1.725, -4.573);

        // create rotation axis
        coord::LineSegment3D rotation_axis( rotation_axis_origin, rotation_axis_end);

        // create moving point
        const linal::Vector3D moving_point( 0.536, 1.029, -4.713);

        // create target point
        const linal::Vector3D target_point( -0.576, 0.674, -3.002);

        // create the pair of target and moving points
        coord::CyclicCoordinateDescent::TargetAndMovingPointPair points( target_point, moving_point);

        const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> point_list( 1, points);

        const coord::CyclicCoordinateDescent ccd;

        const coord::CyclicCoordinateDescent::ResultsAndCoefficients result( ccd.GetOptimalRotationandDistance( rotation_axis, point_list));

        const coord::CyclicCoordinateDescent::ResultsAndCoefficients expected_result( 1.0889, 0.000000052428, 7.99585, 3.70576, 7.08527);
        BCL_Example_Check
        (
          result == expected_result,
          "result.GetOptimalRotation() " + util::Format()( result.GetOptimalRotation()) + " expected_result.GetOptimalRotation() " + util::Format()( expected_result.GetOptimalRotation()) + "\n" +
          "result.GetMinimumDistanceSum() " + util::Format()( result.GetMinimumDistanceSum()) + " expected_result.GetMinimumDistanceSum() " + util::Format()( expected_result.GetMinimumDistanceSum()) + "\n" +
          "result.GetCoefficientA() " + util::Format()( result.GetCoefficientA()) + " expected_result.GetCoefficientA() " + util::Format()( expected_result.GetCoefficientA()) + "\n" +
          "result.GetCoefficientB() " + util::Format()( result.GetCoefficientB()) + " expected_result.GetCoefficientB() " + util::Format()( expected_result.GetCoefficientB()) + "\n" +
          "result.GetCoefficientC() " + util::Format()( result.GetCoefficientC()) + " expected_result.GetCoefficientC() " + util::Format()( expected_result.GetCoefficientC())
        );

        {
          math::TransformationMatrix3D transform( -rotation_axis_origin);
          transform( math::RotationMatrix3D( ( rotation_axis_end - rotation_axis_origin), result.GetOptimalRotation()));
          transform(  rotation_axis_origin);
          linal::Vector3D transformed_moving_point( moving_point);
          transformed_moving_point.Transform( transform);
          const linal::Vector3D expected_transformed_moving_point( -0.575863, 0.674101, -3.00185);
          BCL_MessageDbg( "transformed coordinates of moving_point" + util::Format()( transformed_moving_point));
          BCL_Example_Check
          (
            math::EqualWithinTolerance( transformed_moving_point, expected_transformed_moving_point), "expected " +
            util::Format()( expected_transformed_moving_point) + " but is " +
            util::Format()( transformed_moving_point)
          );
        }
      }

      {
        const linal::Vector3D rotation_axis_origin( 38.972, 24.079, 8.055);
        const linal::Vector3D rotation_axis_end( 37.886, 25.043, 8.022);

        // create rotation axis
        coord::LineSegment3D rotation_axis( rotation_axis_origin, rotation_axis_end);

        // create moving points
        const linal::Vector3D moving_point_a( 40.251, 28.051, 9.993);
        const linal::Vector3D moving_point_b( 40.415, 28.733, 11.251);
        const linal::Vector3D moving_point_c( 41.799,28.692, 11.808);

        // create target points
        const linal::Vector3D target_point_a( 38.258, 25.884, 12.227);
        const linal::Vector3D target_point_b( 37.617, 25.661, 13.498);
        const linal::Vector3D target_point_c( 38.501, 25.072, 14.548);

        const coord::CyclicCoordinateDescent::TargetAndMovingPointPair points_a( target_point_a, moving_point_a);
        const coord::CyclicCoordinateDescent::TargetAndMovingPointPair points_b( target_point_b, moving_point_b);
        const coord::CyclicCoordinateDescent::TargetAndMovingPointPair points_c( target_point_c, moving_point_c);

        storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> point_list;
        point_list.PushBack( points_a);
        point_list.PushBack( points_b);
        point_list.PushBack( points_c);

        const coord::CyclicCoordinateDescent ccd;

        const coord::CyclicCoordinateDescent::ResultsAndCoefficients result( ccd.GetOptimalRotationandDistance( rotation_axis, point_list));
        const coord::CyclicCoordinateDescent::ResultsAndCoefficients expected_result( 0.888471, 0.00000783945, 182.627, 115.165, 141.738);
        BCL_Example_Check
        (
          result == expected_result,
          "result.GetOptimalRotation() " + util::Format()( result.GetOptimalRotation()) + " expected_result.GetOptimalRotation() " + util::Format()( expected_result.GetOptimalRotation()) + "\n" +
          "result.GetMinimumDistanceSum() " + util::Format()( result.GetMinimumDistanceSum()) + " expected_result.GetMinimumDistanceSum() " + util::Format()( expected_result.GetMinimumDistanceSum()) + "\n" +
          "result.GetCoefficientA() " + util::Format()( result.GetCoefficientA()) + " expected_result.GetCoefficientA() " + util::Format()( expected_result.GetCoefficientA()) + "\n" +
          "result.GetCoefficientB() " + util::Format()( result.GetCoefficientB()) + " expected_result.GetCoefficientB() " + util::Format()( expected_result.GetCoefficientB()) + "\n" +
          "result.GetCoefficientC() " + util::Format()( result.GetCoefficientC()) + " expected_result.GetCoefficientC() " + util::Format()( expected_result.GetCoefficientC())
        );

        {
          math::TransformationMatrix3D transform( -rotation_axis_origin);
          transform( math::RotationMatrix3D( ( rotation_axis_end - rotation_axis_origin), result.GetOptimalRotation()));
          transform(  rotation_axis_origin);
          linal::Vector3D transformed_moving_point_a( moving_point_a);
          transformed_moving_point_a.Transform( transform);
          const linal::Vector3D expected_transformed_moving_point_a( 38.2585, 25.8828, 12.227);
          BCL_MessageDbg( "transformed coordinates of moving_point_a" + util::Format()( transformed_moving_point_a));
          BCL_Example_Check
          (
            math::EqualWithinTolerance( transformed_moving_point_a, expected_transformed_moving_point_a), "expected " +
            util::Format()( expected_transformed_moving_point_a) + " but is " +
            util::Format()( transformed_moving_point_a)
          );
        }

        {
          math::TransformationMatrix3D transform( -rotation_axis_origin);
          transform( math::RotationMatrix3D( ( rotation_axis_end - rotation_axis_origin), result.GetOptimalRotation()));
          transform(  rotation_axis_origin);
          linal::Vector3D transformed_moving_point_b( moving_point_b);
          transformed_moving_point_b.Transform( transform);
          const linal::Vector3D expected_transformed_moving_point_b( 37.6186, 25.6596, 13.498);
          BCL_MessageDbg( "transformed coordinates of moving_point_b" + util::Format()( transformed_moving_point_b));
          BCL_Example_Check
          (
            math::EqualWithinTolerance( transformed_moving_point_b, expected_transformed_moving_point_b), "expected " +
            util::Format()( expected_transformed_moving_point_b) + " but is " +
            util::Format()( transformed_moving_point_b)
          );
        }

        {
          math::TransformationMatrix3D transform( -rotation_axis_origin);
          transform( math::RotationMatrix3D( ( rotation_axis_end - rotation_axis_origin), result.GetOptimalRotation()));
          transform(  rotation_axis_origin);
          linal::Vector3D transformed_moving_point_c( moving_point_c);
          transformed_moving_point_c.Transform( transform);
          const linal::Vector3D expected_transformed_moving_point_c( 38.502, 25.0715, 14.5474);
          BCL_MessageDbg( "transformed coordinates of moving_point_a" + util::Format()( transformed_moving_point_c));
          BCL_Example_Check
          (
            math::EqualWithinTolerance( transformed_moving_point_c, expected_transformed_moving_point_c), "expected " +
            util::Format()( expected_transformed_moving_point_c) + " but is " +
            util::Format()( transformed_moving_point_c)
          );
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordCyclicCoordinateDescent

  const ExampleClass::EnumType ExampleCoordCyclicCoordinateDescent::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordCyclicCoordinateDescent())
  );

} // namespace bcl
