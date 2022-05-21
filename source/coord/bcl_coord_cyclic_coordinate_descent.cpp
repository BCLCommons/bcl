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
#include "coord/bcl_coord_cyclic_coordinate_descent.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_line_segment_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::TargetAndMovingPointPair::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent::TargetAndMovingPointPair())
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::ResultsAndCoefficients::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent::ResultsAndCoefficients())
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    CyclicCoordinateDescent::TargetAndMovingPointPair::TargetAndMovingPointPair() :
      m_TargetPoint(),
      m_MovingPoint()
    {
    }

    //! constructor taking a target point and a moving point
    //! @param TARGET_POINT the fixed target point desired to be found
    //! @param MOVING_POINT the point which is being moved in an attempt to reach the target point
    CyclicCoordinateDescent::TargetAndMovingPointPair::TargetAndMovingPointPair
    (
      const linal::Vector3D &TARGET_POINT, const linal::Vector3D &MOVING_POINT
    ) :
      m_TargetPoint( TARGET_POINT),
      m_MovingPoint( MOVING_POINT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TargetAndMovingPointPair
    CyclicCoordinateDescent::TargetAndMovingPointPair *CyclicCoordinateDescent::TargetAndMovingPointPair::Clone() const
    {
      return new TargetAndMovingPointPair( *this);
    }

    //! @brief default constructor
    CyclicCoordinateDescent::ResultsAndCoefficients::ResultsAndCoefficients() :
      m_OptimalRotation(),
      m_MinimumDistanceSum(),
      m_CoefficientA(),
      m_CoefficientB(),
      m_CoefficientC()
    {
    }

    //! @brief constructor taking all of the member variables
    //! @param OPTIMAL_ROTATION the rotation to minimize the distance sum between target and moving points
    //! @param MINIMUM_DISTANCE_SUM the sum of distances between target and moving points that will result after
    //!        the moving points are rotated by "m_OptimalRotation"
    //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
    CyclicCoordinateDescent::ResultsAndCoefficients::ResultsAndCoefficients
    (
      const double OPTIMAL_ROTATION, const double MINIMUM_DISTANCE_SUM, const double COEFFICIENT_A,
      const double COEFFICIENT_B, const double COEFFICIENT_C
    ) :
      m_OptimalRotation( OPTIMAL_ROTATION),
      m_MinimumDistanceSum( MINIMUM_DISTANCE_SUM),
      m_CoefficientA( COEFFICIENT_A),
      m_CoefficientB( COEFFICIENT_B),
      m_CoefficientC( COEFFICIENT_C)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ResultsAndCoefficients
    CyclicCoordinateDescent::ResultsAndCoefficients *CyclicCoordinateDescent::ResultsAndCoefficients::Clone() const
    {
      return new ResultsAndCoefficients( *this);
    }

    //! @brief default constructor
    CyclicCoordinateDescent::CyclicCoordinateDescent()
    {
    }

    //! @brief Clone function
    //! @return pointer to new CyclicCoordinateDescent
    CyclicCoordinateDescent *CyclicCoordinateDescent::Clone() const
    {
      return new CyclicCoordinateDescent( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::TargetAndMovingPointPair::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::ResultsAndCoefficients::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTargetPoint gives the fixed target point desired to be found
    //! @return returns linal::Vector3D which is "m_TargetPoint"
    const linal::Vector3D &CyclicCoordinateDescent::TargetAndMovingPointPair::GetTargetPoint() const
    {
      return m_TargetPoint;
    }

    //! @brief GetMovingPoint gives the point which is being moved in an attempt to reach the target point
    //! @return returns linal::Vector3D which is "m_MovingPoint"
    const linal::Vector3D &CyclicCoordinateDescent::TargetAndMovingPointPair::GetMovingPoint() const
    {
      return m_MovingPoint;
    }

    //! @brief GetOptimalRotation gives the rotation to minimize the distance sum between target and moving points
    //! @return returns "m_OptimalRotation"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetOptimalRotation() const
    {
      return m_OptimalRotation;
    }

    //! @brief GetMinimumDistanceSum gives the sum of distances between target and moving points that will result
    //!        after the moving points are rotated by "m_OptimalRotation"
    //! @return returns "m_MinimumDistanceSum"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetMinimumDistanceSum() const
    {
      return m_MinimumDistanceSum;
    }

    //! @brief GetCoefficientA gives "m_CoefficientA"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "m_CoefficientA"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientA() const
    {
      return m_CoefficientA;
    }

    //! @brief GetCoefficientB gives "GetCoefficientB"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "GetCoefficientB"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientB() const
    {
      return m_CoefficientB;
    }

    //! @brief GetCoefficientC gives "GetCoefficientC"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "GetCoefficientC"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientC() const
    {
      return m_CoefficientC;
    }

    //! @brief GetOptimalRotationandDistance gives rotation to minimize the distance between target and moving points
    //! It calculates the rotation around "ROTATION_AXIS" that needs to be performed on all moving points so that
    //! after the rotation, the sum distance from each of the moving points to corresponding fixed points will be
    //! minimized
    //! @param ROTATION_AXIS this is the axis the rotation will be performed around
    //! @param TARGET_MOVING_POINTS the list of moving points and the corresponding fixed target points the moving
    //!        points are trying to reach
    //! @return returns a ResultsAndCoefficients which holds the results of the calculation
    CyclicCoordinateDescent::ResultsAndCoefficients CyclicCoordinateDescent::GetOptimalRotationandDistance
    (
      const LineSegment3D &ROTATION_AXIS,
      const storage::List< CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_MOVING_POINTS
    ) const
    {
      // create variables to hold the coefficients a, b, and c. These will be summed up over the "TARGET_MOVING_POINTS"
      double coefficient_a_sum( 0);
      double coefficient_b_sum( 0);
      double coefficient_c_sum( 0);

      // iterate through "TARGET_MOVING_POINTS" in order to sum up "coefficient_a_sum", "coefficient_b_sum", and
      // "coefficient_c_sum"
      for
      (
        storage::List< CyclicCoordinateDescent::TargetAndMovingPointPair>::const_iterator
          itr( TARGET_MOVING_POINTS.Begin()), itr_end( TARGET_MOVING_POINTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // calculate the coefficient values for the current target and mvoing point pair
        const storage::VectorND< 3, double> current_a_b_c_coefficient_values
        (
          CalculateCurrentCoefficientValues( ROTATION_AXIS, *itr)
        );

        // add the coefficients of "current_a_b_c_coefficient_values" into "coefficient_a_sum", "coefficient_b_sum",
        // and "coefficient_c_sum"
        coefficient_a_sum += current_a_b_c_coefficient_values.First();
        coefficient_b_sum += current_a_b_c_coefficient_values.Second();
        coefficient_c_sum += current_a_b_c_coefficient_values.Third();
      }

      // message the coefficient sums
      BCL_MessageDbg( "coefficient_a_sum " + util::Format()( coefficient_a_sum));
      BCL_MessageDbg( "coefficient_b_sum " + util::Format()( coefficient_b_sum));
      BCL_MessageDbg( "coefficient_c_sum " + util::Format()( coefficient_c_sum));

      // calculate the cosine of alpha
      const double cos_alpha
      (
        coefficient_b_sum /
          math::Sqrt( math::Sqr( coefficient_b_sum) + math::Sqr( coefficient_c_sum))
      );

      // message the cosine of alpha
      BCL_MessageDbg( "cos_alpha " + util::Format()( cos_alpha));

      // calculate the sin of alpha
      const double sin_alpha
      (
        coefficient_c_sum /
          math::Sqrt( math::Sqr( coefficient_b_sum) + math::Sqr( coefficient_c_sum))
      );

      // message the sin of alpha
      BCL_MessageDbg( "sin_alpha " + util::Format()( sin_alpha));

      // calculate theta, this is the optimal rotation which will minimize the sum square distance between the moving
      // and target points
      const double theta( atan2( sin_alpha, cos_alpha));

      // message theta
      BCL_MessageDbg( "theta " + util::Format()( theta));

      // calculate the distance sum that will result after rotating around "ROTATION_AXIS" by "theta"
      const double distance_sum( CalculateDistanceSum( theta, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum));

      // message the resulting distance sum
      BCL_MessageDbg( "distance_sum_no_rot " + util::Format()( CalculateDistanceSum( 0.0, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum)));
      BCL_MessageDbg( "distance_sum " + util::Format()( distance_sum));

      // create a ResultsAndCoefficients
      const CyclicCoordinateDescent::ResultsAndCoefficients optimal_rotation_and_distance_and_coefficients
      (
        theta, distance_sum, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum
      );

      // return the results and coefficients
      return optimal_rotation_and_distance_and_coefficients;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::TargetAndMovingPointPair::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TargetPoint, ISTREAM);
      io::Serialize::Read( m_MovingPoint, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::TargetAndMovingPointPair::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // Write members
      io::Serialize::Write( m_TargetPoint, OSTREAM, INDENT);
      io::Serialize::Write( m_MovingPoint, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::ResultsAndCoefficients::Read( std::istream &ISTREAM)
    {
      // Read members
      io::Serialize::Read( m_OptimalRotation, ISTREAM);
      io::Serialize::Read( m_MinimumDistanceSum, ISTREAM);
      io::Serialize::Read( m_CoefficientA, ISTREAM);
      io::Serialize::Read( m_CoefficientB, ISTREAM);
      io::Serialize::Read( m_CoefficientC, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::ResultsAndCoefficients::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // Write members
      io::Serialize::Write( m_OptimalRotation, OSTREAM, INDENT);
      io::Serialize::Write( m_MinimumDistanceSum, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientA, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientB, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientC, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief CalculateDistanceSum calculates the sum distance difference between moving and target points
    //! The sum distance difference between the moving points and their corresponding fixed target points will result
    //! if the moving points are rotated by the given rotation angle around the rotation axis used when the
    //! coefficients were calculated.
    //! @param ROTATION the amount of rotation that should be considered
    //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
    double CyclicCoordinateDescent::CalculateDistanceSum
    (
      const double ROTATION, const double COEFFICIENT_A, const double COEFFICIENT_B, const double COEFFICIENT_C
    ) const
    {
      // calculate the distance sum as defined by Canutescu et al. Protein Sci, 2003
      const double distance_sum( COEFFICIENT_A - COEFFICIENT_B * cos( ROTATION) - COEFFICIENT_C * sin( ROTATION));

      // return distance sum
      return distance_sum;
    }

    //! @brief CalculateCurrentCoefficientValues calculates the contribution to the total coefficient values for a
    //! single target and moving point pair
    //! @param ROTATION_AXIS the axis around which rotation will occur
    //! @param CURRENT_COORDINATES the pair of target and moving points
    //! @return a VectorND< 3, double> which holds coefficients a, b, and c, respectively
    storage::VectorND< 3, double> CyclicCoordinateDescent::CalculateCurrentCoefficientValues
    (
      const LineSegment3D &ROTATION_AXIS,
      const CyclicCoordinateDescent::TargetAndMovingPointPair &CURRENT_COORDINATES
    ) const
    {
      // O, the footpoint of the moving atoms along the rotation axis
      const linal::Vector3D footpoint
      (
        linal::CalculateFootpoint
        (
          CURRENT_COORDINATES.GetMovingPoint(), ROTATION_AXIS.GetStartPoint(), ROTATION_AXIS.GetDirection()
        )
      );

      // message the footpoint coordinates
      BCL_MessageDbg( "footprint coordinates are " + util::Format()( footpoint));

      // o^, the unit vector of the rotation axis vector
      linal::Vector3D rotation_axis_unit_vector( ROTATION_AXIS.GetDirection());
      rotation_axis_unit_vector.Normalize();

      // message o^
      BCL_MessageDbg( "rotation_axis_unit_vector " + util::Format()( rotation_axis_unit_vector));

      // f, the distance from O to F
      const double footpoint_to_target_distance( linal::Distance( CURRENT_COORDINATES.GetTargetPoint(), footpoint));

      // message f
      BCL_MessageDbg
      (
        "footpoint_to_target_distance " + util::Format()( footpoint_to_target_distance)
      );

      // f->, the vector from O->F
      const linal::Vector3D footpoint_to_target_vector( CURRENT_COORDINATES.GetTargetPoint() - footpoint);

      // message f->
      BCL_MessageDbg( "footpoint_to_target_vector " + util::Format()( footpoint_to_target_vector));

      // f^, the unit vector of f->
      linal::Vector3D footpoint_to_target_unit_vector( footpoint_to_target_vector);
      footpoint_to_target_unit_vector.Normalize();

      // message f^
      BCL_MessageDbg
      (
        "footpoint_to_target_unit_vector " + util::Format()( footpoint_to_target_unit_vector)
      );

      // r, the distance from O to M
      const double footpoint_to_moving_distance( linal::Distance( CURRENT_COORDINATES.GetMovingPoint(), footpoint));

      // message r
      BCL_MessageDbg
      (
        "footpoint_to_moving_distance " + util::Format()( footpoint_to_moving_distance)
      );

      // r->, the vector from O->M
      const linal::Vector3D footpoint_to_moving_vector( CURRENT_COORDINATES.GetMovingPoint() - footpoint);

      // message r->
      BCL_MessageDbg( "footpoint_to_moving_vector " + util::Format()( footpoint_to_moving_vector));

      // r^, the unit vector of r->
      linal::Vector3D footpoint_to_moving_unit_vector( footpoint_to_moving_vector);
      footpoint_to_moving_unit_vector.Normalize();

      // message r^
      BCL_MessageDbg
      (
        "footpoint_to_moving_unit_vector " + util::Format()( footpoint_to_moving_unit_vector)
      );

      // s, the third orthogonal coordinate axis in the coordinate system defined by and o^, r^
      const linal::Vector3D third_coordinate_axis
      (
        linal::CrossProduct( footpoint_to_moving_unit_vector, rotation_axis_unit_vector)
      );

      // message s
      BCL_MessageDbg( "third_coordinate_axis " + util::Format()( third_coordinate_axis));

      // s^, the unit vector of s
      linal::Vector3D third_coordinate_axis_unit_vector( third_coordinate_axis);
      third_coordinate_axis_unit_vector.Normalize();

      // calculate coefficient a
      const double coefficient_a
      (
        math::Sqr( footpoint_to_moving_distance) + math::Sqr( footpoint_to_target_distance)
      );

      // calculate coefficient b
      const double coefficient_b
      (
        2 * footpoint_to_moving_distance * linal::ScalarProduct
        (
          footpoint_to_target_vector, footpoint_to_moving_unit_vector
        )
      );

      // calculate coefficient c
      const double coefficient_c
      (
        2 * footpoint_to_moving_distance * linal::ScalarProduct
        (
          footpoint_to_target_vector, third_coordinate_axis_unit_vector
        )
      );

      // message the three coefficients
      BCL_MessageDbg
      (
        "coefficient a\t" + util::Format()( coefficient_a) + "\tcoefficient b\t"
        + util::Format()( coefficient_b) + "\tcoefficient c\t" + util::Format()( coefficient_c)
      );

      // return the three coefficients
      return storage::VectorND< 3, double>( coefficient_a, coefficient_b, coefficient_c);
    }

    //! @brief operator== for defining the equality between two ResultsAndCoefficients objects
    //! @param RESULT_A the first ResultsAndCoefficients object
    //! @param RESULT_B the second ResultsAndCoefficients object
    //! @return boolean, true if each component of the two ResultsAndCoefficients objects is equal within tolerance
    bool operator==
    (
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_A,
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_B
    )
    {
      // return true if all values are equal within tolerance, false otherwise
      return math::EqualWithinTolerance( RESULT_A.GetOptimalRotation(), RESULT_B.GetOptimalRotation()) &&
        math::EqualWithinTolerance( RESULT_A.GetMinimumDistanceSum(), RESULT_B.GetMinimumDistanceSum()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientA(), RESULT_B.GetCoefficientA()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientB(), RESULT_B.GetCoefficientB()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientC(), RESULT_B.GetCoefficientC());
    }

    //! @brief operator== for defining the equality between two TargetAndMovingPointPair objects
    //! @param PAIR_A the first TargetAndMovingPointPair object
    //! @param PAIR_B the second TargetAndMovingPointPair object
    //! @return boolean, true if the coordinates of the two TargetAndMovingPointPair objects is equal within tolerance
    bool operator==
    (
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_A,
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_B
    )
    {
      // return true if the coordinates of the target and moving points in each TargetAndMovingPointPair is equal
      // within tolerance
      return math::EqualWithinTolerance( PAIR_A.GetTargetPoint(), PAIR_B.GetTargetPoint()) &&
        math::EqualWithinTolerance( PAIR_A.GetMovingPoint(), PAIR_B.GetMovingPoint());
    }

  } // namespace coord
} // namespace bcl
