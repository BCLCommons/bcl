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

#ifndef BCL_COORD_CYCLIC_COORDINATE_DESCENT_H_
#define BCL_COORD_CYCLIC_COORDINATE_DESCENT_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CyclicCoordinateDescent
    //! @brief Connects two points in space with a series of coordinates
    //! @details Please see Canutescu, A., Dunbrack, R. “Cyclic coordinate descent: A robotics algorithm for protein loop
    //! closure." Protein Science (2003), 12: 963-972
    //!
    //! Given the current location of a set of points, the desired coordinates for those points, and an axis around
    //! which a rotation will occur, this class can provide the optimal rotation around that axis which will minimize
    //! the sum distance difference between the points location after the rotation and the points desired location.
    //!
    //! @see @link example_coord_cyclic_coordinate_descent.cpp @endlink
    //! @author alexanns
    //! @date Aug 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CyclicCoordinateDescent :
      public util::ObjectInterface
    {

    public:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class TargetAndMovingPointPair
      //! @brief Holds a pair of points, the fixed target point and the movable point
      //! The movable point is the point which is being moved in an attempt to reach the target point.
      //!
      //! @see @link example_coord_cyclic_coordinate_descent.cpp @endlink
      //! @author alexanns
      //!
      //! @date Aug 24, 2010
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class BCL_API TargetAndMovingPointPair :
        public util::ObjectInterface
      {

      private:

      //////////
      // data //
      //////////

        //! linal::Vector3D "m_TargetPoint" is the fixed target point desired to be found
        linal::Vector3D m_TargetPoint;

        //! linal::Vector3D "m_MovingPoint" is the point which is being moved in an attempt to reach the target point
        linal::Vector3D m_MovingPoint;

      public:

      //////////
      // data //
      //////////

        //! single instance of that class
        static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! default constructor
        TargetAndMovingPointPair();

        //! constructor taking a target point and a moving point
        //! @param TARGET_POINT the fixed target point desired to be found
        //! @param MOVING_POINT the point which is being moved in an attempt to reach the target point
        TargetAndMovingPointPair( const linal::Vector3D &TARGET_POINT, const linal::Vector3D &MOVING_POINT);

        //! @brief Clone function
        //! @return pointer to new TargetAndMovingPointPair
        TargetAndMovingPointPair *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetClassIdentifier() const;

      ////////////////
      // operations //
      ////////////////

        //! @brief GetTargetPoint gives the fixed target point desired to be found
        //! @return returns linal::Vector3D which is "m_TargetPoint"
        const linal::Vector3D &GetTargetPoint() const;

        //! @brief GetMovingPoint gives the point which is being moved in an attempt to reach the target point
        //! @return returns linal::Vector3D which is "m_MovingPoint"
        const linal::Vector3D &GetMovingPoint() const;

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

      }; // class TargetAndMovingPointPair

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ResultsAndCoefficients
      //! @brief holds the results from CyclicCoordinateDescent::GetOptimalRotationandDistance
      //! The information is :
      //!   a.) the rotation will which minimize the distance between target and moving points
      //!   b.) the distance sum that will results between the target and moving points
      //!   c.) the three coefficients which can be used to calculate the distance sum, S, for any rotation around the
      //!       rotation axis used to compute the coefficients. The function
      //!       CyclicCoordinateDescent::CalculateDistanceSum can be used to calculate S.
      //!       For information about the coefficients, please see the reference given above
      //!       (Canutescu et al. Protein Sci, 2003).
      //!
      //! @see @link example_coord_cyclic_coordinate_descent.cpp @endlink
      //! @author alexanns
      //!
      //! @date Aug 24, 2010
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class BCL_API ResultsAndCoefficients :
        public util::ObjectInterface
      {

      private:

      //////////
      // data //
      //////////

        //! double "m_OptimalRotation" the rotation to minimize the distance sum between target and moving points
        double m_OptimalRotation;

        //! double "m_MinimumDistanceSum" the sum of distances between target and moving points that will result after
        //! the moving points are rotated by "m_OptimalRotation"
        double m_MinimumDistanceSum;

        //! double "m_CoefficientA" please see the reference given above (Canutescu et al. Protein Sci, 2003)
        double m_CoefficientA;

        //! double "m_CoefficientB" please see the reference given above (Canutescu et al. Protein Sci, 2003)
        double m_CoefficientB;

        //! double "m_CoefficientC" please see the reference given above (Canutescu et al. Protein Sci, 2003)
        double m_CoefficientC;

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
        ResultsAndCoefficients();

        //! @brief constructor taking all of the member variables
        //! @param OPTIMAL_ROTATION the rotation to minimize the distance sum between target and moving points
        //! @param MINIMUM_DISTANCE_SUM the sum of distances between target and moving points that will result after
        //!        the moving points are rotated by "m_OptimalRotation"
        //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
        //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
        //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
        ResultsAndCoefficients
        (
          const double OPTIMAL_ROTATION, const double MINIMUM_DISTANCE_SUM, const double COEFFICIENT_A,
          const double COEFFICIENT_B, const double COEFFICIENT_C
        );

        //! @brief Clone function
        //! @return pointer to new ResultsAndCoefficients
        ResultsAndCoefficients *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetClassIdentifier() const;

      ////////////////
      // operations //
      ////////////////

        //! @brief GetOptimalRotation gives the rotation to minimize the distance sum between target and moving points
        //! @return returns "m_OptimalRotation"
        double GetOptimalRotation() const;

        //! @brief GetMinimumDistanceSum gives the sum of distances between target and moving points that will result
        //!        after the moving points are rotated by "m_OptimalRotation"
        //! @return returns "m_MinimumDistanceSum"
        double GetMinimumDistanceSum() const;

        //! @brief GetCoefficientA gives "m_CoefficientA"
        //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
        //! @return returns "m_CoefficientA"
        double GetCoefficientA() const;

        //! @brief GetCoefficientB gives "GetCoefficientB"
        //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
        //! @return returns "GetCoefficientB"
        double GetCoefficientB() const;

        //! @brief GetCoefficientC gives "GetCoefficientC"
        //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
        //! @return returns "GetCoefficientC"
        double GetCoefficientC() const;

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

      }; // class ResultsAndCoefficients

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CyclicCoordinateDescent();

      //! @brief Clone function
      //! @return pointer to new CyclicCoordinateDescent
      CyclicCoordinateDescent *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetOptimalRotationandDistance gives rotation to minimize the distance between target and moving points
      //! It calculates the rotation around "ROTATION_AXIS" that needs to be performed on all moving points so that
      //! after the rotation, the sum distance from each of the moving points to corresponding fixed points will be
      //! minimized
      //! @param ROTATION_AXIS this is the axis the rotation will be performed around
      //! @param TARGET_MOVING_POINTS the list of moving points and the corresponding fixed target points the moving
      //!        points are trying to reach
      //! @return returns a ResultsAndCoefficients which holds the results of the calculation
      ResultsAndCoefficients GetOptimalRotationandDistance
      (
        const LineSegment3D &ROTATION_AXIS,
        const storage::List< TargetAndMovingPointPair> &TARGET_MOVING_POINTS
      ) const;

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

    public:

      //! @brief CalculateDistanceSum calculates the sum distance difference between moving and target points
      //! The sum distance difference between the moving points and their corresponding fixed target points will result
      //! if the moving points are rotated by the given rotation angle around the rotation axis used when the
      //! coefficients were calculated.
      //! @param ROTATION the amount of rotation that should be considered
      //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
      //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
      //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
      double CalculateDistanceSum
      (
        const double ROTATION, const double COEFFICIENT_A, const double COEFFICIENT_B, const double COEFFICIENT_C
      ) const;

    protected:

      //! @brief CalculateCurrentCoefficientValues calculates the contribution to the total coefficient values for a
      //! single target and moving point pair
      //! @param ROTATION_AXIS the axis around which rotation will occur
      //! @param CURRENT_COORDINATES the pair of target and moving points
      //! @return a VectorND< 3, double> which holds coefficients a, b, and c, respectively
      storage::VectorND< 3, double> CalculateCurrentCoefficientValues
      (
        const LineSegment3D &ROTATION_AXIS,
        const CyclicCoordinateDescent::TargetAndMovingPointPair &CURRENT_COORDINATES
      ) const;

    }; // class CyclicCoordinateDescent

    //! @brief operator== for defining the equality between two ResultsAndCoefficients objects
    //! @param RESULT_A the first ResultsAndCoefficients object
    //! @param RESULT_B the second ResultsAndCoefficients object
    //! @return boolean, true if each component of the two ResultsAndCoefficients objects is equal within tolerance
    BCL_API bool operator==
    (
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_A,
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_B
    );

    //! @brief operator== for defining the equality between two TargetAndMovingPointPair objects
    //! @param PAIR_A the first TargetAndMovingPointPair object
    //! @param PAIR_B the second TargetAndMovingPointPair object
    //! @return boolean, true if the coordinates of the two TargetAndMovingPointPair objects is equal within tolerance
    BCL_API bool operator==
    (
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_A,
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_B
    );

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_CYCLIC_COORDINATE_DESCENT_H_ 
