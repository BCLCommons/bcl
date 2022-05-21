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

#ifndef BCL_COORD_H_
#define BCL_COORD_H_

// include the namespace forward header
#include "bcl_coord.fwd.hh"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_coord.h
  //! @brief namespace implementing classes related to coordinate handling - mainly 3D cartesian coordinates
  //! @details biological and chemical entities contain atoms with 3D cartesian coordinates. These can be transformed in
  //!          in space with coordinate transformation which is mainly matrix and vector algebra.
  //!          Additionally, sets of coordinates often need to be superimposed or they need to be compared for similarity.
  //!
  //! @see @link example_coord.cpp @endlink
  //! @author woetzen
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace coord
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //! @brief checks that all coordinates given in the vector COORDINATES are defined
    //! @param COORDINATES vector of coordinates of interest
    //! @return whether all coordinates given in the vector COORDINATES are defined
    BCL_API
    bool AreDefinedCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    );

    //! @brief calculates the Cartesian coordinates from the distance (relative to Zero point) and two angles
    //! @param DISTANCE distance to origin
    //! @param PHI first angle
    //! @param PSI second angle
    //! @return the Cartesian coordinates for given distance and two angles
    linal::Vector3D SphericalCoordinates( const double DISTANCE, const double PHI, const double PSI);

    //! @brief transforms the given vector of coordinates by given TransformationMatrix3D
    //! @param COORDINATES Vector of coordinates of interest
    //! @param TRANSFORMATION_MATRIX TransformationMatrix3D to be used
    BCL_API
    void TransformCoordinates
    (
      util::SiPtrVector< linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
    );

    //! @brief calculate a distance matrix from given vector of coordinates
    //! @param COORDINATES Vector of coordinates of interest
    //! @return distance matrix for given coordinates
    BCL_API
    linal::Matrix< double> CalculateDistanceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    );

    //! @brief calculate difference matrix between two coordinate vectors
    //! @param COORDINATES_A first vector of coordinates of interest
    //! @param COORDINATES_B second vector of coordinates of interest
    //! @return difference matrix between given coordinate vectors
    BCL_API
    linal::Matrix< double> CalculateDifferenceDistanceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    );

    //! computes center of mass of pointset
    BCL_API
    linal::Vector3D CenterOfMass
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const bool SKIP_UNDEFINED = false
    );

    //! compute the square of the radius of gyration for a set of coordinates
    BCL_API
    double SquareRadiusOfGyration( const util::SiPtrVector< const linal::Vector3D> &COORDINATES);

    //! compute the radius of gyration for a set of coordinates
    BCL_API
    double RadiusOfGyration( const util::SiPtrVector< const linal::Vector3D> &COORDINATES);

    //! computes the symmetry of a set of coordinates
    BCL_API
    double SymmetryFactor( const util::ShPtrVector< util::SiPtrVector< const linal::Vector3D> > &COORDINATES);

    //! NeighborWeight - assign each amino amid a weight (0 <= weight <= 1) based on its distance from the amino acid of interest
    BCL_API
    double NeighborWeight( const double DISTANCE, const double LOW_THRESHOLD, const double HIGH_THRESHOLD);

    //! CountNeighbors - calculates the number of neighbors that the current point has
    BCL_API
    double CountNeighbors
    (
      const util::SiPtrVector< const linal::Vector3D> &ALL_POINTS,
      const linal::Vector3D &CURRENT_POINT,
      const double LOW_THRESHOLD,
      const double HIGH_THRESHOLD
    );

    //! NeighborVector - calculates the sum of all vectors between the current point and its neighboring points
    BCL_API
    linal::Vector3D NeighborVector
    (
      const util::SiPtrVector< const linal::Vector3D> &ALL_POINTS,
      const linal::Vector3D &CURRENT_POINT,
      const double LOW_THRESHOLD,
      const double HIGH_THRESHOLD,
      const bool NORMALIZE_BY_NEIGHBOR_COUNT = true
    );

    //! @brief calculate the density of a given set of points
    //! @param COORDINATES the coordinates to be used
    //! @param RESOLUTION the resolution the density is determined at in x, y and z
    //! @return the density in length_unit^3
    BCL_API
    double CalculatePointDensity
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const linal::Vector3D &RESOLUTION
    );

    //! @brief quantize a given set of points to a certain resolution
    //! @param COORDINATES the coordinates to be used
    //! @param RESOLUTION the resolution the points are quantized to in x,y and z
    //! @return a set of 3d indices to a grid of the size RESOLUTION
    BCL_API
    storage::Set< storage::Triplet< int, int, int> > QuantizePoints
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const linal::Vector3D &RESOLUTION
    );

    //! @brief Estimate density of points using 2D convex hulls
    //! @param COORDINATES the coordinates to be used
    //! @param SLICE_WIDTH typical distance between adjacent points
    //! @param MAX_DISTANCE maximum distance between points; if above this, induce a concavity into the convex hull
    //! @param RADIUS radius of each point
    BCL_API double EstimateVolume
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const double &SLICE_WIDTH,
      const double &MAX_DISTANCE,
      const double &RADIUS = 0.0
    );

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_H_
