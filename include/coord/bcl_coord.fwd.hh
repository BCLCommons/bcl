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

#ifndef BCL_COORD_FWD_HH_
#define BCL_COORD_FWD_HH_

// include the dependency file for this header
#include "bcl_coord.depends.fwd.hh"

// This file contains forward declarations for the coord namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace coord
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class Axes;
    class CyclicCoordinateDescent;
    class CylinderCoordinates;
    class GeometricHashStorageClasses;
    class GeometricHashStorageHashMap;
    class GeometricHashStorageInterface;
    class GeometricHashing;
    class GeometryInterface;
    class LineSegment2D;
    class LineSegment3D;
    class MomentOfInertia;
    class MovableEccentric;
    class MovableInterface;
    class MoveCombine;
    class MoveInterface;
    class MoveRotateDefined;
    class MoveRotateRandom;
    class MoveRotateRandomExternalReference;
    class MoveTransformRandom;
    class MoveTranslateDefined;
    class MoveTranslateExternalAxis;
    class MoveTranslateRandom;
    class OrientationInterface;
    class PointCloud;
    class PointToKeyClasses;
    class PointToKeyInterface;
    class PointToKeySphericalRadius;
    class Polygon;
    class Sphere;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class TransformationMatrix3D;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< linal::Vector3D, Axes>                                                    Axis;
    typedef util::Enum< util::ShPtr< GeometricHashStorageInterface>, GeometricHashStorageClasses> GeometricHashStorageClass;
    typedef util::Enum< util::ShPtr< PointToKeyInterface>, PointToKeyClasses>                     PointToKey;

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_FWD_HH_
