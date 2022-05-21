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

#ifndef BCL_COORD_GEOMETRY_INTERFACE_H_
#define BCL_COORD_GEOMETRY_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_line_segment_3d.h"
#include "bcl_coord_movable_interface.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeometryInterface
    //! @brief class that represents shapes in space.  Contains a vector of GeometryInterfaces that
    //!        represent fragments of the complete shape.  In this manner, SSEs can be represented as collections of
    //!        cylinders or boxes to better represent bent SSEs.
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Mar 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeometryInterface :
      public virtual MovableInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns SiPtrVector of GeometryInterfaces that make up this GeometryInterface
      //! @return SiPtrVector of GeometryInterfaces that make up this GeometryInterface
      virtual util::SiPtrVector< const GeometryInterface> GetGeometries() const = 0;

      //! @brief returns the main axis as a LineSegment3D
      //! @return the main axis as a LineSegment3D
      virtual LineSegment3D GetMainAxis() const = 0;

      //! @brief returns the x, y and z extents
      //! @return the requested extents as linal::Vector3D
      virtual linal::Vector3D GetExtents() const;

      //! @brief returns the requested extent
      //! @param AXIS axis of interest
      //! @return the requested extent
      virtual double GetExtent( const Axis &AXIS) const = 0;

      //! @brief returns the radial extent
      //! @return the radial extent
      virtual double GetRadialExtent() const = 0;

      //! @brief returns the orientation of the geometry
      //! @return the orientation of the geometry
      virtual const math::TransformationMatrix3D GetOrientation() const = 0;

    }; // class GeometryInterface

    //! @brief boolean operator GEOMETRY_LHS == GEOMETRY_RHS
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS
    BCL_API bool operator ==( const GeometryInterface &GEOMETRY_LHS, const GeometryInterface &GEOMETRY_RHS);

    //! @brief function to check whether two GEOMETRYs are equal within tolerance
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @param RELATIVE_TOLERANCE relative tolerance
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS within tolerance
    BCL_API bool EqualWithinTolerance
    (
      const GeometryInterface &GEOMETRY_LHS,
      const GeometryInterface &GEOMETRY_RHS,
      const double &RELATIVE_TOLERANCE = 0.001
    );

    //! @brief calculate the shortest connecting linesegment between two geometries
    //! @param GEOMETRY_A first geometry
    //! @param GEOMETRY_B second geometry
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
    //! @return the shortest connecting linesegment between two geometries
    BCL_API storage::Pair< LineSegment3D, bool> ShortestConnectionBetweenGeometries
    (
      const GeometryInterface &GEOMETRY_A, const GeometryInterface &GEOMETRY_B, const double MINIMAL_INTERFACE_LENGTH
    );

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_GEOMETRY_INTERFACE_H_ 
