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
#include "coord/bcl_coord_geometry_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    //! @brief boolean operator GEOMETRY_LHS == GEOMETRY_RHS
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS
    bool operator ==( const GeometryInterface &GEOMETRY_LHS, const GeometryInterface &GEOMETRY_RHS)
    {
      return
      (
        GEOMETRY_LHS.GetOrientation() == GEOMETRY_RHS.GetOrientation() &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetMainAxis().GetLength(), GEOMETRY_RHS.GetMainAxis().GetLength()) &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetExtent( GetAxes().e_X), GEOMETRY_RHS.GetExtent( GetAxes().e_X)) &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetExtent( GetAxes().e_Y), GEOMETRY_RHS.GetExtent( GetAxes().e_Y)) &&
        GEOMETRY_LHS.GetGeometries().GetSize() == GEOMETRY_RHS.GetGeometries().GetSize()
      );
    }

    //! @brief function to check whether two GEOMETRYs are equal within tolerance
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS within tolerance
    bool EqualWithinTolerance
    (
      const GeometryInterface &GEOMETRY_LHS,
      const GeometryInterface &GEOMETRY_RHS,
      const double &RELATIVE_TOLERANCE
    )
    {
      return
      (
        math::SimilarWithinTolerance
        (
          GEOMETRY_LHS.GetOrientation(), GEOMETRY_RHS.GetOrientation(), 0.001, 0.001
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetMainAxis().GetLength(), GEOMETRY_RHS.GetMainAxis().GetLength(), RELATIVE_TOLERANCE
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetExtent( GetAxes().e_X),
          GEOMETRY_RHS.GetExtent( GetAxes().e_X), RELATIVE_TOLERANCE
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetExtent( GetAxes().e_Y),
          GEOMETRY_RHS.GetExtent( GetAxes().e_Y), RELATIVE_TOLERANCE
        ) &&
        GEOMETRY_LHS.GetGeometries().GetSize() == GEOMETRY_RHS.GetGeometries().GetSize()
      );
    }

    //! @brief calculate the shortest connecting linesegment between two geometries
    //! @param GEOMETRY_A first geometry
    //! @param GEOMETRY_B second geometry
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
    //! @return the shortest connecting linesegment between two geometries
    storage::Pair< LineSegment3D, bool> ShortestConnectionBetweenGeometries
    (
      const GeometryInterface &GEOMETRY_A, const GeometryInterface &GEOMETRY_B, const double MINIMAL_INTERFACE_LENGTH
    )
    {
      // make copies of the main axes
      LineSegment3D geometry_a_line( GEOMETRY_A.GetMainAxis());
      LineSegment3D geometry_b_line( GEOMETRY_B.GetMainAxis());

      // shorten them
      geometry_a_line.Shorten( MINIMAL_INTERFACE_LENGTH);
      geometry_b_line.Shorten( MINIMAL_INTERFACE_LENGTH);

      // end
      return ShortestConnectionBetweenLineSegments3D( geometry_a_line, geometry_b_line);
    }

    //! @brief returns the x, y and z extents
    //! @return the requested extents as linal::Vector3D
    linal::Vector3D GeometryInterface::GetExtents() const
    {
      return
        linal::Vector3D
        (
          GetExtent( GetAxes().e_X),
          GetExtent( GetAxes().e_Y),
          GetExtent( GetAxes().e_Z)
        );
    }

  } // namespace coord
} // namespace bcl
