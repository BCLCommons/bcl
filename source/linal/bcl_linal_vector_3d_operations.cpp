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
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //! @brief calculates the absolute difference between two linal::Vector3Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector3D (points)
    Vector3D AbsoluteDifference( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      Vector3D returned_vector( VECTOR_B - VECTOR_A);
      returned_vector.X() = math::Absolute( returned_vector.X());
      returned_vector.Y() = math::Absolute( returned_vector.Y());
      returned_vector.Z() = math::Absolute( returned_vector.Z());
      return returned_vector;
    }

    //! @brief calculates the unit vector starting from one linal::Vector3D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector3D UnitVector( const Vector3D &ORIGIN, const Vector3D &TARGET)
    {
      Vector3D unit_vector( TARGET);
      unit_vector -= ORIGIN;
      unit_vector.Normalize();
      return unit_vector;
    }

    //! @brief dihedral angle between four points (A->B -x-> C->D)
    //! @brief see http://en.wikipedia.org/wiki/Dihedral_angle for reference
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return dihedral angle between four points
    double Dihedral
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D
    )
    {
      // calculate the two cross products (b1xb2 and b2xb3)
      const Vector3D cross_b1_b2( CrossProduct( VECTOR_A, VECTOR_B, VECTOR_B, VECTOR_C));
      const Vector3D cross_b2_b3( CrossProduct( VECTOR_B, VECTOR_C, VECTOR_C, VECTOR_D));

      // calculate the vectors b1 and b2
      const Vector3D b1( VECTOR_B - VECTOR_A);

      // get the distance b -c
      const double distance( Distance( VECTOR_C, VECTOR_B));

      // calculate dihedral angle
      const double dihedral( std::atan2( distance * b1 * cross_b2_b3, cross_b1_b2 * cross_b2_b3));

      return util::IsDefined( dihedral) ? dihedral : 0.0;
    }

    //! @brief calculates coordinates using dihedral angle information (point X in C->B->A->X)
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
    //! @param DIHEDRAL_XABC dihedral angle between X, VECTOR_A, VECTOR_B, and VECTOR_C
    //! @return point X in C->B->A->X
    Vector3D CoordinatesDihedral
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const double DISTANCE_XA,
      const double ANGLE_XAB,
      const double DIHEDRAL_XABC
    )
    {
      const Vector3D a( UnitVector( VECTOR_B, VECTOR_A));
      const Vector3D b( UnitVector( VECTOR_C, VECTOR_B));
      const Vector3D c( CrossProduct( a, b).Normalize());
      const Vector3D d( CrossProduct( a, c).Normalize());
      const Vector3D x
      (
        DISTANCE_XA *
        (
          a * cos( math::g_Pi - ANGLE_XAB) - c * sin( math::g_Pi - ANGLE_XAB) * sin( DIHEDRAL_XABC) +
          d * sin( math::g_Pi - ANGLE_XAB) * cos( DIHEDRAL_XABC)
        )
      );

      return VECTOR_A + x;
    }

    //! @brief calculates coordinates using angle information (point X in B,C->A->X)
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
    //! @param ANGLE_XAC angle between X, VECTOR_A, and VECTOR_C
    //! @param SIDE true if on a side
    //! @param VECTOR_SIDE vector of the side
    //! @return point X in B,C->A->X
    Vector3D CoordinatesAngle
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const double DISTANCE_XA,
      const double ANGLE_XAB,
      const double ANGLE_XAC,
      const bool SIDE,
      const Vector3D &VECTOR_SIDE
    )
    {
      // in plane components and direction
      const Vector3D a( UnitVector( VECTOR_B, VECTOR_A));
      Vector3D b( UnitVector( VECTOR_C, VECTOR_A));

      // determine whether a and b are colinear. In which case, use a different vector C
      while( math::Absolute( ProjAngleCosinus( a, b)) >= 0.99)
      {
        BCL_MessageDbg( "Lines were colinear, CoordinatesAngle returning non-unique vector");
        b.Rotate( math::RotationMatrix3D().SetRand( math::g_Pi));
      }

      const Vector3D c( cos( math::g_Pi - ANGLE_XAB) * a);
      const Vector3D d( cos( math::g_Pi - ANGLE_XAC) * b);

      const double bac_angl( ProjAngle( c, d));
      const double a_dist( c.Norm());
      const double b_dist( d.Norm());
      const double c_dist
      (
        math::Sqrt
        (
          std::max
          (
            0.0,
            ( math::Sqr( a_dist) + math::Sqr( b_dist) - 2 * a_dist * b_dist * cos( bac_angl))
            / math::Sqr( sin( bac_angl)) -
            math::Sqr( a_dist)
          )
        )
      );

      //in a,b plane
      const Vector3D e( CrossProduct( c, d).Normalize());
      const Vector3D f( CrossProduct( e, c).Normalize() * c_dist);
      const Vector3D g( c + f);

      //orthogonal to a,b plane
      Vector3D h( e * math::Sqrt( std::max( 0.0, 1.0 - g.SquareNorm())));

      //side
      if( util::IsDefined( VECTOR_SIDE( 0)))
      {
        const Vector3D i( VECTOR_SIDE - VECTOR_A);
        const double dihe_1( Dihedral( i, Vector3D(), a, b));
        const double dihe_2( Dihedral( e, Vector3D(), a, b));
        if( dihe_2 && ( ( dihe_1 / dihe_2 > 0.0) != SIDE))
        {
          h *= ( -1.0);
        }
      }

      //sum and length
      const Vector3D x( Vector3D( g + h).Normalize() * DISTANCE_XA);

      return VECTOR_A + x;
    }

  } // namespace linal
} // namespace bcl
