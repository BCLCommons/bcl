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

#ifndef BCL_LINAL_VECTOR_3D_OPERATIONS_H_
#define BCL_LINAL_VECTOR_3D_OPERATIONS_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_3d.h"
#include "bcl_linal_vector_operations.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_linal_vector_3d_operations.h
  //! @brief mathematical operators and angle operations for 3d vectors
  //!
  //! @see @link example_linal_vector_3d_operations.cpp @endlink
  //! @author mendenjl
  //! @date Feb 20, 2012
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace linal
  {

  ///////////////
  // operators //
  ///////////////

    //! @brief operator +Vector 3D
    //! @param VECTOR original vector
    //! @return original vector
    inline
    Vector3D operator +( const Vector3D &VECTOR)
    {
      return VECTOR;
    }

    //! @brief operator Vector3D + Vector3D
    //! @param VECTOR_A first vector to be added
    //! @param VECTOR_B second vector to be added
    //! @return sum of VECTOR_A and VECTOR_B
    inline
    Vector3D
    operator +( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      Vector3D sum( VECTOR_A);
      sum += VECTOR_B;
      return sum;
    }

    //! @brief operator Vector3D - Vector3D
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector to be subtracted
    //! @return difference of VECTOR_A and VECTOR_B
    inline
    Vector3D operator -( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      Vector3D diff( VECTOR_A);
      diff -= VECTOR_B;
      return diff;
    }

    //! @brief operator Vector3D + scalar
    //! @param VECTOR vector to be added
    //! @param X scalar to be added
    //! @return sum of VECTOR and X
    inline
    Vector3D
    operator +( const Vector3D &VECTOR, const double &X)
    {
      Vector3D sum( VECTOR);
      sum += X;
      return sum;
    }

    //! @brief operator Vector3D - scalar
    //! @param VECTOR vector
    //! @param X scalar to be subtracted
    //! @return difference of VECTOR and X
    inline
    Vector3D
    operator -( const Vector3D &VECTOR, const double &X)
    {
      Vector3D sum( VECTOR);
      sum -= X;
      return sum;
    }

    //! @brief operator scalar + Vector3D
    //! @param X scalar to be added
    //! @param VECTOR vector to be added
    //! @return sum of X and VECTOR
    inline
    Vector3D operator +( const double &X, const Vector3D &VECTOR)
    {
      Vector3D sum( VECTOR);
      sum += X;
      return sum;
    }

    //! @brief operator * (dot product)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A and VECTOR_B
    inline
    double operator *( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return ( ( VECTOR_A.X() * VECTOR_B.X()) + ( VECTOR_A.Y() * VECTOR_B.Y()) + ( VECTOR_A.Z() * VECTOR_B.Z()));
    }

    //! @brief operator scalar * Vector3D
    //! @param X scalar to multiply
    //! @param VECTOR to multiply
    //! @return X times VECTOR
    inline
    Vector3D operator *( const double &X, const Vector3D &VECTOR)
    {
      Vector3D scaled( VECTOR);
      scaled *= X;
      return scaled;
    }

    //! @brief operator Vector3D * scalar
    //! @param VECTOR to multiply
    //! @param X scalar to multiply
    //! @return VECTOR times X
    inline
    Vector3D operator *( const Vector3D &VECTOR, const double &X)
    {
      Vector3D scaled( VECTOR);
      scaled *= X;
      return scaled;
    }

    //! @brief operator Vector3D / scalar
    //! @param VECTOR to divide
    //! @param X scalar to divide by
    //! @return X times VECTOR
    inline
    Vector3D operator /( const Vector3D &VECTOR, const double &X)
    {
      Vector3D scaled( VECTOR);
      scaled /= X;
      return scaled;
    }

    //! @brief operator -Vector3D
    //! @param VECTOR original vector
    //! @return additive inverse of VECTOR
    inline
    Vector3D operator -( const Vector3D &VECTOR)
    {
      return Vector3D( -VECTOR.X(), -VECTOR.Y(), -VECTOR.Z());
    }

    //! @brief operator scalar - Vector3D
    //! @param X to subtract from
    //! @param VECTOR to subtract by
    //! @return X minus VECTOR
    inline
    Vector3D operator -( const double &X, const Vector3D &VECTOR)
    {
      Vector3D scalar_minus_vector( -VECTOR);
      scalar_minus_vector += X;
      return scalar_minus_vector;
    }

    //! @brief operator ==
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @return true if elements from VECTOR_A == respective elements from VECTOR_B
    inline
    bool operator ==( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return VECTOR_A.X() == VECTOR_B.X() && VECTOR_A.Y() == VECTOR_B.Y() && VECTOR_A.Z() == VECTOR_B.Z();
    }

    //! @brief operator !=
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @return true if any element from VECTOR_A != respective element from VECTOR_B
    inline
    bool operator !=( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return !( VECTOR_A == VECTOR_B);
    }

  //////////////////////
  // Vector functions //
  //////////////////////

    //! @brief calculates the absolute difference between two linal::Vector3Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector3D (points)
    Vector3D AbsoluteDifference( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B);

    //! @brief calculates the unit vector directed from one linal::Vector3D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector3D UnitVector( const Vector3D &ORIGIN, const Vector3D &TARGET);

    //! @brief compares vectors like operator ==, but returns true if difference between both is smaller than EPSILON
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @param RELATIVE_TOLERANCE relative tolerance for comparison
    //! @param ABSOLUTE_TOLERANCE absolute tolerance for comparison
    //! @return true if difference between both is smaller than EPSILON
    inline
    bool EqualWithinTolerance
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const double &RELATIVE_TOLERANCE = 0.001,
      const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon()
    )
    {
      return
        math::EqualWithinTolerance( VECTOR_A.X(), VECTOR_B.X(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
        && math::EqualWithinTolerance( VECTOR_A.Y(), VECTOR_B.Y(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
        && math::EqualWithinTolerance( VECTOR_A.Z(), VECTOR_B.Z(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);
    }

    //! @brief calculates square distance between two linal::Vector3D (points)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return square distance between two linal::Vector3D (points)
    inline
    double SquareDistance( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return math::Sqr( VECTOR_A.X() - VECTOR_B.X())
             + math::Sqr( VECTOR_A.Y() - VECTOR_B.Y())
             + math::Sqr( VECTOR_A.Z() - VECTOR_B.Z());
    }

    //! @brief calculates distance between two linal::Vector3D (points)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return distance between two linal::Vector3D (points)
    inline
    double Distance( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return math::Sqrt( SquareDistance( VECTOR_A, VECTOR_B));
    }

    //! @brief calculate cosine of projection angle
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @return cosine of projection angle between two linal::Vector3D
    inline
    double ProjAngleCosinus( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      const double projection_angle_cosinus( ( VECTOR_A * VECTOR_B) / ( VECTOR_A.Norm() * VECTOR_B.Norm()));
      return projection_angle_cosinus;
    }

    //! @brief calculates cosine of projection angle from three points (A->B and A->C)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @return cosine of projection angle from three points
    inline
    double ProjAngleCosinus( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B, const Vector3D &VECTOR_C)
    {
      return ProjAngleCosinus( VECTOR_B - VECTOR_A, VECTOR_C - VECTOR_A);
    }

    //! @brief calculates cosine of projection angle from four points (A->B and C->D)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return cosine of projection angle from four points
    inline
    double ProjAngleCosinus
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D
    )
    {
      return ProjAngleCosinus( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! @brief calculates projection angle between two linal::Vector3D
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @return projection angle between two linal::Vector3D
    inline
    double ProjAngle( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      const double projection_angle_cosinus( ProjAngleCosinus( VECTOR_A, VECTOR_B));

      // check that the result is defined; this is more robust than checking if the vectors are defined,
      // since even defined vectors may have 0 norms, leading to undefined cosine
      if( !util::IsDefined( projection_angle_cosinus))
      {
        return 0.0; //util::GetUndefined< double>();
      }

      // through numerical drift it could be possible that the value is slightly higher or lower than -1 or 1
      // (VECTOR_A * VECTOR_B) / ( Norm( VECTOR_A) * Norm( VECTOR_B)) is the actual cos angle between vectors ab and cd
      if( projection_angle_cosinus >= double( 1.0))
      {
        // acos( 1) == 0.0
        return 0.0;
      }
      else if( projection_angle_cosinus <= double( -1.0))
      {
        // acos( -1) == pi
        return math::g_Pi;
      }

      // -1 < projection_angle_cosinus < 1, so calling acos( projection_angle_cosinus) yields the proper angle
      return std::acos( projection_angle_cosinus);
    }

    //! @brief calculates projection angle from three points (A->B and A->C)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @return projection angle from three points
    inline
    double ProjAngle( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B, const Vector3D &VECTOR_C)
    {
      return ProjAngle( VECTOR_B - VECTOR_A, VECTOR_C - VECTOR_A);
    }

    //! @brief calculates projection angle from four points (A->B and C->D)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return projection angle from four points
    inline
    double ProjAngle
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D
    )
    {
      return ProjAngle( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! @brief scalar (dot) product of two linal::Vector3Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return scalar (dot) product of VECTOR_A and VECTOR_B
    inline
    double ScalarProduct( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      return ( VECTOR_A * VECTOR_B);
    }

    //! @brief cross product of two linal::Vector3Ds
    //! @param A first vector
    //! @param B second vector
    //! @return cross product of VECTOR_A and VECTOR_B
    inline
    Vector3D CrossProduct( const Vector3D &A, const Vector3D &B)
    {
      return
        Vector3D
        (
          + A.Y() * B.Z() - A.Z() * B.Y(),
          - A.X() * B.Z() + A.Z() * B.X(),
          + A.X() * B.Y() - A.Y() * B.X()
        );
    }

    //! @brief cross product of three points (A->B and A->C)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @return cross product of three points
    inline
    Vector3D CrossProduct
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C
    )
    {
      return CrossProduct( VECTOR_B - VECTOR_A, VECTOR_C - VECTOR_A);
    }

    //! @brief cross product of four points (A->B and C->D)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return cross product of four points
    inline
    Vector3D CrossProduct
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D
    )
    {
      return CrossProduct( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
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
    );

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
    );

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
      const bool SIDE = true,
      const Vector3D &VECTOR_SIDE = Vector3D( util::GetUndefined< double>())
    );

    //the next three methods assume that all atoms/coordinates are in the same simplex: line, plane, tetrahedral
    //the new atom is always connected to the FIRST atom/coordinate

    //! point X in B->A->X where A, B and X are on the same line
    //! @brief calculates point X in B->A->X where A, B and X are on the same line
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @return point X in B->A->X where A, B and X are on the same line
    inline
    Vector3D CoordinatesLinear
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const double DISTANCE_XA
    )
    {
      // get the unit vector directed from B to A
      Vector3D x( UnitVector( VECTOR_B, VECTOR_A));
      // extend it the desired distance
      x *= DISTANCE_XA;
      // add it to point A
      x += VECTOR_A;
      // return the calculated point
      return x;
    }

    //! @brief calculates point X in C,B->A->X
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @return point X in C,B->A->X
    inline
    Vector3D CoordinatesTrigonal
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const double DISTANCE_XA
    )
    {
      // get the average unit vector between B->A and C->A
      Vector3D x( UnitVector( VECTOR_B, VECTOR_A));
      x += UnitVector( VECTOR_C, VECTOR_A);
      x.Normalize();

      // extend the unit vector to the desired length
      x *= DISTANCE_XA;

      // add it to vector A ( the offset)
      x += VECTOR_A;
      return x;
    }

    //! point X in B,C,D->A->X
    //! @brief calculates point X in B,C,D->A->X
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param VECTOR_D fourth point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @return point X in B,C,D->A->X
    inline
    Vector3D CoordinatesTetrahedral
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D,
      const double DISTANCE_XA
    )
    {
      // try the correct way first. This only fails in unusual cases where BC X CD == 0
      Vector3D x( CrossProduct( UnitVector( VECTOR_B, VECTOR_C), UnitVector( VECTOR_D, VECTOR_C)));
      if( x.SquareNorm())
      {
        x.Normalize();
        x *= DISTANCE_XA;
        // next, we need to add to vector A.
        Vector3D ab( VECTOR_B - VECTOR_A);
        ab.Normalize();
        if( ( x * ab) > ( -x * ab))
        {
          return VECTOR_A - x;
        }
        return VECTOR_A + x;
      }

      // get the average unit vector between B->A and C->A, D->A
      Vector3D av( UnitVector( VECTOR_B, VECTOR_A));
      av += UnitVector( VECTOR_C, VECTOR_A);
      av += UnitVector( VECTOR_D, VECTOR_A);
      av.Normalize();

      // extend out to the desired distance from A
      x *= DISTANCE_XA;

      // make position relative to A
      x += VECTOR_A;

      return x;
    }

    //! @brief calculates footpoint of POINT onto the line described by point ORIGIN and vector LINE
    //! @param POINT point of interest
    //! @param ORIGIN origin of line
    //! @param LINE for footpoint calculation
    //! @return footpoint of POINT onto the line described by point ORIGIN and vector LINE
    inline
    Vector3D CalculateFootpoint
    (
      const Vector3D &POINT,
      const Vector3D &ORIGIN,
      const Vector3D &LINE
    )
    {
      // calculate footpoint position in terms of fractional length of LINE vector
      const double footpoint_fraction( ScalarProduct( ( POINT - ORIGIN), LINE) / LINE.SquareNorm());

      // actual footpoint now is ORIGIN + footpoint_fraction * LINE
      const Vector3D footpoint( ORIGIN + footpoint_fraction * LINE);

      return footpoint;
    }

    //! @brief calculates distance of a point POINT from a line described by point ORIGIN and vector LINE
    //! @param POINT point of interest
    //! @param ORIGIN origin of line
    //! @param LINE for footpoint calculation
    //! @return distance of a point POINT from a line described by point ORIGIN and vector LINE
    inline
    double CalculateDistancePointFromLine( const Vector3D &POINT, const Vector3D &ORIGIN, const Vector3D &LINE)
    {
      // calculate footpoint of POINT onto line LINE
      const Vector3D footpoint( CalculateFootpoint( POINT, ORIGIN, LINE));

      // return distance between POINT and footpoint
      return Distance( POINT, footpoint);
    }

  } // namespace linal
} // namespace bcl

#endif //BCL_LINAL_VECTOR_3D_OPERATIONS_H_
