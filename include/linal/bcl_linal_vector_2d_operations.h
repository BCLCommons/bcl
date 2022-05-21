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

#ifndef BCL_LINAL_VECTOR_2D_OPERATIONS_H_
#define BCL_LINAL_VECTOR_2D_OPERATIONS_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_2d.h"
#include "bcl_linal_vector_operations.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_linal_vector_2d_operations.h
  //! @brief mathematical operators and angle operations for 2d vectors
  //!
  //! @see @link example_linal_vector_2d_operations.cpp @endlink
  //! @author mendenjl
  //! @date Dec 05, 2013
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace linal
  {

  ///////////////
  // operators //
  ///////////////

    //! @brief operator +Vector2D
    //! @param VECTOR original vector
    //! @return original vector
    inline
    Vector2D operator +( const Vector2D &VECTOR)
    {
      return VECTOR;
    }

    //! @brief operator Vector2D + Vector2D
    //! @param VECTOR_A first vector to be added
    //! @param VECTOR_B second vector to be added
    //! @return sum of VECTOR_A and VECTOR_B
    inline
    Vector2D
    operator +( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      Vector2D sum( VECTOR_A);
      sum += VECTOR_B;
      return sum;
    }

    //! @brief operator Vector2D - Vector2D
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector to be subtracted
    //! @return difference of VECTOR_A and VECTOR_B
    inline
    Vector2D operator -( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      Vector2D diff( VECTOR_A);
      diff -= VECTOR_B;
      return diff;
    }

    //! @brief operator Vector2D + scalar
    //! @param VECTOR vector to be added
    //! @param X scalar to be added
    //! @return sum of VECTOR and X
    inline
    Vector2D
    operator +( const Vector2D &VECTOR, const double &X)
    {
      Vector2D sum( VECTOR);
      sum += X;
      return sum;
    }

    //! @brief operator Vector2D - scalar
    //! @param VECTOR vector
    //! @param X scalar to be subtracted
    //! @return difference of VECTOR and X
    inline
    Vector2D
    operator -( const Vector2D &VECTOR, const double &X)
    {
      Vector2D sum( VECTOR);
      sum -= X;
      return sum;
    }

    //! @brief operator scalar + Vector2D
    //! @param X scalar to be added
    //! @param VECTOR vector to be added
    //! @return sum of X and VECTOR
    inline
    Vector2D operator +( const double &X, const Vector2D &VECTOR)
    {
      Vector2D sum( VECTOR);
      sum += X;
      return sum;
    }

    //! @brief operator * (dot product)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A and VECTOR_B
    inline
    double operator *( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return ( VECTOR_A.X() * VECTOR_B.X()) + ( VECTOR_A.Y() * VECTOR_B.Y());
    }

    //! @brief operator scalar * Vector2D
    //! @param X scalar to multiply
    //! @param VECTOR to multiply
    //! @return X times VECTOR
    inline
    Vector2D operator *( const double &X, const Vector2D &VECTOR)
    {
      Vector2D scaled( VECTOR);
      scaled *= X;
      return scaled;
    }

    //! @brief operator Vector2D * scalar
    //! @param VECTOR to multiply
    //! @param X scalar to multiply
    //! @return VECTOR times X
    inline
    Vector2D operator *( const Vector2D &VECTOR, const double &X)
    {
      Vector2D scaled( VECTOR);
      scaled *= X;
      return scaled;
    }

    //! @brief operator Vector2D / scalar
    //! @param VECTOR to divide
    //! @param X scalar to divide by
    //! @return X times VECTOR
    inline
    Vector2D operator /( const Vector2D &VECTOR, const double &X)
    {
      Vector2D scaled( VECTOR);
      scaled /= X;
      return scaled;
    }

    //! @brief operator -Vector2D
    //! @param VECTOR original vector
    //! @return additive inverse of VECTOR
    inline
    Vector2D operator -( const Vector2D &VECTOR)
    {
      return Vector2D( -VECTOR.X(), -VECTOR.Y());
    }

    //! @brief operator scalar - Vector2D
    //! @param X to subtract from
    //! @param VECTOR to subtract by
    //! @return X minus VECTOR
    inline
    Vector2D operator -( const double &X, const Vector2D &VECTOR)
    {
      Vector2D scalar_minus_vector( -VECTOR);
      scalar_minus_vector += X;
      return scalar_minus_vector;
    }

    //! @brief operator ==
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @return true if elements from VECTOR_A == respective elements from VECTOR_B
    inline
    bool operator ==( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return VECTOR_A.X() == VECTOR_B.X() && VECTOR_A.Y() == VECTOR_B.Y();
    }

    //! @brief operator !=
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @return true if any element from VECTOR_A != respective element from VECTOR_B
    inline
    bool operator !=( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return !( VECTOR_A == VECTOR_B);
    }

  //////////////////////
  // Vector functions //
  //////////////////////

    //! @brief calculates the absolute difference between two linal::Vector2Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector2D (points)
    Vector2D AbsoluteDifference( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B);

    //! @brief calculates the unit vector directed from one linal::Vector2D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector2D UnitVector( const Vector2D &ORIGIN, const Vector2D &TARGET);

    //! @brief compares vectors like operator ==, but returns true if difference between both is smaller than EPSILON
    //! @param VECTOR_A first vector to be compared
    //! @param VECTOR_B second vector to be compared
    //! @param RELATIVE_TOLERANCE relative tolerance for comparison
    //! @param ABSOLUTE_TOLERANCE absolute tolerance for comparison
    //! @return true if difference between both is smaller than EPSILON
    inline
    bool EqualWithinTolerance
    (
      const Vector2D &VECTOR_A,
      const Vector2D &VECTOR_B,
      const double &RELATIVE_TOLERANCE = 0.001,
      const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon()
    )
    {
      return
        math::EqualWithinTolerance( VECTOR_A.X(), VECTOR_B.X(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
        && math::EqualWithinTolerance( VECTOR_A.Y(), VECTOR_B.Y(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);
    }

    //! @brief calculates square distance between two linal::Vector2D (points)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return square distance between two linal::Vector2D (points)
    inline
    double SquareDistance( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return math::Sqr( VECTOR_A.X() - VECTOR_B.X())
             + math::Sqr( VECTOR_A.Y() - VECTOR_B.Y());
    }

    //! @brief calculates distance between two linal::Vector2D (points)
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return distance between two linal::Vector2D (points)
    inline
    double Distance( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return math::Sqrt( SquareDistance( VECTOR_A, VECTOR_B));
    }

    //! @brief calculate cosine of projection angle
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @return cosine of projection angle between two linal::Vector2D
    inline
    double ProjAngleCosinus( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
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
    double ProjAngleCosinus( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B, const Vector2D &VECTOR_C)
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
      const Vector2D &VECTOR_A,
      const Vector2D &VECTOR_B,
      const Vector2D &VECTOR_C,
      const Vector2D &VECTOR_D
    )
    {
      return ProjAngleCosinus( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! @brief calculates projection angle between two linal::Vector2D
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @return projection angle between two linal::Vector2D
    inline
    double ProjAngle( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      const double projection_angle_cosinus( ProjAngleCosinus( VECTOR_A, VECTOR_B));

      // check that the result is defined; this is more robust than checking if the vectors are defined,
      // since even defined vectors may have 0 norms, leading to undefined cosine
      if( !util::IsDefined( projection_angle_cosinus))
      {
        return util::GetUndefined< double>();
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
    double ProjAngle( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B, const Vector2D &VECTOR_C)
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
      const Vector2D &VECTOR_A,
      const Vector2D &VECTOR_B,
      const Vector2D &VECTOR_C,
      const Vector2D &VECTOR_D
    )
    {
      return ProjAngle( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! @brief scalar (dot) product of two linal::Vector2Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return scalar (dot) product of VECTOR_A and VECTOR_B
    inline
    double ScalarProduct( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      return ( VECTOR_A * VECTOR_B);
    }

    //! @brief cross product of two linal::Vector2Ds
    //! @param A first vector
    //! @param B second vector
    //! @return cross product of VECTOR_A and VECTOR_B (== the determinant of the matrix given by [[VECTOR_A][VECTOR_B]]
    inline
    double CrossProduct( const Vector2D &A, const Vector2D &B)
    {
      return A.X() * B.Y() - A.Y() * B.X();
    }

    //! @brief cross product of three points (A->B and A->C)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @return cross product of three points
    inline
    double CrossProduct
    (
      const Vector2D &VECTOR_A,
      const Vector2D &VECTOR_B,
      const Vector2D &VECTOR_C
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
    double CrossProduct
    (
      const Vector2D &VECTOR_A,
      const Vector2D &VECTOR_B,
      const Vector2D &VECTOR_C,
      const Vector2D &VECTOR_D
    )
    {
      return CrossProduct( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! @brief calculates footpoint of POINT onto the line described by point ORIGIN and vector LINE
    //! @param POINT point of interest
    //! @param ORIGIN origin of line
    //! @param LINE for footpoint calculation
    //! @return footpoint of POINT onto the line described by point ORIGIN and vector LINE
    inline
    Vector2D CalculateFootpoint
    (
      const Vector2D &POINT,
      const Vector2D &ORIGIN,
      const Vector2D &LINE
    )
    {
      // calculate footpoint position in terms of fractional length of LINE vector
      const double footpoint_fraction( ScalarProduct( ( POINT - ORIGIN), LINE) / LINE.SquareNorm());

      // actual footpoint now is ORIGIN + footpoint_fraction * LINE
      const Vector2D footpoint( ORIGIN + footpoint_fraction * LINE);

      return footpoint;
    }

    //! @brief calculates distance of a point POINT from a line described by point ORIGIN and vector LINE
    //! @param POINT point of interest
    //! @param ORIGIN origin of line
    //! @param LINE for footpoint calculation
    //! @return distance of a point POINT from a line described by point ORIGIN and vector LINE
    inline
    double CalculateDistancePointFromLine( const Vector2D &POINT, const Vector2D &ORIGIN, const Vector2D &LINE)
    {
      // calculate footpoint of POINT onto line LINE
      const Vector2D footpoint( CalculateFootpoint( POINT, ORIGIN, LINE));

      // return distance between POINT and footpoint
      return Distance( POINT, footpoint);
    }

  } // namespace linal
} // namespace bcl

#endif //BCL_LINAL_VECTOR_2D_OPERATIONS_H_
