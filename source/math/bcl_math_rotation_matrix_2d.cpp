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
#include "math/bcl_math_rotation_matrix_2d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RotationMatrix2D::s_Instance
    (
      GetObjectInstances().AddInstance( new RotationMatrix2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct empty RotationMatrix2D
    RotationMatrix2D::RotationMatrix2D()
    {
      // set the matrix up as the unitary matrix (no rotation)
      m_RotationMatrix2D( 0, 0) = m_RotationMatrix2D( 1, 1) = 1.0;
    }

    //! construct RotationMatrix2D from RotationAngle
    RotationMatrix2D::RotationMatrix2D( const double &ANGLE)
    {
      SetAngle( ANGLE);
    }

    //! construct RotationMatrix2D from a 2x2 matrix
    RotationMatrix2D::RotationMatrix2D( const linal::Matrix2x2< double> &ROTATION) :
      m_RotationMatrix2D( ROTATION)
    {
    }

    //! virtual copy constructor
    RotationMatrix2D *RotationMatrix2D::Clone() const
    {
      return new RotationMatrix2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RotationMatrix2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return the matrix
    const linal::Matrix2x2< double> &RotationMatrix2D::GetMatrix() const
    {
      return m_RotationMatrix2D;
    }

    //! @brief return the angle that this matrix rotates by
    //! @return rotation angle that this matrix represents (in radians)
    double RotationMatrix2D::GetAngle() const
    {
      return acos( *m_RotationMatrix2D.Begin());
    }

    //! @brief set the angle that this matrix rotates by
    //! @param THETA rotation angle that this matrix represents (in radians)
    void RotationMatrix2D::SetAngle( const double &THETA)
    {
      SetFromCosSin( cos( THETA), sin( THETA));
    }

    //! @brief set given cos(theta) and sin(theta)
    //! @param COS cosine of theta
    //! @param SIN sine theta
    void RotationMatrix2D::SetFromCosSin( const double &COS, const double &SIN)
    {
      m_RotationMatrix2D( 0, 0) = m_RotationMatrix2D( 1, 1) = COS;
      m_RotationMatrix2D( 1, 0) = -( m_RotationMatrix2D( 0, 1) = -SIN);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set random rotation from 0 - MAX_ROTATIONANGLE
    //! @param MAX_ROTATIONANGLE the maximum rotation angle
    RotationMatrix2D &RotationMatrix2D::SetRand( const double MAX_ROTATIONANGLE)
    {
      //random angle between 0 and 2pi
      SetAngle( MAX_ROTATIONANGLE * random::GetGlobalRandom().Double());
      return *this;
    }

    //! @brief set the matrix up to be a givens rotation matrix; which is a numerically stable solution used in tridiagonalization
    //! @param A, B the vector to rotate in the givens rotation matrix problem:
    //! | c -s || A | = | r |
    //! | s  c || B |   | 0 |
    //! @return R (equal to pythag(A,B))
    double RotationMatrix2D::MakeGivens( const double &A, const double &B)
    {
      const double abs_a( Absolute( A)), abs_b( Absolute( B));
      if( abs_a > abs_b)
      {
        if( B == 0.0)
        {
          SetFromCosSin( A > 0.0 ? 1.0 : -1.0, 0.0);
          return abs_a;
        }
        // |A| > |B|, neither B nor A are 0
        const double t( B / A);
        const double u( A > 0.0 ? std::sqrt( 1.0 + t * t) : -std::sqrt( 1.0 + t * t));
        const double c( 1.0 / u);
        SetFromCosSin( c, -c * t);
        return A * u;
      }
      else if( A == 0.0)
      {
        BCL_Assert( A != B, "Givens rotation undefined if A and B are both 0!");
        SetFromCosSin( 0.0, B > 0.0 ? -1.0 : 1.0);
        return abs_b;
      }

      // |A| <= |B|, neither B nor A are 0
      const double t( A / B);
      const double u( B > 0.0 ? std::sqrt( 1.0 + t * t) : -std::sqrt( 1.0 + t * t));
      const double s( -1.0 / u);
      SetFromCosSin( -s * t, s);
      return B * u;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator *= RotationMatrix2D
    RotationMatrix2D &RotationMatrix2D::operator *=( const RotationMatrix2D &MATRIX)
    {
      m_RotationMatrix2D *= ( MATRIX.m_RotationMatrix2D);
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write RotationMatrix2D to std::ostream using the given util::Format
    std::ostream &RotationMatrix2D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_RotationMatrix2D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! write RotationMatrix2D to std::ostream using the given util::Format
    std::ostream &RotationMatrix2D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_RotationMatrix2D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! read RotationMatrix2D from io::IFStream
    std::istream &RotationMatrix2D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_RotationMatrix2D, ISTREAM);
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl

