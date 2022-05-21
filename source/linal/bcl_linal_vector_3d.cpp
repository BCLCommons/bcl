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
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "linal/bcl_linal_vector_3d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Vector3D::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Vector3D::Vector3D() :
      m_X( 0.0),
      m_Y( 0.0),
      m_Z( 0.0)
    {
    }

    //! @brief constructor from a single value
    //! @param VALUE common value for all the elements
    Vector3D::Vector3D( const double &VALUE) :
      m_X( VALUE),
      m_Y( VALUE),
      m_Z( VALUE)
    {
    }

    //! @brief constructor from three values
    //! @param X, Y, Z three elements
    Vector3D::Vector3D( const double &X, const double &Y, const double &Z) :
      m_X( X),
      m_Y( Y),
      m_Z( Z)
    {
    }

    //! @brief constructor from a pointer to three elements
    //! @param PTR_DATA a pointer to three elements
    Vector3D::Vector3D( const double *PTR_DATA) :
      m_X( PTR_DATA[ 0]),
      m_Y( PTR_DATA[ 1]),
      m_Z( PTR_DATA[ 2])
    {
    }

    //! @brief constructor from VectorInterface
    //! @param VECTOR vector
    Vector3D::Vector3D( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "given vector has wrong length");
      std::copy( VECTOR.Begin(), VECTOR.End(), &m_X);
    }

    //! copy constructor
    Vector3D *Vector3D::Clone() const
    {
      return new Vector3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Vector3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns true if all the coordinates are defined
    //! @return boolean true if all coordinates are defined - false otherwise
    bool Vector3D::IsDefined() const
    {
      return util::IsDefined( m_X) && util::IsDefined( m_Y) && util::IsDefined( m_Z);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief set all elements to elements from another vector
    //! @param VECTOR vector, which must contain 3 elements
    //! @return reference to this, after copying vector
    Vector3D &Vector3D::operator =( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "Tried to copy non-3d vector onto vector3d");
      m_X = VECTOR( 0);
      m_Y = VECTOR( 1);
      m_Z = VECTOR( 2);
      return *this;
    }

    //! @brief set all elements to a vector 3d
    //! @param VECTOR the vector to copy into this vector
    //! @return reference to this, after copying vector
    Vector3D &Vector3D::operator =( const Vector3D &VECTOR)
    {
      m_X = VECTOR.m_X;
      m_Y = VECTOR.m_Y;
      m_Z = VECTOR.m_Z;
      return *this;
    }

    //! @brief add a vector to this one
    //! @param VECTOR the vector to add
    //! @return reference to this
    Vector3D &Vector3D::operator +=( const Vector3D &VECTOR)
    {
      m_X += VECTOR.m_X;
      m_Y += VECTOR.m_Y;
      m_Z += VECTOR.m_Z;
      return *this;
    }

    //! @brief add a scalar to this vector
    //! @param SCALAR the scalar to add
    //! @return reference to this
    Vector3D &Vector3D::operator +=( const double &SCALAR)
    {
      m_X += SCALAR;
      m_Y += SCALAR;
      m_Z += SCALAR;
      return *this;
    }

    //! @brief subtract a vector from this one
    //! @param VECTOR the vector to subtract
    //! @return reference to this
    Vector3D &Vector3D::operator -=( const Vector3D &VECTOR)
    {
      m_X -= VECTOR.m_X;
      m_Y -= VECTOR.m_Y;
      m_Z -= VECTOR.m_Z;
      return *this;
    }

    //! @brief subtract a scalar from this vector
    //! @param SCALAR the value to subtract
    //! @return reference to this
    Vector3D &Vector3D::operator -=( const double &SCALAR)
    {
      m_X -= SCALAR;
      m_Y -= SCALAR;
      m_Z -= SCALAR;
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! change all values
    Vector3D &Vector3D::Set( const double &X, const double &Y, const double &Z)
    {
      m_X = X;
      m_Y = Y;
      m_Z = Z;
      return *this;
    }

    //! translate vector with translation
    Vector3D &Vector3D::Translate( const Vector3D &TRANSLATE)
    {
      m_X += TRANSLATE.m_X;
      m_Y += TRANSLATE.m_Y;
      m_Z += TRANSLATE.m_Z;
      return *this;
    }

    //! translate vector with translation
    Vector3D &Vector3D::Translate( const double X, const double Y, const double Z)
    {
      m_X += X;
      m_Y += Y;
      m_Z += Z;
      return *this;
    }

    //! rotate vector with rotation matrix
    Vector3D &Vector3D::Rotate( const math::RotationMatrix3D &ROTATE)
    {
      const double *row_x( ROTATE.m_RotationMatrix3D[ 0]);
      const double *row_y( ROTATE.m_RotationMatrix3D[ 1]);
      const double *row_z( ROTATE.m_RotationMatrix3D[ 2]);

      // this is a simple, unrolled version of operator *= ROTATE.m_RotationMatrix3D
      Set
      (
        m_X * row_x[ 0] + m_Y * row_y[ 0] + m_Z * row_z[ 0],
        m_X * row_x[ 1] + m_Y * row_y[ 1] + m_Z * row_z[ 1],
        m_X * row_x[ 2] + m_Y * row_y[ 2] + m_Z * row_z[ 2]
      );
      return *this;
    }

    //! transform vector with transformation matrix
    Vector3D &Vector3D::Transform( const math::TransformationMatrix3D &TRANSFORM)
    {
      const double *row_x( TRANSFORM.m_TransformationMatrix3D[ 0]);
      const double *row_y( TRANSFORM.m_TransformationMatrix3D[ 1]);
      const double *row_z( TRANSFORM.m_TransformationMatrix3D[ 2]);
      const double *row_t( TRANSFORM.m_TransformationMatrix3D[ 3]); // translation row

      Set
      (
        m_X * row_x[ 0] + m_Y * row_y[ 0] + m_Z * row_z[ 0] + row_t[ 0],
        m_X * row_x[ 1] + m_Y * row_y[ 1] + m_Z * row_z[ 1] + row_t[ 1],
        m_X * row_x[ 2] + m_Y * row_y[ 2] + m_Z * row_z[ 2] + row_t[ 2]
      );
      return *this;
    }

    //! return random translation vector equally distributed in a sphere of RADIUS
    Vector3D &Vector3D::SetRandomTranslation( const double &RADIUS)
    {
      //random angle between 0 and 2pi
      const double phi( 2 * math::g_Pi * random::GetGlobalRandom().Double());

      //random cos angle between -1 and 1
      const double theta( std::acos( random::GetGlobalRandom().Random( double( -1.0), double( 1.0))));

      //random distance
      const double distance( math::Pow( random::GetGlobalRandom().Random( RADIUS * RADIUS * RADIUS), 1.0 / 3.0));

      //product of sin( theta) * distance
      const double sin_theta_times_distance( std::sin( theta) * distance);

      //set the three values
      return Set(
                  std::cos( phi) * sin_theta_times_distance,
                  std::sin( phi) * sin_theta_times_distance,
                  std::cos( theta) * distance
                );
    }

    //! return random translation vector equally distributed in a ellipse of RADII
    Vector3D &Vector3D::SetRandomTranslation( const Vector3D &RADII)
    {
      // rejection sampling method. More efficient than parametric method due to absence of trig function calls
      // On average the loop terminates after 2 pi / 3 = 2.28 iterations (volume of cube / volume of ellipsoid)
      Vector3D radii_sq( math::Sqr( RADII.X()), math::Sqr( RADII.Y()), math::Sqr( RADII.Z()));
      if( RADII.X() && RADII.Y() && RADII.Z())
      {
        do
        {
          m_X = random::GetGlobalRandom().Double( math::Range< double>( -RADII.X(), RADII.X()));
          m_Y = random::GetGlobalRandom().Double( math::Range< double>( -RADII.Y(), RADII.Y()));
          m_Z = random::GetGlobalRandom().Double( math::Range< double>( -RADII.Z(), RADII.Z()));
        } while
          (
            math::Sqr( m_X) / radii_sq.X()
            + math::Sqr( m_Y) / radii_sq.Y()
            + math::Sqr( m_Z) / radii_sq.Z()
            > 1.0
          );
      }
      else
      {
        radii_sq.X() = std::max( radii_sq.X(), 1e-38);
        radii_sq.Y() = std::max( radii_sq.Y(), 1e-38);
        radii_sq.Z() = std::max( radii_sq.Z(), 1e-38);
        do
        {
          m_X = RADII.X() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.X(), RADII.X())) : 0.0;
          m_Y = RADII.Y() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.Y(), RADII.Y())) : 0.0;
          m_Z = RADII.Z() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.Z(), RADII.Z())) : 0.0;
        } while
          (
            math::Sqr( m_X) / radii_sq.X()
            + math::Sqr( m_Y) / radii_sq.Y()
            + math::Sqr( m_Z) / radii_sq.Z()
            > 1.0
          );
      }

      //set the three values
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! print in a user readable format
    //! @return the string containing the information
    std::string Vector3D::ToString() const
    {
      return std::string( util::Format()( m_X) + " " + util::Format()( m_Y) + " " + util::Format()( m_Z));
    }

    //! @brief normalize the vector (such that inner product == 1)
    Vector3D &Vector3D::Normalize()
    {
      const double sqnorm( SquareNorm());
      if( sqnorm != 1.0)
      {
        const double norm( math::Sqrt( sqnorm));
        m_X /= norm;
        m_Y /= norm;
        m_Z /= norm;
      }
      return *this;
    }

    //! @brief normalize the vector (such that inner product == 1)
    //! Overrides function of the same name in VectorInterface for performance reasons
    Vector3D Vector3D::Normalized() const
    {
      Vector3D c( *this);
      return c.Normalize();
    }

    //! @brief norm = length of vector
    //! @return length of vector
    double Vector3D::Norm() const
    {
      return math::Sqrt( SquareNorm());
    }

    //! @brief square norm = square length of vector
    //! @return square length of vector
    double Vector3D::SquareNorm() const
    {
      return m_X * m_X + m_Y * m_Y + m_Z * m_Z;
    }

    //! @brief sum up all elements
    //! @return sum of all elements
    double Vector3D::Sum() const
    {
      return m_X + m_Y + m_Z;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Vector3D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X, ISTREAM);
      io::Serialize::Read( m_Y, ISTREAM);
      io::Serialize::Read( m_Z, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Vector3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Y, OSTREAM, 0) << '\t';
      io::Serialize::Write( m_Z, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
