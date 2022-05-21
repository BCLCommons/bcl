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
#include "math/bcl_math_rotation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_logger_interface.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RotationMatrix3D::s_Instance
    (
      GetObjectInstances().AddInstance( new RotationMatrix3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct empty RotationMatrix3D
    RotationMatrix3D::RotationMatrix3D()
    {
      // set the matrix up as the unitary matrix (no rotation)
      m_RotationMatrix3D( 0, 0) = m_RotationMatrix3D( 1, 1) = m_RotationMatrix3D( 2, 2) = 1.0;
    }

    //! construct RotationMatrix3D from RotationAxis and RotationAngle
    RotationMatrix3D::RotationMatrix3D( const coord::Axis &AXIS, const double &ANGLE)
    {
      Init( AXIS, ANGLE);
    }

    //! construct RotationMatrix3D from any linal::Vector3D around which shall be rotated and the Angle.
    RotationMatrix3D::RotationMatrix3D( const linal::Vector3D &AXIS, const double ANGLE)
    {
      if( !ANGLE)
      {
        // angle is 0, ignore the axis (which does not matter in this case)
        // set the matrix up as the unitary matrix (no rotation)
        m_RotationMatrix3D( 0, 0) = m_RotationMatrix3D( 1, 1) = m_RotationMatrix3D( 2, 2) = 1.0;
      }
      else
      {
        // normalize the provided axis
        linal::Vector3D norm_axis( AXIS);
        norm_axis.Normalize();

        const double rcos( cos( ANGLE));

        // minor improvements in performance and numeric stability can be had by precomputing
        // norm(AXIS) * (1-cos(ANGLE)) and
        // norm(AXIS) * sin(ANGLE)
        linal::Vector3D norm_axis_remain( norm_axis);
        norm_axis_remain *= double( 1.0) - rcos;

        linal::Vector3D norm_axis_rsin( norm_axis);
        norm_axis_rsin *= sin( ANGLE);

        //first col
        m_RotationMatrix3D( 0, 0) =                rcos + norm_axis.X() * norm_axis_remain.X();
        m_RotationMatrix3D( 1, 0) =  norm_axis_rsin.Z() + norm_axis.Y() * norm_axis_remain.X();
        m_RotationMatrix3D( 2, 0) = -norm_axis_rsin.Y() + norm_axis.Z() * norm_axis_remain.X();

        //second col
        m_RotationMatrix3D( 0, 1) = -norm_axis_rsin.Z() + norm_axis.X() * norm_axis_remain.Y();
        m_RotationMatrix3D( 1, 1) =                rcos + norm_axis.Y() * norm_axis_remain.Y();
        m_RotationMatrix3D( 2, 1) =  norm_axis_rsin.X() + norm_axis.Z() * norm_axis_remain.Y();

        //third col
        m_RotationMatrix3D( 0, 2) =  norm_axis_rsin.Y() + norm_axis.X() * norm_axis_remain.Z();
        m_RotationMatrix3D( 1, 2) = -norm_axis_rsin.X() + norm_axis.Y() * norm_axis_remain.Z();
        m_RotationMatrix3D( 2, 2) =                rcos + norm_axis.Z() * norm_axis_remain.Z();
      }
    }

    //! construct RotationMatrix3D from three RotationAxis and three RotationAngles
    RotationMatrix3D::RotationMatrix3D( const coord::Axis AXES[ 3], const double ANGLE[ 3])
    {
      Init( AXES[ 0], ANGLE[ 0]);
      RotationMatrix3D b( AXES[ 1], ANGLE[ 1]), c( AXES[ 2], ANGLE[ 2]);
      operator *= ( b);
      operator *= ( c);
    }

    // assume rotation of PHI around x-axis, THETA around y-axis and PSI around z-axis
    //! construct RotationMatrix3D from three Euler angles
    RotationMatrix3D::RotationMatrix3D( const double PHI, const double THETA, const double PSI)
    {
      // pre calculate some values for efficiency
      const double c_x( cos( PHI));         // A
      const double s_x( sin( PHI));         // B
      const double c_y( cos( THETA));       // C
      const double s_y( sin( THETA));       // D
      const double c_z( cos( PSI));         // E
      const double s_z( sin( PSI));         // F
      const double c_x_s_y( c_x * s_y);     // AD
      const double s_x_s_y( s_x * s_y);     // BD

      m_RotationMatrix3D( 0, 0) =  c_y * c_z ;
      m_RotationMatrix3D( 0, 1) = -c_y * s_z;
      m_RotationMatrix3D( 0, 2) =  s_y;

      m_RotationMatrix3D( 1, 0) =  s_x_s_y * c_z + c_x * s_z;
      m_RotationMatrix3D( 1, 1) = -s_x_s_y * s_z + c_x * c_z;
      m_RotationMatrix3D( 1, 2) = -s_x * c_y;

      m_RotationMatrix3D( 2, 0) = -c_x_s_y * c_z + s_x * s_z;
      m_RotationMatrix3D( 2, 1) =  c_x_s_y * s_z + s_x * c_z;
      m_RotationMatrix3D( 2, 2) =  c_x * c_y;
    }

    //! construct RotationMatrix3D from a 3x3 matrix
    RotationMatrix3D::RotationMatrix3D( const linal::Matrix3x3< double> &ROTATION) :
      m_RotationMatrix3D( ROTATION)
    {
    }

    //! virtual copy constructor
    RotationMatrix3D *RotationMatrix3D::Clone() const
    {
      return new RotationMatrix3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RotationMatrix3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return the matrix
    const linal::Matrix3x3< double> &RotationMatrix3D::GetMatrix() const
    {
      return m_RotationMatrix3D;
    }

    //! @brief return the X Y or Z axis in respect to their reference coordinate system x, y and z (euclidean unit vectors)
    //! @param AXIS axis that will be returned
    //! @return axis in coordinates of reference coordinate system
    linal::Vector3D RotationMatrix3D::GetAxis( const coord::Axis &AXIS) const
    {
      // axis is stored in the columns
      return linal::Vector3D( m_RotationMatrix3D( 0, AXIS), m_RotationMatrix3D( 1, AXIS), m_RotationMatrix3D( 2, AXIS));
    }

    //! @brief return the axis of rotation; result is a unit vector
    //! @return axis of rotation; if this rotation matrix is applied to it, the vector will be unchanged
    linal::Vector3D RotationMatrix3D::GetRotationAxis() const
    {
      linal::Vector3D ret
      (
        m_RotationMatrix3D( 2, 1) - m_RotationMatrix3D( 1, 2),
        m_RotationMatrix3D( 0, 2) - m_RotationMatrix3D( 2, 0),
        m_RotationMatrix3D( 1, 0) - m_RotationMatrix3D( 0, 1)
      );
      const double norm( ret.Norm());
      if( norm == 0.0)
      {
        ret += 1.0 / Sqrt( 3.0);
      }
      else
      {
        ret /= norm;
      }
      return ret;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return the euler angles for that rotation matrix
    //! http://www.gamedev.net/reference/articles/article1691.asp
    //! @return linal::Vector3D with 1st to 3rd element being alpha, beta, gamma - rotations around x y and z
    linal::Vector3D RotationMatrix3D::EulerAnglesXYZ() const
    {
      linal::Vector3D euler_angles_x_y_z( 0.0, 0.0, 0.0);

      euler_angles_x_y_z( 1) = asin( m_RotationMatrix3D( 0,2));
      const double cos_b     = cos( euler_angles_x_y_z( 1))   ;
      // check for no gimbal lock
      if( Absolute( cos_b) > 0.005)
      {
        // x axis angle
        euler_angles_x_y_z( 0) = atan2( -m_RotationMatrix3D( 1, 2) / cos_b, m_RotationMatrix3D( 2, 2) / cos_b);
        // z axis angle
        euler_angles_x_y_z( 2) = atan2( -m_RotationMatrix3D( 0, 1) / cos_b, m_RotationMatrix3D( 0, 0) / cos_b);
      }
      // there was a gimbal lock
      else
      {
        // x axis angle
        euler_angles_x_y_z( 0) = 0;
        // z axis angle
        euler_angles_x_y_z( 2) = atan2( m_RotationMatrix3D( 1, 0), m_RotationMatrix3D( 1, 1));
      }

      // add 2 pi to all negative angles
      for( double *ptr( euler_angles_x_y_z.Begin()), *ptr_end( euler_angles_x_y_z.End()); ptr != ptr_end; ++ptr)
      {
        if( *ptr < double( 0.0))
        {
          *ptr += 2 * g_Pi;
        }
      }

      // end
      return euler_angles_x_y_z;
    }

    //! @brief return the Euler angles according to z-x'-z'' convention
    //! @return Euler angles according to z-x'-z'' convention
    storage::VectorND< 3, double> RotationMatrix3D::EulerAnglesZXZ() const
    {
      const double beta( std::acos( m_RotationMatrix3D( 2, 2)));
      const double gamma
      (
        std::acos
        (
          m_RotationMatrix3D( 2, 0) / beta == 0.0
          ? std::numeric_limits< double>::epsilon()
          : std::sin( beta)
        )
      );
      const double alpha
      (
        std::asin
        (
          m_RotationMatrix3D( 0, 2) / beta == 0.0
          ? std::numeric_limits< double>::epsilon()
          : std::sin( beta)
        )
      );
      const storage::VectorND< 3, double> euler_angles( alpha, beta, gamma);

      return euler_angles;
    }

    //http://journals.iucr.org/j/issues/2002/05/00/vi0166/index.html
    //! calculate effective rotation angle
    double RotationMatrix3D::EffectiveRotationAngle() const
    {
      //make sure that initial argument for std::acos is between [-1,1]
      return std::acos( std::max( std::min( ( m_RotationMatrix3D.Trace() - 1.0) / 2.0, 1.0), -1.0));
    }

    // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
    //! set random rotation from 0 - Pi
    RotationMatrix3D &RotationMatrix3D::SetRand( const double MAX_ROTATIONANGLE)
    {
      do
      {
        //random angle between -pi and pi
        const double a1( MAX_ROTATIONANGLE * ( 2.0 * random::GetGlobalRandom().Double() - 1));
        const double a3( MAX_ROTATIONANGLE * ( 2.0 * random::GetGlobalRandom().Double() - 1));
        //random cos angle between -pi/2 and 3pi/2
        double a2
        (
          std::acos
          (
            random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0))
          ) + g_Pi * 0.5
        );
        if( random::GetGlobalRandom().Boolean())
        {
          a2 += ( a2 < g_Pi ? g_Pi : -g_Pi);
        }
        a2 *= MAX_ROTATIONANGLE / g_Pi;

        ( *this = RotationMatrix3D( a1, a2, a3));
      } while( this->EffectiveRotationAngle() > MAX_ROTATIONANGLE);
      return *this;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct from alpha, beta and gamma angles according to z-x-z convention
    //! @param ALPHA rotation angle around LCS z-axis
    //! @param BETA rotation angle around LCS x'-axis
    //! @param GAMMA rotation angle around LCS z''-axis
    RotationMatrix3D RotationMatrix3D::CreateZXZ( const double ALPHA, const double BETA, const double GAMMA)
    {
      // 3x3 matrix with 0.0 as elements
      linal::Matrix3x3< double> zxz_rotation;

      zxz_rotation( 0, 0) =  cos( ALPHA) * cos( GAMMA) - sin( ALPHA) * cos( BETA) * sin( GAMMA);
      zxz_rotation( 0, 1) = -cos( ALPHA) * sin( GAMMA) - sin( ALPHA) * cos( BETA) * cos( GAMMA);
      zxz_rotation( 0, 2) =  sin(  BETA) * sin( ALPHA);

      zxz_rotation( 1, 0) =  sin( ALPHA) * cos( GAMMA) + cos( ALPHA) * cos( BETA) * sin( GAMMA);
      zxz_rotation( 1, 1) = -sin( ALPHA) * sin( GAMMA) + cos( ALPHA) * cos( BETA) * cos( GAMMA);
      zxz_rotation( 1, 2) = -sin(  BETA) * cos( ALPHA);

      zxz_rotation( 2, 0) =  sin( BETA) * sin( GAMMA);
      zxz_rotation( 2, 1) =  sin( BETA) * cos( GAMMA);
      zxz_rotation( 2, 2) =  cos( BETA);

      return RotationMatrix3D( zxz_rotation);
    }

    //! private helper function to set rotation matrix data from RotationAxis and RotationAngles
    void RotationMatrix3D::Init( const coord::Axis &AXIS, const double &ANGLE)
    {
      m_RotationMatrix3D( AXIS, AXIS) = 1;
      m_RotationMatrix3D( ( AXIS + 1) % 3, ( AXIS + 1) % 3) =  cos( ANGLE);
      m_RotationMatrix3D( ( AXIS + 2) % 3, ( AXIS + 2) % 3) =  cos( ANGLE);
      m_RotationMatrix3D( ( AXIS + 1) % 3, ( AXIS + 2) % 3) = -sin( ANGLE);
      m_RotationMatrix3D( ( AXIS + 2) % 3, ( AXIS + 1) % 3) =  sin( ANGLE);
    }

  ///////////////
  // operators //
  ///////////////

    //! operator *= RotationMatrix3D
    RotationMatrix3D &RotationMatrix3D::operator *=( const RotationMatrix3D &MATRIX)
    {
      m_RotationMatrix3D *= ( MATRIX.m_RotationMatrix3D);
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write RotationMatrix3D to std::ostream using the given util::Format
    std::ostream &RotationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_RotationMatrix3D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! write RotationMatrix3D to std::ostream using the given util::Format
    std::ostream &RotationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_RotationMatrix3D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! read RotationMatrix3D from io::IFStream
    std::istream &RotationMatrix3D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_RotationMatrix3D, ISTREAM);
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
