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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "math/bcl_math_quaternion.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_angle.h"

// external includes - sorted alphabetically

namespace bcl
{
  static double DistanceFromLine
  (
    const linal::Vector3D &ORIGIN, const linal::Vector3D &AXIS, const linal::Vector3D &POINT
  )
  {
    //! A line is defined by y = ORIGIN + x * AXIS
    const linal::Vector3D distance
    (
      ORIGIN - ( ( ORIGIN - POINT) * AXIS) / ( AXIS * AXIS) * AXIS - POINT
    );

    return distance.Norm();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_quaternion.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathQuaternion :
    public ExampleInterface
  {
  public:

    ExampleMathQuaternion *Clone() const
    {
      return new ExampleMathQuaternion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      math::Quaternion qa, qb, qc;
      const linal::Vector3D v( 1.7, 1.2, 1.8), origin( 0.0, 0.0, 0.0);
      linal::Vector3D axis( -1.0, -1.0, 0.0), vx( 1.0, 0.0, 0.0);
      const double alpha( math::Angle::Radian( 32.0));

      BCL_MessageStd( "Rotation of point x about a UNIT axis.");
      BCL_MessageStd( "The axis: " + util::Format()( axis));
      BCL_MessageStd( "The normalized axis: " + util::Format()( axis.Normalize()));
      BCL_MessageStd( "The angle alpha in degree: " + util::Format()( math::Angle::Degree( alpha)));

      BCL_MessageStd( "===== The point x:");
      BCL_MessageStd( util::Format()( v));
      BCL_MessageStd
      (
        "Test: norm of the vector x:  " + util::Format()( v.Norm())
      );
      BCL_MessageStd( "Test: distance of x from axis");
      BCL_MessageStd( util::Format()( DistanceFromLine( origin, axis, v)));

      qa.Equal( cos( alpha / 2.), sin( alpha / 2.) * axis);
      BCL_MessageStd
      (
        "The rotation is described by the quaternion q = ( cos(alpha/2) , sin(alpha/2) * a ):"
      );
      BCL_MessageStd( util::Format()( qa));
      BCL_MessageStd( "Test: norm of Quaternion q: " + util::Format()( qa.Norm()));

      qb.VectorToQuaternion( v);
      BCL_MessageStd( "Write vector x to Quaternion p: ");
      BCL_MessageStd( util::Format()( qb));

      BCL_MessageStd( "Test: norm of Quaternion p:   " + util::Format()( qb.Norm()));
      BCL_MessageStd( "Just another test: the normalized Quaternion p:  ");
      BCL_MessageStd( util::Format()( qb.Normalized()));

      //  BCL_Assert(  qa.Normalized().Norm() - 1. < LIMIT );

      BCL_MessageStd( "To rotate x we need also the conjugate of q, here written as q*:");
      BCL_MessageStd( util::Format()( math::Conjugate( qa)));

      BCL_MessageStd( "The rotated point p' is given by: p' = q p q*:");
      qc = qa * qb * math::Conjugate( qa);
      BCL_MessageStd( util::Format()( qc));
      BCL_MessageStd( "Test: norm of rotated point p':   " + util::Format()( qc.Norm()));

      const linal::Vector3D tmp( qc[ 1], qc[ 2], qc[ 3]);

      BCL_MessageStd
      (
        "Test: distance of p' to axis a:  " + util::Format()( DistanceFromLine( origin, axis, tmp))
      );

      BCL_MessageStd( "Comparison to Rotation using Matrix");

      BCL_MessageStd( util::Format()( v));

      BCL_MessageStd( "Test: distance from axis:");
      BCL_MessageStd( util::Format()( DistanceFromLine( origin, axis, v)));

      BCL_MessageStd( util::Format()( v * v));

      BCL_MessageStd
      (
        "A hidden usage of Quaternions in Constructor RotationMatrix3D( Vector3D , double )"
        " for building a rotation matrix for Rotation around arbitrary axes"
      );

      const linal::Vector3D uni( 1, 1, 1);
      const math::RotationMatrix3D rotate( uni, math::g_Pi / 2.0);

      BCL_MessageStd( "Before Rotation:\nPoint to rotate: " + util::Format()( vx));
      BCL_MessageStd
      (
        "Test1: distance from axis: " + util::Format()( DistanceFromLine( origin, uni, vx))
      );
      BCL_MessageStd
      (
        "Test2: angle: " + util::Format()( math::Angle::Degree( linal::ProjAngle( uni, vx)))
      );

      // the actual rotation
      vx.Rotate( rotate);

      BCL_MessageStd( "After Rotation:\n Rotated point:" + util::Format()( vx));
      BCL_MessageStd
      (
        "Test1: distance from axis: " + util::Format()( DistanceFromLine( origin, uni, vx))
      );
      BCL_MessageStd
      (
        "Test2: angle:  " + util::Format()( math::Angle::Degree( linal::ProjAngle( uni, vx)))
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathQuaternion

  const ExampleClass::EnumType ExampleMathQuaternion::s_Instance
  (
    GetExamples().AddEnum( ExampleMathQuaternion())
  );

} // namespace bcl
