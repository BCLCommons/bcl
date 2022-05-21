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
#include "math/bcl_math_rotation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_rotation_matrix_3d.cpp
  //!
  //! @author butkiem1, fischea
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRotationMatrix3D :
    public ExampleInterface
  {
  public:

    ExampleMathRotationMatrix3D *Clone() const
    { return new ExampleMathRotationMatrix3D( *this);}

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
    /////////////////
    // constructor //
    /////////////////

      //default constructor
      const math::RotationMatrix3D rotationmatrix3D_default;

      //construct from axes and angle
      const math::RotationMatrix3D rotationmatrix3D_x( coord::GetAxes().e_X, math::g_Pi / 2.5);
      const math::RotationMatrix3D rotationmatrix3D_y( coord::GetAxes().e_Y, math::g_Pi / 2.5);
      const math::RotationMatrix3D rotationmatrix3D_z( coord::GetAxes().e_Z, math::g_Pi / 2.5);

      //construct from rotation axis (vector3D) and angle
      const math::RotationMatrix3D rotationmatrix3D_123( linal::Vector3D( 1.0, 2.0, 3.0), math::g_Pi / 2.0);

      //construct from three axes and three angles
      const coord::Axis axes[]   = { coord::GetAxes().e_Y, coord::GetAxes().e_Z, coord::GetAxes().e_X};
      const double             angles[] = { 1.0, 2.0, math::g_Pi};
      math::RotationMatrix3D rotationmatrix3D_yzx( axes, angles);

      //construct from phi, psi, theta (three euler angles)
      math::RotationMatrix3D rotationmatrix3D_euler( math::g_Pi / 1.0, math::g_Pi / 2.0, math::g_Pi / 3.0);

      BCL_MessageStd( util::Format()( rotationmatrix3D_euler.EulerAnglesXYZ()));

      // construct from alpha beta and gamma - according to zxz convention
      math::RotationMatrix3D rotationmatrix3D_zxz( math::RotationMatrix3D::CreateZXZ( 0.5, 1.2, 1.67));

      // construct from different euler angles to check if euler angle extraction and reusage works
      const double angle_step_size( math::g_Pi / 1.5);
      for( double alpha( -2 * math::g_Pi); alpha < 2 * math::g_Pi; alpha += angle_step_size)
      {
        for( double beta( -2 * math::g_Pi); beta < 2 * math::g_Pi; beta += angle_step_size)
        {
          for( double gamma( -2 * math::g_Pi); gamma < 2 * math::g_Pi; gamma += angle_step_size)
          {
            // construct from angles
            const math::RotationMatrix3D rotate_xyz_orig( math::RotationMatrix3D( alpha, beta, gamma));
            // extract angles
            const linal::Vector3D euler_angles( rotate_xyz_orig.EulerAnglesXYZ());
            // construct from extracted angles
            const math::RotationMatrix3D rotate_xyz_reconstructed( math::RotationMatrix3D( euler_angles( 0), euler_angles( 1), euler_angles( 2)));

            // write original and extracted angles
            BCL_MessageDbg
            (
              "a, b, y:\n" +
              util::Format()( alpha) + "\t" + util::Format()( euler_angles( 0)) + "\n" +
              util::Format()( beta) + "\t" + util::Format()( euler_angles( 1)) + "\n" +
              util::Format()( gamma) + "\t" + util::Format()( euler_angles( 2))
            );

            // identical test vectors are rotate with original and reconstructed rotation matrix
            linal::Vector3D test_a( 1.1, 1.5, 1.7);
            linal::Vector3D test_b( 1.1, 1.5, 1.7);
            test_a.Rotate( rotate_xyz_orig);
            test_b.Rotate( rotate_xyz_reconstructed);

            // check the resulting points are equal
            BCL_Example_Check
            (
              math::EqualWithinTolerance( test_a, test_b),
              "vectors should be the same:\n" + util::Format()( test_a) + "\n" + util::Format()( test_b) +
              "\n" + util::Format()( rotate_xyz_orig) + "\n" + util::Format()( rotate_xyz_reconstructed)
            );
          }
        }
      }

      //construct from 3x3 matrix
      math::RotationMatrix3D rotationmatrix3D_frommatrix( linal::Matrix3x3< double>( 1.0));

      //randomize rotation angle to a maximum of math::g_Pi / 2
      math::RotationMatrix3D rotationmatrix3D_rand;
      rotationmatrix3D_rand.SetRand( math::g_Pi / double( 10));

      //calculate effective rotation angles
      BCL_MessageStd
      (
        "effective rotation angles"
          "\nrotationmatrix3D_default:     " + util::Format()( rotationmatrix3D_default.EffectiveRotationAngle())
        + "\nrotationmatrix3D_x:           " + util::Format()( rotationmatrix3D_x.EffectiveRotationAngle())
        + "\nrotationmatrix3D_y:           " + util::Format()( rotationmatrix3D_y.EffectiveRotationAngle())
        + "\nrotationmatrix3D_z:           " + util::Format()( rotationmatrix3D_z.EffectiveRotationAngle())
        + "\nrotationmatrix3D_123:         " + util::Format()( rotationmatrix3D_123.EffectiveRotationAngle())
        + "\nrotationmatrix3D_yzx:         " + util::Format()( rotationmatrix3D_yzx.EffectiveRotationAngle())
        + "\nrotationmatrix3D_euler:       " + util::Format()( rotationmatrix3D_euler.EffectiveRotationAngle())
        + "\nrotationmatrix3D_frommatrix:  " + util::Format()( rotationmatrix3D_frommatrix.EffectiveRotationAngle())
        + "\nrotationmatrix3D_rand(pi/10): " + util::Format()( rotationmatrix3D_rand.EffectiveRotationAngle())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathRotationMatrix3D

  const ExampleClass::EnumType ExampleMathRotationMatrix3D::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRotationMatrix3D())
  );

} // namespace bcl
