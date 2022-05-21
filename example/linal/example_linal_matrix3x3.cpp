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
#include "linal/bcl_linal_matrix3x3.h"

// includes from bcl - sorted alphabetically
#include "example_check_macros.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix3x3.cpp
  //! @brief this example tests the implementation of the 3x3 matrix.
  //!
  //! @author fischea
  //! @date Oct 11, 2015
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrix3x3 :
    public ExampleInterface
  {

  public:

  //////////
  // data //
  //////////

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to a new ExampleLinalMatrix3x3
    ExampleLinalMatrix3x3 *Clone() const
    {
      return new ExampleLinalMatrix3x3( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create three test matrices
      linal::Matrix3x3< double> mat_1( 0.0);
      mat_1( 0, 0) = 1.0;
      mat_1( 1, 1) = 1.0;
      mat_1( 2, 2) = 1.0;
      linal::Matrix3x3< double> mat_2( 0.0);
      mat_2( 0, 0) = 0.5 / std::sqrt( 0.5);
      mat_2( 0, 1) = 0.5 / std::sqrt( 0.5);
      mat_2( 1, 0) = -0.5 / std::sqrt( 0.5);
      mat_2( 1, 1) = 0.5 / std::sqrt( 0.5);
      mat_2( 2, 2) = 1.0;
      linal::Matrix3x3< double> mat_3( 0.0);
      mat_3( 0, 1) = -1.0;
      mat_3( 1, 0) = 1.0;
      mat_3( 2, 2) = 1.0;
      linal::Matrix3x3< double> mat_4( 0.0);
      mat_4( 0, 0) = -0.5 / std::sqrt( 0.5);
      mat_4( 0, 1) = -0.5 / std::sqrt( 0.5);
      mat_4( 1, 0) = 0.5 / std::sqrt( 0.5);
      mat_4( 1, 1) = -0.5 / std::sqrt( 0.5);
      mat_4( 2, 2) = -1.0;

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( mat_1.GetClassIdentifier(), ( GetStaticClassName< linal::Matrix3x3< double> >()));

    ////////////////
    // operations //
    ////////////////

      // computation of Euler angles between mat_1 and mat_2
      const linal::Vector3D euler_1( linal::Matrix3x3< double>::ComputeEulerAngles( mat_1, mat_2));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_1( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_1( 1)) &&
        math::EqualWithinTolerance( math::g_Pi / 4.0, euler_1( 2)),
        true,
        "Incorrect computation of the euler angles between mat_1 and mat_2."
      );
      const linal::Vector3D euler_1_xyz
      (
        linal::Matrix3x3< double>::ComputeEulerAnglesXYZ( mat_1, mat_2)
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_1_xyz( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_1_xyz( 1)) &&
        math::EqualWithinTolerance( 1.75 * math::g_Pi, euler_1_xyz( 2)),
        true,
        "Incorrect computation of the euler angles between mat_1 and mat_2."
      );

      // computation of Euler angles between mat_1 and mat_3
      const linal::Vector3D euler_2( linal::Matrix3x3< double>::ComputeEulerAngles( mat_1, mat_3));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_2( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_2( 1)) &&
        math::EqualWithinTolerance( -math::g_Pi / 2.0, euler_2( 2)),
        true,
        "Incorrect computation of the euler angles between mat_1 and mat_3."
      );
      const linal::Vector3D euler_2_xyz
      (
        linal::Matrix3x3< double>::ComputeEulerAnglesXYZ( mat_1, mat_3)
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_2_xyz( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_2_xyz( 1)) &&
        math::EqualWithinTolerance( math::g_Pi / 2.0, euler_2_xyz( 2)),
        true,
        "Incorrect computation of the euler angles between mat_1 and mat_3."
      );

      // computation of euler angles between mat_2 and mat_3
      const linal::Vector3D euler_3( linal::Matrix3x3< double>::ComputeEulerAngles( mat_2, mat_3));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_3( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_3( 1)) &&
        math::EqualWithinTolerance( -0.75 * math::g_Pi, euler_3( 2)),
        true,
        "Incorrect computation of the euler angles between mat_2 and mat_3."
      );
      const linal::Vector3D euler_3_xyz
      (
        linal::Matrix3x3< double>::ComputeEulerAnglesXYZ( mat_2, mat_3)
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, euler_3_xyz( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_3_xyz( 1)) &&
        math::EqualWithinTolerance( 0.75 * math::g_Pi, euler_3_xyz( 2)),
        true,
        "Incorrect computation of the euler angles between mat_2 and mat_3."
      );

      // computation of euler angles between mat_1 and mat_4
      const linal::Vector3D euler_4_xyz
      (
        linal::Matrix3x3< double>::ComputeEulerAnglesXYZ( mat_1, mat_4)
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( math::g_Pi, euler_4_xyz( 0)) &&
        math::EqualWithinTolerance( 0.0, euler_4_xyz( 1)) &&
        math::EqualWithinTolerance( 0.75 * math::g_Pi, euler_4_xyz( 2)),
        true,
        "Incorrect computation of the euler angles between mat_1 and mat_4."
      );

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleLinalMatrix3x3

  //! single instance of this class
  const ExampleClass::EnumType ExampleLinalMatrix3x3::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrix3x3())
  );

} // namespace bcl
