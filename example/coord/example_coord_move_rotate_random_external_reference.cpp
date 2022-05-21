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
#include "coord/bcl_coord_move_rotate_random_external_reference.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_rotate_random_external_reference.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveRotateRandomExternalReference :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveRotateRandomExternalReference *Clone() const
    {
      return new ExampleCoordMoveRotateRandomExternalReference( *this);
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

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      assemble::LocatorSSE locator( 'B', 91, 109);
      const util::SiPtr< const assemble::SSE> located_sse( locator.Locate( protein_model));
      BCL_ExampleIndirectAssert
      (
        located_sse.IsDefined(),
        true,
        "Finding sse" + util::Format()( locator) + " in " + pdb_filename
      );
      util::ShPtr< assemble::SSE> sse_orig( located_sse->Clone());

      util::ShPtr< assemble::SSE> sse;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      coord::MoveRotateRandomExternalReference def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetMethod(),
        coord::MoveRotateRandomExternalReference::e_Internal,
        "default constructor"
      );

      // create other moves
      const double rotation_angle( math::g_Pi / 6);
      const double max_rotation_angle( math::Sqrt( 3.0) * rotation_angle);
      coord::MoveRotateRandomExternalReference move_rotate_xyz_internal( rotation_angle);
      coord::MoveRotateRandomExternalReference move_rotate_xyz_internal_reference
      (
        rotation_angle,
        math::TransformationMatrix3D( 1.0, 1.0, 1.0)
      );
      coord::MoveRotateRandomExternalReference move_rotate_xyz_external
      (
        rotation_angle,
        math::TransformationMatrix3D(),
        coord::MoveRotateRandomExternalReference::e_External
      );
      coord::MoveRotateRandomExternalReference move_rotate_xy_internal_translate
      (
        linal::Vector3D( rotation_angle, rotation_angle, 0.0),
        math::TransformationMatrix3D(),
        coord::MoveRotateRandomExternalReference::e_InternalTranslate
      );
      coord::MoveRotateRandomExternalReference move_rotate_xy_internal_rotate
      (
        linal::Vector3D( rotation_angle, rotation_angle, 0.0),
        math::TransformationMatrix3D(),
        coord::MoveRotateRandomExternalReference::e_InternalRotate
      );

      // test Clone
      util::ShPtr< coord::MoveRotateRandomExternalReference> clone_construct
      (
        move_rotate_xy_internal_translate.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< coord::MoveRotateRandomExternalReference>(),
        clone_construct->GetClassIdentifier()
      );

      // check GetMethod
      BCL_ExampleCheck( move_rotate_xyz_external.GetMethod(), coord::MoveRotateRandomExternalReference::e_External);

      // check SetMethod
      def_construct.SetMethod( coord::MoveRotateRandomExternalReference::e_InternalTranslate);
      BCL_ExampleIndirectCheck
      (
        move_rotate_xyz_external.GetMethod(),
        coord::MoveRotateRandomExternalReference::e_External,
        "SetMethod"
      );

    ///////////////
    // operators //
    ///////////////

      // test internal move with reference at origin
      sse = sse_orig.HardCopy();
      move_rotate_xyz_internal.Move( *sse);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          ( sse->GetCenter() - sse_orig->GetCenter()).Norm(),
          double( 0),
          double( 0.0001)
        ),
        true,
        "Internal XYZ rotation"
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          sse->GetRotation().EffectiveRotationAngle(),
          sse_orig->GetRotation().EffectiveRotationAngle(),
          max_rotation_angle
        ),
        true,
        "Internal XYZ rotation"
      );

      // test internal move with reference at (1.0, 1.0, 1.0)
      sse = sse_orig.HardCopy();
      move_rotate_xyz_internal.Move( *sse);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          ( sse->GetCenter() - linal::Vector3D( 1.0, 1.0, 1.0)).Norm(),
          ( sse_orig->GetCenter() - linal::Vector3D( 1.0, 1.0, 1.0)).Norm(),
          double( 0.0001)
        ),
        true,
        "Internal XYZ rotation w/ reference orientation"
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          sse->GetRotation().EffectiveRotationAngle(),
          sse_orig->GetRotation().EffectiveRotationAngle(),
          max_rotation_angle
        ),
        true,
        "Internal XYZ rotation w/ reference orientation"
      );

      // test external move with reference at origin
      sse = sse_orig.HardCopy();
      move_rotate_xyz_external.Move( *sse);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          sse->GetRotation().EffectiveRotationAngle(),
          sse_orig->GetRotation().EffectiveRotationAngle(),
          max_rotation_angle
        ),
        true,
        "External XYZ rotation"
      );

      // test internal translate move with reference at origin
      sse = sse_orig.HardCopy();
      move_rotate_xy_internal_translate.Move( *sse);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          ( sse->GetCenter() - sse_orig->GetCenter()).Norm(),
          double( 0),
          double( 0.0001)
        ),
        true,
        "Internal translation, XYZ rotation"
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          sse->GetRotation().EffectiveRotationAngle(),
          sse_orig->GetRotation().EffectiveRotationAngle(),
          max_rotation_angle
        ),
        true,
        "Internal translation, XYZ rotation"
      );

      // test internal rotate move with reference at origin
      sse = sse_orig.HardCopy();
      move_rotate_xy_internal_rotate.Move( *sse);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          sse->GetRotation().EffectiveRotationAngle(),
          sse_orig->GetRotation().EffectiveRotationAngle(),
          max_rotation_angle
        ),
        true,
        "Internal rotation, XYZ rotation"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( move_rotate_xy_internal_translate);

      // read the object back in
      coord::MoveRotateRandomExternalReference move_read;
      ReadBCLObject( move_read);

      BCL_ExampleIndirectCheck
      (
        move_read.GetMethod(),
        move_rotate_xy_internal_translate.GetMethod(),
        "read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveRotateRandomExternalReference

  const ExampleClass::EnumType ExampleCoordMoveRotateRandomExternalReference::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveRotateRandomExternalReference())
  );

} // namespace bcl
