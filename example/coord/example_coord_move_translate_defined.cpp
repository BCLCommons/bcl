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
#include "coord/bcl_coord_move_translate_defined.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_translate_defined.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveTranslateDefined :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveTranslateDefined *Clone() const
    {
      return new ExampleCoordMoveTranslateDefined( *this);
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

      // create LocatorSSE to an SSE "sse_orig" and initialize with Clone of the SSE from "protein_model"
      assemble::LocatorSSE locator( 'B', 91, 109);
      util::ShPtr< assemble::SSE> sse_orig( locator.Locate( protein_model)->Clone());

      // create ShPtr to an SSE
      util::ShPtr< assemble::SSE> sse;

      const linal::Vector3D translation( 10.0, 15.0, 20.0);
      const linal::Vector3D zero_translation( 0.0, 0.0, 0.0);

      // create moves
      coord::MoveTranslateDefined move_translate_xyz_internal
      (
        translation, true
      );
      coord::MoveTranslateDefined move_translate_xyz_external
      (
        translation, false
      );
      coord::MoveTranslateDefined move_translate_zero_internal
      (
        zero_translation, true
      );
      coord::MoveTranslateDefined move_translate_zero_external
      (
        zero_translation, false
      );

      // set "sse" to hard copy of "sse_orig"
      sse = sse_orig.HardCopy();

      // create Vector3D "center_start" which is initialized to the center of "sse"
      const linal::Vector3D center_start( sse->GetCenter());

      // create double "rotation_start" and initialize to rotation angle of "sse"
      const double rotation_start( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "starting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd
      (
        "starting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle())
      );

      // test move with defined translation using internal coordinates of "sse"
      BCL_MessageStd( "test move with defined translation using internal coordinates");
      sse = sse_orig.HardCopy();
      move_translate_xyz_internal.Move( *sse);
      // create Vector3D "center_xyz_internal" and initialize to the center of the moved "sse"
      const linal::Vector3D center_xyz_internal( sse->GetCenter());
      // create double "rotation_xyz_internal" and set to rotation angle of moved "sse"
      const double rotation_xyz_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd
      (
        "resulting rotation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle())
      );
      // put the moved "sse" into "protein_model"
      protein_model.Replace( sse);
      // check that no rotation occured
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_xyz_internal),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_xyz_internal)
      );
      // check that the correct translation occured
      // (center_start - center_xyz_internal).Norm() == translation.Norm()
      BCL_Example_Check
      (
        math::EqualWithinTolerance( center_start + linal::Vector3D( translation).Rotate( sse_orig->GetRotation()), center_xyz_internal),
        "resulting should be: " + util::Format()( center_start + linal::Vector3D( translation).Rotate( sse_orig->GetRotation())) + " but is: " + util::Format()( center_xyz_internal)
      );

      // test move with defined translation using external coordinates
      BCL_MessageStd( "test move with defined translation using external coordinates");
      // set "sse" to copy of "sse_orig"
      sse = sse_orig.HardCopy();
      // do the move
      move_translate_xyz_external.Move( *sse);
      // create Vector3D "center_xyz_external" and initialize with center of moved "sse"
      const linal::Vector3D center_xyz_external( sse->GetCenter());
      // create double "rotation_xyz_external" and initialize with rotation of moved "sse"
      const double rotation_xyz_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd
      (
        "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle())
      );
      // put "sse" into "protein_model"
      protein_model.Replace( sse);
      // make sure no rotation occured
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_xyz_external),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_xyz_external)
      );
      // make sure translation occured correctly
      BCL_Example_Check
      (
        math::EqualWithinTolerance( ( center_start - center_xyz_external).Norm(), translation.Norm()),
        "translation should be longer: " + util::Format()( center_start - center_xyz_external) +
        "\nexpected but result:\n" + util::Format()( translation) +
        "\nlength: " + util::Format()( ( center_start - center_xyz_external).Norm())
        + " != " + util::Format()( translation.Norm())
      );

      // test that no translation occurs if all values of translation vector are zero for internal coordinates
      BCL_MessageStd
      (
        "test that no translation occurs if all values of translation vector are zero for internal coordinates"
      );
      // set sse to hard copy of "sse_orig"
      sse = sse_orig.HardCopy();
      // move "sse"
      move_translate_zero_internal.Move( *sse);
      // create Vector3D "center_zero_internal" and initialize with center of moved "sse"
      const linal::Vector3D center_zero_internal( sse->GetCenter());
      // create double "rotation_zero_internal" and initialize with rotation of moved "sse"
      const double rotation_zero_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd
      (
        "starting position in space" + util::Format()( center_start) +
        "resulting position in space: " + util::Format()( center_zero_internal)
      );
      BCL_MessageStd
      (
        "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle())
      );
      // put "sse" into "protein_model"
      protein_model.Replace( sse);
      // check that no rotation occurs
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_zero_internal),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_zero_internal)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( center_start, center_zero_internal, 0.01),
        "translation should not have moved : " + util::Format()( center_zero_internal)
      );

      // test that no translation occurs if all values of translation vector are zero for external coordinates
      BCL_MessageStd
      (
        "test that no translation occurs if all values of translation vector are zero for external coordinates"
      );
      // move "sse"
      move_translate_zero_external.Move( *sse);
      // create Vector3D "center_zero_external" an initialize to center of "sse"
      const linal::Vector3D center_zero_external( sse->GetCenter());
      // create double "rotation_zero_external" and initialize to rotation angle of "sse"
      const double rotation_zero_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd
      (
        "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle())
      );
      // put "sse" into "protein_model"
      protein_model.Replace( sse);
      // test that no rotation occurs
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_zero_external),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_zero_external)
      );
      // test that no translation occurs
      BCL_Example_Check
      (
        math::EqualWithinTolerance( center_start, center_zero_external),
        "tranalation should not have changed: " + util::Format()( center_zero_external)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveTranslateDefined

  const ExampleClass::EnumType ExampleCoordMoveTranslateDefined::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveTranslateDefined())
  );

} // namespace bcl
