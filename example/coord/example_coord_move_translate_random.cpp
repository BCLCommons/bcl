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
#include "coord/bcl_coord_move_translate_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_translate_random.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveTranslateRandom :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveTranslateRandom *Clone() const
    {
      return new ExampleCoordMoveTranslateRandom( *this);
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
      util::ShPtr< assemble::SSE> sse_orig( locator.Locate( protein_model)->Clone());

      util::ShPtr< assemble::SSE> sse;

      const double max_translation( 3.0);
      coord::MoveTranslateRandom move_translate_xyz_internal( max_translation, true);
      coord::MoveTranslateRandom move_translate_xyz_external( max_translation, false);
      coord::MoveTranslateRandom move_translate_y_internal( linal::Vector3D( 0.0, max_translation, 0.0), true);
      coord::MoveTranslateRandom move_translate_y_external( linal::Vector3D( 0.0, max_translation, 0.0), false);

      const double max_translation_norm( math::Sqrt( 3 * math::Sqr( max_translation)));
      sse = sse_orig.HardCopy();
      const linal::Vector3D center_start( sse->GetCenter());
      const double rotation_start( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "starting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "starting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));

      BCL_MessageStd( "applying move_translate_xyz_internal");
      sse = sse_orig.HardCopy();
      move_translate_xyz_internal.Move( *sse);
      const linal::Vector3D center_xyz_internal( sse->GetCenter());
      const double rotation_xyz_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_xyz_internal),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_xyz_internal)
      );
      BCL_Example_Check
      (
        ( center_start - center_xyz_internal).Norm() < max_translation_norm,
        "translation should not be longer than: " + util::Format()( max_translation_norm)
      );

      BCL_MessageStd( "applying move_translate_xyz_external");
      sse = sse_orig.HardCopy();
      move_translate_xyz_external.Move( *sse);
      const linal::Vector3D center_xyz_external( sse->GetCenter());
      const double rotation_xyz_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_xyz_external), "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_xyz_external)
      );
      BCL_Example_Check
      (
        ( center_start - center_xyz_external).Norm() < max_translation_norm, "translation should not be longer than: " + util::Format()( max_translation_norm)
      );

      BCL_MessageStd( "applying move_translate_y_internal");
      sse = sse_orig.HardCopy();
      move_translate_y_internal.Move( *sse);
      const linal::Vector3D center_y_internal( sse->GetCenter());
      const double rotation_y_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_y_internal),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_y_internal)
      );
      BCL_Example_Check
      (
        ( center_start - center_y_internal).Norm() < max_translation_norm,
        "translation should not be longer than: " + util::Format()( max_translation_norm)
      );

      BCL_MessageStd( "applying move_translate_y_external");
      move_translate_y_external.Move( *sse);
      const linal::Vector3D center_y_external( sse->GetCenter());
      const double rotation_y_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rotation_start, rotation_y_external),
        "there should be no change in the rotation: " + util::Format()( rotation_start - rotation_y_external)
      );
      BCL_Example_Check
      (
        ( center_start - center_y_external).Norm() < max_translation_norm,
        "translation should not be longer than: " + util::Format()( max_translation_norm) + " but is: " +
        util::Format()( ( center_start - center_y_external).Norm())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveTranslateRandom

  const ExampleClass::EnumType ExampleCoordMoveTranslateRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveTranslateRandom())
  );

} // namespace bcl
