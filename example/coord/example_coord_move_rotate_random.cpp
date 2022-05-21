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
#include "coord/bcl_coord_move_rotate_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_rotate_random.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveRotateRandom :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveRotateRandom *Clone() const
    {
      return new ExampleCoordMoveRotateRandom( *this);
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

      const double rot_angle( math::g_Pi / 18.0);
      const double max_rot_angle( math::Sqrt( 3.0) * rot_angle);
      coord::MoveRotateRandom move_rotate_xyz_internal( rot_angle, true);
      coord::MoveRotateRandom move_rotate_xyz_external( rot_angle, false);
      coord::MoveRotateRandom move_rotate_y_internal( linal::Vector3D( 0.0, rot_angle, 0.0), true);
      coord::MoveRotateRandom move_rotate_y_external( linal::Vector3D( 0.0, rot_angle, 0.0), false);

      sse = sse_orig.HardCopy();
      const linal::Vector3D center_start( sse->GetCenter());
      const double rotation_start( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "starting position in space: " + util::Format()( center_start));
      BCL_MessageStd( "starting orientation in space: " + util::Format()( rotation_start));

      BCL_MessageStd( "applying move_rotate_xyz_internal");
      sse = sse_orig.HardCopy();
      move_rotate_xyz_internal.Move( *sse);
      const linal::Vector3D center_xyz_internal( sse->GetCenter());
      const double rotation_xyz_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( center_xyz_internal));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( rotation_xyz_internal));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( ( center_start - center_xyz_internal).Norm(), double( 0), double( 0.0001)),
        "the center shall stay the same after applying internal rotation"
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( rotation_xyz_internal, rotation_start, max_rot_angle),
        "the maximal change in the rotation shall not exceed: " + util::Format()( max_rot_angle)
      );

      BCL_MessageStd( "applying move_rotate_xyz_external");
      sse = sse_orig.HardCopy();
      move_rotate_xyz_external.Move( *sse);
      const linal::Vector3D center_xyz_external( sse->GetCenter());
      const double rotation_xyz_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( center_xyz_external));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( rotation_xyz_external));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( rotation_xyz_external, rotation_start, max_rot_angle),
        "the maximal change in the rotation shall not exceed: " + util::Format()( max_rot_angle) + " but is " +
        util::Format()( math::Absolute( rotation_xyz_external - rotation_start))
      );

      BCL_MessageStd( "applying move_rotate_y_internal");
      sse = sse_orig.HardCopy();
      move_rotate_y_internal.Move( *sse);
      const linal::Vector3D center_y_internal( sse->GetCenter());
      const double rotation_y_internal( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( center_y_internal));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( rotation_y_internal));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( ( center_start - center_y_internal).Norm(), double( 0), double( 0.0001)),
        "the center shall stay the same after applying internal rotation"
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( rotation_y_internal, rotation_start, max_rot_angle),
        "the maximal change in the rotation shall not exceed: " + util::Format()( max_rot_angle)
      );

      BCL_MessageStd( "applying move_rotate_y_external");
      sse = sse_orig.HardCopy();
      move_rotate_y_external.Move( *sse);
      const linal::Vector3D center_y_external( sse->GetCenter());
      const double rotation_y_external( sse->GetRotation().EffectiveRotationAngle());
      BCL_MessageStd( "resulting position in space: " + util::Format()( center_y_external));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( rotation_y_external));
      protein_model.Replace( sse);
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( rotation_y_external, rotation_start, max_rot_angle),
        "the maximal change in the rotation shall not exceed: " + util::Format()( max_rot_angle)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveRotateRandom

  const ExampleClass::EnumType ExampleCoordMoveRotateRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveRotateRandom())
  );

} // namespace bcl
