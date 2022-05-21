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
#include "coord/bcl_coord_move_transform_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_transform_random.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveTransformRandom :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveTransformRandom *Clone() const
    {
      return new ExampleCoordMoveTransformRandom( *this);
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

      coord::MoveTransformRandom move_transform_xyz_internal( 3.0, math::g_Pi, true);
      coord::MoveTransformRandom move_transform_xyz_external( 3.0, math::g_Pi, false);
      coord::MoveTransformRandom move_transform_y_internal( linal::Vector3D( 0.0, 3.0, 0.0), linal::Vector3D( 0.0, math::g_Pi, 0.0), true);
      coord::MoveTransformRandom move_transform_y_external( linal::Vector3D( 0.0, 3.0, 0.0), linal::Vector3D( 0.0, math::g_Pi, 0.0), false);

      sse = sse_orig.HardCopy();
      BCL_MessageStd( "starting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "starting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));

      BCL_MessageStd( "applying move_transform_xyz_internal");
      sse = sse_orig.HardCopy();
      move_transform_xyz_internal.Move( *sse);
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);

      BCL_MessageStd( "applying move_transform_xyz_external");
      sse = sse_orig.HardCopy();
      move_transform_xyz_external.Move( *sse);
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);

      BCL_MessageStd( "applying move_transform_y_internal");
      sse = sse_orig.HardCopy();
      move_transform_y_internal.Move( *sse);
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);

      BCL_MessageStd( "applying move_transform_y_external");
      move_transform_y_external.Move( *sse);
      BCL_MessageStd( "resulting position in space: " + util::Format()( sse->GetCenter()));
      BCL_MessageStd( "resulting orientation in space: " + util::Format()( sse->GetRotation().EffectiveRotationAngle()));
      protein_model.Replace( sse);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveTransformRandom

  const ExampleClass::EnumType ExampleCoordMoveTransformRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveTransformRandom())
  );

} // namespace bcl
