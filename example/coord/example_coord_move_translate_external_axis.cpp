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
#include "coord/bcl_coord_move_translate_external_axis.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_translate_external_axis.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveTranslateExternalAxis :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveTranslateExternalAxis *Clone() const
    {
      return new ExampleCoordMoveTranslateExternalAxis( *this);
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

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      coord::MoveTranslateExternalAxis def_construct;

      // test constructor from range and object
      coord::MoveTranslateExternalAxis translate_construct( 2.0, 10.0, coord::GetAxes().e_Z);

      // test Clone
      util::ShPtr< coord::MoveTranslateExternalAxis> clone_construct
      (
        translate_construct.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< coord::MoveTranslateExternalAxis>(),
        clone_construct->GetClassIdentifier()
      );

    ////////////////
    // operations //
    ////////////////

      // test Move
      protein_model.Transform( math::Inverse( protein_model.GetOrientation()));
      protein_model.Translate( linal::Vector3D( 5.0, 5.0, 5.0));
      translate_construct.Move( protein_model);
      BCL_ExampleIndirectCheck
      (
        !math::EqualWithinTolerance( protein_model.GetCenter().X(), 5.0) &&
          !math::EqualWithinTolerance( protein_model.GetCenter().Y(), 5.0) &&
          math::EqualWithinTolerance( protein_model.GetCenter().Z(), 5.0),
        true,
        "Move"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( translate_construct);

      // read the object back in
      coord::MoveTranslateExternalAxis move_read;
      ReadBCLObject( move_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveTranslateExternalAxis

  const ExampleClass::EnumType ExampleCoordMoveTranslateExternalAxis::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveTranslateExternalAxis())
  );
  
} // namespace bcl
