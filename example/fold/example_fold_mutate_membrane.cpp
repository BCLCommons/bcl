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
#include "fold/bcl_fold_mutate_membrane.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "coord/bcl_coord_move_translate_random.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_membrane.cpp
  //!
  //! @author weinerbe
  //! @date Nov 28, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateMembrane :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateMembrane *Clone() const
    {
      return new ExampleFoldMutateMembrane( *this);
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
      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );

      // add a membrane
      util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, util::CloneToShPtr( biol::Membrane()));

      // build move
      const coord::MoveTranslateRandom translate_move
      (
        linal::Vector3D( 0.0, 0.0, 2.0),
        linal::Vector3D( 0.0, 0.0, 10.0),
        false
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateMembrane def_construct;
      BCL_ExampleIndirectCheck( def_construct( protein_model).GetArgument().IsDefined(), false, "default constructor");

      // test constructor from move
      const std::string scheme( "membrane_mutate_example");
      const fold::MutateMembrane move_construct( translate_move, scheme);

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      BCL_ExampleCheck( move_construct.GetScheme(), scheme);

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtr< assemble::ProteinModel> mutated_model( move_construct( protein_model).GetArgument());
      BCL_ExampleAssert( mutated_model.IsDefined(), true);
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        mutated_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );
      BCL_ExampleCheck( sp_membrane->GetCenter().Z() >= 2.0 || sp_membrane->GetCenter().Z() <= -2.0, true);

    //////////////////////
    // input and output //
    //////////////////////

      // test write and read
      WriteBCLObject( move_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleCheck( def_construct.GetScheme(), move_construct.GetScheme());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateMembrane

  const ExampleClass::EnumType ExampleFoldMutateMembrane::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateMembrane())
  );

} // namespace bcl
