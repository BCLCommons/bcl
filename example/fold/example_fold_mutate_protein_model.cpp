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
#include "fold/bcl_fold_mutate_protein_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_move_translate_random.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model.cpp
  //!
  //! @author weinerbe
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModel :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModel *Clone() const
    {
      return new ExampleFoldMutateProteinModel( *this);
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
      // get a protein model and place it at the origin
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );
      protein_model.Translate( -protein_model.GetCenter());

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
      fold::MutateProteinModel def_construct;
      BCL_ExampleIndirectCheck( def_construct( protein_model).GetArgument().IsDefined(), false, "default constructor");

      // test constructor from move
      const std::string scheme( "protein_mutate_example");
      const fold::MutateProteinModel move_construct( translate_move, scheme);

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
      BCL_ExampleCheck( mutated_model->GetCenter().Z() >= 2.0 || mutated_model->GetCenter().Z() <= -2.0, true);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModel

  const ExampleClass::EnumType ExampleFoldMutateProteinModel::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModel())
  );

} // namespace bcl
