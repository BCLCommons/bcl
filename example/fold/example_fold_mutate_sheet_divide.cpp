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
#include "fold/bcl_fold_mutate_sheet_divide.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sheet_divide.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetDivide :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetDivide *Clone() const
    {
      return new ExampleFoldMutateSheetDivide( *this);
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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2yv8_sheet_ideal.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));

      // join the adjacent strand if any exists
      model.Join( biol::GetSSTypes().STRAND, false);

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain( model.GetChain( 'A'));

      // collect the sheets and create pointer to second one
      util::ShPtrVector< assemble::Domain> sheets( assemble::CollectorSheet().Collect( model));
      util::ShPtr< assemble::Domain> this_sheet( sheets( 0));

      // initialize min and max translations
      linal::Vector3D min_translations_x( 10.0,  0.0,  0.0);
      linal::Vector3D max_translations_x( 10.0,  0.0,  0.0);
      linal::Vector3D min_translations_y(  0.0, 10.0,  0.0);
      linal::Vector3D max_translations_y(  0.0, 10.0,  0.0);
      linal::Vector3D min_translations_z(  0.0,  0.0, 10.0);
      linal::Vector3D max_translations_z(  0.0,  0.0, 10.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor from x translations
      BCL_MessageStd( "test constructor");
      fold::MutateSheetDivide mutate_x( 5, 2, false, min_translations_x, max_translations_x);
      BCL_Example_Check
      (
        mutate_x.GetMinSheetSize()        == 5 &&
        mutate_x.GetMinDividedSheetSize() == 2 &&
        mutate_x.GetFormBetaSandwich()    == false &&
        mutate_x.GetMinTranslations()     == min_translations_x &&
        mutate_x.GetMaxTranslations()     == max_translations_x,
        "The constructor for mutate_x failed " + util::Format()( mutate_x)
      );

      // constructor from y translations
      fold::MutateSheetDivide mutate_y( 5, 2, false, min_translations_y, max_translations_y);
      BCL_Example_Check
      (
        mutate_y.GetMinSheetSize()        == 5 &&
        mutate_y.GetMinDividedSheetSize() == 2 &&
        mutate_y.GetFormBetaSandwich()    == false &&
        mutate_y.GetMinTranslations()     == min_translations_y &&
        mutate_y.GetMaxTranslations()     == max_translations_y,
        "The constructor for mutate_y failed " + util::Format()( mutate_y)
      );

      // constructor from z translations
      fold::MutateSheetDivide mutate_z( 4, 1, false, min_translations_z, max_translations_z);
      BCL_Example_Check
      (
        mutate_z.GetMinSheetSize()        == 4 &&
        mutate_z.GetMinDividedSheetSize() == 1 &&
        mutate_z.GetFormBetaSandwich()    == false &&
        mutate_z.GetMinTranslations()     == min_translations_z &&
        mutate_z.GetMaxTranslations()     == max_translations_z,
        "The constructor for mutate_z failed " + util::Format()( mutate_z)
      );

      // constructor for forming beta sandwiches
      fold::MutateSheetDivide mutate_beta_sandwich( 5, 2, true);
      BCL_Example_Check
      (
        mutate_beta_sandwich.GetMinSheetSize()        == 5 &&
        mutate_beta_sandwich.GetMinDividedSheetSize() == 2 &&
        mutate_beta_sandwich.GetFormBetaSandwich()    == true,
        "The constructor for mutate_beta_sandwich " + util::Format()( mutate_beta_sandwich)
      );

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< fold::MutateSheetDivide> sp_mutate( mutate_x.Clone());
      BCL_Example_Check
      (
        sp_mutate->GetMinSheetSize()        == mutate_x.GetMinSheetSize() &&
        sp_mutate->GetMinDividedSheetSize() == mutate_x.GetMinDividedSheetSize() &&
        sp_mutate->GetMinTranslations()     == mutate_x.GetMinTranslations() &&
        sp_mutate->GetMaxTranslations()     == mutate_x.GetMaxTranslations(),
        "The cloned mutate is different!\n" + util::Format()( *sp_mutate) + "\nvs\n" + util::Format()( mutate_x)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( mutate_x.GetClassIdentifier(), GetStaticClassName( mutate_x));

      // test GetMinSheetSize()
      BCL_ExampleCheck( mutate_x.GetMinSheetSize(), 5);

      // test GetMinDividedSheetSize()
      BCL_ExampleCheck( mutate_x.GetMinDividedSheetSize(), 2);

      // test GetMinTranslations()
      BCL_MessageStd( "test GetMinTranslations()");
      BCL_Example_Check
      (
        mutate_x.GetMinTranslations() == min_translations_x,
        "The GetMinTranslations() for mutate_x should return\n" + util::Format()( min_translations_x) +
        "\nnot\n" + util::Format()( mutate_x.GetMinTranslations())
      );

      // test GetMaxTranslations()
      BCL_MessageStd( "test GetMaxTranslations()");
      BCL_Example_Check
      (
        mutate_x.GetMaxTranslations() == min_translations_x,
        "The GetMaxTranslations() for mutate_x should return\n" + util::Format()( min_translations_x) +
        "\nnot\n" + util::Format()( mutate_x.GetMaxTranslations())
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator()");

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_x");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_x( mutate_x( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_x.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_x.GetArgument()),
        AddExampleOutputPathToFilename( mutate_x, "mutate_sheet_divide_x.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_y");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_y( mutate_y( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_y.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_y.GetArgument()),
        AddExampleOutputPathToFilename( mutate_y, "mutate_sheet_divide_y.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_z");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_z( mutate_z( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_z.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_z.GetArgument()),
        AddExampleOutputPathToFilename( mutate_z, "mutate_sheet_divide_z.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_beta_sandwich");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_sandwich( mutate_beta_sandwich( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_sandwich.GetArgument().IsDefined(), true,
        "The mutate_result_sandwich returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_sandwich.GetArgument()),
        AddExampleOutputPathToFilename( mutate_beta_sandwich, "mutate_sheet_divide_sandwich.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( mutate_z);
      fold::MutateSheetDivide mutate_read( 4, 1, false);
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_z.GetMinSheetSize()        == mutate_read.GetMinSheetSize() &&
        mutate_z.GetMinDividedSheetSize() == mutate_read.GetMinDividedSheetSize() &&
        mutate_z.GetMinTranslations()     == mutate_read.GetMinTranslations() &&
        mutate_z.GetMaxTranslations()     == mutate_read.GetMaxTranslations(),
        "The read mutate is different\n" + util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate_z)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetDivide

  const ExampleClass::EnumType ExampleFoldMutateSheetDivide::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetDivide())
  );

} // namespace bcl
