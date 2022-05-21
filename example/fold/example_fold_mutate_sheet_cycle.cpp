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
#include "fold/bcl_fold_mutate_sheet_cycle.h"

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
  //! @example example_fold_mutate_sheet_cycle.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetCycle :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetCycle *Clone() const
    {
      return new ExampleFoldMutateSheetCycle( *this);
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

      // build strand counts
      const math::Range< size_t> nd_1_1( 1, 1);
      const math::Range< size_t> nd_1_2( 1, 2);
      const math::Range< size_t> nd_2_2( 2, 2);
      const math::Range< size_t> nd_2_3( 2, 3);
      const math::Range< size_t> nd_4_4( 4, 4);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      fold::MutateSheetCycle mutate_def;

      // constructor from boolean and strand counts
      BCL_MessageStd( "test constructor");
      fold::MutateSheetCycle mutate_a( true, true, false, nd_1_1);
      BCL_Example_Check
      (
        mutate_a.GetPreserveSheetGeometry() &&
        mutate_a.GetPreserveStrandOrientations() &&
        !mutate_a.GetRotateSubset() &&
        mutate_a.GetNumberRotations() == nd_1_1 &&
        mutate_a.GetSubsetSize() == nd_2_3,
        "The default constructor failed " + util::Format()( mutate_a)
      );

      // constructor from boolean and strand counts
      fold::MutateSheetCycle mutate_b( true, false, false, nd_2_2);

      // constructor from boolean and strand counts
      fold::MutateSheetCycle mutate_c( true, true, false, nd_2_2);

      // constructor from boolean and strand counts
      fold::MutateSheetCycle mutate_d( true, true, true, nd_1_1, nd_4_4);

      // constructor from boolean and strand counts
      fold::MutateSheetCycle mutate_e( false, false, false, nd_2_2);

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< fold::MutateSheetCycle> sp_mutate( mutate_a.Clone());
      BCL_Example_Check
      (
        sp_mutate->GetPreserveSheetGeometry()      == mutate_a.GetPreserveSheetGeometry() &&
        sp_mutate->GetPreserveStrandOrientations() == mutate_a.GetPreserveStrandOrientations() &&
        sp_mutate->GetRotateSubset()               == mutate_a.GetRotateSubset() &&
        sp_mutate->GetNumberRotations()            == mutate_a.GetNumberRotations() &&
        sp_mutate->GetSubsetSize()                 == mutate_a.GetSubsetSize(),
        "The cloned mutate is different!\n" + util::Format()( *sp_mutate) + "\nvs\n" + util::Format()( mutate_a)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( GetStaticClassName( mutate_a), mutate_a.GetClassIdentifier());

      // test GetPreserveSheetGeometry()
      BCL_MessageStd( "test GetPreserveSheetGeometry()");
      BCL_Example_Check
      (
        !mutate_e.GetPreserveSheetGeometry(),
        "The GetPreserveSheetGeometry() should return false for mutate_e not true!"
      );

      // test GetPreserveStrandOrientations()
      BCL_MessageStd( "test GetPreserveStrandOrientations()");
      BCL_Example_Check
      (
        !mutate_e.GetPreserveStrandOrientations(),
        "The GetPreserveStrandOrientations() should return false for mutate_e not true!"
      );

      // test GetRotateSubset()
      BCL_MessageStd( "test GetRotateSubset()");
      BCL_Example_Check
      (
        !mutate_e.GetRotateSubset(),
        "The GetRotateSubset() should return false for mutate_e not true!"
      );

      // test GetNumberRotations()
      BCL_MessageStd( "test GetNumberRotations()");
      BCL_Example_Check
      (
        mutate_c.GetNumberRotations() == nd_2_2,
        "The GetNumberRotations() for mutate_c should return\n" + util::Format()( nd_2_2)
        + "\nnot\n" + util::Format()( mutate_c.GetNumberRotations())
      );

      // test GetSubsetSize()
      BCL_MessageStd( "test GetSubsetSize()");
      BCL_Example_Check
      (
        mutate_c.GetSubsetSize() == nd_2_3,
        "The GetSubsetSize() for mutate_c should return\n" + util::Format()( nd_2_3)
        + "\nnot\n" + util::Format()( mutate_c.GetSubsetSize())
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator()");

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_a");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_a( mutate_a( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_a.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_a.GetArgument()),
        AddExampleOutputPathToFilename( mutate_a, "mutate_sheet_cycle_a.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_b");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_b( mutate_b( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_b.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_b.GetArgument()),
        AddExampleOutputPathToFilename( mutate_b, "mutate_sheet_cycle_b.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_c");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_c( mutate_c( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_c.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_c.GetArgument()),
        AddExampleOutputPathToFilename( mutate_c, "mutate_sheet_cycle_c.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_d");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_d( mutate_d( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_d.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_d.GetArgument()),
        AddExampleOutputPathToFilename( mutate_d, "mutate_sheet_cycle_d.pdb")
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_e");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_e( mutate_e( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_e.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
      );

      // write out the model
      Proteins::WriteChainToPDB
      (
        assemble::Chain( chain->GetSequence(), *mutate_result_e.GetArgument()),
        AddExampleOutputPathToFilename( mutate_e, "mutate_sheet_cycle_e.pdb")
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( mutate_a);
      fold::MutateSheetCycle mutate_read;
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_a.GetPreserveSheetGeometry()      == mutate_read.GetPreserveSheetGeometry() &&
        mutate_a.GetPreserveStrandOrientations() == mutate_read.GetPreserveStrandOrientations() &&
        mutate_a.GetRotateSubset()               == mutate_read.GetRotateSubset() &&
        mutate_a.GetNumberRotations()            == mutate_read.GetNumberRotations() &&
        mutate_a.GetSubsetSize()                 == mutate_read.GetSubsetSize(),
        "The read mutate is different\n" + util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate_a)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetCycle

  const ExampleClass::EnumType ExampleFoldMutateSheetCycle::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetCycle())
  );

} // namespace bcl
