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
#include "fold/bcl_fold_mutate_sheet_register_shift.h"

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
  //! @example example_fold_mutate_sheet_register_shift.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetRegisterShift :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetRegisterShift *Clone() const
    {
      return new ExampleFoldMutateSheetRegisterShift( *this);
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
      // construct min_sse_sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 0;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 999;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;

      // initialize probabilities
      storage::VectorND< 2, double> regular_prob( 1.0, 0.0);
      storage::VectorND< 2, double> flip_prob( 0.0, 1.0);

    ///////////////
    // read 2yv8 //
    ///////////////

      BCL_MessageStd( "reading disordered sheet from 2yv8");
      // initialize pdb filename
      const std::string pdb_filename_a( AddExampleInputPathToFilename( e_Biology, "2yv8_sheet.pdb"));

      // get the protein model
      assemble::ProteinModel model_a
      (
        Proteins::GetModel( pdb_filename_a, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain_a( model_a.GetChain( 'A'));

      // get the sheet
      util::ShPtr< assemble::Domain> sheet_a( assemble::CollectorSheet().Collect( model_a)( 0));
      BCL_MessageStd( "Sheet in 2yv8\n" + sheet_a->GetIdentification());

    ///////////////
    // read 2qv3 //
    ///////////////

      BCL_MessageStd( "reading disordered sheet from 2qv3");
      // initialize pdb filename
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "2qv3_sheet.pdb"));

      // get the protein model
      assemble::ProteinModel model_b
      (
        Proteins::GetModel( pdb_filename_b, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain_b( model_b.GetChain( 'A'));

      // get the sheet
      util::ShPtr< assemble::Domain> sheet_b( assemble::CollectorSheet().Collect( model_b)( 0));
      BCL_MessageStd( "Sheet in 2qv3\n" + sheet_b->GetIdentification());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      fold::MutateSheetRegisterShift mutate_def;

      // test default constructor
      BCL_MessageStd( "testing constructor from shift probabilities - only regular shifts");
      fold::MutateSheetRegisterShift mutate( regular_prob);

      // test default constructor
      BCL_MessageStd( "testing constructor from shift probabilities - only flip shifts");
      fold::MutateSheetRegisterShift mutate_flip( flip_prob);

      // test clone
      BCL_MessageStd( "testing clone()");
      util::ShPtr< fold::MutateSheetRegisterShift> sp_mutate( mutate.Clone());
      BCL_Example_Check
      (
        sp_mutate->GetShiftProbabilities() == mutate.GetShiftProbabilities(),
        "The cloned mutate is different\n" + util::Format()( *sp_mutate) + "\nvs\n" + util::Format()( mutate)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_MessageStd
      (
        "This class has the following identifier " + mutate.GetClassIdentifier()
      );

      // test the GetStaticClassName
      BCL_ExampleCheck( sp_mutate->GetClassIdentifier(), "bcl::fold::MutateSheetRegisterShift");
      BCL_ExampleCheck( GetStaticClassName( mutate), "bcl::fold::MutateSheetRegisterShift");

      // test GetShiftProbabilities()
      BCL_MessageStd( "testing GetShiftProbabilities()");
      BCL_ExampleCheck( mutate.GetShiftProbabilities(), regular_prob);

    ///////////////
    // operators //
    ///////////////

    ///////////////////
    // regular shift //
    ///////////////////

      {
        BCL_MessageStd( "testing operator() with regular shift");

        // test operator()
        BCL_MessageStd( "testing operator() with anti-parallel sheet from 2yv8");
        math::MutateResult< assemble::Domain> mutate_result_a( mutate( *sheet_a));

        // make sure the return result is not empty
        BCL_ExampleIndirectAssert
        (
          mutate_result_a.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
        );

        // construct a chain from the sheet
        assemble::Chain new_chain_a( chain_a->GetSequence(), *mutate_result_a.GetArgument());

        // write out the model
        Proteins::WriteChainToPDB
        (
          new_chain_a, AddExampleOutputPathToFilename( mutate, "mutate_sheet_register_shift_a.pdb")
        );

        // test operator()
        BCL_MessageStd( "testing operator() with parallel sheet from 2qv3");
        math::MutateResult< assemble::Domain> mutate_result_b( mutate( *sheet_b));

        // make sure the return result is not empty
        BCL_ExampleIndirectAssert
        (
          mutate_result_b.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
        );

        // construct a chain from the sheet
        assemble::Chain new_chain_b( chain_b->GetSequence(), *mutate_result_b.GetArgument());

        // write out the model
        Proteins::WriteChainToPDB
        (
          new_chain_b, AddExampleOutputPathToFilename( mutate, "mutate_sheet_register_shift_b.pdb")
        );
      }

    ////////////////
    // flip shift //
    ////////////////

      {
        BCL_MessageStd( "testing operator() with flip shift");

        // test operator()
        BCL_MessageStd( "testing operator() with anti-parallel sheet from 2yv8");
        math::MutateResult< assemble::Domain> mutate_result_a( mutate_flip( *sheet_a));

        // make sure the return result is not empty
        BCL_ExampleIndirectAssert
        (
          mutate_result_a.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
        );

        // construct a chain from the sheet
        assemble::Chain new_chain_a( chain_a->GetSequence(), *mutate_result_a.GetArgument());

        // write out the model
        Proteins::WriteChainToPDB
        (
          new_chain_a, AddExampleOutputPathToFilename( mutate_flip, "mutate_sheet_register_shift_flip_a.pdb")
        );

        // test operator()
        BCL_MessageStd( "testing operator() with parallel sheet from 2qv3");
        math::MutateResult< assemble::Domain> mutate_result_b( mutate_flip( *sheet_b));

        // make sure the return result is not empty
        BCL_ExampleIndirectAssert
        (
          mutate_result_b.GetArgument().IsDefined(), true, "The mutate returned an empty argument"
        );

        // construct a chain from the sheet
        assemble::Chain new_chain_b( chain_b->GetSequence(), *mutate_result_b.GetArgument());

        // write out the model
        Proteins::WriteChainToPDB
        (
          new_chain_b, AddExampleOutputPathToFilename( mutate_flip, "mutate_sheet_register_shift_flip_b.pdb")
        );
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "testing read/write");
      WriteBCLObject( mutate);
      fold::MutateSheetRegisterShift mutate_read;
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_read.GetShiftProbabilities() == mutate.GetShiftProbabilities(),
        "The read mutate is different\n" +
        util::Format()( mutate_read.GetShiftProbabilities()) + "\nvs\n" + util::Format()( mutate.GetShiftProbabilities())
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetRegisterShift

  const ExampleClass::EnumType ExampleFoldMutateSheetRegisterShift::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetRegisterShift())
  );
  
} // namespace bcl
