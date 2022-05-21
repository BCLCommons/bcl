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
#include "fold/bcl_fold_mutate_domain_shuffle.h"

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
  //! @example example_fold_mutate_domain_shuffle.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateDomainShuffle :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateDomainShuffle *Clone() const
    {
      return new ExampleFoldMutateDomainShuffle( *this);
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
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain( model.GetChain( 'A'));

      // get the sheet
      util::ShPtrVector< assemble::Domain> sheets_collected( assemble::CollectorSheet().Collect( model));
      util::ShPtr< assemble::Domain> this_sheet( sheets_collected( 1));
      BCL_MessageStd( "Sheet in 2yv8\n" + this_sheet->GetTopology()->GetOrderedIdentification());
      const size_t nr_strands( this_sheet->GetNumberSSEs());

      // initialize variables
      const size_t nr_swaps( 2);
      const std::string scheme( "this_mutate");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      fold::MutateDomainShuffle mutate_def;

      BCL_MessageStd( "test constructor from a number of swaps and a scheme");
      fold::MutateDomainShuffle mutate( nr_swaps, false, scheme);

      BCL_MessageStd( "test constructor from a number of swaps and a scheme");
      fold::MutateDomainShuffle mutate_one( 1, false, scheme);

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< fold::MutateDomainShuffle> sp_mutate( mutate.Clone());
      BCL_Example_Check
      (
        sp_mutate->GetMaxNumberSwaps() == mutate.GetMaxNumberSwaps() &&
        sp_mutate->GetScheme() == mutate.GetScheme(),
        "The cloned mutate is different!\n" + util::Format()( *sp_mutate) + "\nvs\n" + util::Format()( mutate)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( GetStaticClassName( mutate), mutate.GetClassIdentifier());

      // test GetMaxNumberSwaps
      BCL_MessageStd( "test GetMaxNumberSwaps()");
      BCL_ExampleCheck( mutate.GetMaxNumberSwaps(), nr_swaps);

      // test GetScheme
      BCL_MessageStd( "test GetScheme()");
      BCL_ExampleCheck( mutate.GetScheme(), scheme);

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with only one swap");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_a( mutate_one( *this_sheet));
      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_a.GetArgument().IsDefined(), true, "The mutate_one returned an empty argument!!"
      );

      // construct a chain
      util::ShPtr< assemble::Chain>
        sp_new_chain_a( new assemble::Chain( chain->GetSequence(), *mutate_result_a.GetArgument()));
      // construct a model
      assemble::ProteinModel new_model_a( sp_new_chain_a);
      // write out the model
      Proteins::WriteModelToPDB( new_model_a, AddExampleOutputPathToFilename( mutate, "mutate_domain_shuffle_a.pdb"));

      // collect the sheet back
      util::ShPtrVector< assemble::Domain> sheets_a( assemble::CollectorSheet().Collect( new_model_a));

      // make sure there is only one sheet
      BCL_ExampleAssert( sheets_a.GetSize(), 1);

      // print out the sheet found
      BCL_MessageStd
      (
        "Sheet after mutate_one\n" + sheets_a( 0)->GetTopology()->GetOrderedIdentification()
      );

      // make sure the sheet is of the same size
      BCL_ExampleCheck( sheets_a( 0)->GetNumberSSEs(), nr_strands);

      // test with a two swaps
      BCL_MessageStd( "testing operator() with two swaps");
      math::MutateResult< assemble::Domain> mutate_result_b( mutate( *this_sheet));
      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_b.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
      );

      // construct a new chain
      util::ShPtr< assemble::Chain>
        sp_new_chain_b( new assemble::Chain( chain->GetSequence(), *mutate_result_b.GetArgument()));
      // construct a model
      assemble::ProteinModel new_model_b( sp_new_chain_b);
      // print out the model
      Proteins::WriteModelToPDB( new_model_b, AddExampleOutputPathToFilename( mutate, "mutate_domain_shuffle_b.pdb"));

      // collect the sheet back
      util::ShPtrVector< assemble::Domain> sheets_b( assemble::CollectorSheet().Collect( new_model_b));

      // make sure there is only one sheet
      BCL_ExampleAssert( sheets_b.GetSize(), 1);

      // print out the sheet found
      BCL_MessageStd( "Sheet after mutate_one\n" + sheets_b( 0)->GetTopology()->GetOrderedIdentification());
      // make sure the sheet is of the same size
      BCL_Example_Check
      (
        sheets_b( 0)->GetNumberSSEs() == nr_strands,
        "The sheet should have " + util::Format()( nr_strands) + " strands but instead it has " +
        util::Format()( sheets_b( 0)->GetNumberSSEs())
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( mutate);
      fold::MutateDomainShuffle mutate_read;
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_read.GetMaxNumberSwaps() == mutate.GetMaxNumberSwaps() &&
        mutate_read.GetScheme() == mutate.GetScheme(),
        "The read mutate is different\n" + util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateDomainShuffle

  const ExampleClass::EnumType ExampleFoldMutateDomainShuffle::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateDomainShuffle())
  );

} // namespace bcl
