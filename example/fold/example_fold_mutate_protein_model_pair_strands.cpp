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
#include "fold/bcl_fold_mutate_protein_model_pair_strands.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_pair_strands.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelPairStrands :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelPairStrands *Clone() const
    {
      return new ExampleFoldMutateProteinModelPairStrands( *this);
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
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 6;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

      // join the adjacent strand if any exists
      model.Join( biol::GetSSTypes().STRAND, false);

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain( model.GetChain( 'A'));

      // create new chain
      util::ShPtr< assemble::Chain> new_chain( new assemble::Chain( chain->GetSequence()));

      // create ShPtrs on new SSEs
      util::ShPtr< assemble::SSE> strand_68_76  ( assemble::LocatorSSE( 'A',  68,  76).Locate( model)->HardCopy());
      util::ShPtr< assemble::SSE> strand_82_89  ( assemble::LocatorSSE( 'A',  82,  89).Locate( model)->HardCopy());
      util::ShPtr< assemble::SSE> strand_120_125( assemble::LocatorSSE( 'A', 120, 125).Locate( model)->HardCopy());
      util::ShPtr< assemble::SSE> strand_152_159( assemble::LocatorSSE( 'A', 152, 159).Locate( model)->HardCopy());

      // create varying models
      // create model with two unpaired strands
      assemble::ProteinModel model_with_two_unpaired_strands( new_chain);
      model_with_two_unpaired_strands.Insert( strand_120_125);
      model_with_two_unpaired_strands.Insert( strand_152_159);

      // create model with one unpaired strand and a sheet of two strands
      assemble::ProteinModel model_with_one_sheet_one_unpaired_strand( new_chain);
      model_with_one_sheet_one_unpaired_strand.Insert( strand_68_76);
      model_with_one_sheet_one_unpaired_strand.Insert( strand_82_89);
      model_with_one_sheet_one_unpaired_strand.Insert( strand_120_125);

      // create model with two unpaired strands and a sheet of two strands
      assemble::ProteinModel model_with_one_sheet_two_unpaired_strands( new_chain);
      model_with_one_sheet_two_unpaired_strands.Insert( strand_68_76);
      model_with_one_sheet_two_unpaired_strands.Insert( strand_82_89);
      model_with_one_sheet_two_unpaired_strands.Insert( strand_120_125);
      model_with_one_sheet_two_unpaired_strands.Insert( strand_152_159);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      fold::MutateProteinModelPairStrands mutate_def;

      BCL_MessageStd( "test constructor from a scheme");
      fold::MutateProteinModelPairStrands mutate( "test_scheme");

      BCL_MessageStd( "test copy constructor");
      fold::MutateProteinModelPairStrands mutate_copy( mutate);

      BCL_MessageStd( "test clone");
      util::ShPtr< fold::MutateProteinModelPairStrands> sp_mutate( mutate.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( mutate.GetClassIdentifier(), GetStaticClassName( mutate));

      // test GetScheme
      BCL_ExampleCheck( mutate.GetScheme(), "test_scheme");

    ///////////////
    // operators //
    ///////////////

      // model with two unpaired strands
      {
        // test the operator()
        BCL_MessageStd( "testing operator() with model with two unpaired strands");

        // write the model before
        Proteins::WriteModelToPDB
        (
          model_with_two_unpaired_strands, AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_a_0.pdb")
        );

        // apply the mutate and store the result
        math::MutateResult< assemble::ProteinModel> result_a( mutate( model_with_two_unpaired_strands));

        // make sure the result is successful
        BCL_ExampleIndirectAssert
        (
          result_a.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
        );

        // write the mutate model
        Proteins::WriteModelToPDB
        (
          *result_a.GetArgument(), AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_a_1.pdb")
        );
      }

      // model with one unpaired strand and a sheet of two strands
      {
        // test the operator()
        BCL_MessageStd
        (
          "testing operator() with model with one unpaired strand and a sheet of two strands"
        );

        // write the model before
        Proteins::WriteModelToPDB
        (
          model_with_one_sheet_one_unpaired_strand,
          AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_b_0.pdb")
        );

        // apply the mutate and store the result
        math::MutateResult< assemble::ProteinModel> result_b( mutate( model_with_one_sheet_one_unpaired_strand));

        // make sure the result is successful
        BCL_ExampleIndirectAssert
        (
          result_b.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
        );

        // write the mutate model
        Proteins::WriteModelToPDB
        (
          *result_b.GetArgument(), AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_b_1.pdb")
        );
      }

      // model with two unpaired strands and a sheet of two strands
      {
        // test the operator()
        BCL_MessageStd
        (
          "testing operator() with model with two unpaired strands and a sheet of two strands"
        );

        // write the model before
        Proteins::WriteModelToPDB
        (
          model_with_one_sheet_two_unpaired_strands,
          AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_c_0.pdb")
        );

        // apply the mutate and store the result
        math::MutateResult< assemble::ProteinModel> result_c( mutate( model_with_one_sheet_two_unpaired_strands));

        // make sure the result is successful
        BCL_ExampleIndirectAssert
        (
          result_c.GetArgument().IsDefined(), true, "the mutate_a returned an empty argument!!"
        );

        // write the mutate model
        Proteins::WriteModelToPDB
        (
          *result_c.GetArgument(), AddExampleOutputPathToFilename( mutate, "mutate_pair_strand_c_1.pdb")
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "testing read write")
      WriteBCLObject( mutate);
      fold::MutateProteinModelPairStrands mutate_read;
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_read.GetScheme() == mutate.GetScheme(),
        "The read mutate is different\n" + util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelPairStrands

  const ExampleClass::EnumType ExampleFoldMutateProteinModelPairStrands::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelPairStrands())
  );

} // namespace bcl
