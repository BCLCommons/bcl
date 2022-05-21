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
#include "fold/bcl_fold_mutate_sheet_register_fix.h"

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
  //! @example example_fold_mutate_sheet_register_fix.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetRegisterFix :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetRegisterFix *Clone() const
    {
      return new ExampleFoldMutateSheetRegisterFix( *this);
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

    ///////////////
    // read 2yv8 //
    ///////////////

      BCL_MessageStd( "reading disordered sheet from 2yv8");
      // initialize pdb filename
      const std::string pdb_filename_a( AddExampleInputPathToFilename( e_Biology, "2yv8_sheet_disordered.pdb"));

      // get the protein model
      assemble::ProteinModel model_a( Proteins::GetModel( pdb_filename_a, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

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
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "2qv3_sheet_disordered.pdb"));

      // get the protein model
      assemble::ProteinModel model_b( Proteins::GetModel( pdb_filename_b, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

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
      fold::MutateSheetRegisterFix mutate;

      // test clone
      BCL_MessageStd( "testing clone()");
      util::ShPtr< math::MutateInterface< assemble::Domain> > sp_mutate( mutate.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_MessageStd
      (
        "This class has the following identifier " + mutate.GetClassIdentifier()
      );
      BCL_ExampleCheck( sp_mutate->GetClassIdentifier(), "bcl::fold::MutateSheetRegisterFix");

      // test the GetStaticClassName
      BCL_ExampleCheck( GetStaticClassName( mutate), "bcl::fold::MutateSheetRegisterFix");

    ///////////////
    // operators //
    ///////////////

      // test operator()
      BCL_MessageStd( "testing operator() with anti-parallel sheet from 2yv8");
      math::MutateResult< assemble::Domain> mutate_result_a( mutate( *sheet_a));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_a.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
      );

      // construct a chain from the sheet
      assemble::Chain new_chain_a( chain_a->GetSequence(), *mutate_result_a.GetArgument());

      // write out the model
      Proteins::WriteChainToPDB( new_chain_a, AddExampleOutputPathToFilename( mutate, "mutate_sheet_register_fix_a.pdb"));

      // test operator()
      BCL_MessageStd( "testing operator() with parallel sheet from 2qv3");
      math::MutateResult< assemble::Domain> mutate_result_b( mutate( *sheet_b));

      // make sure the return result is not empty
       BCL_ExampleIndirectAssert
       (
         mutate_result_b.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
       );

      // construct a chain from the sheet
      assemble::Chain new_chain_b( chain_b->GetSequence(), *mutate_result_b.GetArgument());

      // write out the model
      Proteins::WriteChainToPDB( new_chain_b, AddExampleOutputPathToFilename( mutate, "mutate_sheet_register_fix_b.pdb"));

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetRegisterFix

  const ExampleClass::EnumType ExampleFoldMutateSheetRegisterFix::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetRegisterFix())
  );
  
} // namespace bcl
