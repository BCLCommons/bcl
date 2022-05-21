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
#include "fold/bcl_fold_mutate_sheet_twist.h"

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
  //! @example example_fold_mutate_sheet_twist.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetTwist :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetTwist *Clone() const
    {
      return new ExampleFoldMutateSheetTwist( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2CRT.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));

      // join the adjacent strand if any exists
      model.Join( biol::GetSSTypes().STRAND, false);

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain( model.GetChain( 'A'));

      // collect the sheets and create pointer to second one
      util::ShPtrVector< assemble::Domain> sheets( assemble::CollectorSheet().Collect( model));
      util::ShPtr< assemble::Domain> this_sheet( sheets( 1));

      // set expected angle deviation
      const math::Range< double> angle_range( 0.0, math::Angle::Radian( 20.0));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      fold::MutateSheetTwist mutate_def;

      // constructor from an angle deviation
      BCL_MessageStd( "test constructor from an angle deviation");
      fold::MutateSheetTwist mutate( angle_range);

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< fold::MutateSheetTwist> sp_mutate( mutate.Clone());
      BCL_ExampleIndirectCheck
      (
        sp_mutate->GetAngleRange(), mutate.GetAngleRange(), "Clone"
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( mutate.GetClassIdentifier(), GetStaticClassName( mutate));

      // test GetMaxNumberSwaps
      BCL_ExampleCheck( mutate.GetAngleRange(), angle_range);

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator()");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result( mutate( *this_sheet));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert( mutate_result.GetArgument().IsDefined(), true, "mutate operator()");

      // construct a chain
      util::ShPtr< assemble::Chain> sp_new_chain
      (
        new assemble::Chain( chain->GetSequence(), *mutate_result.GetArgument())
      );

      // construct a model
      assemble::ProteinModel new_model( sp_new_chain);

      // write out the model
      Proteins::WriteModelToPDB( new_model, AddExampleOutputPathToFilename( mutate, "mutate_sheet_twist.pdb"));

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( mutate);
      fold::MutateSheetTwist mutate_read;
      ReadBCLObject( mutate_read);
      BCL_ExampleIndirectCheck( mutate.GetAngleRange().GetString(), mutate_read.GetAngleRange().GetString(), "I/O");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetTwist

  const ExampleClass::EnumType ExampleFoldMutateSheetTwist::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetTwist())
  );

} // namespace bcl
