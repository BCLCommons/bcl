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
#include "fold/bcl_fold_mutate_sse_bend_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sse_bend_random.cpp
  //!
  //! @author karakam
  //! @date Jan 24, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSSEBendRandom :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSSEBendRandom *Clone() const
    {
      return new ExampleFoldMutateSSEBendRandom( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // get the SSEs
      util::ShPtr< assemble::SSE>
        sp_helix( Proteins::GetSSE( pdb_filename, 'A', 23, 34, biol::GetAAClasses().e_AABackBone)),
        sp_strand( Proteins::GetSSE( pdb_filename, 'A', 1, 7, biol::GetAAClasses().e_AABackBone));

      // create the random phi psi change pair
      math::Range< double> range_default( math::Angle::Radian( -20.0), math::Angle::Radian( 20.0));
      math::Range< double> range( math::Angle::Radian( -30.0), math::Angle::Radian( 30.0));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::MutateSSEBendRandom default_mutate;
      BCL_ExampleCheck
      (
        default_mutate.GetPhiChangeRange() == range_default &&
        default_mutate.GetPsiChangeRange() == range_default,
        true
      );

      // construct mutates from locator and phi psi ranges
      fold::MutateSSEBendRandom mutate( range, range, biol::AASequenceFlexibility::e_Bidirectional, "bend");

    /////////////////
    // data access //
    /////////////////

      // test GetPhiChangeRange()
      BCL_ExampleCheck( mutate.GetPhiChangeRange(), range);

      // test GetPsiChangeRange()
      BCL_ExampleCheck( mutate.GetPsiChangeRange(), range);

      // test GetBendingDirection()
      BCL_ExampleCheck( mutate.GetBendingDirection(), biol::AASequenceFlexibility::e_Bidirectional);

      // test GetScheme()
      BCL_ExampleCheck( mutate.GetScheme(), "bend");

    ///////////////
    // operators //
    ///////////////

      // test operator() helix
      math::MutateResult< assemble::SSE> result_a( mutate( *sp_helix));
      BCL_ExampleCheck( result_a.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_a.GetArgument(), *sp_helix),
        1,
        "there should be only 1 phi/psi change"
      );
      Proteins::WriteSSEToPDB
      (
        result_a.GetArgument(), AddExampleOutputPathToFilename( mutate, "1ubi_bent_model_a.pdb")
      );

      // test operator() with strand
      math::MutateResult< assemble::SSE> result_b( mutate( *sp_strand));
      BCL_ExampleCheck( result_b.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_b.GetArgument(), *sp_strand),
        1,
        "there should be only 1 phi/psi change"
      );
      Proteins::WriteSSEToPDB
      (
        result_b.GetArgument(), AddExampleOutputPathToFilename( mutate, "1ubi_bent_model_b.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // check read write
      WriteBCLObject( mutate);
      fold::MutateSSEBendRandom mutate_read;
      ReadBCLObject( mutate_read);
      BCL_ExampleIndirectCheck
      (
        mutate_read.GetScheme() == mutate.GetScheme() &&
        mutate_read.GetPhiChangeRange().GetString() == mutate.GetPhiChangeRange().GetString() &&
        mutate_read.GetPsiChangeRange().GetString() == mutate.GetPsiChangeRange().GetString(),
        true,
        "Read() function failed"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSSEBendRandom

  const ExampleClass::EnumType ExampleFoldMutateSSEBendRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSSEBendRandom())
  );

} // namespace bcl
