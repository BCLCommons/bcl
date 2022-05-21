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
#include "fold/bcl_fold_mutate_sse_bend_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sse_bend_ramachandran.cpp
  //!
  //! @author karakam
  //! @date Jan 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSSEBendRamachandran :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSSEBendRamachandran *Clone() const
    {
      return new ExampleFoldMutateSSEBendRamachandran( *this);
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
      math::Range< size_t> one_change( 1, 1);
      math::Range< size_t> three_changes( 3, 3);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::MutateSSEBendRamachandran default_mutate;
      BCL_ExampleCheck
      (
        default_mutate.GetScheme() == GetStaticClassName< fold::MutateSSEBendRamachandran>() &&
        default_mutate.GetNrResiduesChangeRange() == one_change,
        true
      );

      // construct mutates from locator and phi psi ranges
      fold::MutateSSEBendRamachandran mutate_one
      (
        one_change, biol::AASequenceFlexibility::e_Bidirectional, "bend_one"
      );
      fold::MutateSSEBendRamachandran mutate_three
      (
        three_changes, biol::AASequenceFlexibility::e_Bidirectional, "bend_three"
      );

    /////////////////
    // data access //
    /////////////////

      // test GetNrResiduesChangeRange()
      BCL_ExampleCheck( mutate_one.GetNrResiduesChangeRange(), one_change);
      BCL_ExampleCheck( mutate_three.GetNrResiduesChangeRange(), three_changes);

      // test GetBendingDirection()
      BCL_ExampleCheck( mutate_one.GetBendingDirection(), biol::AASequenceFlexibility::e_Bidirectional);
      BCL_ExampleCheck( mutate_three.GetBendingDirection(), biol::AASequenceFlexibility::e_Bidirectional);

      // test GetScheme()
      BCL_ExampleCheck( mutate_one.GetScheme(), "bend_one");
      BCL_ExampleCheck( mutate_three.GetScheme(), "bend_three");

    ///////////////
    // operators //
    ///////////////

      // test operator() helix one change and write model to pdb
      BCL_MessageStd( "test operator() with one phi/psi change with helix");
      math::MutateResult< assemble::SSE> result_a( mutate_one( *sp_helix));
      BCL_ExampleCheck( result_a.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_a.GetArgument(), *sp_helix),
        1,
        "helix should have 1 phi/psi change"
      );
      Proteins::WriteSSEToPDB
      (
        result_a.GetArgument(), AddExampleOutputPathToFilename( mutate_one, "1ubi_bent_rama_a.pdb")
      );

      // test operator() helix three change
      BCL_MessageStd( "test operator() with three phi/psi changes with helix");
      math::MutateResult< assemble::SSE> result_b( mutate_three( *sp_helix));
      BCL_ExampleCheck( result_b.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_b.GetArgument(), *sp_helix),
        3,
        "helix should have only 3 phi/psi changes"
      );
      Proteins::WriteSSEToPDB
      (
        result_b.GetArgument(), AddExampleOutputPathToFilename( mutate_one, "1ubi_bent_rama_b.pdb")
      );

      // test operator() helix one change
      BCL_MessageStd( "test operator() with one phi/psi changes with strand");
      math::MutateResult< assemble::SSE> result_c( mutate_one( *sp_strand));
      BCL_ExampleCheck( result_c.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_c.GetArgument(), *sp_strand),
        1,
        "strand should have only 1 phi/psi change"
      );
      Proteins::WriteSSEToPDB
      (
        result_c.GetArgument(), AddExampleOutputPathToFilename( mutate_one, "1ubi_bent_rama_c.pdb")
      );

      // test operator() helix one change
      BCL_MessageStd( "test operator() with three phi/psi changes with strand");
      math::MutateResult< assemble::SSE> result_d( mutate_three( *sp_strand));
      BCL_ExampleCheck( result_d.GetArgument().IsDefined(), true);
      BCL_ExampleIndirectCheck
      (
        biol::AASequenceFlexibility::GetNumberDifferentPhiPsi( *result_d.GetArgument(), *sp_strand),
        3,
        "strand should have only 3 phi/psi changes"
      );
      Proteins::WriteSSEToPDB
      (
        result_d.GetArgument(), AddExampleOutputPathToFilename( mutate_one, "1ubi_bent_rama_d.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSSEBendRamachandran

  const ExampleClass::EnumType ExampleFoldMutateSSEBendRamachandran::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSSEBendRamachandran())
  );

} // namespace bcl
