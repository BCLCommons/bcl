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
#include "fold/bcl_fold_mutate_protein_model_sse_resize.h"

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
  //! @example example_fold_mutate_protein_model_sse_resize.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEResize :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEResize *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSEResize( *this);
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
      // minimum sse sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;

      // initialize pdb filename to be read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // read the model
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // initialize locator
      const assemble::LocatorSSE locator_helix_23_34( 'A', 23, 34);

      // get the original helix
      const util::SiPtr< const assemble::SSE>
        sp_original_helix_23_34( locator_helix_23_34.Locate( native_model));

      // extend shrink probabilities
      const double prob_shrink( 0.0);
      const double prob_extend( 1.0);

      // length increment
      const math::Range< size_t> increment_shrink_1( 1, 1);
      const math::Range< size_t> increment_extend_2( 2, 2);
      const math::Range< size_t> increment_extend_10( 10, 10);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      fold::MutateProteinModelSSEResize mutate_default;

      // construct one that shrinks
      BCL_MessageStd( "construct shrink_1");
      fold::MutateProteinModelSSEResize
        mutate_shrink_1
        (
          util::CloneToShPtr( locator_helix_23_34),
          prob_shrink,
          increment_shrink_1,
          biol::AASequenceFlexibility::e_Bidirectional,
          false,
          min_sse_sizes
        );

      // construct one that extends
      BCL_MessageStd( "construct extend_2");
      fold::MutateProteinModelSSEResize
        mutate_extend_2
        (
          util::CloneToShPtr( locator_helix_23_34),
          prob_extend,
          increment_extend_2,
          biol::AASequenceFlexibility::e_Bidirectional,
          false,
          min_sse_sizes
        );

      // construct one that extend too much
      BCL_MessageStd( "construct extend_10");
      fold::MutateProteinModelSSEResize
        mutate_extend_10
        (
          util::CloneToShPtr( locator_helix_23_34),
          prob_extend,
          increment_extend_10,
          biol::AASequenceFlexibility::e_Bidirectional,
          false,
          min_sse_sizes
        );

      BCL_Debug( mutate_extend_10);

      // construct one that extends on only one side (by 2 residues), but recenters the sse afterwards
      BCL_MessageStd( "construct extend_2_one_side_recenter");
      fold::MutateProteinModelSSEResize mutate_extend_2_one_side_recenter
      (
        util::CloneToShPtr( locator_helix_23_34),
        prob_extend,
        increment_extend_2,
        biol::AASequenceFlexibility::e_CTerminal,
        true,
        min_sse_sizes
      );

      // clone
      BCL_MessageStd( "test clone");
      util::ShPtr< fold::MutateProteinModelSSEResize> sp_mutate_cloned( mutate_shrink_1.Clone());
      BCL_Example_Check
      (
        sp_mutate_cloned->GetExtendProbability() == mutate_shrink_1.GetExtendProbability() &&
        sp_mutate_cloned->GetLengthChangeRange() == mutate_shrink_1.GetLengthChangeRange() &&
        sp_mutate_cloned->GetSide() == mutate_shrink_1.GetSide() &&
        sp_mutate_cloned->GetMinSSESizes() == mutate_shrink_1.GetMinSSESizes(),
        "the cloned mutate differs from mutate_shrink_1" +
        util::Format()( *sp_mutate_cloned) + "\nvs\n" + util::Format()( mutate_shrink_1)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_MessageStd( "test GetStatisClassName");
      const std::string correct_static_class_name( "bcl::fold::MutateProteinModelSSEResize");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelSSEResize>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< fold::MutateProteinModelSSEResize>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelSSEResize>() == mutate_shrink_1.GetClassIdentifier(),
        "GetClassIdentifier gives " + mutate_shrink_1.GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetExtendShrinkProbability()
      BCL_MessageStd( "test GetExtendShrinkProbability");
      BCL_Example_Check
      (
        mutate_shrink_1.GetExtendProbability() == prob_shrink,
        "GetExtendShrinkProbability gives " + util::Format()( mutate_shrink_1.GetExtendProbability()) +
        " instead of " + util::Format()( prob_shrink)
      );

      // check GetLengthChangeRange()
      BCL_MessageStd( "test GetLengthChangeRange");
      BCL_Example_Check
      (
        mutate_shrink_1.GetLengthChangeRange() == increment_shrink_1,
        "GetLengthChangeRange gives " + util::Format()( mutate_shrink_1.GetLengthChangeRange()) +
        " instead of " + util::Format()( increment_shrink_1)
      );

      // check GetSide()
      BCL_MessageStd( "test GetChangeBothEnds");
      BCL_Example_Check
      (
        mutate_shrink_1.GetSide() == biol::AASequenceFlexibility::e_Bidirectional,
        "GetSide gives " + util::Format()( mutate_shrink_1.GetSide()) +
        " instead of bidirectional"
      );

      // check GetMinSSESizes()
      BCL_MessageStd( "test GetMinSSESizes");
      BCL_Example_Check
      (
        mutate_shrink_1.GetMinSSESizes() == min_sse_sizes,
        "GetMinSSESizes gives " + util::Format()( mutate_shrink_1.GetMinSSESizes()) +
        " instead of " + util::Format()( min_sse_sizes)
      );

    ///////////////
    // operators //
    ///////////////

      // initialize locators to locate SSE after they have been mutated
      assemble::LocatorSSE locator_helix_shrink_1( 'A', 24, 33);
      assemble::LocatorSSE locator_helix_extend_2( 'A', 21, 36);
      assemble::LocatorSSE locator_helix_extend_2_left( 'A', 21, 34);
      assemble::LocatorSSE locator_helix_extend_2_right( 'A', 23, 36);
      assemble::LocatorSSE locator_helix_extend_10( 'A', 13, 44);

      // start testing the operator
      BCL_MessageStd( "testing operator()");

      // first shrink the helix
      BCL_MessageStd( "test shrink_1 mutate");
      math::MutateResult< assemble::ProteinModel> result_shrink_1( mutate_shrink_1( native_model));

      // check that a non-empty model is returned
      BCL_Example_Check
      (
        result_shrink_1.GetArgument().IsDefined(),
        "The mutate failed when it should have succeeded"
      );

      // check the model returned for the SSE size
      util::SiPtr< const assemble::SSE>
        sp_helix_shrink_1( locator_helix_shrink_1.Locate( *result_shrink_1.GetArgument()));
      BCL_Example_Check
      (
        sp_helix_shrink_1.IsDefined(),
        "The shrunk helix 24-33 cannot be located in the returned model from mutate"
      );

      // extend the helix
      BCL_MessageStd( "test extend_2 mutate");
      math::MutateResult< assemble::ProteinModel> result_extend_2( mutate_extend_2( native_model));

      // check that a non-empty model is returned
      BCL_Example_Check
      (
        result_extend_2.GetArgument().IsDefined(),
        "The mutate failed when it should have succeeded"
      );

      // check the model returned for the SSE size
      util::SiPtr< const assemble::SSE>
        sp_helix_extend_2( locator_helix_extend_2.Locate( *result_extend_2.GetArgument()));
      BCL_Example_Check
      (
        sp_helix_extend_2.IsDefined(),
        "The extended helix 21-36 cannot be located in the returned model from mutate"
      );

      // extend the helix too much
      BCL_MessageStd( "test extend_10");
      math::MutateResult< assemble::ProteinModel> result_extend_10( mutate_extend_10( native_model));

      // check that a non-empty model is returned
      BCL_Example_Check
      (
        !result_extend_10.GetArgument().IsDefined(),
        "The mutate should have failed but instead it succeeded"
      );

      // extend the helix by two residues, but only on one side, recenter after resizing
      BCL_MessageStd( "test extend_2_one_side_recenter mutate");
      math::MutateResult< assemble::ProteinModel> result_extend_2_one_side_recenter
      (
        mutate_extend_2_one_side_recenter( native_model)
      );

      // check that a non-empty model is returned
      BCL_Example_Check
      (
        result_extend_2_one_side_recenter.GetArgument().IsDefined(),
        "The mutate failed when it should have succeeded"
      );

      // check the model returned for the SSE size
      util::SiPtr< const assemble::SSE>
        sp_helix_extend_2_left( locator_helix_extend_2_left.Locate( *result_extend_2_one_side_recenter.GetArgument()));
      util::SiPtr< const assemble::SSE>
        sp_helix_extend_2_right( locator_helix_extend_2_right.Locate( *result_extend_2_one_side_recenter.GetArgument()));
      BCL_Example_Check
      (
        sp_helix_extend_2_left.IsDefined() || sp_helix_extend_2_right.IsDefined(),
        "Either helix 21-34 or 23-36 should have been located in the returned model from mutate"
      );

      // check that center of extended and recentered helix agrees with original center
      BCL_Example_Check
      (
        (
          sp_helix_extend_2_left.IsDefined()
          && sp_helix_extend_2_left->GetCenter() == sp_original_helix_23_34->GetCenter()
        )
        ||
        (
          sp_helix_extend_2_right.IsDefined()
          && sp_helix_extend_2_right->GetCenter() == sp_original_helix_23_34->GetCenter()
        ),
        "The center of the helix should have been reset after the extension!"
      );

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

  }; //end ExampleFoldMutateProteinModelSSEResize

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEResize::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEResize())
  );

} // namespace bcl
