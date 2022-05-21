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
#include "fold/bcl_fold_mutate_protein_model_sse_move.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_locator_criteria_wrapper.h"
#include "fold/bcl_fold_placement_sse_next_to_sse.h"
#include "fold/bcl_fold_placement_sse_short_loop.h"
#include "fold/bcl_fold_placement_strand_next_to_sheet.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_fold_mutate_protein_model_sse_move.cpp
  //!
  //! @author karakam
  //! @date Aug 5, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEMove :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEMove *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSEMove( *this);
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

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1J27_bad_model.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // construct locator
      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
        locator_helix_63_81( assemble::LocatorSSE( 'A', 63, 81));

      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
        locator_strand_2_13( assemble::LocatorSSE( 'A', 2, 13));

      // construct random locator
      const assemble::LocatorSSERandom locator_sse_random;
      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
        locator_sse_random_crit_wrap( locator_sse_random);

      // placement of an SSE next to another SSE
      const fold::PlacementSSENextToSSE placement_next( locator_sse_random_crit_wrap);

      // placement of an SSE connected with a short loop
      const fold::PlacementSSEShortLoop placement_short( 7);

      // placement for a strand next to a sheet
      const fold::PlacementStrandNextToSheet placement_sheet( 0.5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct helix next to sse
      fold::MutateProteinModelSSEMove mutate_helix_next( locator_helix_63_81, placement_next, "mutate_helix_next");
      fold::MutateProteinModelSSEMove mutate_strand_next( locator_strand_2_13, placement_next);
      fold::MutateProteinModelSSEMove mutate_helix_short( locator_helix_63_81, placement_short);
      fold::MutateProteinModelSSEMove mutate_strand_short( locator_strand_2_13, placement_short);
      fold::MutateProteinModelSSEMove mutate_strand_sheet( locator_strand_2_13, placement_sheet);

    /////////////////
    // data access //
    /////////////////

      // check GetScheme() function
      BCL_ExampleCheck( mutate_helix_next.GetScheme(), "mutate_helix_next");

    ///////////////
    // operators //
    ///////////////

      {
        math::MutateResult< assemble::ProteinModel> result( mutate_helix_next( model));
        BCL_ExampleCheck( result.GetArgument().IsDefined(), true);
        Proteins::WriteModelToPDB
        (
          *result.GetArgument(), AddExampleOutputPathToFilename( mutate_helix_next, "mutate_sse_move_helix_next.pdb")
        );
      }

      {
        math::MutateResult< assemble::ProteinModel> result( mutate_strand_next( model));
        BCL_ExampleCheck( result.GetArgument().IsDefined(), true);
        Proteins::WriteModelToPDB
        (
          *result.GetArgument(), AddExampleOutputPathToFilename( mutate_strand_next, "mutate_sse_move_strand_next.pdb")
        );
      }

      {
        math::MutateResult< assemble::ProteinModel> result( mutate_helix_short( model));
        BCL_ExampleCheck( result.GetArgument().IsDefined(), true);
        Proteins::WriteModelToPDB
        (
          *result.GetArgument(), AddExampleOutputPathToFilename( mutate_helix_short, "mutate_sse_move_helix_short.pdb")
        );
      }

      {
        math::MutateResult< assemble::ProteinModel> result( mutate_strand_short( model));
        BCL_ExampleCheck( result.GetArgument().IsDefined(), true);
        Proteins::WriteModelToPDB
        (
          *result.GetArgument(), AddExampleOutputPathToFilename( mutate_strand_short, "mutate_sse_move_strand_short.pdb")
        );
      }

      {
        math::MutateResult< assemble::ProteinModel> result( mutate_strand_sheet( model));
        BCL_ExampleCheck( result.GetArgument().IsDefined(), true);
        Proteins::WriteModelToPDB
        (
          *result.GetArgument(), AddExampleOutputPathToFilename( mutate_helix_next, "mutate_sse_move_strand_sheet.pdb")
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEMove

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEMove::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEMove())
  );
  
} // namespace bcl
