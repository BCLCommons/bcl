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
#include "fold/bcl_fold_mutate_sheet_sort.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sheet_sort.cpp
  //!
  //! @author karakam
  //! @date Mar 18, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetSort :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetSort *Clone() const
    {
      return new ExampleFoldMutateSheetSort( *this);
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

      // pdb filename for 1JL1
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1JL1A.pdb"));
      storage::Map< biol::SSType, size_t> sse_sizes;
      sse_sizes[ biol::GetSSTypes().HELIX] = 7;
      sse_sizes[ biol::GetSSTypes().STRAND] = 4;
      sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // read the model and create reference on chain
      const assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_sizes));
      const assemble::Chain &chain( *model.GetChains().FirstElement());

      // collect sheet
      const util::ShPtrVector< assemble::Domain> sp_sheets( assemble::CollectorSheet().Collect( model));
      const util::ShPtr< assemble::Domain> sp_sheet( sp_sheets.FirstElement());

      const util::SiPtr< const assemble::SSE> sp_strand_5_13( assemble::LocatorSSE( 'A', 5, 13).Locate( model));
      const util::SiPtr< const assemble::SSE> sp_strand_18_28( assemble::LocatorSSE( 'A', 18, 28).Locate( model));
      const util::SiPtr< const assemble::SSE> sp_strand_31_42( assemble::LocatorSSE( 'A', 31, 42).Locate( model));
      const util::SiPtr< const assemble::SSE> sp_strand_64_69( assemble::LocatorSSE( 'A', 64, 69).Locate( model));
      const util::SiPtr< const assemble::SSE> sp_strand_115_120( assemble::LocatorSSE( 'A', 115, 120).Locate( model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "testing default constructor");
      const fold::MutateSheetSort mutate_default;
      BCL_ExampleCheck( mutate_default.GetReverseSortProbability(), 0.0);

      // construct from a range
      BCL_MessageStd( "testing constructor a");
      const fold::MutateSheetSort mutate_a( 0.0);
      BCL_ExampleCheck( mutate_a.GetReverseSortProbability(), 0.0);

      // construct from range and consecutive
      BCL_MessageStd( "testing constructor b");
      const fold::MutateSheetSort mutate_b( 1.0);

    /////////////////
    // data access //
    /////////////////

      // test GetReverseSortProbability()
      BCL_MessageStd( "GetReverseSortProbability()");
      BCL_ExampleCheck( mutate_b.GetReverseSortProbability(), 1.0);

    ////////////////
    // operations //
    ////////////////

      // testing operator() with mutate_a
      const math::MutateResult< assemble::Domain> result_a( mutate_a( *sp_sheet));
      // make sure the result is defined
      BCL_ExampleCheck( result_a.GetArgument().IsDefined(), true);
      const util::ShPtr< assemble::Domain> sp_sheet_a( result_a.GetArgument());
      // make sure the size is 5
      BCL_ExampleCheck( sp_sheet_a->GetNumberSSEs(), 5);
      // make sure the order of strands are correct
      BCL_ExampleCheck( sp_sheet_a->GetTopology().IsDefined(), true);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 0), sp_strand_5_13);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 1), sp_strand_18_28);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 2), sp_strand_31_42);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 3), sp_strand_64_69);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 4), sp_strand_115_120);
      // construct a chain from the sheet and write it to a pdb file
      const assemble::Chain chain_a( chain.GetSequence(), *sp_sheet_a);
      Proteins::WriteChainToPDB( chain_a, AddExampleOutputPathToFilename( mutate_a, "mutate_sheet_sort_a.pdb"));

      // testing operator() with mutate_b
      BCL_MessageStd( "testing operator() with mutate_b");
      const math::MutateResult< assemble::Domain> result_b( mutate_b( *sp_sheet));
      // make sure the result is defined
      BCL_ExampleCheck( result_b.GetArgument().IsDefined(), true);
      const util::ShPtr< assemble::Domain> sp_sheet_b( result_b.GetArgument());
      // make sure the size is 3
      BCL_ExampleCheck( sp_sheet_b->GetNumberSSEs(), 5);
      // make sure the order of strands are correct
      BCL_ExampleCheck( sp_sheet_b->GetTopology().IsDefined(), true);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 0), sp_strand_115_120);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 1), sp_strand_64_69);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 2), sp_strand_31_42);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 3), sp_strand_18_28);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 4), sp_strand_5_13);
      // construct a chain from the sheet and write it to a pdb file
      const assemble::Chain chain_b( chain.GetSequence(), *sp_sheet_b);
      Proteins::WriteChainToPDB( chain_b, AddExampleOutputPathToFilename( mutate_b, "mutate_sheet_sort_b.pdb"));

    //////////////////////
    // input and output //
    //////////////////////

      // test read/write
      BCL_MessageStd( "testing Read/Write");
      WriteBCLObject( mutate_b);
      fold::MutateSheetSort mutate_read;
      ReadBCLObject( mutate_read);
      BCL_ExampleCheck( mutate_read.GetReverseSortProbability(), mutate_b.GetReverseSortProbability());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetSort

  const ExampleClass::EnumType ExampleFoldMutateSheetSort::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetSort())
  );

} // namespace bcl
