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
#include "fold/bcl_fold_mutate_sheet_order.h"

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
  //! @example example_fold_mutate_sheet_order.cpp
  //!
  //! @author karakam
  //! @date Mar 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSheetOrder :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSheetOrder *Clone() const
    {
      return new ExampleFoldMutateSheetOrder( *this);
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
      const fold::MutateSheetOrder mutate_default;
      BCL_ExampleCheck( mutate_default.GetParallelProbability(), 0.0);

      // construct from a range
      BCL_MessageStd( "testing constructor a");
      const fold::MutateSheetOrder mutate_a( 0.0);
      BCL_ExampleCheck( mutate_a.GetParallelProbability(), 0.0);

      // construct from range and consecutive
      BCL_MessageStd( "testing constructor b");
      const fold::MutateSheetOrder mutate_b( 1.0);

    /////////////////
    // data access //
    /////////////////

      // test GetReverseSortProbability()
      BCL_MessageStd( "GetParallelProbability()");
      BCL_ExampleCheck( mutate_b.GetParallelProbability(), 1.0);

    ////////////////
    // operations //
    ////////////////

      // testing operator() with mutate_a
      const math::MutateResult< assemble::Domain> result_a( mutate_a( *sp_sheet));

      // make sure the result is defined
      BCL_ExampleCheck( result_a.GetArgument().IsDefined(), true);
      const util::ShPtr< assemble::Domain> sp_sheet_a( result_a.GetArgument());

      // make sure it has 5 SSEs
      const util::SiPtrVector< const assemble::SSEGeometryInterface> sses_a( sp_sheet_a->GetSSEs());
      BCL_ExampleCheck( sses_a.GetSize(), 5);

      // calculate packings
      const assemble::SSEGeometryPacking packing_a_1( *sses_a( 2), *sses_a( 1));
      const assemble::SSEGeometryPacking packing_a_2( *sses_a( 1), *sses_a( 0));
      const assemble::SSEGeometryPacking packing_a_3( *sses_a( 0), *sses_a( 3));
      const assemble::SSEGeometryPacking packing_a_4( *sses_a( 3), *sses_a( 4));

      // make sure all packings are anti-parallel
      const assemble::SSEGeometryPacking::Orientation antiparallel( assemble::SSEGeometryPacking::e_AntiParallel);
      BCL_ExampleCheck( packing_a_1.GetOrientation(), antiparallel);
      BCL_ExampleCheck( packing_a_2.GetOrientation(), antiparallel);
      BCL_ExampleCheck( packing_a_3.GetOrientation(), antiparallel);
      BCL_ExampleCheck( packing_a_4.GetOrientation(), antiparallel);

      // construct a chain from the sheet and write it to a pdb file
      const assemble::Chain chain_a( chain.GetSequence(), *sp_sheet_a);
      Proteins::WriteChainToPDB( chain_a, AddExampleOutputPathToFilename( mutate_a, "mutate_sheet_order_a.pdb"));

      // testing operator() with mutate_b
      BCL_MessageStd( "testing operator() with mutate_b");
      const math::MutateResult< assemble::Domain> result_b( mutate_b( *sp_sheet));

      // make sure the result is defined
      BCL_ExampleCheck( result_b.GetArgument().IsDefined(), true);
      const util::ShPtr< assemble::Domain> sp_sheet_b( result_b.GetArgument());

      // make sure the size is 5
      const util::SiPtrVector< const assemble::SSEGeometryInterface> sses_b( sp_sheet_b->GetSSEs());
      BCL_ExampleCheck( sses_b.GetSize(), 5);

      // calculate packings
      const assemble::SSEGeometryPacking packing_b_1( *sses_b( 2), *sses_b( 1));
      const assemble::SSEGeometryPacking packing_b_2( *sses_b( 1), *sses_b( 0));
      const assemble::SSEGeometryPacking packing_b_3( *sses_b( 0), *sses_b( 3));
      const assemble::SSEGeometryPacking packing_b_4( *sses_b( 3), *sses_b( 4));

      // make sure all packings are anti-parallel
      const assemble::SSEGeometryPacking::Orientation parallel( assemble::SSEGeometryPacking::e_Parallel);
      BCL_ExampleCheck( packing_b_1.GetOrientation(), parallel);
      BCL_ExampleCheck( packing_b_2.GetOrientation(), parallel);
      BCL_ExampleCheck( packing_b_3.GetOrientation(), parallel);
      BCL_ExampleCheck( packing_b_4.GetOrientation(), parallel);

      // construct a chain from the sheet and write it to a pdb file
      const assemble::Chain chain_b( chain.GetSequence(), *sp_sheet_b);
      Proteins::WriteChainToPDB( chain_b, AddExampleOutputPathToFilename( mutate_b, "mutate_sheet_order_b.pdb"));

    //////////////////////
    // input and output //
    //////////////////////

      // test read/write
      BCL_MessageStd( "testing Read/Write");
      WriteBCLObject( mutate_b);
      fold::MutateSheetOrder mutate_read;
      ReadBCLObject( mutate_read);
      BCL_ExampleCheck( mutate_read.GetParallelProbability(), mutate_b.GetParallelProbability());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSheetOrder

  const ExampleClass::EnumType ExampleFoldMutateSheetOrder::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSheetOrder())
  );

} // namespace bcl
