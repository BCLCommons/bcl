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
#include "assemble/bcl_assemble_sse_compare.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_compare.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSECompare :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSECompare *Clone() const
    {
      return new ExampleAssembleSSECompare( *this);
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
      // initialize read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //instantiate sequences
      BCL_MessageStd( "building sequences from pdb chains");
      // initialize sequence
      biol::AASequence seq( *( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone).GetSequences()( 0)));

      BCL_MessageStd( "sequences are built");
      BCL_MessageStd( "Creating sses");

      util::ShPtr< assemble::SSE> sse_15_24( new assemble::SSE( seq.SubSequence( 14, 10),biol::GetSSTypes().HELIX));
      util::ShPtr< assemble::SSE> sse_17_29( new assemble::SSE( seq.SubSequence( 16, 13),biol::GetSSTypes().HELIX));
      util::ShPtr< assemble::SSE> sse_30_35( new assemble::SSE( seq.SubSequence( 29, 6),biol::GetSSTypes().HELIX));
      util::ShPtr< assemble::SSE> sse_40_45( new assemble::SSE( seq.SubSequence( 39, 6),biol::GetSSTypes().HELIX));
      util::ShPtr< assemble::SSE> sse_40_48( new assemble::SSE( seq.SubSequence( 39, 9),biol::GetSSTypes().HELIX));

    /////////////////
    // SSELessThan //
    /////////////////

      // test SSELessThan
      BCL_MessageStd( "testing SSELessThan, inserting sses");
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan> set_sse_less_than;
      set_sse_less_than.Insert( sse_40_45);
      set_sse_less_than.Insert( sse_15_24);
      set_sse_less_than.Insert( sse_30_35);
      set_sse_less_than.Insert( sse_17_29);
      set_sse_less_than.Insert( sse_40_48);

      // print out the inserted sse
      BCL_MessageStd( "The inserted sses are:");
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::const_iterator
          sse_itr( set_sse_less_than.Begin()), sse_itr_end( set_sse_less_than.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        BCL_MessageStd
        (
          util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " to " +
          util::Format()( ( *sse_itr)->GetLastAA()->GetSeqID())
        );
      }

      // check that all the sses have been inserted correctly
      BCL_ExampleIndirectCheck( set_sse_less_than.GetSize(), 5, "Set Insert using sse_less_than");

      BCL_MessageStd( "all sses have been inserted");

    //////////////////////////
    // SSELessThanNoOverlap //
    //////////////////////////

      // test SSELessThan
      BCL_MessageStd( "testing SSELessThanNoOverlap, inserting sses");
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> set_sse_less_than_no_overlap;
      set_sse_less_than_no_overlap.Insert( sse_40_45);
      set_sse_less_than_no_overlap.Insert( sse_15_24);
      set_sse_less_than_no_overlap.Insert( sse_30_35);
      set_sse_less_than_no_overlap.Insert( sse_17_29);
      set_sse_less_than_no_overlap.Insert( sse_40_48);

      // print out the inserted sse
      BCL_MessageStd( "The inserted sses are:");
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( set_sse_less_than_no_overlap.Begin()), sse_itr_end( set_sse_less_than_no_overlap.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        BCL_MessageStd
        (
          util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " to " +
          util::Format()( ( *sse_itr)->GetLastAA()->GetSeqID())
        );
      }

      // check that all the sses have been inserted correctly
      BCL_ExampleIndirectCheck( set_sse_less_than_no_overlap.GetSize(), 3, "Set Insert using sse_less_than_no_overlap");

      BCL_MessageStd( "all sses have been inserted");

    ///////////////////////
    // SSELessThanBySize //
    ///////////////////////

      // test SSELessThanBySize
      BCL_MessageStd( "testing SSELessThanBySize, inserting sses");
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanBySize> set_sse_less_than_by_size;
      set_sse_less_than_by_size.Insert( sse_40_45);
      set_sse_less_than_by_size.Insert( sse_15_24);
      set_sse_less_than_by_size.Insert( sse_30_35);
      set_sse_less_than_by_size.Insert( sse_17_29);
      set_sse_less_than_by_size.Insert( sse_40_48);

      // print out the inserted sse
      BCL_MessageStd( "The inserted sses are:");
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanBySize>::const_iterator
          sse_itr( set_sse_less_than_by_size.Begin()), sse_itr_end( set_sse_less_than_by_size.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        BCL_MessageStd
        (
          util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " to " +
          util::Format()( ( *sse_itr)->GetLastAA()->GetSeqID())
        );
      }

      // check that all the sses have been inserted correctly
      BCL_ExampleIndirectCheck( set_sse_less_than_by_size.GetSize(), 5, "Set Insert using SSELessThanBySize");

      BCL_MessageStd( "all sses have been inserted");

    ////////////////
    // SSECompare //
    ////////////////

      BCL_MessageStd( "testing SSECompare");

      // search for sse_40_48
      BCL_MessageStd( "searching for existing sse");
      const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::const_iterator itr_compare
      (
        std::find_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompare( *sse_40_48))
      );

      // assert it is found and found correctly indeed
      BCL_ExampleIndirectAssert( itr_compare != set_sse_less_than.End(), true, "Find with SSECompare");
      BCL_ExampleIndirectCheck( **itr_compare, *sse_40_48, "Find with SSECompare");

      // search for sse_40_48
      BCL_MessageStd( "counting matches to existing sse");

      // assert that no matches are found
      BCL_ExampleCheck
      (
        std::count_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompare( *sse_40_48)),
        1
      );

      // create an overlapping sse
      util::ShPtr< assemble::SSE> sse_14_29( new assemble::SSE( seq.SubSequence( 13, 16), biol::GetSSTypes().HELIX));

      // search for sse_14_29
      BCL_MessageStd( "searching for not-inserted sse");

      // assert that it is not found
      BCL_ExampleCheck
      (
        std::find_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompare( *sse_14_29))
        == set_sse_less_than.End(),
        true
      );

      BCL_MessageStd( "counting machtes to not-inserted sse");

      // assert that no matches are found
      BCL_ExampleCheck
      (
        std::count_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompare( *sse_14_29)),
        0
      );

    ///////////////////////
    // SSECompareOverlap //
    ///////////////////////

      BCL_MessageStd( "testing SSECompareOverlap");

      // search for sse_40_45
      BCL_MessageStd( "searching for existing sse, should find 40_45 first");

      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::iterator
      itr_compare_overlap
      (
        std::find_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompareOverlap( *sse_40_45))
      );

      // assert that match is found and indeed found correctly
      BCL_ExampleIndirectAssert( itr_compare_overlap != set_sse_less_than.End(), true, "Find with SSECompareOverlap");
      BCL_ExampleIndirectCheck( **itr_compare_overlap, *sse_40_45, "Find with SSECompareOverlap");

      // search for sse_40_45
      BCL_MessageStd( "continuing search from the previous match, should find 40_48 now");

      // increment itr
      ++itr_compare_overlap;

      // continue searching from the last position
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::const_iterator
      itr_compare_overlap_2
      (
        std::find_if( itr_compare_overlap, set_sse_less_than.End(), assemble::SSECompareOverlap( *sse_40_45))
      );

      // assert that match is found and indeed found correctly
      BCL_ExampleIndirectAssert( itr_compare_overlap_2 != set_sse_less_than.End(), true, "Find with SSECompareOverlap");
      BCL_ExampleIndirectCheck( **itr_compare_overlap_2, *sse_40_48, "Find with SSECompareOverlap");

      BCL_MessageStd( "counting overlapping sses with 40_45");

      // count equal ones
      BCL_ExampleCheck
      (
        std::count_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompareOverlap( *sse_40_45)),
        2
      );

      // search for sse_40_45
      BCL_MessageStd( "searching for non-inserted but overlapping sse, ");
      const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
      itr_compare_overlap_3
      (
        std::find_if
        (
          set_sse_less_than_no_overlap.Begin(),
          set_sse_less_than_no_overlap.End(),
          assemble::SSECompareOverlap( *sse_40_48)
        )
      );

      // assert that match is found and indeed found correctly
      BCL_ExampleIndirectAssert
      (
        itr_compare_overlap_3 != set_sse_less_than_no_overlap.End(),
        true,
        "find with SSECompareOverlap"
      );
      BCL_ExampleIndirectCheck( **itr_compare_overlap_3, *sse_40_45, "Find with SSECompareOverlap");

      // search for sse_40_48
      BCL_MessageStd( "counting matches to non-inserted but overlapping sse");

      // count equal ones
      BCL_ExampleCheck
      (
        std::count_if
        (
          set_sse_less_than_no_overlap.Begin(),
          set_sse_less_than_no_overlap.End(),
          assemble::SSECompareOverlap( *sse_40_48)
        ),
        1
      );

      // create a non-overlapping sse
      util::ShPtr< assemble::SSE> sse_60_70( new assemble::SSE( seq.SubSequence( 59, 11), biol::GetSSTypes().HELIX));

      // search for sse_14_29
      BCL_MessageStd( "searching for not-inserted non-overlapping sse");

      // assert that it is not found
      BCL_ExampleCheck
      (
        std::find_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompareOverlap( *sse_60_70))
        == set_sse_less_than.End(),
        true
      );

      BCL_MessageStd( "counting matches to not-inserted non-overlapping sse");

      // assert that no matches are found
      BCL_ExampleCheck
      (
        std::count_if( set_sse_less_than.Begin(), set_sse_less_than.End(), assemble::SSECompareOverlap( *sse_60_70)),
        0
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSECompare

  const ExampleClass::EnumType ExampleAssembleSSECompare::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSECompare())
  );

} // namespace bcl
