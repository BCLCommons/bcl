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
#include "assemble/bcl_assemble_collector_sheet.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_sheet.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorSheet :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorSheet *Clone() const
    {
      return new ExampleAssembleCollectorSheet( *this);
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

      // expected number of sheets
      const size_t expected_nr_sheets( 2);
      // locate the SSEs
      const assemble::SSE &sse_1_5( *assemble::LocatorSSE( 'A', 1, 5).Locate( model));
      const assemble::SSE& sse_10_14( *assemble::LocatorSSE( 'A', 10, 14).Locate( model));
      const assemble::SSE& sse_20_26( *assemble::LocatorSSE( 'A', 20, 26).Locate( model));
      const assemble::SSE& sse_34_39( *assemble::LocatorSSE( 'A', 34, 39).Locate( model));
      const assemble::SSE& sse_50_55( *assemble::LocatorSSE( 'A', 50, 55).Locate( model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a collector
      assemble::CollectorSheet collector;

      // test copy constructor
      BCL_MessageStd( "Testing copy constructor");
      assemble::CollectorSheet collector_copy( collector);

      // test clone
      BCL_MessageStd( "Testing clone function");
      util::ShPtr< assemble::CollectorSheet> sp_collector( collector.Clone());

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( collector.GetClassIdentifier(), GetStaticClassName( collector));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // collect sheets
      BCL_MessageStd( "Collecting sheets")
      util::ShPtrVector< assemble::Domain> sheets( collector.Collect( model));

      BCL_MessageStd( "number of collected sheets " + util::Format()( sheets.GetSize()));
      // check the number of sheets found
      BCL_ExampleIndirectAssert( sheets.GetSize(), expected_nr_sheets, "Collector sheet count");

      BCL_MessageStd( "Following sheets were found:")
      size_t sheet_ctr( 0);
      // iterate over sheets
      for
      (
        util::ShPtrVector< assemble::Domain>::const_iterator sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
        sheet_itr != sheet_itr_end; ++sheet_itr
      )
      {
        // increment sheet counter
        ++sheet_ctr;
        // output the sheet
        BCL_MessageStd
        (
          "The sheet # " + util::Format()( sheet_ctr) + "\n" + ( *sheet_itr)->GetTopology()->GetOrderedIdentification()
        );
      }

      // create reference to the first and second sheet
      const util::ShPtr< assemble::Domain> sp_first_sheet( sheets( 0));
      const util::ShPtr< assemble::Domain> sp_second_sheet( sheets( 1));

      // check the SSEs in the first sheet
      BCL_MessageStd( "Checking the SSEs in the first sheet:");
      BCL_ExampleIndirectCheck
      (
        sp_first_sheet->DoesContain( sse_1_5) && sp_first_sheet->DoesContain( sse_10_14),
        true,
        "Collect"
      );

      // check the order vector for the first sheet
      BCL_MessageStd( "Checking the order of SSEs in the first sheet:");
      BCL_ExampleIndirectCheck
      (
        sp_first_sheet->GetTopology()->GetElements(),
        util::SiPtrVector< const assemble::SSEGeometryInterface>::Create( sse_10_14, sse_1_5),
        "Collect"
      );

      // check the SSEs in the second sheet
      BCL_MessageStd( "Checking the SSEs in the second sheet:");
      BCL_ExampleIndirectCheck
      (
        sp_second_sheet->DoesContain( sse_20_26) &&
        sp_second_sheet->DoesContain( sse_34_39) &&
        sp_second_sheet->DoesContain( sse_50_55),
        true,
        "Collect"
      );

      // check the order vector for the second sheet
      BCL_MessageStd( "Checking the order of SSEs in the second sheet:");
      BCL_ExampleIndirectCheck
      (
        sp_second_sheet->GetTopology()->GetElements(),
        util::SiPtrVector< const assemble::SSEGeometryInterface>::Create( sse_34_39, sse_20_26, sse_50_55),
        "The second sheet's order vector should be 34-39, 20-26 and 50-55"
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorSheet

  const ExampleClass::EnumType ExampleAssembleCollectorSheet::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorSheet())
  );

} // namespace bcl
