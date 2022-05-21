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
#include "assemble/bcl_assemble_locator_sub_domain_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_sub_domain_random.cpp
  //!
  //! @author karakam
  //! @date Mar 16, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorSubDomainRandom :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorSubDomainRandom *Clone() const
    {
      return new ExampleAssembleLocatorSubDomainRandom( *this);
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

      // read the model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_sizes));

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
      const assemble::LocatorSubDomainRandom locator_default;
      BCL_ExampleCheck( locator_default.GetSizeRange(), math::Range< size_t>( 0, 0));
      BCL_ExampleCheck( locator_default.GetLocateConsecutive(), false);
      BCL_ExampleCheck( locator_default.GetUseTopologyOrder(), false);

      // construct from a range
      BCL_MessageStd( "testing constructor a");
      const math::Range< size_t> range( 3, 3);
      const assemble::LocatorSubDomainRandom locator_a( range);
      BCL_ExampleCheck( locator_a.GetSizeRange(), range);
      BCL_ExampleCheck( locator_a.GetLocateConsecutive(), false);
      BCL_ExampleCheck( locator_a.GetUseTopologyOrder(), false);

      // construct from range and consecutive
      BCL_MessageStd( "testing constructor b");
      const assemble::LocatorSubDomainRandom locator_b( range, true, false);
      BCL_ExampleCheck( locator_b.GetSizeRange(), range);
      BCL_ExampleCheck( locator_b.GetLocateConsecutive(), true);
      BCL_ExampleCheck( locator_b.GetUseTopologyOrder(), false);

      // construct from range and consecutive and use topology order
      BCL_MessageStd( "testing constructor c");
      const assemble::LocatorSubDomainRandom locator_c( range, true, true);

    /////////////////
    // data access //
    /////////////////

      // test GetSizeRange()
      BCL_MessageStd( "testing GetSizeRange()");
      BCL_ExampleCheck( locator_c.GetSizeRange(), range);

      // test GetLocateConsecutive()
      BCL_MessageStd( "testing GetLocateConsecutive()");
      BCL_ExampleCheck( locator_c.GetLocateConsecutive(), true);

      // test GetUseTopologyOrder()
      BCL_MessageStd( "testing GetUseTopologyOrder()");
      BCL_ExampleCheck( locator_c.GetUseTopologyOrder(), true);

    ////////////////
    // operations //
    ////////////////

      // testing locate with locator_a
      BCL_MessageStd( "testing Locate() with locator_a");
      util::ShPtr< assemble::Domain> sp_sheet_a( locator_a.Locate( sp_sheet));
      // make sure it is defined
      BCL_ExampleCheck( sp_sheet_a.IsDefined(), true);
      // make sure the size is 3
      BCL_ExampleCheck( sp_sheet_a->GetNumberSSEs(), 3);
      // make sure the order of strands are correct
      BCL_ExampleCheck( sp_sheet_a->GetTopology().IsDefined(), true);

      // handle 32 and 64 bit versions separately
      util::SiPtr< const assemble::SSE> expected_strands[ 3];
      expected_strands[ 0] = sp_strand_5_13;
      expected_strands[ 1] = sp_strand_31_42;
      expected_strands[ 2] = sp_strand_115_120;
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 0), expected_strands[ 0]);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 1), expected_strands[ 1]);
      BCL_ExampleCheck( sp_sheet_a->GetTopology()->GetElements()( 2), expected_strands[ 2]);

      // testing locate with locator_a
      BCL_MessageStd( "testing Locate() with locator_b");
      util::ShPtr< assemble::Domain> sp_sheet_b( locator_b.Locate( sp_sheet));
      // make sure it is defined
      BCL_ExampleCheck( sp_sheet_b.IsDefined(), true);
      // make sure the size is 3
      BCL_ExampleCheck( sp_sheet_b->GetNumberSSEs(), 3);
      // make sure the order of strands are correct
      BCL_ExampleCheck( sp_sheet_b->GetTopology().IsDefined(), true);

      expected_strands[ 0] = sp_strand_5_13;
      expected_strands[ 1] = sp_strand_18_28;
      expected_strands[ 2] = sp_strand_31_42;
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 0), expected_strands[ 0]);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 1), expected_strands[ 1]);
      BCL_ExampleCheck( sp_sheet_b->GetTopology()->GetElements()( 2), expected_strands[ 2]);

      // testing locate with locator_c
      BCL_MessageStd( "testing Locate() with locator_c");
      util::ShPtr< assemble::Domain> sp_sheet_c( locator_c.Locate( sp_sheet));
      // make sure it is defined
      BCL_ExampleCheck( sp_sheet_c.IsDefined(), true);
      // make sure the size is 3
      BCL_ExampleCheck( sp_sheet_c->GetNumberSSEs(), 3);
      // make sure the order of strands are correct
      expected_strands[ 0] = sp_strand_31_42;
      expected_strands[ 1] = sp_strand_18_28;
      expected_strands[ 2] = sp_strand_5_13;
      BCL_ExampleCheck( sp_sheet_c->GetTopology().IsDefined(), true);
      BCL_ExampleCheck( sp_sheet_c->GetTopology()->GetElements()( 0), expected_strands[ 0]);
      BCL_ExampleCheck( sp_sheet_c->GetTopology()->GetElements()( 1), expected_strands[ 1]);
      BCL_ExampleCheck( sp_sheet_c->GetTopology()->GetElements()( 2), expected_strands[ 2]);

    //////////////////////
    // input and output //
    //////////////////////

      // test read/write
      BCL_MessageStd( "testing Read/Write");
      WriteBCLObject( locator_b);
      assemble::LocatorSubDomainRandom locator_read;
      ReadBCLObject( locator_read);
      BCL_ExampleCheck( locator_read.GetSizeRange(), locator_b.GetSizeRange());
      BCL_ExampleCheck( locator_read.GetLocateConsecutive(), locator_b.GetLocateConsecutive());
      BCL_ExampleCheck( locator_read.GetUseTopologyOrder(), locator_b.GetUseTopologyOrder());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorSubDomainRandom

  const ExampleClass::EnumType ExampleAssembleLocatorSubDomainRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorSubDomainRandom())
  );

} // namespace bcl
