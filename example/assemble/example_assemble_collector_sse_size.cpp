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
#include "assemble/bcl_assemble_collector_sse_size.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_assemble_collector_sse_size.cpp
  //! @brief this example demonstrates how the collector collects all sses of a certian size from a domain
  //!
  //! @author woetzen, karakam
  //! @date Jun 21, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorSSESize :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorSSESize *Clone() const
    {
      return new ExampleAssembleCollectorSSESize( *this);
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

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      const assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      const math::Range< size_t> range( 5, 8);
      storage::Map< biol::SSType, math::Range< size_t> > range_map_a;
      range_map_a[ biol::GetSSTypes().HELIX] = range;
      range_map_a[ biol::GetSSTypes().STRAND] = range;
      range_map_a[ biol::GetSSTypes().COIL] = range;

      storage::Map< biol::SSType, math::Range< size_t> > range_map_b;
      range_map_b[ biol::GetSSTypes().HELIX] = math::Range< size_t>( 2, 6);
      range_map_b[ biol::GetSSTypes().STRAND] = math::Range< size_t>( 3, 8);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from range
      assemble::CollectorSSESize collector( range);

      // construct from range map
      assemble::CollectorSSESize collector_a( range_map_a);
      assemble::CollectorSSESize collector_b( range_map_b);

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( collector.GetClassIdentifier(), GetStaticClassName< assemble::CollectorSSESize>());

      // size range
      BCL_ExampleCheck( collector.GetRanges(), range_map_a);
      BCL_ExampleCheck( collector_a.GetRanges(), range_map_a);
      BCL_ExampleCheck( collector_b.GetRanges(), range_map_b);

    ////////////////
    // operations //
    ////////////////

      // collect
      util::SiPtrList< const assemble::SSE> collected_sses( collector.Collect( protein_model));
      BCL_ExampleIndirectCheck( collected_sses.GetSize(), 5, "collect check if number of sses matches");

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( collector);
      assemble::CollectorSSESize collector_read( math::Range< size_t>( 0, 0));
      ReadBCLObject( collector_read);
      BCL_ExampleIndirectCheck( collector_read.GetRanges(), range_map_a, "write and read of bcl object");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorSSESize

  const ExampleClass::EnumType ExampleAssembleCollectorSSESize::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorSSESize())
  );

} // namespace bcl
