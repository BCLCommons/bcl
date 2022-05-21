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
#include "assemble/bcl_assemble_collector_topology_combined.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_topology.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_topology_combined.cpp
  //!
  //! @author weinerbe
  //! @date Oct 26, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorTopologyCombined :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorTopologyCombined *Clone() const
    {
      return new ExampleAssembleCollectorTopologyCombined( *this);
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

      // create protein model from pdb
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 3;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // locate 2 strands
      const assemble::LocatorSSE locator_1_7( 'A', 1, 7);
      const assemble::LocatorSSE locator_10_17( 'A', 10, 17);
      util::ShPtr< assemble::SSE> sp_strand_1_7( locator_1_7.Locate( protein_model)->Clone());
      util::ShPtr< assemble::SSE> sp_strand_10_17( locator_10_17.Locate( protein_model)->Clone());

      // move the SSEs far away so they are no longer part of the topology
      const linal::Vector3D translation( 0.0, 0.0, 50.0);
      sp_strand_1_7->Translate( translation);
      sp_strand_10_17->Translate( translation);

      // update the protein model
      protein_model.Replace( sp_strand_1_7);
      protein_model.Replace( sp_strand_10_17);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const assemble::CollectorTopologyCombined collector;

    ////////////////
    // operations //
    ////////////////

      // test Collect
      const util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());
      const util::ShPtrVector< assemble::Topology> all_topologies( collector.Collect( sses));
      BCL_ExampleCheck( all_topologies.GetSize(), 2);

      // test CalculateTopology
      const assemble::Topology complete_topology( collector.CalculateTopology( sses));
      BCL_ExampleCheck( complete_topology.GetGraph().GetNumberVertices(), 6);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorTopologyCombined

  const ExampleClass::EnumType ExampleAssembleCollectorTopologyCombined::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorTopologyCombined())
  );

} // namespace bcl
