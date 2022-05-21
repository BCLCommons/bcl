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
#include "graph/bcl_graph_edge_cover_ring_perception.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_edge_cover_ring_perception.cpp
  //!
  //! @author mendenjl
  //! @date September 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphEdgeCoverRingPerception :
    public ExampleInterface
  {
  public:

    ExampleGraphEdgeCoverRingPerception *Clone() const
    {
      return new ExampleGraphEdgeCoverRingPerception( *this);
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
      // insert vertices
      storage::Vector< size_t> vertices( 9, 0);

      linal::Matrix< size_t> edges( size_t( 9), size_t( 9), size_t( 0));

      // insert edges
      edges( 0, 1) = 1;
      edges( 1, 2) = 1;
      edges( 1, 3) = 1;
      edges( 2, 4) = 1;
      edges( 2, 8) = 1;
      edges( 3, 4) = 1;
      edges( 4, 5) = 1;
      edges( 5, 6) = 1;
      edges( 5, 7) = 1;
      edges( 0, 3) = 1;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      graph::ConstGraph< size_t, size_t> test_graph_a( vertices, edges, 0);

      // test constructor
      graph::EdgeCoverRingPerception all_rings( test_graph_a);

      BCL_ExampleCheck( graph::EdgeCoverRingPerception( test_graph_a).GetRings().GetSize(), 2);

      // test clone, first by making sure that it returns a defined object
      if
      (
        BCL_ExampleCheck
        (
          util::ShPtr< graph::EdgeCoverRingPerception>( graph::EdgeCoverRingPerception().Clone()).IsDefined(),
          true
        )
      )
      {
        // test that clone actually copies the object
        BCL_ExampleCheck
        (
          util::ShPtr< graph::EdgeCoverRingPerception>( all_rings.Clone())->GetRings().GetSize(),
          2
        );
      }

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      BCL_ExampleCheck( graph::EdgeCoverRingPerception().GetClassIdentifier(), "bcl::graph::EdgeCoverRingPerception");

      // test a harder example: buckminsterfullerene.  There should be 32 rings (20 6-membered, and 12 5-membered)
      // setup input stream
      io::IFStream input_sdf;
      // read in molecule
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "c60_fullerene.sdf"));
      chemistry::FragmentComplete small_mol( sdf::FragmentFactory::MakeFragment( input_sdf, sdf::e_Remove));
      const graph::ConstGraph< size_t, size_t> const_graph_c60( chemistry::ConformationGraphConverter()( small_mol));
      const graph::EdgeCoverRingPerception bucky_rings( const_graph_c60);

      storage::Vector< size_t> c60_fullerene_ring_sizes;
      for
      (
        storage::List< graph::Ring>::const_iterator itr( bucky_rings.GetRings().Begin()),
          itr_end( bucky_rings.GetRings().End());
        itr != itr_end;
        ++itr
      )
      {
        c60_fullerene_ring_sizes.PushBack( itr->GetSize());
      }

      //  There should be 32 rings (20 6-membered, and 12 5-membered)
      BCL_ExampleIndirectCheck( c60_fullerene_ring_sizes.GetSize(), 32, "graph::EdgeCoverRingPerception()");
      BCL_ExampleIndirectCheck
      (
        size_t( std::count( c60_fullerene_ring_sizes.Begin(), c60_fullerene_ring_sizes.End(), size_t( 5))),
        12,
        "graph::EdgeCoverRingPerception()"
      );
      BCL_ExampleIndirectCheck
      (
        size_t( std::count( c60_fullerene_ring_sizes.Begin(), c60_fullerene_ring_sizes.End(), size_t( 6))),
        20,
        "graph::EdgeCoverRingPerception()"
      );

    ///////////////
    // operators //
    ///////////////

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphEdgeCoverRingPerception

  const ExampleClass::EnumType ExampleGraphEdgeCoverRingPerception::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphEdgeCoverRingPerception())
  );

} // namespace bcl
