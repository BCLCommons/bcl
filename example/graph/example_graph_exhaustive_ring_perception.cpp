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
#include "graph/bcl_graph_exhaustive_ring_perception.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_graph_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_exhaustive_ring_perception.cpp
  //!
  //! @author mueller
  //! @date 01/01/10
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphExhaustiveRingPerception :
    public ExampleInterface
  {
  public:

    ExampleGraphExhaustiveRingPerception *Clone() const
    {
      return new ExampleGraphExhaustiveRingPerception( *this);
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
      // create test graphs
      graph::GraphWithData< size_t, size_t> test_graph_a;

      // insert vertices
      test_graph_a.AddVertex( 1);
      test_graph_a.AddVertex( 2);
      test_graph_a.AddVertex( 3);
      test_graph_a.AddVertex( 4);
      test_graph_a.AddVertex( 5);
      test_graph_a.AddVertex( 6);
      test_graph_a.AddVertex( 7);
      test_graph_a.AddVertex( 8);
      test_graph_a.AddVertex( 9);

      // insert edges
      test_graph_a.AddEdge( 1, 2, 1);
      test_graph_a.AddEdge( 2, 3, 1);
      test_graph_a.AddEdge( 2, 4, 1);
      test_graph_a.AddEdge( 3, 5, 1);
      test_graph_a.AddEdge( 3, 9, 1);
      test_graph_a.AddEdge( 4, 5, 1);
      test_graph_a.AddEdge( 5, 6, 1);
      test_graph_a.AddEdge( 6, 7, 1);
      test_graph_a.AddEdge( 6, 8, 1);

      test_graph_a.AddEdge( 1, 4, 1);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor
      graph::ExhaustiveRingPerception all_rings( test_graph_a);

      BCL_MessageCrt( "all rings: " + util::Format()( all_rings.GetRings()));

      BCL_ExampleCheck
      (
        graph::ExhaustiveRingPerception( test_graph_a).GetRings().GetSize(), 3
      );

      // test clone, first by making sure that it returns a defined object
      if
      (
        BCL_ExampleCheck
        (
          util::ShPtr< graph::ExhaustiveRingPerception>( graph::ExhaustiveRingPerception().Clone()).IsDefined(),
          true
        )
      )
      {
        // test that clone actually copies the object
        BCL_ExampleCheck
        (
          util::ShPtr< graph::ExhaustiveRingPerception>( all_rings.Clone())->GetRings().GetSize(),
          3
        );
      }

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      BCL_ExampleCheck( graph::ExhaustiveRingPerception().GetClassIdentifier(), "bcl::graph::ExhaustiveRingPerception");

      // test GetPaths
      BCL_ExampleCheck( all_rings.GetPaths().IsEmpty(), true);

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

  }; //end ExampleGraphExhaustiveRingPerception

  const ExampleClass::EnumType ExampleGraphExhaustiveRingPerception::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphExhaustiveRingPerception())
  );

} // namespace bcl
