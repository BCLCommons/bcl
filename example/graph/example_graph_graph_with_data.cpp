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
#include "graph/bcl_graph_graph_with_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    // explicit template instantiation
    template class VertexWithData< size_t, size_t>;
    template class EdgeWithData< size_t, size_t>;
  } // namespace graph

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_graph_with_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphGraphWithData :
    public ExampleInterface
  {
  public:

    ExampleGraphGraphWithData *Clone() const
    {
      return new ExampleGraphGraphWithData( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a test graph
      graph::GraphWithData< size_t, size_t> test_graph;

      // test copy constructor
      graph::GraphWithData< size_t, size_t> test_graph_a( test_graph);

      // test Clone
      util::ShPtr< graph::GraphWithData< size_t, size_t> > test_graph_b( test_graph_a.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // temporary variables with the data contained in each of the vertices
      const size_t vertex_one( 32), vertex_two( 23), vertex_three( 66), vertex_four( 99), vertex_five( 100);

      // add vertices to the test graph
      test_graph.AddVertex( vertex_one);
      test_graph.AddVertex( vertex_two);
      test_graph.AddVertex( vertex_three);
      test_graph.AddVertex( vertex_four);
      test_graph.AddVertex( vertex_five);

      // connect the vertices
      test_graph.AddEdge( vertex_one, vertex_two, 1);
      test_graph.AddEdge( vertex_two, vertex_three, 1);
      test_graph.AddEdge( vertex_two, vertex_four, 1);
      test_graph.AddEdge( vertex_one, vertex_four, 1);
      test_graph.AddEdge( vertex_five, vertex_two, 1);
      test_graph.AddEdge( vertex_five, vertex_three, 1);
      test_graph.AddEdge( vertex_five, vertex_four, 1);

      // there should be five vertices and 14 edges (e.g. 7 undirected edges)
      BCL_ExampleCheck( test_graph.GetVertices().GetSize(), 5);
      BCL_ExampleCheck( test_graph.GetNumberEdges(), 14);

      // output the Graph
      BCL_MessageStd( "Graph: " + util::Format()( test_graph) + "\n");

      // check if there is a vertex containing the number 23
      BCL_ExampleCheck( test_graph.HasVertex( 23), true);

      // delete Vertex( 23) and output the new Graph
      test_graph.DeleteVertex( test_graph.FindVertex( 23));
      BCL_MessageStd
      (
        "Graph after deletion of vertex containing 23: " + util::Format()( test_graph) + "\n"
      );

      // check if there is a vertex containing the number 23 (it should no longer exist)
      BCL_ExampleIndirectCheck
      (
        test_graph.HasVertex( 23),
        false,
        "test_graph.DeleteVertex( test_graph.FindVertex( 23))"
      );

      // check if Vertex( 32) and Vertex( 99) are connected. They should be
      BCL_ExampleCheck( test_graph.AreVerticesConnected( 32, 99), true);

      // perform the same check on Vertex( 32) and Vertex( 100), which are not linked
      BCL_ExampleCheck( test_graph.AreVerticesConnected( 32, 100), false);

      // delete the edge between Vertex( 32) and Vertex( 99)
      test_graph.DeleteEdge( 32, 99);

      BCL_MessageStd
      (
        "Graph after deletion of Edge( 32, 99): " + util::Format()( test_graph) + "\n"
      );

      // check if there is still a connection between Vertex( 32) and Vertex( 99)
      BCL_ExampleIndirectCheck
      (
        test_graph.AreVerticesConnected( 32, 99), false, "test_graph.DeleteEdge( 32, 99);"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphGraphWithData

  const ExampleClass::EnumType ExampleGraphGraphWithData::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphGraphWithData())
  );

} // namespace bcl
