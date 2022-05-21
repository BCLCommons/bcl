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
#include "graph/bcl_graph_subgraph.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_subgraph.cpp
  //!
  //! @author mendenjl
  //! @date Jan 18, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphSubgraph :
    public ExampleInterface
  {
  public:

    ExampleGraphSubgraph *Clone() const
    {
      return new ExampleGraphSubgraph( *this);
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

      graph::Subgraph< size_t, size_t> default_subgraph;

    /////////////////
    // data access //
    /////////////////

      // test that the default constructed object was constructed properly
      BCL_ExampleCheck( default_subgraph.GetSize(), 0);
      BCL_ExampleCheck( default_subgraph.ToGraph().GetSize(), 0);

      // make a graph containing 5 vertices connected in a chain as follows:
      // 4-3-1-0-2
      linal::Matrix< size_t> edges( 5, 5, size_t( 0));
      edges( 4, 3) = edges( 3, 1) = edges( 1, 0) = edges( 0, 2) = 1;

      // make the vertices of the graph be the vertex id * 2
      storage::Vector< size_t> vertex_id( 5);
      vertex_id( 0) = 0;
      vertex_id( 1) = 2;
      vertex_id( 2) = 4;
      vertex_id( 3) = 6;
      vertex_id( 4) = 8;

      // make the graph and a non-owning pointer to it
      graph::ConstGraph< size_t, size_t> graph( vertex_id, edges);
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > graph_ptr( &graph, false);

      // now make a few different sets of vertices to make subgraphs with
      storage::Vector< size_t> empty;
      storage::Vector< size_t> zero( 1, size_t( 0));
      storage::Vector< size_t> three_one( storage::Vector< size_t>::Create( 3, 1));
      storage::Vector< size_t> three_one_four( storage::Vector< size_t>::Create( 3, 1, 4));

      // now create the subgraphs
      graph::Subgraph< size_t, size_t> subgraph_empty( graph_ptr, empty);
      graph::Subgraph< size_t, size_t> subgraph_0( graph_ptr, zero);
      graph::Subgraph< size_t, size_t> subgraph_3_1( graph_ptr, three_one);
      graph::Subgraph< size_t, size_t> subgraph_3_1_4( graph_ptr, three_one_four);

      // test get size
      BCL_ExampleCheck( subgraph_empty.GetSize(), 0);
      BCL_ExampleCheck( subgraph_0.GetSize(), 1);
      BCL_ExampleCheck( subgraph_3_1.GetSize(), 2);
      BCL_ExampleCheck( subgraph_3_1_4.GetSize(), 3);

      // test GetParentGraphPtr
      BCL_ExampleCheck( subgraph_empty.GetParentGraphPtr(), graph_ptr);

      // test get vertex indices
      BCL_ExampleCheck( subgraph_empty.GetVertexIndices().GetSize(), 0);
      BCL_ExampleCheck( subgraph_0.GetVertexIndices()( 0), 0);
      BCL_ExampleCheck( subgraph_3_1.GetVertexIndices()( 0), 3);
      BCL_ExampleCheck( subgraph_3_1.GetVertexIndices()( 1), 1);

    ////////////////
    // operations //
    ////////////////

      // test to graph
      BCL_ExampleCheck( subgraph_3_1.ToGraph().GetSize(), 2);
      BCL_ExampleCheck( subgraph_3_1.ToGraph().NumEdges(), 2); // 1-3, 3-1
      BCL_ExampleCheck( subgraph_3_1_4.ToGraph().NumEdges(), 4); // 1-3, 3-1, 1-4, 4-1

      // test GetComplement()
      storage::Vector< size_t> complement_3_1( subgraph_3_1.GetComplement().GetVertexIndices());
      if( BCL_ExampleIndirectCheck( complement_3_1.GetSize(), 3, "GetComplement"))
      {
        BCL_ExampleIndirectCheck( complement_3_1( 0), 0, "GetComplement");
        BCL_ExampleIndirectCheck( complement_3_1( 1), 2, "GetComplement");
        BCL_ExampleIndirectCheck( complement_3_1( 2), 4, "GetComplement");
      }

      // test GetEdgeIndices()
      storage::List< storage::Pair< size_t, size_t> > edge_indices( subgraph_3_1_4.GetEdgeIndices());
      storage::Vector< storage::Pair< size_t, size_t> > edge_indices_vec( edge_indices.Begin(), edge_indices.End());
      if( BCL_ExampleIndirectCheck( edge_indices_vec.GetSize(), 4, "GetEdgeIndices"))
      {
        BCL_ExampleIndirectCheck
        (
          edge_indices_vec( 0),
          ( storage::Pair< size_t, size_t>( 3, 1)),
          "GetEdgeIndices"
        );

        BCL_ExampleIndirectCheck
        (
          edge_indices_vec( 1),
          ( storage::Pair< size_t, size_t>( 3, 4)),
          "GetEdgeIndices"
        );

        BCL_ExampleIndirectCheck
        (
          edge_indices_vec( 2),
          ( storage::Pair< size_t, size_t>( 1, 3)),
          "GetEdgeIndices"
        );

        BCL_ExampleIndirectCheck
        (
          edge_indices_vec( 3),
          ( storage::Pair< size_t, size_t>( 4, 3)),
          "GetEdgeIndices"
        );
      }

      // test GetIdsOfInteriorVertices.  For 3-1-4, 3 and 4 are interior vertices because neither has edges
      // that lead outside the subgraph
      storage::Vector< size_t> interior_vertices( subgraph_3_1_4.GetIdsOfInteriorVertices());
      if( BCL_ExampleIndirectCheck( interior_vertices.GetSize(), 2, "GetIdsOfInteriorVertices"))
      {
        BCL_ExampleIndirectCheck
        (
          interior_vertices.FirstElement(),
          3,
          "GetIdsOfInteriorVertices"
        );

        BCL_ExampleIndirectCheck
        (
          interior_vertices.LastElement(),
          4,
          "GetIdsOfInteriorVertices"
        );
      }

      // for 3-1, there should be no interior vertices, since both 3 and 1 have connections to vertices outside the
      // subgraph
      BCL_ExampleCheck( subgraph_3_1.GetIdsOfInteriorVertices().GetSize(), 0);

      // test GetAdjacentEdgeIndices
      // there should be only 1 edge adjacent to 4-3-1 (1-0)
      storage::List< storage::Pair< size_t, size_t> > adjacent_edge_indices( subgraph_3_1_4.GetAdjacentEdgeIndices());
      if( BCL_ExampleIndirectCheck( adjacent_edge_indices.GetSize(), 1, "GetAdjacentEdgeIndices"))
      {
        BCL_ExampleIndirectCheck
        (
          adjacent_edge_indices.FirstElement(),
          ( storage::Pair< size_t, size_t>( 1, 0)),
          "GetAdjacentEdgeIndices"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphSubgraph

  const ExampleClass::EnumType ExampleGraphSubgraph::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphSubgraph())
  );

} // namespace bcl
