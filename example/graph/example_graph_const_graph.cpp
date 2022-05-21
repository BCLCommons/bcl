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
#include "graph/bcl_graph_const_graph.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_const_graph.cpp
  //!
  //! @author mendenjl
  //! @date July 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphConstGraph :
    public ExampleInterface
  {
  public:

    ExampleGraphConstGraph *Clone() const
    {
      return new ExampleGraphConstGraph( *this);
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

      // create a const graph with string for vertex data and size_t's for edge data
      graph::ConstGraph< std::string, size_t> default_const_graph;

      // check that the default const graph does is empty and has no edges
      BCL_ExampleCheck( default_const_graph.GetSize(), 0);
      BCL_ExampleCheck( default_const_graph.NumEdges(), 0);

      // create a vector of strings to be the vertex data
      storage::Vector< std::string> subject_names( 3);
      subject_names( 0) = "politics";
      subject_names( 1) = "physics";
      subject_names( 2) = "chemistry";

      // create a matrix containing 5 for subjects that are closely related
      linal::Matrix< size_t> subject_relations_directed( 3, 3, util::GetUndefined< size_t>());
      subject_relations_directed( 2, 1) = 5; // Chemical physics

      // also make an undirected version of the matrix
      linal::Matrix< size_t> subject_relations_undirected( subject_relations_directed);
      subject_relations_undirected( 1, 2) = 5; // Physical chemistry

      // construct a graph of the relation
      // by default, if the directness is unspecified, the graph will automatically be made undirected
      graph::ConstGraph< std::string, size_t> subject_graph_undirected
      (
        subject_names,
        subject_relations_directed,
        util::GetUndefined< size_t>()
      );
      BCL_ExampleCheck( subject_graph_undirected.GetSize(), subject_names.GetSize());
      BCL_ExampleCheck( subject_graph_undirected.NumEdges(), 2);
      BCL_ExampleCheck( subject_graph_undirected.GetVertices(), subject_names);
      BCL_ExampleCheck( subject_graph_undirected.IsDirected(), false);
      BCL_ExampleCheck( subject_graph_undirected.IsUndirected(), true);
      BCL_ExampleCheck( subject_graph_undirected.GetVertexData( 1), "physics");
      BCL_ExampleCheck( subject_graph_undirected.GetEdgeDataMatrix(), subject_relations_undirected);

      // We can force a graph to be directed by passing two bools
      // The first bool indicates that the
      // The second bool indicates that the matrix should not be forced into undirectedness
      graph::ConstGraph< std::string, size_t> subject_graph_directed
      (
        subject_names,
        subject_relations_directed,
        util::GetUndefined< size_t>(),
        true, // graph shall be directed
        false // do not make the ininitial connectivity undirected
      );
      BCL_ExampleCheck( subject_graph_directed.GetSize(), subject_names.GetSize());
      BCL_ExampleCheck( subject_graph_directed.NumEdges(), 1);
      BCL_ExampleCheck( subject_graph_directed.GetVertices(), subject_names);
      BCL_ExampleCheck( subject_graph_directed.IsDirected(), true);
      BCL_ExampleCheck( subject_graph_directed.IsUndirected(), false);
      BCL_ExampleCheck( subject_graph_directed.GetVertexData( 1), "physics");
      BCL_ExampleCheck
      (
        std::equal
        (
          subject_graph_directed.GetEdgeDataMatrix().Begin(),
          subject_graph_directed.GetEdgeDataMatrix().End(),
          subject_relations_directed.Begin()
        ),
        true
      );

      // test out clone
      util::ShPtr< graph::ConstGraph< std::string, size_t> > shptr_subject_graph_directed
      (
        subject_graph_directed.Clone()
      );
      BCL_ExampleAssert( shptr_subject_graph_directed.IsDefined(), true);
      BCL_ExampleCheck( shptr_subject_graph_directed->GetSize(), subject_names.GetSize());
      BCL_ExampleCheck( shptr_subject_graph_directed->NumEdges(), 1);
      BCL_ExampleCheck( shptr_subject_graph_directed->GetVertices(), subject_names);
      BCL_ExampleCheck( shptr_subject_graph_directed->IsDirected(), true);
      BCL_ExampleCheck( shptr_subject_graph_directed->IsUndirected(), false);
      BCL_ExampleCheck( shptr_subject_graph_directed->GetVertexData( 1), "physics");

    /////////////////
    // data access //
    /////////////////

      // make sure that the neighbor data was determined properly
      BCL_ExampleCheck( subject_graph_directed.GetNeighborData( 0).GetSize(), 0);
      BCL_ExampleCheck( subject_graph_directed.GetNeighborData( 1).GetSize(), 0);
      BCL_ExampleCheck( subject_graph_directed.GetNeighborData( 2).GetSize(), 1);
      BCL_ExampleCheck( subject_graph_directed.GetNeighborData( 2)( 0), 5);
      BCL_ExampleCheck( subject_graph_directed.GetNeighborIndices( 2).GetSize(), 1);
      BCL_ExampleCheck( subject_graph_directed.GetNeighborIndices( 2)( 0), 1);
      BCL_ExampleCheck( subject_graph_directed.GetEdgeData( 2, 1), 5);
      BCL_ExampleCheck( subject_graph_directed.GetEdgeData( 1, 2), util::GetUndefined< size_t>());
      BCL_ExampleCheck( subject_graph_directed.AreConnected( 2, 1), true);
      BCL_ExampleCheck( subject_graph_directed.AreConnected( 1, 2), false);
      BCL_ExampleCheck( subject_graph_undirected.AreConnected( 1, 2), true);
      BCL_ExampleCheck( subject_graph_directed.AreVerticesConnected( "chemistry", "physics"), true);
      BCL_ExampleCheck( subject_graph_directed.AreVerticesConnected( "physics", "chemistry"), false);
      BCL_ExampleCheck( subject_graph_undirected.AreVerticesConnected( "physics", "chemistry"), true);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // get the subgraph with just physics and chemistry
      storage::Vector< size_t> chemistry_and_physics( 2);
      chemistry_and_physics( 0) = 2; // put chemistry first this time.  Now chemistry should be at index 0
      chemistry_and_physics( 1) = 1;
      graph::ConstGraph< std::string, size_t> chemistry_and_physics_undirected_subgraph
      (
        subject_graph_undirected.GetSubgraph( chemistry_and_physics)
      );
      storage::Vector< size_t> chem_and_phys_mapping;
      util::SiPtr< storage::Vector< size_t> > sip_chem_and_phys_mapping( &chem_and_phys_mapping);
      graph::ConstGraph< std::string, size_t> chemistry_and_physics_undirected_subgraph_w_mapping
      (
        subject_graph_undirected.GetSubgraph( chemistry_and_physics, sip_chem_and_phys_mapping)
      );
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.AreVerticesConnected( "chemistry", "physics"), true);
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.AreVerticesConnected( "physics", "chemistry"), true);
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.AreVerticesConnected( "physics", "physics"), false);
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.GetSize(), 2);
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.IsUndirected(), true);
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.GetVertexData( 0), "chemistry");
      BCL_ExampleCheck( chemistry_and_physics_undirected_subgraph.GetVertexData( 1), "physics");

      size_t number_mapped( 0);
      for( size_t i( 0); i < chem_and_phys_mapping.GetSize(); ++i)
      {
        if( util::IsDefined( chem_and_phys_mapping( i)))
        {
          ++number_mapped;
        }
      }

      BCL_ExampleIndirectCheck
      (
        number_mapped,
        chemistry_and_physics_undirected_subgraph_w_mapping.GetSize(),
        "subgraph has correct number of vertices mapped"
      );

      // test basic connectivity
      BCL_MessageStd
      (
        "Here is the basic connectivity for the undirected subject const graph:\n"
        + subject_graph_undirected.GetBasicConnectivity()
      );

      // get a map from vertex datum to count
      storage::Map< std::string, size_t> vertex_counts
      (
        subject_graph_undirected.GetVertexTypeCountMap()
      );
      BCL_ExampleIndirectCheck( vertex_counts.GetSize(), 3, "GetVertexTypeCountMap");
      BCL_ExampleIndirectCheck( vertex_counts[ "physics"], 1, "GetVertexTypeCountMap");
      BCL_ExampleIndirectCheck( vertex_counts[ "music"], 0, "GetVertexTypeCountMap");

      // get a map from edges to counts
      storage::Map< size_t, size_t> edge_counts
      (
        subject_graph_undirected.GetEdgeTypeCountMap()
      );
      BCL_ExampleIndirectCheck( edge_counts.GetSize(), 1, "GetEdgeTypeCountMap");

      // check that the first element is 5, for which there should be a count of 2
      BCL_ExampleIndirectCheck( edge_counts.Begin()->first, 5, "GetEdgeTypeCountMap");
      BCL_ExampleIndirectCheck( edge_counts.Begin()->second, 2, "GetEdgeTypeCountMap");

      // check that editing an edge for an undirected edge changes both edges
      subject_graph_undirected.ChangeEdge( 2, 1, util::GetUndefined< size_t>());
      BCL_ExampleIndirectCheck( subject_graph_undirected.NumEdges(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborData( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborIndices( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborData( 2).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborIndices( 2).GetSize(), 0, "ChangeEdge");

      // restore the edge to its original condition
      subject_graph_undirected.ChangeEdge( 2, 1, 5);
      BCL_ExampleIndirectCheck( subject_graph_undirected.NumEdges(), 2, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborData( 1).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborIndices( 1).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborData( 2).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborIndices( 2).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborData( 2)( 0), 5, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_undirected.GetNeighborIndices( 2)( 0), 1, "ChangeEdge");

      // on a directed graph, removing an edge should remove only that edge
      subject_graph_directed.ChangeEdge( 2, 1, util::GetUndefined< size_t>());
      BCL_ExampleIndirectCheck( subject_graph_directed.NumEdges(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborData( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborIndices( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborData( 2).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborIndices( 2).GetSize(), 0, "ChangeEdge");

      // restore the edge to its original condition
      subject_graph_directed.ChangeEdge( 2, 1, 5);
      BCL_ExampleIndirectCheck( subject_graph_directed.NumEdges(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborData( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborIndices( 1).GetSize(), 0, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborData( 2).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborIndices( 2).GetSize(), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborData( 2)( 0), 5, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.GetNeighborIndices( 2)( 0), 1, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.AreVerticesConnected( "chemistry", "physics"), true, "ChangeEdge");
      BCL_ExampleIndirectCheck( subject_graph_directed.AreVerticesConnected( "physics", "chemistry"), false, "ChangeEdge");

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( chemistry_and_physics_undirected_subgraph, default_const_graph),
        true,
        "ConstGraph I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphConstGraph

  const ExampleClass::EnumType ExampleGraphConstGraph::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphConstGraph())
  );

} // namespace bcl
