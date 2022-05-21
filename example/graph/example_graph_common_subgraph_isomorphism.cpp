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
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_connectivity.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @class Egalitarian
  //! @brief a functor that always returns true no matter what object is compared
  //!
  //! @author mendenjl
  //! @date Aug 09, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class t_DataType>
  class Egalitarian :
    public util::BinaryFunctionInterface< t_DataType, t_DataType, bool>
  {
  public:

    //! single instance of that class
    static const util::SiPtr< const util::ObjectInterface> s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new class_name< t_typename_a>
    Egalitarian *Clone() const
    {
      return new Egalitarian( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @return true
    bool operator ()( const t_DataType &, const t_DataType &) const
    {
      // always return true
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  protected:

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  };

  // instantiate s_Instance
  template< typename t_DataType>
  const util::SiPtr< const util::ObjectInterface> Egalitarian< t_DataType>::s_Instance
  (
    GetObjectInstances().AddInstance( new Egalitarian< t_DataType>())
  );

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_common_subgraph_isomorphism.cpp
  //!
  //! @author mendenjl
  //! @date Aug 09, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphCommonSubgraphIsomorphism :
    public ExampleInterface
  {
  public:

    //! @brief Clone
    ExampleGraphCommonSubgraphIsomorphism *Clone() const
    {
      return new ExampleGraphCommonSubgraphIsomorphism( *this);
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

    //! @brief Run tests the common subgraph isomorphism class
    //! @note  This class tests connected and unconnected implementations of the common subgraph isomorphism, as well
    //! @note  as the isomorphisms ability to use user-provided comparison objects
    int Run() const
    {
      // These parameters control how large of isomorphisms we test.  If you want to test larger examples, or examples
      // with more data, change these parameters
      const size_t number_of_test_graphs = 5;
      const size_t number_of_edge_data   = 4; // Edge data are assigned randomly from 1 to number_of_edge_data
      const size_t number_of_vertex_data = 4;
      const size_t min_graph_size        = 8;
      const size_t max_graph_size        = 15;
      const double average_connections   = 5;

      BCL_Assert( average_connections < min_graph_size, " too many connections for some graphs");
      BCL_Assert( average_connections >= 2.5, " too few connections -> high probability of unconnected graph");
      BCL_Assert( number_of_edge_data > 1, " not enough edge data to indicate connected/unconnected!");
      BCL_Assert( number_of_vertex_data > 0, " must allow at least one vertex datum!");

      // create a new egalitarian object to compare vertex edges
      util::ShPtr< util::BinaryFunctionInterface< size_t, size_t, bool> > do_not_compare
      (
        new Egalitarian< size_t>()
      );

      // create a new egalitarian object to compare vertex edges
      util::ShPtr< util::BinaryFunctionInterface< size_t, size_t, bool> > do_compare
      (
        *math::Comparisons< size_t>::GetEnums().e_Equal
      );

      graph::CommonSubgraphIsomorphism< size_t, size_t> connected_subgraph;

      // this test generates a randomly connected graph_a, copies it to graph_b, and then randomly permutes the
      // vertices in test_graph_b, and then runs the isomorphism
      // the result should be the same as the size of the graph.
      for( size_t graph_number = 0; graph_number < number_of_test_graphs; graph_number++)
      {
        const size_t graph_size( random::GetGlobalRandom().Random< size_t>( min_graph_size, max_graph_size));

        // make test_graph_a be a random graph given the min/max graph size, number of edge data, and so forth,

        // Construct vectors for edge and vertex datum ratios
        storage::Vector< double> vertex_datum_ratios( number_of_vertex_data);
        storage::Vector< double> edge_datum_ratios( number_of_edge_data);

        // use random numbers in for the vertex datum ratios
        for( size_t a = 0; a < number_of_vertex_data; a++)
        {
          vertex_datum_ratios( a) = random::GetGlobalRandom().Double();
        }

        // do roughly the same thing for edge_datum ratios, except datum 0 (unconnected edges)
        edge_datum_ratios( 0) = 0.0;
        for( size_t a = 1; a < number_of_edge_data; a++)
        {
          edge_datum_ratios( a) = random::GetGlobalRandom().Double();
        }

        // set the ratio of the first edge_datum (which represents unconnected edges)
        // to make the average_connections work out right
        const double connected_probability
        (
          double( average_connections * graph_size)
          /
          double( graph_size * ( graph_size - 1))
        );
        const double unconnection_probability( 1.0 - connected_probability);
        const double sum_edge_proportions
        (
          std::accumulate( edge_datum_ratios.Begin(), edge_datum_ratios.End(), 0.0)
        );

        edge_datum_ratios( 0) = unconnection_probability * ( sum_edge_proportions / connected_probability);
        BCL_MessageDbg( " edge datum ratios: " + util::Format()( edge_datum_ratios));
        BCL_MessageDbg( " vertex datum ratios: " + util::Format()( vertex_datum_ratios));

        // create a test graph
        graph::ConstGraph< size_t, size_t> test_graph_a
        (
          graph::ConstGraph< size_t, size_t>::MakeRandomUndirectedGraph
          (
            graph_size,
            edge_datum_ratios,
            vertex_datum_ratios
          )
        );

        // now we will create an isomorphism of test_graph_a

        // make the const-graph with the reordered vertices and edges
        graph::ConstGraph< size_t, size_t> test_graph_b( test_graph_a);

        // reorder the rows-columns of the matrix of edge data, store the permutation used
        storage::Vector< size_t> actual_permutation( test_graph_b.Shuffle());

        // check to make sure we got a different permutation, otherwise we are just testing equality of the graphs
        // (which could happen with a really odd choice of random numbers)
        if( actual_permutation( 0) == 0)
        {
          size_t last_in_order = 0;
          while( last_in_order < graph_size && actual_permutation( last_in_order) == last_in_order)
          {
            last_in_order++;
          }

          if( last_in_order == graph_size)
          {
            // find a new graph
            continue;
          }
        }

        BCL_MessageDbg( "graph a:\n" + test_graph_a.GetBasicConnectivity());
        BCL_MessageDbg( "graph b:\n" + test_graph_b.GetBasicConnectivity());
        BCL_MessageDbg( "actual permutation:\n" + util::Format()( actual_permutation));

        // test whether the graph can find the largest unconnected substructure (i.e. the whole graph)
        util::Stopwatch unconnected_timer;
        unconnected_timer.Reset();
        graph::CommonSubgraphIsomorphism< size_t, size_t>
        unconnected_subgraph( graph::CommonSubgraphIsomorphismBase::e_Unconnected),
        unconnected_subgraph_larger_smaller( graph::CommonSubgraphIsomorphismBase::e_Unconnected);

        unconnected_subgraph.FindIsomorphism
        (
          test_graph_a,
          test_graph_b,
          graph::CommonSubgraphIsomorphismBase::EstimateUpperBounds( test_graph_a, test_graph_b)
        );

        BCL_MessageStd
        (
          "unconnected algorithm time: "
          + util::Format()( unconnected_timer.GetProcessDuration())
        );

        BCL_ExampleIndirectCheck
        (
          unconnected_subgraph.GetIsomorphism().GetSize(), graph_size,
          "Finding the largest unconnected common substructure"
        );

        BCL_MessageStd
        (
          "    detected isomorphism size: "
          + util::Format()( unconnected_subgraph.GetIsomorphism().GetSize())
          + " actual isomorphism size: " + util::Format()( graph_size)
        );

        const bool graph_is_connected( graph::Connectivity::IsConnected( test_graph_a));

        // ensure that the graph is connected before trying the connected algorithm
        if( graph_is_connected) // unconnected graph.  Try again
        {
          BCL_MessageDbg
          (
            " Generated graph (is connected):\n" + util::Format()( test_graph_a)
          );
        }
        else
        {
          BCL_MessageDbg
          (
            " Generated graph (is unconnected):\n" + util::Format()( test_graph_a)
          );
        }

        // test the connected version of the algorithm
        util::Stopwatch connected_timer;
        connected_timer.Reset();

        connected_subgraph.SetGraphs( test_graph_a, test_graph_b);

        connected_subgraph.FindIsomorphism
        (
          test_graph_a,
          test_graph_b,
          graph::CommonSubgraphIsomorphismBase::EstimateUpperBounds( test_graph_a, test_graph_b)
        );

        BCL_MessageStd
        (
          "connected algorithm time: " + util::Format()( connected_timer.GetProcessDuration())
        );

        BCL_MessageStd
        (
          "    isomorphism size: " + util::Format()( connected_subgraph.GetIsomorphism().GetSize())
        );

        if( graph_is_connected) // connected graph, should find an isomorphism of graph_size
        {
          BCL_ExampleIndirectCheck
          (
            connected_subgraph.GetIsomorphism().GetSize(), graph_size,
            "Finding largest connected common substructure on a connected graph"
          );
        }
        else
        {
          BCL_ExampleIndirectCheck
          (
            connected_subgraph.GetIsomorphism().GetSize() < graph_size, true,
            "Finding largest connected common substructure on an unconnected graph"
          );
        }

        if // the vertex data were not equal under the permutation, then test whether vertex data are checked
        (
          !std::equal
          (
            test_graph_a.GetVertices().Begin(),
            test_graph_a.GetVertices().End(),
            test_graph_b.GetVertices().Begin()
          )
        )
        {
          // this tests that the vertex data are checked
          // make a graph containing the vertex data in the same order as graph a, but the edge data in the same
          // order as graph b
          graph::ConstGraph< size_t, size_t> test_graph_c
          (
            test_graph_a.GetVertices(),
            test_graph_b.GetEdgeDataMatrix(),
            0
          );

          graph::CommonSubgraphIsomorphism< size_t, size_t>
            unconnected_subgraph_mismatched_verts( graph::CommonSubgraphIsomorphismBase::e_Unconnected);

          unconnected_subgraph_mismatched_verts.FindIsomorphism
          (
            test_graph_a,
            test_graph_c,
            graph::CommonSubgraphIsomorphismBase::EstimateUpperBounds( test_graph_a, test_graph_c)
          );

          BCL_ExampleIndirectCheck
          (
            unconnected_subgraph_mismatched_verts.GetIsomorphism().GetSize() < graph_size, true,
            "Finding the largest unconnected common substructure, considering vertex data"
          );

          graph::CommonSubgraphIsomorphism< size_t, size_t>
          unconnected_subgraph_ignore_vertices
          (
            graph::CommonSubgraphIsomorphismBase::e_Unconnected,
            do_compare,
            do_not_compare
          );

          unconnected_subgraph_ignore_vertices.FindIsomorphism
          (
            test_graph_a,
            test_graph_c,
            std::min( test_graph_a.GetSize(), test_graph_c.GetSize())
          );

          BCL_ExampleIndirectCheck
          (
            unconnected_subgraph_ignore_vertices.GetIsomorphism().GetSize(), graph_size,
            "Finding the largest unconnected common substructure while ignoring vertex data"
          );

          graph::CommonSubgraphIsomorphism< size_t, size_t>
          unconnected_subgraph_ignore_vertices_and_edges
          (
            graph::CommonSubgraphIsomorphismBase::e_Unconnected,
            do_not_compare,
            do_not_compare
          );

          unconnected_subgraph_ignore_vertices_and_edges.FindIsomorphism
          (
            test_graph_a,
            test_graph_c,
            std::min( test_graph_a.GetSize(), test_graph_c.GetSize())
          );

          BCL_ExampleIndirectCheck
          (
            unconnected_subgraph_ignore_vertices_and_edges.GetIsomorphism().GetSize(), graph_size,
            "Finding the largest unconnected common substructure while ignoring vertex & edge data"
          );

          // now construct a larger test graph
          graph::ConstGraph< size_t, size_t> test_graph_large
          (
            graph::ConstGraph< size_t, size_t>::MakeRandomUndirectedGraph
            (
              graph_size + 4,
              edge_datum_ratios,
              vertex_datum_ratios
            )
          );

          graph::CommonSubgraphIsomorphism< size_t, size_t>
            unconnected_subgraph_larger_smaller( graph::CommonSubgraphIsomorphismBase::e_Unconnected);

          // find the isomorphism, giving the largest graph second
          unconnected_subgraph.FindIsomorphism
          (
            test_graph_a,
            test_graph_large,
            graph::CommonSubgraphIsomorphismBase::EstimateUpperBounds( test_graph_a, test_graph_large)
          );

          // find the isomorphism, giving the largest graph first
          unconnected_subgraph_larger_smaller.FindIsomorphism
          (
            test_graph_large,
            test_graph_a,
            graph::CommonSubgraphIsomorphismBase::EstimateUpperBounds( test_graph_large, test_graph_a)
          );

          // the resulting isomorphisms should be of the same size
          BCL_ExampleIndirectCheck
          (
            unconnected_subgraph.GetIsomorphism().GetSize(),
            unconnected_subgraph_larger_smaller.GetIsomorphism().GetSize(),
            "whether the order of the graphs changes the isomorphism size"
          );

          // determine whether the isorphisms are identical by constructing the inverse map and checking for
          // equality
          storage::Map< size_t, size_t> inverse_isomorphism_of_inverse_subgraph_isomorphism;
          for
          (
            storage::Map< size_t, size_t>::const_iterator
              iso( unconnected_subgraph_larger_smaller.GetIsomorphism().Begin()),
              iso_end( unconnected_subgraph_larger_smaller.GetIsomorphism().End());
            iso != iso_end;
            ++iso
          )
          {
            inverse_isomorphism_of_inverse_subgraph_isomorphism[ iso->second] = iso->first;
          }

          if
          (
            BCL_ExampleIndirectCheck
            (
              unconnected_subgraph.GetIsomorphism().GetSize(),
              unconnected_subgraph_larger_smaller.GetIsomorphism().GetSize(),
              "whether the order of the graphs changes the isomorphism size"
            )
          )
          {
            // only do the full check for equality if the sizes were equal, otherwise could walk off the end of the
            // smaller map
            BCL_ExampleIndirectCheck
            (
              std::equal
              (
                unconnected_subgraph.GetIsomorphism().Begin(),
                unconnected_subgraph.GetIsomorphism().End(),
                inverse_isomorphism_of_inverse_subgraph_isomorphism.Begin()
              ),
              true,
              "whether the order of the graphs changes the isomorphism size"
            );
          }
        }
      }

      // check that finding multiple isomorphisms works
      linal::Matrix< size_t> three_membered_ring( 3, 3, size_t( 1));
      three_membered_ring( 0, 0) = three_membered_ring( 1, 1) = three_membered_ring( 2, 2) = 0;
      storage::Vector< size_t> vertex_data( 3, size_t( 0));

      graph::ConstGraph< size_t, size_t> graph( vertex_data, three_membered_ring);

      // there should be 6 isomorphisms of  1-2-3
      //                                     \-/
      // onto itself:
      // 1,2,3
      // 1,3,2
      // 2,1,3
      // 2,3,1
      // 3,1,2
      // 3,2,1
      connected_subgraph.FindIsomorphism
      (
        graph, graph, 3, 3, storage::Vector< storage::Vector< size_t> >(), true
      );
      BCL_ExampleCheck( connected_subgraph.GetIsomorphisms().GetSize(), 6);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleGraphCommonSubgraphIsomorphism

  const ExampleClass::EnumType ExampleGraphCommonSubgraphIsomorphism::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphCommonSubgraphIsomorphism())
  );
} // namespace bcl
