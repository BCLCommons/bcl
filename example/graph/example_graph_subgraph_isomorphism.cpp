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
#include "graph/bcl_graph_subgraph_isomorphism.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_subgraph_isomorphism.cpp
  //!
  //! @author mendenjl
  //! @date Jun 07, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphSubgraphIsomorphism :
    public ExampleInterface
  {
  public:

    //! @brief Clone
    ExampleGraphSubgraphIsomorphism *Clone() const
    {
      return new ExampleGraphSubgraphIsomorphism( *this);
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
      const size_t number_of_test_graphs = 25;
      const size_t number_of_edge_data   = 4; // Edge data are assigned randomly from 1 to number_of_edge_data
      const size_t number_of_vertex_data = 4;
      const size_t min_graph_size        = 50; // these graphs are at the upper limit of what the normal csi can handle
      const size_t max_graph_size        = 60;
      const double average_connections   = 5;

      BCL_Assert( average_connections < min_graph_size, " too many connections for some graphs");
      BCL_Assert( average_connections >= 2.5, " too few connections -> high probability of unconnected graph");
      BCL_Assert( number_of_edge_data > 1, " not enough edge data to indicate connected/unconnected!");
      BCL_Assert( number_of_vertex_data > 0, " must allow at least one vertex datum!");

      graph::SubgraphIsomorphism< size_t, size_t> connected_subgraph;

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

        // test the connected version of the algorithm
        util::Stopwatch connected_timer( false);
        connected_timer.Reset();
        connected_timer.Start();

        connected_subgraph.SetSubgraphExternalOwnership( test_graph_a);
        connected_subgraph.SetGraphExternalOwnership( test_graph_b);
        BCL_ExampleCheck( connected_subgraph.FindIsomorphism(), true);

        BCL_MessageStd
        (
          "Algorithm time: " + util::Format()( connected_timer.GetProcessDuration())
        );

        BCL_MessageStd
        (
          "    isomorphism size: " + util::Format()( connected_subgraph.GetIsomorphism().GetSize())
        );

        BCL_ExampleIndirectCheck
        (
          connected_subgraph.GetIsomorphism().GetSize(), graph_size,
          "Finding fragment isomorphism on a graph"
        );

        connected_subgraph.SetSubgraphExternalOwnership( test_graph_b);
        connected_subgraph.SetGraphExternalOwnership( test_graph_a);
        connected_subgraph.FindIsomorphism();
        BCL_ExampleIndirectCheck
        (
          connected_subgraph.GetIsomorphism().GetSize(),
          graph_size,
          "order independence of subgraph isomorphism"
        );
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
      connected_subgraph.SetGraphExternalOwnership( graph);
      connected_subgraph.SetSubgraphExternalOwnership( graph);
      connected_subgraph.FindAllIsomorphisms();
      BCL_ExampleCheck( connected_subgraph.GetIsomorphisms().GetSize(), 6);
      connected_subgraph.FindDisparateIsomorphisms();
      BCL_ExampleCheck( connected_subgraph.GetIsomorphisms().GetSize(), 1);

      io::IFStream input;

      util::Stopwatch timer_si( false);
      util::Stopwatch timer_csi( false);
      timer_si.Reset();
      timer_csi.Reset();
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "113mer.sdf"));
      chemistry::FragmentEnsemble big_chain( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);
      const chemistry::FragmentComplete &mol( big_chain.GetMolecules().FirstElement());
      chemistry::AtomVector< chemistry::AtomComplete> atoms_reordered( mol.GetAtomInfo(), mol.GetBondInfo());
      const size_t size( mol.GetNumberAtoms());
      storage::Vector< size_t> new_ordering( size);
      for( size_t i( 0); i < size; ++i)
      {
        new_ordering( i) = i;
      }
      // randomly reorder the vector
      new_ordering.Shuffle();
      // reorder the atom vector accordingly
      atoms_reordered.Reorder( new_ordering);

      // create a new fragment complete with the reordered vector
      const chemistry::FragmentComplete reordered_mol( atoms_reordered, "test");

      // create graphs
      graph::ConstGraph< size_t, size_t> graph_a( chemistry::ConformationGraphConverter()( mol));
      graph::ConstGraph< size_t, size_t> graph_b( chemistry::ConformationGraphConverter()( reordered_mol));

      // find the isomorphism using SubgraphIsomorphism
      timer_si.Start();
      graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
      isomorphism.SetGraphExternalOwnership( graph_a);
      isomorphism.SetSubgraphExternalOwnership( graph_b);
      BCL_ExampleIndirectCheck( isomorphism.FindIsomorphism(), true, "subgraph isomorphism works");
      timer_si.Stop();

      // the tolerance is necessary here because the molecule has one plane of symmetry, and fortunately
      // this only changes the position of two atoms, which are consecutive in the file
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        isomorphism.GetIsomorphism(),
        new_ordering,
        size_t( 2),
        "subgraph isomorphism found the actual isomorphism"
      );

      timer_csi.Start();
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_isomorphism;
      csi_isomorphism.FindIsomorphism( graph_a, graph_b, graph_a.GetSize(), graph_a.GetSize());
      BCL_ExampleIndirectCheck( csi_isomorphism.GetIsomorphism().GetSize(), size, "csi isomorphism works");
      timer_csi.Stop();
      BCL_MessageStd
      (
        "Finding isomorphisms of 113mer randomly reordered "
        " with SubgraphIsomorphism took " + timer_si.GetTotalTime().GetTimeAsHourMinuteSecondMilliSeconds()
      );
      BCL_MessageStd
      (
        "Finding isomorphisms of 113mer molecules randomly reordered "
        " with CommonSubgraphIsomorphism took " + timer_csi.GetTotalTime().GetTimeAsHourMinuteSecondMilliSeconds()
      );
      BCL_ExampleIndirectCheck
      (
        timer_si.GetTotalTime() < timer_csi.GetTotalTime(),
        true,
        "Subgraph isomorphism should be faster than finding largest common substructure!"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleGraphSubgraphIsomorphism

  const ExampleClass::EnumType ExampleGraphSubgraphIsomorphism::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphSubgraphIsomorphism())
  );
} // namespace bcl
