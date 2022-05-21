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

#ifndef BCL_GRAPH_EDGE_COVER_RING_PERCEPTION_H_
#define BCL_GRAPH_EDGE_COVER_RING_PERCEPTION_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_connectivity.h"
#include "bcl_graph_path.h"
#include "bcl_graph_ring.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

    namespace
    {
      //! @brief count the number of comparisons required to reach an integer in a binary search on an integral array of size SIZE
      //! @param INTEGER the integer
      //! @param SIZE the size of the search
      //! @return number of binary comparisons required
      size_t BinarySearchDistance( const size_t &INTEGER, const size_t &SIZE)
      {
        size_t position( SIZE >> 1);
        size_t lower_bound( 0), upper_bound( SIZE - 1);
        size_t n_comparisons( 0);
        while( position != INTEGER)
        {
          ++n_comparisons;
          if( position > INTEGER)
          {
            upper_bound = position - 1;
          }
          else
          {
            lower_bound = position + 1;
          }
          position = ( lower_bound + upper_bound) >> 1;
        }
        return n_comparisons;
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EdgeCoverRingPerception
    //! @brief determines the set of rings for a given graph such that every cyclic edge is included in at least one ring
    //! @details Only works with undirected graphs currently, but could be modified to support directed graphs when there is need
    //!
    //! @see @link example_graph_edge_cover_ring_perception.cpp @endlink
    //! @author mendenjl
    //! @date 09/01/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EdgeCoverRingPerception :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! rings in the graph
      storage::List< Ring> m_Rings;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EdgeCoverRingPerception();

      //! @brief constructor from a ConstGraph
      //! @param GRAPH graph for ring detection
      template< typename t_GraphType>
        EdgeCoverRingPerception
      (
        const t_GraphType &GRAPH
      )
      {
        const size_t graph_size( GRAPH.GetSize());

        // there can be no rings in a graph with size < 3
        if( graph_size < 3)
        {
          return;
        }

        // make sure we didn't get a directed graph
        BCL_Assert( GRAPH.IsUndirected(), "EdgeCoverRingPerception only works on undirected graphs");

        // if the graph has more than one component, time and memory can be saved splitting it into
        // components first
        storage::List< storage::Vector< size_t> > components( Connectivity::GetComponents( GRAPH));
        if( components.GetSize() > 1)
        {
          for
          (
            storage::List< storage::Vector< size_t> >::const_iterator
              itr_component( components.Begin()), itr_component_end( components.End());
            itr_component != itr_component_end;
            ++itr_component
          )
          {
            const storage::Vector< size_t> &component( *itr_component);
            EdgeCoverRingPerception perception_component( GRAPH.GetSubgraph( component));
            // get the rings, translate them back to the old indices
            for
            (
              storage::List< Ring>::const_iterator
                itr( perception_component.GetRings().Begin()), itr_end( perception_component.GetRings().End());
              itr != itr_end;
              ++itr
            )
            {
              // create the translated vertices
              storage::Vector< size_t> vertices( itr->GetVertices());
              for
              (
                storage::Vector< size_t>::iterator itr_vertices( vertices.Begin()), itr_vertices_end( vertices.End());
                itr_vertices != itr_vertices_end;
                ++itr_vertices
              )
              {
                *itr_vertices = component( *itr_vertices);
              }
              // append the ring
              m_Rings.PushBack( Ring( graph_size, vertices));
            }
          }
          return;
        }

        storage::Vector< size_t> shortest_cycles_for_vertex( graph_size, graph_size + 1);
        storage::Vector< size_t> cyclic_vertices;
        cyclic_vertices.AllocateMemory( graph_size);
        storage::Vector< size_t> vertex_is_cyclic( graph_size, size_t( 1));

        // find the shortest cycle lengths for each vertex and edge
        {
          // visit vertices out-of-order, this algorithm is a degree of complexity faster on weakly connected graphs
          // like those from proteins than if they are visited in order, given that atoms are often arranged sequentially
          // while shuffling would be easiest, it has the side effect of making the state of the random
          // number generator dependent on the number of molecules read in, which makes reproducibility in unit tests
          // much more difficult. Therefore, we instead use the strategy of conducting a pseudo-random binary search
          // across the set of indices by sorting by depth along the binary search path
          storage::Vector< size_t> order;
          order.AllocateMemory( graph_size);
          if( graph_size > size_t( 1))
          {
            size_t max_search_len( 1);

            {
              size_t graph_size_depl( graph_size);
              while( graph_size_depl)
              {
                ++max_search_len;
                graph_size_depl >>= 1;
              }
            }

            storage::Vector< storage::List< size_t> > search_points( max_search_len);
            for( size_t i( 0); i < graph_size; ++i)
            {
              search_points( BinarySearchDistance( i, graph_size)).PushBack( i);
            }
            for( size_t i( 0); i < max_search_len; ++i)
            {
              order.InsertElements( order.End(), search_points( i).Begin(), search_points( i).End());
            }
          }
          for( size_t a( 0); a < graph_size; ++a)
          {
            const size_t vertex_a( order( a));
            shortest_cycles_for_vertex( vertex_a) = Connectivity::LengthOfSmallestCycleWithVertex( GRAPH, vertex_a, vertex_is_cyclic);
            if( shortest_cycles_for_vertex( vertex_a) > graph_size)
            {
              vertex_is_cyclic( vertex_a) = 0;
            }
            else
            {
              cyclic_vertices.PushBack( vertex_a);
            }
          }
          cyclic_vertices.Sort( std::less< size_t>());
        }

        // test whether there were any cyclic vertices, if not, no need to proceed since we just have chains
        if( cyclic_vertices.IsEmpty())
        {
          return;
        }

        // if fewer than N / Sqrt(2) of the vertices are cyclic, then it will be faster and require less memory to just
        // create a new graph with the cylical vertices
        if( double( graph_size) / math::Sqrt( 2.0) > cyclic_vertices.GetSize())
        {
          EdgeCoverRingPerception perception_cyclic_vertices( GRAPH.GetSubgraph( cyclic_vertices));
          // get the rings, translate them back to the old indices
          for
          (
            storage::List< Ring>::const_iterator
              itr( perception_cyclic_vertices.GetRings().Begin()), itr_end( perception_cyclic_vertices.GetRings().End());
            itr != itr_end;
            ++itr
          )
          {
            // create the translated vertices
            storage::Vector< size_t> vertices( itr->GetVertices());
            for
            (
              storage::Vector< size_t>::iterator itr_vertices( vertices.Begin()), itr_vertices_end( vertices.End());
              itr_vertices != itr_vertices_end;
              ++itr_vertices
            )
            {
              *itr_vertices = cyclic_vertices( *itr_vertices);
            }
            // append the ring
            m_Rings.PushBack( Ring( graph_size, vertices));
          }
          return;
        }

        linal::Matrix< size_t>   shortest_cycles_for_edge( graph_size, graph_size, graph_size + 1);
        for( size_t vertex_a( 0); vertex_a < graph_size; ++vertex_a)
        {
          if( shortest_cycles_for_vertex( vertex_a) > graph_size)
          {
            // skip vertices that are not part of rings
            continue;
          }
          shortest_cycles_for_edge( vertex_a, vertex_a) = shortest_cycles_for_vertex( vertex_a);

          if( GRAPH.GetNeighborIndices( vertex_a).GetSize() == size_t( 2))
          {
            // both edges have the same smallest cycle edge length
            const size_t neighbor_a( GRAPH.GetNeighborIndices( vertex_a)( 0));
            const size_t neighbor_b( GRAPH.GetNeighborIndices( vertex_a)( 1));
            const size_t cycle_length( shortest_cycles_for_vertex( vertex_a));
            shortest_cycles_for_edge( vertex_a, neighbor_a) = cycle_length;
            shortest_cycles_for_edge( vertex_a, neighbor_b) = cycle_length;
            shortest_cycles_for_edge( neighbor_a, vertex_a) = cycle_length;
            shortest_cycles_for_edge( neighbor_b, vertex_a) = cycle_length;
            continue;
          }

          for( size_t vertex_b( 0); vertex_b < vertex_a; ++vertex_b)
          {
            if( shortest_cycles_for_vertex( vertex_b) <= graph_size)
            {
              shortest_cycles_for_edge( vertex_a, vertex_b) =
                Connectivity::LengthOfSmallestCycleWithEdge( GRAPH, vertex_a, vertex_b);
              shortest_cycles_for_edge( vertex_b, vertex_a) = shortest_cycles_for_edge( vertex_a, vertex_b);
            }
          }
        }

        // for each vertex in the graph
        for( size_t source_vertex_number( 0); source_vertex_number < graph_size; ++source_vertex_number)
        {
          // check that the vertex is part of a cycle
          if( shortest_cycles_for_vertex( source_vertex_number) > graph_size)
          {
            continue;
          }

          // get the distances from the vertex to all other vertices
          util::ShPtr< storage::Vector< size_t> > sp_distances
          (
            Connectivity::DistancesToOtherVertices( GRAPH, source_vertex_number)
          );

          const storage::Vector< size_t> &distances( *sp_distances);

          // get the indices of all connected vertices (aka: the neighborhood)
          const typename t_GraphType::t_EdgeTargetsOfVertex &neighborhood
          (
            GRAPH.GetNeighborIndices( source_vertex_number)
          );

          // for each neighbor
          for
          (
            typename t_GraphType::t_EdgeTargetsOfVertex::const_iterator
              itr_neighbor( neighborhood.Begin()),
              itr_neighbor_end( neighborhood.End());
            itr_neighbor != itr_neighbor_end;
            ++itr_neighbor
          )
          {
            if( source_vertex_number < *itr_neighbor)
            {
              // for each neighbor with vertex >
              m_Rings.Append
              (
                FindCyclesOfMinimalLengthBetween
                (
                  GRAPH,
                  source_vertex_number,
                  *itr_neighbor,
                  distances,
                  shortest_cycles_for_edge
                )
              );
            }
          }
        }
      }

      //! clone the object
      EdgeCoverRingPerception *Clone() const
      {
        return new EdgeCoverRingPerception( *this);
      }

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief Get all rings in the graph
      //! @return list of rings (all closed paths)
      const storage::List< Ring> &GetRings() const
      {
        return m_Rings;
      }

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    /////////////
    // methods //
    /////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief find the shortest cycles involving a given edge in the graph
      //! @param GRAPH the graph that is being examined
      //! @param SOURCE the lowest index vertex that will be in the resulting cycles
      //! @param TARGET another vertex to include in the path
      //! @param DISTANCES a vector containing distances from SOURCE to all vertices
      //! @param SHORTEST_CYCLES a matrix with lengths of the shortest cycle involving each edge
      template< typename t_GraphType>
      storage::List< Ring> FindCyclesOfMinimalLengthBetween
      (
        const t_GraphType &GRAPH,
        const size_t &SOURCE,
        const size_t &TARGET,
        const storage::Vector< size_t> &DISTANCES,
        const linal::Matrix< size_t> &SHORTEST_CYCLES
      )
      {
          
        // This function finds the *unique* cycles of the graph of minimal length between SOURCE and TARGET that originate
        // at SOURCE.
        // For a cycle to be unique, there must be no other choice of SOURCE and TARGET that would yield the same ring,
        // since we will be iterating over all other vertices, and sorting the rings would be painfully slow; moreover,
        // it could cause potentially the paths to be explored N^2 times too many

        // To ensure the uniqueness of a ring, all rings will be ordered in the following way:
        // the first vertex will be lowest-index vertex on the ring which has an edge on the ring with >= the cyclic
        // distance between SOURCE and TARGET

        // record the size of the graph
        const size_t graph_size( GRAPH.GetSize());

        // find the length of the cycle
        const size_t cycle_size( SHORTEST_CYCLES( SOURCE, TARGET));

        storage::List< Path> new_paths;

        storage::List< Ring> new_rings;

        if( cycle_size > graph_size)
        {
          return new_rings;
        }

        // make an edge object that will be the starting point for all cycles that we create
        storage::List< Path> old_paths( 1, Path( graph_size, storage::Vector< size_t>::Create( SOURCE, TARGET), false));

        // since the path already contains the source and target, we should now be looking for vertices at distance 2 from
        // the source (unless the cycle size is 3, in which case we will only need to consider the last vertex to add)
        size_t distance_from_source( 2);

        // the apex is the turn-around point, where we start looking for vertices with a shorter distance from the
        // source vertex
        const size_t apex( cycle_size / 2);

        // expand the paths out until they are all of length apex, pruning bad paths along the way
        while( distance_from_source <= apex)
        {
          // try to expand every path in old_paths by adding a new vertex
          for
          (
            storage::List< Path>::const_iterator itr_old_paths( old_paths.Begin()), itr_old_paths_end( old_paths.End());
            itr_old_paths != itr_old_paths_end;
            ++itr_old_paths
          )
          {
            // get the last element from the old path, which we will try to extend
            const size_t end_vertex( itr_old_paths->LastElement());
            for
            (
              storage::Vector< size_t>::const_iterator
                itr_neighbors( GRAPH.GetNeighborIndices( end_vertex).Begin()),
                itr_neighbors_end( GRAPH.GetNeighborIndices( end_vertex).End());
              itr_neighbors != itr_neighbors_end;
              ++itr_neighbors
            )
            { // consider adding each neighboring vertex to form a new path

              // don't follow paths that would double back towards the vertex prematurely, in which case their distance
              // will not be == distance from source.
              // Since we know the cycle's length, we know that such paths does not form a cycle with SOURCE and TARGET on it
              if( DISTANCES( *itr_neighbors) != distance_from_source)
              {
                continue;
              }

              // make sure that the path does not retrace it's steps ( adding something it already contains)
              if( itr_old_paths->Contains( *itr_neighbors))
              {
                continue;
              }

              // if the edge between end_vertex and *itr_neighbors has the same cycle size, then both vertices must be
              // greater than SOURCE, otherwise calling this function with end_vertex and *itr_neighbors would produce
              // an equivalent ring
              if
              (
                ( *itr_neighbors > SOURCE && end_vertex > SOURCE)
                || SHORTEST_CYCLES( end_vertex, *itr_neighbors) != cycle_size
              )
              {
                // the new vertex can be added to this path to form a new one

                // make a vector with the vertices of the old path
                storage::Vector< size_t> new_vertices( itr_old_paths->GetVertices());

                // add the new vertex
                new_vertices.PushBack( *itr_neighbors);

                // make a new path from new_vertices and add it to the new paths list
                new_paths.PushBack( Path( graph_size, new_vertices, false));
              }
            }
          }

          // finished expanding all paths, move the new paths into the old paths list so they can be expanded further
          old_paths.InternalData().swap( new_paths.InternalData());

          // and reset the new paths list
          new_paths.Reset();

          // increase the distance
          ++distance_from_source;
        }

        // the apex was reached, reset the distance from source according to whether the graph was of even or odd size
        if( cycle_size & size_t( 1)) // odd cycles have two vertices at the apex distance
        {
          distance_from_source = apex;
        }
        else // even cycles only have the apex
        {
          distance_from_source = apex - 1;
        }

        // now return towards the source vertex, stopping when all paths need only one more vertex to join up
        // to the source vertex (this is because there is some special ordering necessary for adding the final vertex)
        while( distance_from_source > 1)
        {
          // try to expand every path in old_paths by adding a new vertex
          for
          (
            storage::List< Path>::const_iterator itr_old_paths( old_paths.Begin()), itr_old_paths_end( old_paths.End());
            itr_old_paths != itr_old_paths_end;
            ++itr_old_paths
          )
          {
            // get the last element from the old path, which we will try to extend
            const size_t end_vertex( itr_old_paths->LastElement());
            for
            (
              storage::Vector< size_t>::const_iterator
                itr_neighbors( GRAPH.GetNeighborIndices( end_vertex).Begin()),
                itr_neighbors_end( GRAPH.GetNeighborIndices( end_vertex).End());
              itr_neighbors != itr_neighbors_end;
              ++itr_neighbors
            ) // consider adding each neighbor of end_vertex to form a new path
            {
              // don't follow paths that would lead further away from the vertex
              // Because we are past the apex, all paths need to return towards SOURCE to form the shortest cycles
              if( DISTANCES( *itr_neighbors) != distance_from_source)
              {
                continue;
              }

              // make sure that the path does not retrace it's steps ( adding something it already contains)
              if( itr_old_paths->Contains( *itr_neighbors))
              {
                continue;
              }

              // if the edge between end_vertex and *itr_neighbors has the same cycle size, then both vertices must be
              // greater than SOURCE, otherwise calling this function with end_vertex and *itr_neighbors would produce
              // an equivalent ring
              if
              (
                ( *itr_neighbors > SOURCE && end_vertex > SOURCE)
                || SHORTEST_CYCLES( end_vertex, *itr_neighbors) != cycle_size
              )
              {
                // the new vertex can be added to this path to form a new one

                // make a vector with the vertices of the old path
                storage::Vector< size_t> new_vertices( itr_old_paths->GetVertices());

                // add the new vertex
                new_vertices.PushBack( *itr_neighbors);

                // make a new path from new_vertices and add it to the new paths list
                new_paths.PushBack( Path( graph_size, new_vertices, false));
              }
            }
          }

          // finished expanding all paths, move the new paths into the old paths list so they can be expanded further
          old_paths.InternalData().swap( new_paths.InternalData());

          // and reset the new paths list
          new_paths.Reset();

          // decrease the distance
          --distance_from_source;
        }

        // now distance from source == 1, walk through the list of neighbors, looking for anything with distance == 1 and
        // index > target
        for
        (
          storage::List< Path>::const_iterator itr_old_paths( old_paths.Begin()), itr_old_paths_end( old_paths.End());
          itr_old_paths != itr_old_paths_end;
          ++itr_old_paths
        )
        {
          // get the last element from the old path, which we will try to extend 1 more
          const size_t end_vertex( itr_old_paths->LastElement());
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( GRAPH.GetNeighborIndices( end_vertex).Begin()),
              itr_neighbors_end( GRAPH.GetNeighborIndices( end_vertex).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          ) // consider adding each neighbor of end_vertex to form a cycle path
          {
            // don't follow paths that would lead further away from the vertex
            // Because we are past the apex, all paths need to return towards SOURCE to form the shortest cycles
            if( DISTANCES( *itr_neighbors) != distance_from_source)
            {
              continue;
            }

            // make sure that the path does not retrace it's steps ( adding something it already contains)
            if( itr_old_paths->Contains( *itr_neighbors))
            {
              continue;
            }

            // if the edge between end_vertex and *itr_neighbors has the same cycle size, then both vertices must be
            // greater than SOURCE, otherwise calling this function with end_vertex and *itr_neighbors would produce
            // an equivalent ring
            if
            (
              ( *itr_neighbors > SOURCE && end_vertex > SOURCE)
              || SHORTEST_CYCLES( end_vertex, *itr_neighbors) != cycle_size
            )
            {
              // if the edge between SOURCE and *itr_neighbors has the same cycle size, then *itr_neighbors must be
              // greater than SOURCE, and additionally, greater than TARGET, which otherwise calling this function with
              // SOURCE and *itr_neighbors would yield an equivalent ring.
              if( *itr_neighbors > TARGET || SHORTEST_CYCLES( SOURCE, *itr_neighbors) != cycle_size)
              {
                storage::Vector< size_t> new_vertices( itr_old_paths->GetVertices());
                new_vertices.PushBack( *itr_neighbors);
                new_rings.PushBack( Ring( graph_size, new_vertices));
              }
            }
          }
        }

        return new_rings;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class EdgeCoverRingPerception

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_EDGE_COVER_RING_PERCEPTION_H_

