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

#ifndef BCL_GRAPH_CONNECTIVITY_H_
#define BCL_GRAPH_CONNECTIVITY_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_path.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <numeric>
#include <string>

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connectivity
    //! @brief Determines connectivity information from graphs (currently just ConstGraph)
    //!
    //! @see @link example_graph_common_subgraph_isomorphism.cpp @endlink
    //! @author mendenjl
    //! @date August 23, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Connectivity :
      public util::ObjectInterface
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief DistancesToOtherVertices
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertices to find distances starting from
      //! @return a shared pointer to the vector containing distances starting from VERTEX.
      //! @note any vertices that are unreachable from the vertex at VERTEX are skipped
      //! @note this works appropriately on both directed and undirected graphs
      template< typename t_Graph>
      static util::ShPtr< storage::Vector< size_t> > DistancesToOtherVertices
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX,
        const size_t &MAX_DISTANCE_OF_INTEREST = util::GetUndefinedSize_t()
      )
      {
        const size_t unseen_flag( util::GetUndefinedSize_t());
        const size_t size( GRAPH.GetSize());
        util::ShPtr< storage::Vector< size_t> > distances_shptr
        (
          new storage::Vector< size_t>( size, unseen_flag)
        );

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< size_t> &distances( *distances_shptr);

        // initialize seen_vertices with the first vertex
        storage::Vector< size_t> seen_vertices_queue( 1, VERTEX);
        // allocate enough memory for the queue to contain all vertices in the graph
        seen_vertices_queue.AllocateMemory( size);

        size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search
        size_t distance( 1);    // initial distance will be 0

        distances( VERTEX) = 0; // VERTEX is where the search starts, so its distance is 0

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( vertex_queue_position < seen_vertices_queue.GetSize() && seen_vertices_queue.GetSize() < size)
        {
          // loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
          // not already in the queue, then add them to the queue and record their distance.
          // Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
          for
          (
            const size_t last_vertex_at_distance( seen_vertices_queue.GetSize());
            vertex_queue_position < last_vertex_at_distance && seen_vertices_queue.GetSize() < size;
            vertex_queue_position++
          )
          {
            const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));
            for
            (
              size_t i( 0), number_seen( target_row.GetSize());
              i < number_seen;
              i++
            )
            {
              const size_t new_vertex( target_row( i));
              if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                seen_vertices_queue.PushBack( new_vertex);
                distances( new_vertex) = distance;

                if( seen_vertices_queue.GetSize() == size) // the entire graph has been reached
                {
                  break;
                }
              }
            }
          }

          distance++;
          if( distance > MAX_DISTANCE_OF_INTEREST)
          {
            break;
          }
        }
        return distances_shptr;
      } // DistancesToOtherVertices

      //! @brief DistancesToOtherVertices
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertices to find distances starting from
      //! @param VERTEX_TO_IGNORE index of the vertex which should not be followed
      //! @return a shared pointer to the vector containing distances starting from VERTEX.
      //! @note any vertices that are unreachable from the vertex at VERTEX are skipped
      //! @note this works appropriately on both directed and undirected graphs
      template< typename t_Graph>
      static util::ShPtr< storage::Vector< size_t> > DirectedDistancesToOtherVertices
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX,
        const size_t &VERTEX_TO_IGNORE
      )
      {
        const size_t unseen_flag( util::GetUndefinedSize_t());
        const size_t size( GRAPH.GetSize());
        util::ShPtr< storage::Vector< size_t> > distances_shptr
        (
          new storage::Vector< size_t>( size, unseen_flag)
        );

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< size_t> &distances( *distances_shptr);

        // initialize seen_vertices with the first vertex
        storage::Vector< size_t> seen_vertices_queue( 1, VERTEX);
        // allocate enough memory for the queue to contain all vertices in the graph
        seen_vertices_queue.AllocateMemory( size);

        size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search
        size_t distance( 1);    // initial distance will be 0

        distances( VERTEX) = 0; // VERTEX is where the search starts, so its distance is 0
        distances( VERTEX_TO_IGNORE) = 0;

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( vertex_queue_position < seen_vertices_queue.GetSize() && seen_vertices_queue.GetSize() < size)
        {
          // loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
          // not already in the queue, then add them to the queue and record their distance.
          // Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
          for
          (
            const size_t last_vertex_at_distance( seen_vertices_queue.GetSize());
            vertex_queue_position < last_vertex_at_distance && seen_vertices_queue.GetSize() < size;
            vertex_queue_position++
          )
          {
            const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));
            for
            (
              size_t i( 0), number_seen( target_row.GetSize());
              i < number_seen;
              i++
            )
            {
              const size_t new_vertex( target_row( i));
              if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                seen_vertices_queue.PushBack( new_vertex);
                distances( new_vertex) = distance;

                if( seen_vertices_queue.GetSize() == size) // the entire graph has been reached
                {
                  break;
                }
              }
            }
          }

          distance++;
        }
        distances( VERTEX_TO_IGNORE) = unseen_flag;
        return distances_shptr;
      } // DistancesToOtherVertices

      //! @brief BreadthFirstSearchDirectedEdge
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertices to find distances starting from
      //! @param VERTEX_TO_IGNORE index of the vertex which should not be followed
      //! @return a shared pointer to a vector containing breadth-first search order of vertices visited
      template< typename t_Graph>
      static util::ShPtr< storage::Vector< size_t> > BreadthFirstSearchDirectedEdge
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX,
        const size_t &VERTEX_TO_IGNORE
      )
      {
        const size_t unseen_flag( util::GetUndefinedSize_t());
        const size_t size( GRAPH.GetSize());

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< size_t> distances( size, unseen_flag);

        // initialize seen_vertices with the first vertex
        util::ShPtr< storage::Vector< size_t> > seen_vertices_queue_sp( new storage::Vector< size_t>( 1, VERTEX));
        storage::Vector< size_t> &seen_vertices_queue( *seen_vertices_queue_sp);
        // allocate enough memory for the queue to contain all vertices in the graph
        seen_vertices_queue.AllocateMemory( size);

        size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search
        size_t distance( 1);    // initial distance will be 0

        distances( VERTEX) = 0; // VERTEX is where the search starts, so its distance is 0
        distances( VERTEX_TO_IGNORE) = 0;

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( vertex_queue_position < seen_vertices_queue.GetSize() && seen_vertices_queue.GetSize() < size)
        {
          // loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
          // not already in the queue, then add them to the queue and record their distance.
          // Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
          for
          (
            const size_t last_vertex_at_distance( seen_vertices_queue.GetSize());
            vertex_queue_position < last_vertex_at_distance && seen_vertices_queue.GetSize() < size;
            vertex_queue_position++
          )
          {
            const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));
            for
            (
              size_t i( 0), number_seen( target_row.GetSize());
              i < number_seen;
              i++
            )
            {
              const size_t new_vertex( target_row( i));
              if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                seen_vertices_queue.PushBack( new_vertex);
                distances( new_vertex) = distance;

                if( seen_vertices_queue.GetSize() == size) // the entire graph has been reached
                {
                  break;
                }
              }
            }
          }

          distance++;
        }
        distances( VERTEX_TO_IGNORE) = unseen_flag;
        return seen_vertices_queue_sp;
      } // DistancesToOtherVertices

      //! @brief LengthOfSmallestCycleWithVertex finds the length of the shortest cycle beginning at VERTEX
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertex to find the girth from
      //! @param CAN_VISIT 0 if the associated vertex cannot be visited (normally vertices that are already known to be non-cyclical), non-zero if it can
      //! @return length of the shortest cycle beginning at VERTEX (may be undefined)
      //! @note any vertices that are unreachable from the vertex at VERTEX are skipped
      //! @note this works appropriately on both directed and undirected graphs
      template< typename t_Graph>
      static size_t LengthOfSmallestCycleWithVertex
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX,
        const storage::Vector< size_t> CAN_VISIT = storage::Vector< size_t>()
      )
      {
        size_t shortest_branch( std::numeric_limits< size_t>::max());

        // if there are only 0 or 1 vertices, then there can be no cycle, so return immediately
        if( !GRAPH.IsDirected() && GRAPH.GetNeighborIndices( VERTEX).GetSize() <= 1)
        {
          return shortest_branch;
        }

        const size_t unseen_flag( util::GetUndefinedSize_t());
        const size_t size( GRAPH.GetSize());

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< size_t> distances( size, unseen_flag);

        // initialize seen_vertices with the first vertex
        storage::Vector< size_t> seen_vertices_queue( 1, VERTEX);

        // store the first branch # (index of the vertex connected to VERTEX which first reached this particular vertex
        // If multiple branches reached a particular vertex simultaneously, then they are pushed back into the list
        storage::Vector< size_t> branch_number( size, unseen_flag);

        // allocate enough memory for the queue to contain all vertices in the graph
        seen_vertices_queue.AllocateMemory( size);

        size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search

        distances( VERTEX) = 0; // VERTEX is where the search starts, so its distance is 0
        branch_number( VERTEX) = VERTEX;

        bool cycle_found( false);

        size_t distance( 1);    // initial distance will be 0

        storage::Vector< size_t> can_visit;
        if( CAN_VISIT.GetSize() != size)
        {
          can_visit = storage::Vector< size_t>( size, size_t( 1));
        }
        const storage::Vector< size_t> &can_visit_ref( CAN_VISIT.GetSize() != size ? can_visit : CAN_VISIT);

        {
          // assign the 1st-degree neighbors a branch number = their index
          const typename t_Graph::t_EdgeTargetsOfVertex &row( GRAPH.GetNeighborIndices( VERTEX));
          for( size_t i( 0), number_seen( row.GetSize()); i < number_seen; i++)
          {
            if( can_visit_ref( row( i)))
            {
              distances( row( i)) = distance;
              branch_number( row( i)) = row( i);
              seen_vertices_queue.PushBack( row( i));
            }
          }
          ++vertex_queue_position;
          ++distance;
        }

        // if there are only 0 or 1 visitable vertices, then there can be no cycle, so return immediately
        if( !GRAPH.IsDirected() && seen_vertices_queue.GetSize() <= 1)
        {
          return shortest_branch;
        }

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( vertex_queue_position < seen_vertices_queue.GetSize() && !cycle_found)
        {
          // loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
          // not already in the queue, then add them to the queue and record their distance.
          // Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
          for
          (
            const size_t last_vertex_at_distance( seen_vertices_queue.GetSize());
            vertex_queue_position < last_vertex_at_distance;
            vertex_queue_position++
          )
          {
            const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));

            for( size_t i( 0), number_seen( target_row.GetSize()); i < number_seen; i++)
            {
              const size_t new_vertex( target_row( i));
              if( !can_visit_ref( new_vertex))
              {
                // do nothing; can't visit the associated vertex
              }
              else if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                seen_vertices_queue.PushBack( new_vertex);
                distances( new_vertex) = distance;
                branch_number( new_vertex) = branch_number( current_vertex);
              }
              else if // check for cycles
              (
                new_vertex != VERTEX
                && branch_number( current_vertex) != branch_number( new_vertex)
              )
              {
                cycle_found = true;
                shortest_branch = std::min( shortest_branch, distance + distances( new_vertex));
              }
            }
          }

          distance++;
        }

        return shortest_branch;
      } // LengthOfSmallestCycleWithVertex

      //! @brief IsInThreeMemberedCycle fast check for being in a 3-membered cycle
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertex to find the girth from
      //! @return true if VERTEX is in a 3-membered, undirected cycle
      template< typename t_Graph>
      static bool IsInThreeMemberedCycle
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX
      )
      {
        // if there are only 0 or 1 vertices, then there can be no cycle, so return immediately
        const typename t_Graph::t_EdgeTargetsOfVertex &neighbors( GRAPH.GetNeighborIndices( VERTEX));
        const size_t number_neighbors( neighbors.GetSize());
        if( !number_neighbors || ( !GRAPH.IsDirected() && number_neighbors <= 1))
        {
          return false;
        }

        if( GRAPH.IsDirected())
        {
          // directed graph; fall back to slow algorithm
          return LengthOfSmallestCycleWithVertex( GRAPH, VERTEX) == size_t( 3);
        }

        // make an incidence vector using a string with 0 to denote non-connections, 1 to denote connections
        const char unseen_flag( '\0');
        std::string incidence( GRAPH.GetSize(), unseen_flag);
        for( size_t i( 0); i < number_neighbors; ++i)
        {
          incidence[ neighbors( i)] = '1';
        }

        // ignore connections to self
        incidence[ VERTEX] = '\0';

        // walk over connections
        for( size_t i( 0); i < number_neighbors; ++i)
        {
          const size_t neighbor( neighbors( i));
          const typename t_Graph::t_EdgeTargetsOfVertex &neighbors_neighbors( GRAPH.GetNeighborIndices( neighbor));
          for
          (
            typename t_Graph::t_EdgeTargetsOfVertex::const_iterator
              itr( neighbors_neighbors.Begin()), itr_end( neighbors_neighbors.End());
            itr != itr_end;
            ++itr
          )
          {
            if( incidence[ *itr] == '1' && *itr != neighbor)
            {
              return true;
            }
          }
        }

        return false;
      } // IsInThreeMemberedCycle

      //! @brief LengthOfSmallestCycleWithEdge finds the length of the shortest cycle containing the edge
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX_A index of either vertex for the edge
      //! @param VERTEX_B index of the other vertex for the edge
      //! @param CAN_VISIT 0 if the associated vertex cannot be visited (normally vertices that are already known to be non-cyclical), non-zero if it can
      //! @return length of the shortest cycle traveling the edge between VERTEX_A and VERTEX_B
      //! @note any vertices that are unreachable from the vertex at VERTEX are skipped
      //! @note this works appropriately on both directed and undirected graphs
      template< typename t_Graph>
      static size_t LengthOfSmallestCycleWithEdge
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX_A,
        const size_t &VERTEX_B,
        const storage::Vector< size_t> CAN_VISIT = storage::Vector< size_t>()
      )
      {
        size_t shortest_branch( std::numeric_limits< size_t>::max());

        // if there are only 0 or 1 vertices, then there can be no cycle, so return immediately
        if
        (
          !GRAPH.IsDirected()
          &&
          ( GRAPH.GetNeighborIndices( VERTEX_A).GetSize() <= 1 || GRAPH.GetNeighborIndices( VERTEX_B).GetSize() <= 1)
        )
        {
          return shortest_branch;
        }

        // make sure the vertices are really connected
        if( !GRAPH.AreConnected( VERTEX_A, VERTEX_B))
        {
          return shortest_branch;
        }

        const size_t unseen_flag( util::GetUndefinedSize_t());
        const size_t size( GRAPH.GetSize());

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< size_t> distances( size, unseen_flag);

        // initialize seen_vertices with the first vertex
        storage::Vector< size_t> seen_vertices_queue( 1, VERTEX_A);

        // store the first branch # (index of the vertex connected to VERTEX_A which first reached this particular vertex
        // If multiple branches reached a particular vertex simultaneously, then they are pushed back into the list
        storage::Vector< size_t> branch_number( size, unseen_flag);

        // allocate enough memory for the queue to contain all vertices in the graph
        seen_vertices_queue.AllocateMemory( size);

        size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search

        distances( VERTEX_A) = 0; // VERTEX is where the search starts, so its distance is 0
        branch_number( VERTEX_A) = VERTEX_A;

        bool cycle_found( false);

        storage::Vector< size_t> can_visit;
        if( CAN_VISIT.GetSize() != size)
        {
          can_visit = storage::Vector< size_t>( size, size_t( 1));
        }
        const storage::Vector< size_t> &can_visit_ref( CAN_VISIT.GetSize() != size ? can_visit : CAN_VISIT);

        size_t distance( 1);    // initial distance will be 0
        {
          // assign the 1st-degree neighbors a branch number = their index
          const typename t_Graph::t_EdgeTargetsOfVertex &row( GRAPH.GetNeighborIndices( VERTEX_A));
          for( size_t i( 0), number_seen( row.GetSize()); i < number_seen; i++)
          {
            if( can_visit_ref( row( i)))
            {
              distances( row( i)) = distance;
              branch_number( row( i)) = row( i);
              seen_vertices_queue.PushBack( row( i));
            }
          }
          ++vertex_queue_position;
          ++distance;
        }
        if( seen_vertices_queue.GetSize() == size_t( 1))
        {
          return shortest_branch;
        }

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( vertex_queue_position < seen_vertices_queue.GetSize() && !cycle_found)
        {
          // loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
          // not already in the queue, then add them to the queue and record their distance.
          // Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
          for
          (
            const size_t last_vertex_at_distance( seen_vertices_queue.GetSize());
            vertex_queue_position < last_vertex_at_distance;
            vertex_queue_position++
          )
          {
            const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));

            for( size_t i( 0), number_seen( target_row.GetSize()); i < number_seen; i++)
            {
              const size_t new_vertex( target_row( i));
              if( !can_visit_ref( new_vertex))
              {
              }
              else if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                seen_vertices_queue.PushBack( new_vertex);
                distances( new_vertex) = distance;
                branch_number( new_vertex) = branch_number( current_vertex);
              }
              else if // check for cycles that include the target edge
              (
                new_vertex != VERTEX_A
                && branch_number( current_vertex) != branch_number( new_vertex)
                &&
                (
                  branch_number( current_vertex) == VERTEX_B
                  || branch_number( new_vertex) == VERTEX_B
                )
              )
              {
                cycle_found = true;
                shortest_branch = std::min( shortest_branch, distance + distances( new_vertex));
              }
            }
          }

          distance++;
        }

        return shortest_branch;
      } // LengthOfSmallestCycleWithEdge

      //! @brief GetVerticesReachableFrom
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertices to find the connected component starting from
      //! @return a shared pointer to the vector containing the indices that are reachable from VERTEX.
      template< typename t_Graph>
      static util::ShPtr< storage::Vector< size_t> > GetVerticesReachableFrom
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX
      )
      {
        util::ShPtr< storage::Vector< size_t> > distances( DistancesToOtherVertices( GRAPH, VERTEX));
        const size_t unseen_flag( util::GetUndefinedSize_t());

        // the distances vector will have the undefined value on all unreachable vertices
        // The connected component is the indices of the defined positions in the vector, so
        // replace the defined distances with their indices
        const size_t size( distances->GetSize());
        size_t       number_unreachable( 0);
        for( size_t i( 0); i < size; i++)
        {
          if( ( *distances)( i) != unseen_flag)
          {
            ( *distances)( i) = i;
          }
          else
          {
            number_unreachable++;
          }
        }

        std::remove( distances->Begin(), distances->End(), unseen_flag); // and remove the unseen-flags from the vector
        distances->Resize( GRAPH.GetSize() - number_unreachable);
        return distances;
      }

      //! @brief GetVerticesReachableFrom
      //! @param GRAPH the graph in which the vertex resides
      //! @param VERTEX index of the vertices to find the connected component starting from
      //! @param VERTEX_TO_IGNORE index of the vertex which should not be followed
      //! @return a shared pointer to the vector containing the indices that are reachable from VERTEX.
      template< typename t_Graph>
      static util::ShPtr< storage::Vector< size_t> > GetVerticesReachableFromDirectedEdge
      (
        const t_Graph &GRAPH,
        const size_t &VERTEX,
        const size_t &VERTEX_TO_IGNORE
      )
      {
        util::ShPtr< storage::Vector< size_t> > distances( DirectedDistancesToOtherVertices( GRAPH, VERTEX, VERTEX_TO_IGNORE));
        const size_t unseen_flag( util::GetUndefinedSize_t());

        // the distances vector will have the undefined value on all unreachable vertices
        // The connected component is the indices of the defined positions in the vector, so
        // replace the defined distances with their indices
        const size_t size( distances->GetSize());
        size_t       number_unreachable( 0);
        for( size_t i( 0); i < size; i++)
        {
          if( ( *distances)( i) != unseen_flag)
          {
            ( *distances)( i) = i;
          }
          else
          {
            number_unreachable++;
          }
        }

        std::remove( distances->Begin(), distances->End(), unseen_flag); // and remove the unseen-flags from the vector
        distances->Resize( GRAPH.GetSize() - number_unreachable);
        return distances;
      }

      //! @brief Get the connected components of the graph
      //! @param GRAPH the graph for which to find components
      //! @return a list of indices of the connected components of the graph
      template< typename t_Graph>
      static storage::List< storage::Vector< size_t> > GetComponents( const t_Graph &GRAPH)
      {
        BCL_Assert
        (
          !GRAPH.IsDirected() || !GRAPH.GetSize(),
          "GetComponents is currently only supported on undirected graphs"
        );

        std::vector< bool> have_seen( GRAPH.GetSize(), false);
        storage::List< storage::Vector< size_t> > components;

        for( size_t i( 0), size( GRAPH.GetSize()); i < size; i++)
        {
          if( !have_seen[ i])
          {
            components.PushBack( *GetVerticesReachableFrom( GRAPH, i));

            for
            (
              storage::Vector< size_t>::const_iterator
                itr( components.LastElement().Begin()),
                itr_end( components.LastElement().End());
              itr != itr_end;
              ++itr
            )
            {
              have_seen[ *itr] = true;
            }
          }
        }

        return components;
      }

      //! @brief Get the connected components of the graph as new graphs
      //! @param GRAPH the graph whose components to find
      //! @return a list of graphs containing the connected components from this graph
      template< typename t_Graph>
      static storage::List< t_Graph> GetComponentsAsNewGraphs( const t_Graph &GRAPH)
      {
        storage::List< storage::Vector< size_t> > component_indices( GetComponents( GRAPH));
        storage::List< t_Graph> components;
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator
            itr( component_indices.Begin()),
            itr_end( component_indices.End());
          itr != itr_end;
          ++itr
        )
        {
          components.PushBack( GRAPH.GetSubgraph( *itr));
        }

        return components;
      }

      //! @brief GetVerticesReachableFrom
      //! @param GRAPH the graph whose components to find
      //! @return true if all vertices of the graph are mutually reachable ( e.g. strong connectivity)
      template< typename t_Graph>
      static bool IsConnected( const t_Graph &GRAPH)
      {
        const size_t size( GRAPH.GetSize());
        if( size == 0)
        {
          return true;
        }
        else if( GetVerticesReachableFrom( GRAPH, 0)->GetSize() == size) // The first vertex can reach all other vertices
        {
          if( GRAPH.IsDirected())
          {
            // The harder, directed, case.  The fastest way O(N^2) is:
            // 1. construct the transpose graph and
            // 2. test the transpose graph and the original graph with GetVerticesReachableFrom(0)on vertex 0 as well
            //    but this O(N^3) algorithm should suffice for now
            for( size_t i( 1); i < size; i++)
            {
              if( GetVerticesReachableFrom( GRAPH, i)->GetSize() != size)
              {
                return false;
              }
            }
          }

          return true;
        }

        return false;
      }

      //! @brief Find a minimal path between two vertices, using edges as weights
      //! @param GRAPH the graph for which to find the minimal path
      //! @param VERTEX_A the index of the first vertex to find a path between
      //! @param VERTEX_B the index of the second vertex to find a path between
      //! @return the shortest (or minimal) path between A and B, provided that one exists
      template< typename t_Graph>
      static Path FindMinimalPath( const t_Graph &GRAPH, const size_t &VERTEX_A, const size_t &VERTEX_B)
      {
        const double unseen_flag( std::numeric_limits< double>::max());
        const size_t size( GRAPH.GetSize());

        BCL_Assert
        (
          VERTEX_A < size,
          "Vertex A " + util::Format()( VERTEX_A) + " was beyond the graph size: " + util::Format()( size)
        );
        BCL_Assert
        (
          VERTEX_B < size,
          "Vertex B " + util::Format()( VERTEX_B) + " was beyond the graph size: " + util::Format()( size)
        );

        // maintain a list of shortest paths to each vertex
        storage::Vector< storage::Vector< size_t> > paths( size);

        // distances will hold the distance of each vertex to the current vertex
        storage::Vector< double> distances( size, unseen_flag);

        distances( VERTEX_A) = 0; // VERTEX is where the search starts, so its distance is 0
        paths( VERTEX_A).PushBack( VERTEX_A);

        // Construct a list of all vertices that still have neighbors that are not yet in the minimum spanning tree
        storage::List< size_t> available_vertices;
        available_vertices.PushBack( VERTEX_A);

        // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
        // unless all vertices are put into the graph.
        while( !available_vertices.IsEmpty() && distances( VERTEX_B) == unseen_flag)
        {
          size_t best_vertex_start( 0), best_vertex_terminal( 0);
          double best_distance( std::numeric_limits< double>::max());
          // walk through all available vertices
          for
          (
            typename storage::List< size_t>::iterator
              itr_available( available_vertices.Begin()), itr_available_end( available_vertices.End());
            itr_available != itr_available_end;
            // iteration in loop
          )
          {
            // current vertex
            const size_t current_vertex( *itr_available);

            // current distance
            const double current_distance( distances( current_vertex));

            // keep track of whether there were any available nodes still for this vertex
            bool found_available_node( false);

            // scan all connections of the vertex
            // walk over remaining edges
            // target row is a reference to the edges reachable from the current vertex
            const typename t_Graph::t_EdgeTargetsOfVertex &target_row( GRAPH.GetNeighborIndices( current_vertex));
            const typename t_Graph::t_EdgeDataOfVertex &target_row_data( GRAPH.GetNeighborData( current_vertex));
            for
            (
              size_t i( 0), number_seen( target_row.GetSize());
              i < number_seen;
              ++i
            )
            {
              const size_t new_vertex( target_row( i));
              if( distances( new_vertex) == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
              {
                found_available_node = true;
                // compute the new distance to this node
                const double new_distance( current_distance + target_row_data( i));
                if( new_distance < best_distance)
                {
                  best_vertex_start = current_vertex;
                  best_vertex_terminal = new_vertex;
                  best_distance = new_distance;
                }
              }
            }
            if( !found_available_node)
            {
              // no more edges that lead to other vertices from this node, remove it
              storage::List< size_t>::iterator itr_old( itr_available);
              ++itr_available;
              available_vertices.Remove( itr_old);
            }
            else
            {
              // continue to the next node
              ++itr_available;
            }
          }
          if( best_distance < unseen_flag)
          {
            distances( best_vertex_terminal) = best_distance;
            paths( best_vertex_terminal) = paths( best_vertex_start);
            paths( best_vertex_terminal).PushBack( best_vertex_terminal);
            available_vertices.PushBack( best_vertex_terminal);
          }
        }
        return Path( size, paths( VERTEX_B), GRAPH.IsUndirected());
      }

    }; // class GraphConnectivity

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_CONNECTIVITY_H_

