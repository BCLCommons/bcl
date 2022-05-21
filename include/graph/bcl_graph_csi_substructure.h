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

#ifndef BCL_GRAPH_CSI_SUBSTRUCTURE_H_
#define BCL_GRAPH_CSI_SUBSTRUCTURE_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_common_subgraph_isomorphism.h"
#include "bcl_graph_connectivity.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically
#include <utility>

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CSISubstructure
    //! @brief Class which houses static routines related to substructure and isomorphism search
    //!
    //! @see @link example_graph_csi_substructure.cpp @endlink
    //! @author mendenjl
    //! @date Mar 21, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CSISubstructure
    {

    public:

      //! @brief get the vertices from a subgraph that have isomorphic edges and vertex types to those in another graph
      //! @param SUBGRAPH The fragment to search for
      //! @param GRAPH The graph to search for the subgraph within
      //! @return vector of matching vertices from GRAPH for each vertex in SUBGRAPH
      template< typename t_VertexType, typename t_EdgeType>
      static storage::Vector< storage::Vector< size_t> > GetVertexMatchingMatrixForSubgraph
      (
        const ConstGraph< t_VertexType, t_EdgeType> &SUBGRAPH,
        const ConstGraph< t_VertexType, t_EdgeType> &GRAPH
      )
      {
        static util::Stopwatch s_vmm_creation( "VMM creation", util::Time( 1, 0), util::Message::e_Verbose, true, false);
        s_vmm_creation.Start();
        // store the sizes of the graphs for easy access
        const size_t subgraph_size( SUBGRAPH.GetSize());
        const size_t graph_size( GRAPH.GetSize());

        // setup storage for the matching vertices
        storage::Vector< storage::Vector< size_t> > vertex_matching_matrix( subgraph_size);

        // if graphs are equal sized, employ additional optimizations
        if( subgraph_size == graph_size)
        {
          storage::List< storage::VectorND< 2, storage::Vector< size_t> > > initial_vmm;

          // create a vector with all the indices
          storage::Vector< size_t> index_vector( graph_size);
          for( size_t i( 0); i < graph_size; ++i)
          {
            index_vector( i) = i;
          }

          storage::VectorND< 2, storage::Vector< size_t> > initial_group( index_vector, index_vector);

          if
          (
            !RefineMatchingVerticesForIsomorphism
            (
              SUBGRAPH,
              GRAPH,
              initial_group,
              initial_vmm,
              &GetVertexType< t_VertexType, t_EdgeType>
            )
          )
          {
            s_vmm_creation.Stop();
            return vertex_matching_matrix;
          }

          storage::List< storage::VectorND< 2, storage::Vector< size_t> > > vmm_post_environment;

          // track the number of possible permutations given the current group size
          uint64_t permutations( 1);

          // based on permutations, determine whether to filter by distance as well
          bool filter_by_distance( false);
          for
          (
            storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::const_iterator
              itr( initial_vmm.Begin()), itr_end( initial_vmm.End());
            itr != itr_end;
            ++itr
          )
          {
            const size_t old_size( vmm_post_environment.GetSize());
            if( itr->First().GetSize() == size_t( 1))
            {
              vmm_post_environment.PushBack( *itr);
              continue;
            }
            if
            (
              !RefineMatchingVerticesForIsomorphism
              (
                SUBGRAPH,
                GRAPH,
                *itr,
                vmm_post_environment,
                &GetVertexEnvironment< t_VertexType, t_EdgeType>
              )
            )
            {
              s_vmm_creation.Stop();
              return vertex_matching_matrix;
            }
            if( filter_by_distance)
            {
              // if filtering by distance, there is no further need to track permutations
              continue;
            }

            const uint64_t group_size( vmm_post_environment.GetSize() - old_size);
            // get the size of this group factorial
            const uint64_t permutations_of_group( math::Factorial( group_size));
            if( permutations_of_group > graph_size)
            {
              filter_by_distance = true;
              continue;
            }
            permutations *= permutations_of_group;
            if( permutations > graph_size)
            {
              filter_by_distance = true;
            }
          }

          storage::List< storage::VectorND< 2, storage::Vector< size_t> > > vmm_final;
          if( filter_by_distance)
          {
            // Filter by distance to further reduce the number of options to consider
            for
            (
              storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::const_iterator
                itr( vmm_post_environment.Begin()), itr_end( vmm_post_environment.End());
              itr != itr_end;
              ++itr
            )
            {
              if( itr->First().GetSize() == size_t( 1))
              {
                vmm_final.PushBack( *itr);
                continue;
              }
              if
              (
                !RefineMatchingVerticesForIsomorphism
                (
                  SUBGRAPH,
                  GRAPH,
                  *itr,
                  vmm_final,
                  &GetVertexDistanceCounts< t_VertexType, t_EdgeType>
                )
              )
              {
                s_vmm_creation.Stop();
                return vertex_matching_matrix;
              }
            }
          }
          else
          {
            // there are very few possibilities, so it would likely be quicker to check them individually
            // swap the internal data so that vmm_final contains the right groups
            vmm_final.InternalData().swap( vmm_post_environment.InternalData());
          }

          for
          (
            storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::const_iterator
              itr( vmm_final.Begin()), itr_end( vmm_final.End());
            itr != itr_end;
            ++itr
          )
          {
            for
            (
              storage::Vector< size_t>::const_iterator
                itr_group_a( itr->First().Begin()), itr_group_a_end( itr->First().End());
              itr_group_a != itr_group_a_end;
              ++itr_group_a
            )
            {
              vertex_matching_matrix( *itr_group_a) = itr->Second();
            }
          }

          s_vmm_creation.Stop();
          return vertex_matching_matrix;
        }
        else
        {
          // make a vector for each graph containing maps whose keys contain bond-type, atom-type pairs
          // and whose values are the count of that in the graph and subgraph

          storage::Vector< storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t> >
            vertex_data_graph( GetVertexEnvironments( GRAPH));

          // make a vector for each graph containing maps whose keys contain bond-type, atom-type pairs
          // and whose values are the count of that in the graph and subgraph
          storage::Vector< storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t> >
            vertex_data_subgraph( GetVertexEnvironments( SUBGRAPH));

          for( size_t index_in_subgraph( 0); index_in_subgraph < subgraph_size; ++index_in_subgraph)
          {
            // get the vertex of interest
            const size_t &subgraph_datum( SUBGRAPH.GetVertexData( index_in_subgraph));

            // get the environment around this vertex
            const storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t>
              &vertex_datum_subgraph( vertex_data_subgraph( index_in_subgraph));

            // get the row of the vertex matching vectors that this data will belong to
            storage::Vector< size_t> &vmv_row( vertex_matching_matrix( index_in_subgraph));

            // allocate sufficient memory to hold all the vertices if necessary
            vmv_row.AllocateMemory( graph_size);
            for( size_t index_in_graph( 0); index_in_graph < graph_size; ++index_in_graph)
            {
              // add the vertices in GRAPH that are comparable to this vertex in SUBGRAPH

              // ignore vertices with values != the current vertex in the subgraph
              if( subgraph_datum != GRAPH.GetVertexData( index_in_graph))
              {
                continue;
              }

              // ignore vertices in GRAPH that have fewer connections than the vertex from the subgraph
              if
              (
                SUBGRAPH.GetNeighborData( index_in_subgraph).GetSize()
                > GRAPH.GetNeighborData( index_in_graph).GetSize()
              )
              {
                continue;
              }

              // check that the environment around this vertex in the subgraph is contained in the environment for GRAPH
              if( !CSISubstructure::IsContainedIn( vertex_datum_subgraph, vertex_data_graph( index_in_graph)))
              {
                continue;
              }

              // the vertex appears to be a possible match.
              vmv_row.PushBack( index_in_graph);
            }

            // if there were no compatible vertices, then clearly SUBGRAPH is not a subgraph of GRAPH, so stop
            if( vmv_row.GetSize() == 0)
            {
              break;
            }
          }
        }
        s_vmm_creation.Stop();
        return vertex_matching_matrix;
      }

      //! @brief get the vertices from a graph that have isomorphic edges and vertex types to those in another graph
      //! @param GRAPH_A, GRAPH_B the graphs of the isomorphism
      //! @param MATCHING_VERTEX_GROUP indices of graph a (1st vector) that thus far appear isomorphic to indices in graph b (2nd vector)
      //! @param MATCHING_VERTICES list of matching vertices to append results to
      //! @param FUNCTION function used to generate a key for each vertex in a graph
      //! @return true if there are still matching group sizes -> if an isomorphism is still possible
      template< typename t_KeyType, typename t_ArgType>
      static bool RefineMatchingVerticesForIsomorphism
      (
        const t_ArgType &GRAPH_A,
        const t_ArgType &GRAPH_B,
        const storage::VectorND< 2, storage::Vector< size_t> > &MATCHING_VERTEX_GROUP,
        storage::List< storage::VectorND< 2, storage::Vector< size_t> > > &MATCHING_VERTICES,
        t_KeyType ( *FUNCTION)( const t_ArgType &, const size_t &)
      )
      {
        // store the sizes of the graph for easy access
        const storage::Vector< size_t> &group_a( MATCHING_VERTEX_GROUP.First());
        const storage::Vector< size_t> &group_b( MATCHING_VERTEX_GROUP.Second());

        const size_t group_size( MATCHING_VERTEX_GROUP.First().GetSize());
        storage::List< t_KeyType> graph_b_environments;
        storage::List< storage::VectorND< 2, storage::Vector< size_t> > > new_groups;
        for( size_t index_in_group_b( 0); index_in_group_b < group_size; ++index_in_group_b)
        {
          const size_t vertex_b_id( group_b( index_in_group_b));
          t_KeyType env( FUNCTION( GRAPH_B, vertex_b_id));
          storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::iterator itr_id( new_groups.Begin());
          for
          (
            typename storage::List< t_KeyType>::const_iterator
              itr_env( graph_b_environments.Begin()), itr_env_end( graph_b_environments.End());
            itr_env != itr_env_end;
            ++itr_env, ++itr_id
          )
          {
            if( env == *itr_env)
            {
              itr_id->Second().PushBack( vertex_b_id);
              break;
            }
          }
          if( itr_id == new_groups.End())
          {
            new_groups.PushBack();
            new_groups.LastElement().Second().PushBack( vertex_b_id);
            graph_b_environments.PushBack( env);
          }
        }

        // allocate memory for the corresponding members of graph A
        for
        (
          storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::iterator
            itr( new_groups.Begin()), itr_end( new_groups.End());
          itr != itr_end;
          ++itr
        )
        {
          itr->First().AllocateMemory( itr->Second().GetSize());
        }

        for( size_t index_in_group_a( 0); index_in_group_a < group_size; ++index_in_group_a)
        {
          const size_t vertex_a_id( group_a( index_in_group_a));
          t_KeyType env( FUNCTION( GRAPH_A, vertex_a_id));
          storage::List< storage::VectorND< 2, storage::Vector< size_t> > >::iterator itr_id( new_groups.Begin());
          typename storage::List< t_KeyType>::const_iterator
              itr_env( graph_b_environments.Begin()), itr_env_end( graph_b_environments.End());
          for( ; itr_env != itr_env_end; ++itr_env, ++itr_id)
          {
            if( itr_id->First().GetSize() != itr_id->Second().GetSize() && *itr_env == env)
            {
              itr_id->First().PushBack( vertex_a_id);
              break;
            }
          }
          if( itr_env == itr_env_end)
          {
            return false;
          }
        }
        // add the new groups to matching vertices; faster/better than calling Append because this involves
        // no new allocation
        MATCHING_VERTICES.InternalData().splice( MATCHING_VERTICES.End(), new_groups.InternalData(), new_groups.Begin(), new_groups.End());
        return true;
      }

      //! @brief determine whether all the elements in one vector are also in another
      //! @param SUBSET, FULLSET vectors whose data to consider the subset of
      //! @return whether SUBSET is fully contained in FULLSET
      template< typename t_Data>
      static bool IsContainedIn
      (
        const storage::Vector< t_Data> &SUBSET,
        const storage::Vector< t_Data> &FULLSET
      )
      {
        // an empty subset can always be found in the fullset
        if( SUBSET.GetSize() == 0)
        {
          return true;
        }

        // a subset cannot be larger than the fullset
        if( SUBSET.GetSize() > FULLSET.GetSize())
        {
          return false;
        }

        // get unique datum counts from the subset
        storage::Map< t_Data, size_t> subset_data_counts;
        for( size_t index( 0), size( SUBSET.GetSize()); index != size; ++index)
        {
          subset_data_counts[ SUBSET( index)]++;
        }

        // walk through the fullset vector
        for
        (
          typename storage::Vector< t_Data>::const_iterator itr_full( FULLSET.Begin()), itr_full_end( FULLSET.End());
          itr_full != itr_full_end;
          ++itr_full
        )
        {
          // see if the data from the fullset vector exists in the part of the subset not already found in the fullset
          typename storage::Map< t_Data, size_t>::iterator itr_count( subset_data_counts.Find( *itr_full));

          if( itr_count == subset_data_counts.End()) // the current datum is not in the subset, so go on
          {
            continue;
          }

          // the current datum was in the subset, so decrease the number of this datum that we must still find in the
          // fullset to establish that it contains the subset
          itr_count->second--;

          // if the fullset has enough of this datum to satisfy the subset, then remove this element from the subset
          // data counting map
          if( itr_count->second == 0)
          {
            subset_data_counts.RemoveElement( itr_count);

            // if no more elements remain in the subset counting map, then the fullset contains all elements
            // in the subset
            if( subset_data_counts.IsEmpty())
            {
              return true;
            }
          }
        }

        // the elements of FULLSET are not contained in the subset
        return false;
      }

      //! @brief determine the number of elements in one vector are also in another
      //! @param VECTOR_A a vector of t_Data
      //! @param VECTOR_B another vector
      //! @return the number of elements in A that could be matched to something in B
      template< typename t_Data>
      static size_t GetOverlap
      (
        const storage::Vector< t_Data> &VECTOR_A,
        const storage::Vector< t_Data> &VECTOR_B
      )
      {
        // nothing to be found if either vector is empty
        if( VECTOR_A.IsEmpty() || VECTOR_B.IsEmpty())
        {
          return 0;
        }

        // get unique datum counts from the subset
        storage::Map< t_Data, size_t> subset_data_counts;
        for( size_t index( 0), size( VECTOR_A.GetSize()); index != size; ++index)
        {
          subset_data_counts[ VECTOR_A( index)]++;
        }

        // track overlap
        size_t overlap( 0);

        // walk through the fullset vector
        for
        (
          typename storage::Vector< t_Data>::const_iterator itr_full( VECTOR_B.Begin()), itr_full_end( VECTOR_B.End());
          itr_full != itr_full_end;
          ++itr_full
        )
        {
          // see if the data from the fullset vector exists in the part of the subset not already found in the fullset
          typename storage::Map< t_Data, size_t>::iterator itr_count( subset_data_counts.Find( *itr_full));

          if( itr_count == subset_data_counts.End()) // the current datum is not in the subset, so go on
          {
            continue;
          }

          // the current datum was in the subset, so decrease the number of this datum that we must still find in the
          // fullset to establish that it contains the subset
          itr_count->second--;
          ++overlap;

          // if the fullset has enough of this datum to satisfy the subset, then remove this element from the subset
          // data counting map
          if( itr_count->second == 0)
          {
            subset_data_counts.RemoveElement( itr_count);

            // if no more elements remain in the subset counting map, then the fullset contains all elements
            // in the subset
            if( subset_data_counts.IsEmpty())
            {
              return overlap;
            }
          }
        }

        // the elements of FULLSET are not contained in the subset
        return overlap;
      }

      //! @brief determine whether all the elements in a map from element to count are in another
      //! @param SUBSET, FULLSET maps from data to counts of that data inside another container
      //! @return whether SUBSET is fully contained in FULLSET
      template< typename t_Data>
      static bool IsContainedIn
      (
        const storage::Map< t_Data, size_t> &SUBSET,
        const storage::Map< t_Data, size_t> &FULLSET
      )
      {
        // an empty subset can always be found in the fullset
        if( SUBSET.IsEmpty())
        {
          return true;
        }

        // a subset cannot be larger than the fullset
        if( SUBSET.GetSize() > FULLSET.GetSize())
        {
          return false;
        }

        // walk through the subset data count map
        auto itr_full( FULLSET.Begin()), itr_full_end( FULLSET.End());
        for
        (
          typename storage::Map< t_Data, size_t>::const_iterator
            itr( SUBSET.Begin()), itr_end( SUBSET.End());
          itr != itr_end;
          ++itr, ++itr_full
        )
        {
          while( itr_full != itr_full_end && itr_full->first < itr->first)
          {
            ++itr_full;
          }
          if( itr_full == itr_full_end || itr_full->first != itr->first || itr->second > itr_full->second)
          {
            return false;
          }
        }

        // all members of SUBSET were found in FULLSET
        return true;
      }

      //! @brief get a map whose keys are edge-types, vertex-types around a vertex
      //! @details the map maps the number of edges of that bond & atom type for a vertex
      //! @param GRAPH the graph to make the map for
      //! @param VERTEX the vertex index to consider
      //! @return a map from edge/connected-vertex type to number of occurances in the edges connected to VERTEX
      template< typename t_VertexType, typename t_EdgeType>
      static storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t>
        GetVertexEnvironment( const ConstGraph< t_VertexType, t_EdgeType> &GRAPH, const size_t &VERTEX)
      {
        // make the map that will store this vertex's environment
        storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t> environment;

        // get the vector of indices that this vertex is connected to
        const storage::Vector< size_t> &neighbor_indices( GRAPH.GetNeighborIndices( VERTEX));

        // get the vector of neighbor data
        const storage::Vector< t_EdgeType> &neighbor_data( GRAPH.GetNeighborData( VERTEX));

        // collect the bond-type atom-type pairs
        for
        (
          size_t neighbor_number( 0), number_neighbors( neighbor_indices.GetSize());
          neighbor_number < number_neighbors;
          ++neighbor_number
        )
        {
          // get the index of this neighbor in the graph
          const size_t neighbor_id( neighbor_indices( neighbor_number));

          // increment the count of that edge in the map
          ++environment[ std::make_pair( GRAPH.GetVertexData( neighbor_id), neighbor_data( neighbor_number))];
        }

        return environment;
      }

      //! @brief get a map whose keys are edge-types, vertex-types around a vertex
      //! @details the map maps the number of edges of that bond & atom type for a vertex
      //! @param GRAPH the graph to make the map for
      //! @param VERTEX the vertex index to consider
      //! @return a map from edge/connected-vertex type to number of occurances in the edges connected to VERTEX
      template< typename t_VertexType, typename t_EdgeType>
      static storage::Vector< storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t> >
        GetVertexEnvironments( const ConstGraph< t_VertexType, t_EdgeType> &GRAPH)
      {
        const size_t graph_size( GRAPH.GetSize());
        // make the vector that will hold all the environments
        storage::Vector< storage::Map< std::pair< t_VertexType, t_EdgeType>, size_t> > environments( graph_size);
        // walk over the graph and collect bond-type atom-type pairs
        for( size_t index_in_graph( 0); index_in_graph < graph_size; ++index_in_graph)
        {
          // get the map that will store this node's information
          environments( index_in_graph) = GetVertexEnvironment( GRAPH, index_in_graph);
        }

        return environments;
      }

      template< typename t_VertexType, typename t_EdgeType>
      static storage::Map< std::pair< t_VertexType, size_t>, size_t> GetVertexDistanceCounts
      (
        const ConstGraph< t_VertexType, t_EdgeType> &GRAPH,
        const size_t &VERTEX
      )
      {
        util::ShPtr< storage::Vector< size_t> > distances( Connectivity::DistancesToOtherVertices( GRAPH, VERTEX));

        storage::Map< std::pair< t_VertexType, size_t>, size_t> counts;
        for( size_t i( 0), graph_size( GRAPH.GetSize()); i < graph_size; ++i)
        {
          ++counts[ std::make_pair( GRAPH.GetVertexData( i), ( *distances)( i))];
        }
        return counts;
      }

      template< typename t_VertexType, typename t_EdgeType>
      static std::pair< t_VertexType, size_t> GetVertexType
      (
        const ConstGraph< t_VertexType, t_EdgeType> &GRAPH,
        const size_t &VERTEX
      )
      {
        return std::make_pair( GRAPH.GetVertexData( VERTEX), GRAPH.GetNeighborData( VERTEX).GetSize());
      }

    }; // class CSISubstructure

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_CSI_SUBSTRUCTURE_H_
