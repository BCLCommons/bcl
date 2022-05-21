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

#ifndef BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_H_
#define BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_common_subgraph_isomorphism_base.h"
#include "bcl_graph_subgraph.h"
#include "math/bcl_math_comparisons.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CommonSubgraphIsomorphism
    //! @brief Solve the subgraph isomorphism problem on a ConstGraph with arbitrary data types
    //! @details allows user to choose either connected or unconnected solutions
    //!   The algorithm for unconnected csi is given in see Krissinel, Evgeny B., Henrick, Kim,
    //!   "Common subgraph isomorphism detection by backtracking search",
    //!   European Bioinformatics Institute, Genome Campus, Hinxton, Cambridge CB10 ISD, U.K.
    //!   The connected version is a simple refinement of the unconnected version wherein the vertex that is picked
    //!   is required to be connected to the existing isomorphism, if there is one
    //! @see @link example_graph_common_subgraph_isomorphism.cpp @endlink
    //! @author mendenjl
    //! @date Aug 09, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class CommonSubgraphIsomorphism :
      public CommonSubgraphIsomorphismBase
    {

    private:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class IntermediateSolution
      //! @brief Class used internally to find common substructure(s)
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Aug 09, 2011
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class IntermediateSolution
      {

      protected:

        //! a protected typedef for how the data is stored internally
        //! by putting this in the protected section, it can only be accessed by this class and deriving classes
        typedef storage::Vector< storage::Vector< size_t> > t_MatchingVertexVector;

        //! holder of the class that does the edge comparison
        typedef util::BinaryFunctionInterface< t_EdgeData, t_EdgeData, bool> t_EdgeComparison;

      //////////
      // data //
      //////////

        size_t m_IsomorphismSizeUpperBound; //!< Upper bound on the size of the isomorphism
        size_t m_GraphASize;                //!< Size of graph A, cached for convenience
        size_t m_GraphBSize;                //!< Size of graph B, cached for convenience
        const size_t m_UndefinedVertexID;   //!< Undefined size_t for this function
        CommonSubgraphIsomorphismBase::SolutionType m_SolutionType; //!< the solution type to seek

        //! minimum size acceptable for the isomorphism
        //! this effectively prunes the branches of the isomorphism and can greatly speed up the algorithm
        size_t m_MinIsomorphismSize;

        size_t m_BestIsomorphismSize;    //!< number of vertices in the largest sub-isomorphism found so far
        size_t m_CurrentIsomorphismSize; //!< number of vertices in the current isomorphism

        //! the two graphs being studied
        util::SiPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_GraphA;
        util::SiPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_GraphB;

        //! the largest isomorphism(s) found
        //! m_GraphA(i) <-> m_GraphB( m_BestIsomorphism(i)),
        //! if m_BestIsomorphism(i) != m_UndefinedVertexID
        storage::List< storage::Vector< size_t> > m_BestIsomorphisms;

        //! the current isomorphism
        //! m_GraphA(i) <-> m_GraphB( m_CurrentIsomorphism(i)),
        //! if m_CurrentIsomorphism(i) != m_UndefinedVertexID
        storage::Vector< size_t> m_CurrentIsomorphism;

        //! m_CouldAddToIsomorphism(i) == 1 iff vertex i in m_GraphA is not in the current isomorphism but could be added
        //! m_CouldAddToIsomorphism(i) == 0 iff vertex i cannot be added to the isomorphism for any reason
        storage::Vector< size_t> m_CouldAddToIsomorphism;

        //! vector of the maximum number of MatchingVertexVectors that might be needed
        //! while finding the isomorphism
        //! Storing these vectors means that memory will be allocated only O(m) times, vs. O(e^m)
        //!   if memory were allocated/deallocated within the recursive functions
        //! this is effectively a 3-dimensional vector, initially of size
        //! m_IsomorphismSizeUpperBound X m_GraphASize X 0.
        storage::Vector< t_MatchingVertexVector> m_MatchingVertexStorage;

        //! method for comparing edges
        const t_EdgeComparison &m_EdgeComparer;

        //! a vector containing the number of neighbors of a vertex in m_GraphA/B that are
        //! currently in the isomorphism (only used/tracked if searching for a connected solution)
        storage::Vector< size_t> m_NeighborsInIsomorphismA;

        //! a vector containing the number of neighbors of a vertex in m_GraphA/B that are
        //! currently in the isomorphism (only used/tracked if searching for a connected solution)
        storage::Vector< size_t> m_NeighborsInIsomorphismB;

        //! each t_MatchingVertexVector in m_ReachableVertexStorage is a subset of
        //! of a corresponding t_MatchingVertexVector in m_MatchingVertexStorage
        //! Only used/tracked if searching for a connected solution
        storage::Vector< t_MatchingVertexVector> m_ReachableVertexStorage;

        //! Reference to the storage containing potential vertex matchings
        //! m_ReachableVertexStorage for connected solutions,
        //! m_MatchingVertexStorage  for unconnected solutions
        storage::Vector< t_MatchingVertexVector> &m_VertexMatchingStorage;

      protected:

        //! @brief ensure that ID cannot be added to the isomorphism again until AllowInIsomorphism is called
        //! @param ID of vertex most recently added or otherwise to be excluded
        void ForbidFromIsomorphism( const size_t &ID)
        {
          m_CouldAddToIsomorphism( ID) = 0;
        }

        //! @brief mark vertex at index ID as eligible for inclusion in the isomorphism
        void AllowInIsomorphism( const size_t &ID)
        {
          m_CouldAddToIsomorphism( ID) = 1;
        }

        //! @brief Remove a vertex from graph a and b to the current isomorphism
        //! @param ID_IN_GRAPHA The vertex (index) from graph A to remove from the current isomorphism
        //! @param ID_IN_GRAPHB The vertex (index) from graph B to remove from the isomorphism.  Used in derived class
        void RemoveFromIsomorphism( const size_t &ID_IN_GRAPHA, const size_t &ID_IN_GRAPHB)
        {
          // make sure that the vertex is in the isomorphism.  If it is not, something must have
          // went wrong in the algorithm
          BCL_Assert
          (
            m_CouldAddToIsomorphism( ID_IN_GRAPHA) == 0,
            "Cannot remove a vertex from the isomorphism that was not already added to it"
          );
          AllowInIsomorphism( ID_IN_GRAPHA); // allow ID_IN_GRAPHA to be added to the isomorphism again
          m_CurrentIsomorphism( ID_IN_GRAPHA) = m_UndefinedVertexID; // Map of ID_IN_GRAPHA to the undefined vertex
          m_CurrentIsomorphismSize--; // record the decrement of the current isomorphism size

          if( m_SolutionType != e_Unconnected)
          {
            const storage::Vector< size_t> &neighbors_a( m_GraphA->GetNeighborIndices( ID_IN_GRAPHA));

            // decrement count of all neighbors of ID_IN_GRAPHA because they are now connected to by 1 less vertex in A
            for
            (
              size_t target_index( 0), neighbors_a_size( neighbors_a.GetSize());
              target_index != neighbors_a_size;
              ++target_index
            )
            {
              m_NeighborsInIsomorphismA( neighbors_a( target_index))--;
            }

            // decrement count of all neighbors of ID_IN_GRAPHB because they are now connected to by 1 less vertex in B
            const storage::Vector< size_t> &neighbors_b( m_GraphB->GetNeighborIndices( ID_IN_GRAPHB));
            for
            (
              size_t target_index( 0), neighbors_b_size( neighbors_b.GetSize());
              target_index != neighbors_b_size;
              ++target_index
            )
            {
              m_NeighborsInIsomorphismB( neighbors_b( target_index))--;
            }
          }
        }

        //! @brief Add a vertex from graph a and b to the current isomorphism
        //! @param ID_IN_GRAPHA The vertex (index) from graph A to add to the current isomorphism
        //! @param ID_IN_GRAPHB The vertex (index) from graph B to add to the current isomorphism
        void AddToIsomorphism( const size_t &ID_IN_GRAPHA, const size_t &ID_IN_GRAPHB)
        {
          // make sure that the vertex isn't already in the isomorphism.  If it is, something must have
          // went wrong in the algorithm
          BCL_Assert
          (
            m_CouldAddToIsomorphism( ID_IN_GRAPHA) == 1,
            "Cannot add a vertex to the isomorphism that was already in it"
          );
          ForbidFromIsomorphism( ID_IN_GRAPHA); // forbid ID_IN_GRAPHA from being added to the isomorphism twice
          m_CurrentIsomorphism( ID_IN_GRAPHA) = ID_IN_GRAPHB; // Map of ID_IN_GRAPHA to ID_IN_GRAPHB
          m_CurrentIsomorphismSize++; // record the increment of the current isomorphism size

          if( m_SolutionType != e_Unconnected)
          {
            // increment count of all neighbors of ID_IN_GRAPHA because they are now connected to by 1 more vertex in A
            const storage::Vector< size_t> &neighbors_a( m_GraphA->GetNeighborIndices( ID_IN_GRAPHA));
            for
            (
              size_t target_index( 0), neighbors_a_size( neighbors_a.GetSize());
              target_index != neighbors_a_size;
              ++target_index
            )
            {
              m_NeighborsInIsomorphismA( neighbors_a( target_index))++;
            }

            // increment count of all neighbors of ID_IN_GRAPHB because they are now connected to by 1 more vertex in B
            const storage::Vector< size_t> &neighbors_b( m_GraphB->GetNeighborIndices( ID_IN_GRAPHB));
            for
            (
              size_t target_index( 0), neighbors_b_size( neighbors_b.GetSize());
              target_index != neighbors_b_size;
              ++target_index
            )
            {
              m_NeighborsInIsomorphismB( neighbors_b( target_index))++;
            }
          }
        }

        //! @brief The recursive part of the algorithm
        void Backtrack()
        {
          if( m_BestIsomorphismSize < m_IsomorphismSizeUpperBound)
          {
            if( Extendable()) // can the current isomorphism be extended?
            {
              // yep, so pick a vertex in graph a and see what it could be isomorphic to in graph b
              const size_t picked_vertex_index( PickVertex());

              // get the potential matches in graph b to the picked_vertex_index from graph A (from the vertex matching)
              const storage::Vector< size_t> &mappable_vertices
              (
                m_VertexMatchingStorage( m_CurrentIsomorphismSize)( picked_vertex_index)
              );

              for // for each vertex in graph b that might be isomorphic to the picked vertex (from graph a)
              (
                storage::Vector< size_t>::const_iterator
                  itr_mappable_vertices( mappable_vertices.Begin()),
                  itr_mappable_vertices_end( mappable_vertices.End());
                itr_mappable_vertices != itr_mappable_vertices_end
                  && m_BestIsomorphismSize < m_IsomorphismSizeUpperBound;
                ++itr_mappable_vertices
              )
              {
                // add picked_vertex_index <-> *itr_mappable_vertices to the isomorphism and see how far that gets us
                AddToIsomorphism( picked_vertex_index, *itr_mappable_vertices);
                Refine( picked_vertex_index, *itr_mappable_vertices);
                Backtrack();
                RemoveFromIsomorphism( picked_vertex_index, *itr_mappable_vertices);
                // try connecting the picked_vertex to something else if we can
              }

              if( m_BestIsomorphismSize < m_IsomorphismSizeUpperBound) // if we already reached the maximum size, we can stop already
              {
                // otherwise, remove the vertex most recently picked from consideration
                ForbidFromIsomorphism( picked_vertex_index);

                // see if that makes the isomorphism any bigger.  Effectively, this just makes us choose the next
                // potential vertex in PickVertex
                Backtrack();

                // add back the most recently picked vertex so that later isomorphisms can consider it
                AllowInIsomorphism( picked_vertex_index);
              }
            }
            else if( m_BestIsomorphismSize < m_CurrentIsomorphismSize && m_CurrentIsomorphismSize >= m_MinIsomorphismSize)
            {
              // cannot extend the current isomorphism anymore, and we have found a new best isomorphism
              m_BestIsomorphismSize = m_CurrentIsomorphismSize;
              m_BestIsomorphisms = storage::List< storage::Vector< size_t> >( 1, m_CurrentIsomorphism);
            }
          }
        }

        //! @brief finds all isomorphisms larger than the min isomorphism size (rather than just one as backtrack does)
        void FindAll()
        {
          if( Extendable()) // can the current isomorphism be extended?
          {
            // yep, so pick a vertex in graph a and see what it could be isomorphic to in graph b
            const size_t picked_vertex_index( PickVertex());

            // get the potential matches in graph b to the picked_vertex_index from graph A (from the vertex matching)
            const storage::Vector< size_t> &mappable_vertices
            (
              m_VertexMatchingStorage( m_CurrentIsomorphismSize)( picked_vertex_index)
            );

            for // for each vertex in graph b that might be isomorphic to the picked vertex (from graph a)
            (
              storage::Vector< size_t>::const_iterator
                itr_mappable_vertices( mappable_vertices.Begin()),
                itr_mappable_vertices_end( mappable_vertices.End());
              itr_mappable_vertices != itr_mappable_vertices_end;
              ++itr_mappable_vertices
            )
            {
              // add picked_vertex_index <-> *itr_mappable_vertices to the isomorphism and see how far that gets us
              AddToIsomorphism( picked_vertex_index, *itr_mappable_vertices);
              Refine( picked_vertex_index, *itr_mappable_vertices);
              FindAll();
              RemoveFromIsomorphism( picked_vertex_index, *itr_mappable_vertices);
              // try connecting the picked_vertex to something else if we can
            }

            // otherwise, remove the vertex most recently picked from consideration
            ForbidFromIsomorphism( picked_vertex_index);

            // see if that makes the isomorphism any bigger.  Effectively, this just makes us choose the next
            // potential vertex in PickVertex
            FindAll();

            // add back the most recently picked vertex so that later isomorphisms can consider it
            AllowInIsomorphism( picked_vertex_index);
          }
          else if( m_CurrentIsomorphismSize == m_IsomorphismSizeUpperBound)
          {
            // cannot extend the current isomorphism anymore, and the isomorphism is at least of the minimum size, so
            // add it to the list
            m_BestIsomorphismSize = m_IsomorphismSizeUpperBound - 1;
            m_BestIsomorphisms.PushBack( m_CurrentIsomorphism);
          }
        }

        //! @brief Check if there it is possible that the isomorphism can be extended to an isomorphism that is larger
        //! @brief than the current isomorphism and m_BestIsomorphism
        //! @return true if a larger isomorphism is possible
        bool Extendable() const
        {
          // get a reference to the current vmv
          const t_MatchingVertexVector &vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize));

          // how large does the isomorphism have to get to be relevant?
          // larger than the current isomorphism, the best isomorphism, and at least as big as the min isomorphism size
          const size_t goal_size
          (
            std::max( m_CurrentIsomorphismSize, std::max( m_BestIsomorphismSize, m_MinIsomorphismSize - 1))
          );

          // uncomment this to see the progression of isomorphism sizes
          //BCL_Message
          //(
          //  util::Message::e_Verbose,
          //  util::Format()( m_CurrentIsomorphismSize)
          //  + "<=" + util::Format()( m_BestIsomorphismSize)
          //  + "<=" + util::Format()( goal_size)
          //  + "<=" + util::Format()( m_IsomorphismSizeUpperBound)
          //);

          // track whether the size is sufficient
          bool size_is_sufficient( false);
          for( size_t index( 0), max_potential_size( m_CurrentIsomorphismSize); index < m_GraphASize; ++index)
          {
            // if a vertex in graph A is eligible for inclusion
            if
            (
              m_CouldAddToIsomorphism( index)
                && vmv( index).GetSize() != 0 // and has at least one remaining potential match in graph b
            )
            {
              // then the isomorphism might include this vertex, so increment the potential size this isomorphism may
              // grow to
              if( ++max_potential_size > goal_size)
              {
                // return when it is clear that the current isomorphism could grow larger than the goal size
                size_is_sufficient = true;
              }
            }
          }

          if( m_SolutionType == e_Unconnected)
          {
            // for unconnected solutions, it is enough to know whether the size is sufficient
            return size_is_sufficient;
          }

          // if the size was insufficient, then clearly the isomorphism cannot be extended
          if( !size_is_sufficient)
          {
            return false;
          }

          // connected solution, check to see that there is at least one vertex eligible for inclusion in the isomorphism
          // that can be reached from the existing isomorphism
          const t_MatchingVertexVector &connected_vmv( m_ReachableVertexStorage( m_CurrentIsomorphismSize));
          for( size_t index( 0); index != m_GraphASize; ++index)
          {
            if( m_CouldAddToIsomorphism( index) && connected_vmv( index).GetSize() != 0)
            {
              return true;
            }
          }

          return false;
        }

        //! @brief Pick a vertex to extend the existing isomorphism with
        //! @return index of the vertex with the smallest non-zero number of mappable vertices
        size_t PickVertex() const
        {
          const t_MatchingVertexVector &vmv( m_VertexMatchingStorage( m_CurrentIsomorphismSize));

          size_t best_index( 0);
          for
          (
            size_t index( 0), fewest_mappable_vertices( std::numeric_limits< size_t>::max());
            index != m_GraphASize;
            ++index
          )
          {
            // first check whether the vertex at this index is eligible for addition to the isomorphism
            if( m_CouldAddToIsomorphism( index))
            {
              // vertex is eligible.

              const size_t mappable_vertices( vmv( index).GetSize()); // number of vertices that index could map to in B
              if( mappable_vertices != size_t( 0) && mappable_vertices < fewest_mappable_vertices)
              {
                // choose the vertex with the least, non-zero number of potential matches in B
                fewest_mappable_vertices = mappable_vertices;
                best_index = index;
                if( mappable_vertices == size_t( 1))
                {
                  // stop if we find a vertex with only 1 potential match.
                  break;
                }
              }
            }
          }

          return best_index;
        }

        //! @brief Refine Update the current matching vertex structure with the most newly-added match (from graph a/b)
        //! @param LAST_INDEX_FROM_GRAPH_A index of the most-recently added vertex from graph a
        //! @param LAST_INDEX_FROM_GRAPH_B index of the most-recently added vertex from graph b
        void Refine
        (
          const size_t &LAST_INDEX_FROM_GRAPH_A,
          const size_t &LAST_INDEX_FROM_GRAPH_B
        )
        {
          // remove any data from the current vertex matching vector
          ClearVertexMatchingVector( m_MatchingVertexStorage( m_CurrentIsomorphismSize));

          // get a reference to the vertex matching vector that resulted in the current vertex being added to the isomorphism
          // which is the last one on the stack of vertex matching vectors
          t_MatchingVertexVector &old_vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize - 1));

          // get a reference to the current vertex matching vector
          t_MatchingVertexVector &new_vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize));

          // iterate over all the vertices in graph a
          for
          (
            typename ConstGraph< t_VertexData, t_EdgeData>::AdjacencyMatrixIterator
              itr_a( m_GraphA->GetAdjacencyMatrixIterator( LAST_INDEX_FROM_GRAPH_A));
            itr_a.GetTarget() < m_GraphASize;
            ++itr_a
          )
          {
            const size_t index_in_graph_a( itr_a.GetTarget());

            // check whether index_in_graph_a is eligible for addition to the isomorphism
            if( !m_CouldAddToIsomorphism( index_in_graph_a))
            {
              continue;
            }

            const t_EdgeData &edge_a_data( itr_a.GetData());

            // BCL_Assert
            // (
            //   itr_a.GetData() == m_GraphA->GetEdgeData( LAST_INDEX_FROM_GRAPH_A, itr_a.GetTarget()),
            //   "Wrong edge data: " + util::Format()( LAST_INDEX_FROM_GRAPH_A)
            //   + " " + util::Format()( itr_a.GetTarget()) + " Graph connectivity: " + m_GraphA->GetBasicConnectivity()
            // );

            // Get a reference to the VMV row that we'll be recording the potential matches in graph b in
            storage::Vector< size_t> &new_row( new_vmv( index_in_graph_a));

            // Get a reference to the equivalent row of the old vmv VMV row
            // this row contains a superset of the vertices that will go in the new row
            const storage::Vector< size_t> &old_row( old_vmv( index_in_graph_a));

            // Keep the size of the old row since we'll be iterating over it
            const size_t old_row_size( old_row.GetSize());

            // make sure the new row has enough space to hold all the elements in the old row
            new_row.AllocateMemory( old_row_size);

            for( size_t index_in_vmv( 0); index_in_vmv < old_row_size; ++index_in_vmv)
            {
              // if the new edge in graph a is equivalent to the potential new edge in graph b
              // then add the new vertex in graph b to the VMV.  Ignore LAST_INDEX_FROM_GRAPH_B,
              // since it would already be added
              if
              (
                m_EdgeComparer
                (
                  m_GraphB->GetEdgeData( LAST_INDEX_FROM_GRAPH_B, old_row( index_in_vmv)),
                  edge_a_data
                )
                && old_row( index_in_vmv) != LAST_INDEX_FROM_GRAPH_B
              )
              {
                // could still add old_row( index_in_vmv) to the isomorphism
                new_row.PushBack( old_row( index_in_vmv));
              }
            }
          }

          if( m_SolutionType != e_Unconnected)
          {
            // get a reference to the reachable vmv
            t_MatchingVertexVector &reachable_vmv( m_ReachableVertexStorage( m_CurrentIsomorphismSize));

            ClearVertexMatchingVector( reachable_vmv);

            for( size_t index_in_graph_a( 0); index_in_graph_a != m_GraphASize; ++index_in_graph_a)
            {
              if( m_NeighborsInIsomorphismA( index_in_graph_a))
              {
                // get a reference to the vector that will store the reachable potential matches for vertex at
                // index_in_graph_a in graph B
                storage::Vector< size_t> &reachable_row( reachable_vmv( index_in_graph_a));

                // allocate enough memory to hold all the connectivity-independent matches
                reachable_row.AllocateMemory( new_vmv( index_in_graph_a).GetSize());

                for
                (
                  storage::Vector< size_t>::const_iterator
                    itr_vertex_b( new_vmv( index_in_graph_a).Begin()),
                    itr_vertex_b_end( new_vmv( index_in_graph_a).End());
                  itr_vertex_b != itr_vertex_b_end;
                  ++itr_vertex_b
                )
                {
                  // add the subset of the normal vmv which is reachable
                  if( m_NeighborsInIsomorphismB( *itr_vertex_b) != 0)
                  {
                    reachable_row.PushBack( *itr_vertex_b);
                  }
                }
              }
            }
          }
        }

        //! @brief Resize (to 0) each vector in VMM(0...m_GraphASize).  This effectively clears the vectors, without
        //!        actually deallocating any memory
        void ClearVertexMatchingVector( t_MatchingVertexVector &VMM)
        {
          for( size_t a = 0; a < m_GraphASize; a++)
          {
            VMM( a).Resize( 0);
          }
        }

        //! @brief setup the intermediate solution for a new search
        //! @param GRAPH_A should be the smaller of the two graphs for maximum speed
        //! @param GRAPH_B should be the larger of the two graphs for maximum speed
        //! @param MAX_SIZE_OF_ISOMORPHISM the maximum size of the isomorphism, based on whatever criteria are available
        //! @param MIN_SIZE_OF_ISOMORPHISM do not calculate the maximum isomorphism if it would be smaller than this
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        void Setup
        (
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B,
          const size_t MAX_SIZE_OF_ISOMORPHISM,
          const size_t MIN_SIZE_OF_ISOMORPHISM,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          m_GraphA = GRAPH_A;
          m_GraphB = GRAPH_B;
          m_GraphASize = m_GraphA->GetSize();
          m_GraphBSize = m_GraphB->GetSize();
          BCL_Assert
          (
            m_GraphA->GetUnconnectedEdgeData() == m_GraphB->GetUnconnectedEdgeData(),
            "CSI requires the graphs to have the same unconnected edge data value"
          );
          m_IsomorphismSizeUpperBound = std::min( MAX_SIZE_OF_ISOMORPHISM, std::min( m_GraphASize, m_GraphBSize));
          m_MinIsomorphismSize = std::max( MIN_SIZE_OF_ISOMORPHISM, size_t( 1));
          m_BestIsomorphismSize = 0;
          m_CurrentIsomorphismSize = 0;
          m_BestIsomorphisms = storage::List< storage::Vector< size_t> >
                               (
                                 1,
                                 storage::Vector< size_t>( m_GraphASize, m_UndefinedVertexID)
                               );
          m_CurrentIsomorphism.Resize( m_GraphASize);
          m_CouldAddToIsomorphism.Resize( m_GraphASize);
          m_CurrentIsomorphism.SetAllElements( m_UndefinedVertexID);
          m_CouldAddToIsomorphism.SetAllElements( 1);

          m_MatchingVertexStorage.Resize( m_IsomorphismSizeUpperBound + 1); // allocate room for the vertex matching
          m_MatchingVertexStorage( 0) = MATCHING_VERTICES;
          for( size_t index( 1); index <= m_IsomorphismSizeUpperBound; ++index)
          {
            m_MatchingVertexStorage( index).Resize( m_GraphASize);
          }
          if( m_SolutionType != e_Unconnected)
          {
            m_NeighborsInIsomorphismA.Resize( m_GraphASize);
            m_NeighborsInIsomorphismB.Resize( m_GraphBSize);
            m_NeighborsInIsomorphismA.SetAllElements( 0);
            m_NeighborsInIsomorphismB.SetAllElements( 0);
            m_ReachableVertexStorage.Resize( m_IsomorphismSizeUpperBound + 1); // allocate room for the vertex matching
            m_ReachableVertexStorage( 0) = MATCHING_VERTICES;
            for( size_t index( 1); index <= m_IsomorphismSizeUpperBound; ++index)
            {
              m_ReachableVertexStorage( index).Resize( m_GraphASize);
            }
          }
        }

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor ( optional in this case)
        IntermediateSolution()
        {
        }

        //! @brief the primary constructor.  Sets up the intermediate solution on graph_a and graph_b
        //! @param MAX_SIZE_OF_ISOMORPHISM should be the maximum size that the isomorphism might become,
        //!        based on any information available.  If the calling class can determine that the isomorphism must
        //!        be even smaller, then set this number accordingly
        //! minimum size acceptable for the isomorphism
        //! @param IGNORE_BRANCHES_SMALLER_THAN do not calculate the maximum isomorphism if it would be smaller than this number
        //! @param EDGE_COMPARER how to compare edges ( must override bool operator()(t_EdgeData, t_EdgeData) const;
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @param SOLUTION_TYPE true -> narrow search to only find connected subgraph isomorphisms
        IntermediateSolution
        (
          const SolutionType &SOLUTION_TYPE,
          const t_EdgeComparison &EDGE_COMPARER
        ) :
          m_IsomorphismSizeUpperBound( 0),
          m_GraphASize( 0),
          m_GraphBSize( 0),
          m_UndefinedVertexID( util::GetUndefinedSize_t()),
          m_SolutionType( SOLUTION_TYPE),
          m_BestIsomorphismSize( 0),
          m_CurrentIsomorphismSize( 0),
          m_GraphA(),
          m_GraphB(),
          m_EdgeComparer( EDGE_COMPARER),
          m_VertexMatchingStorage( SOLUTION_TYPE != e_Unconnected ? m_ReachableVertexStorage : m_MatchingVertexStorage)
        {
        }

        //! @brief Runs the common substructure algorithm
        //! @param GRAPH_A should be the smaller of the two graphs for maximum speed
        //! @param GRAPH_B should be the larger of the two graphs for maximum speed
        //! @param MAX_SIZE_OF_ISOMORPHISM the maximum size of the isomorphism, based on whatever criteria are available
        //! @param MIN_SIZE_OF_ISOMORPHISM do not calculate the maximum isomorphism if it would be smaller than this
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @return a map of graph a indices to graph b indices for the best subisomorphism
        storage::Map< size_t, size_t> FindIsomorphism
        (
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B,
          const size_t MAX_SIZE_OF_ISOMORPHISM,
          const size_t MIN_SIZE_OF_ISOMORPHISM,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          //BCL_Debug( MATCHING_VERTICES);
          const SolutionType original_solution_type( m_SolutionType);
          size_t effective_min_size( MIN_SIZE_OF_ISOMORPHISM);
          if( m_SolutionType == e_GreedyUnconnected)
          {
            m_SolutionType = e_Connected;
            effective_min_size = 1;
          }
          Setup( GRAPH_A, GRAPH_B, MAX_SIZE_OF_ISOMORPHISM, effective_min_size, MATCHING_VERTICES);
          if( m_IsomorphismSizeUpperBound > 0)
          {
            Backtrack(); // find the largest isomorphism
          }

          // if no isomorphism was found, return an empty map
          if( m_BestIsomorphisms.IsEmpty())
          {
            return storage::Map< size_t, size_t>();
          }

          // get the first (and only) isomorphism from the map
          const storage::Vector< size_t> &best_isomorphism( m_BestIsomorphisms.FirstElement());
          storage::Map< size_t, size_t> isomorphism;
          for( size_t index( 0); index < m_GraphASize; index++)
          {
            // walk through the vector,  store any vertices that are defined
            if( best_isomorphism( index) != m_UndefinedVertexID)
            {
              isomorphism[ index] = best_isomorphism( index);
            }
          }

          m_SolutionType = original_solution_type;
          if
          (
            original_solution_type == e_GreedyUnconnected
            && isomorphism.GetSize() < MAX_SIZE_OF_ISOMORPHISM
          )
          {
            storage::Vector< size_t> is_available_graph_b( m_GraphBSize, size_t( 1));
            for( size_t index( 0); index < m_GraphBSize; ++index)
            {
              if( m_NeighborsInIsomorphismB( index))
              {
                is_available_graph_b( index) = size_t( 0);
              }
            }
            for( size_t index( 0); index < m_GraphASize; ++index)
            {
              if( best_isomorphism( index) != m_UndefinedVertexID)
              {
                is_available_graph_b( best_isomorphism( index)) = size_t( 0);
              }
            }
            storage::Vector< storage::Vector< size_t> > matching_vertices_reduced( m_GraphASize);
            size_t n_could_match( 0);
            for( size_t index( 0); index < m_GraphASize; index++)
            {
              // update the vertex matching matrix, removing all vertices already in the isomorphism as well as their
              // neighbors
              if( best_isomorphism( index) != m_UndefinedVertexID || m_NeighborsInIsomorphismA( index))
              {
                continue;
              }
              for
              (
                storage::Vector< size_t>::const_iterator
                  itr( MATCHING_VERTICES( index).Begin()), itr_end( MATCHING_VERTICES( index).End());
                itr != itr_end;
                ++itr
              )
              {
                if( is_available_graph_b( *itr))
                {
                  matching_vertices_reduced( index).PushBack( *itr);
                }
              }
              if( !matching_vertices_reduced( index).IsEmpty())
              {
                ++n_could_match;
              }
            }
            size_t best_isomorphism_connected_size( m_BestIsomorphismSize);
            if( n_could_match)
            {
              isomorphism.InsertElements
              (
                FindIsomorphism
                (
                  GRAPH_A,
                  GRAPH_B,
                  n_could_match,
                  std::max( best_isomorphism_connected_size, MIN_SIZE_OF_ISOMORPHISM) - best_isomorphism_connected_size,
                  matching_vertices_reduced
                )
              );
              best_isomorphism_connected_size = isomorphism.GetSize();
            }
            if( best_isomorphism_connected_size < MIN_SIZE_OF_ISOMORPHISM)
            {
              return storage::Map< size_t, size_t>();
            }
          }

          return isomorphism;
        }

        //! @brief Runs the common substructure algorithm until all isomorphisms of at least MIN_SIZE_OF_ISOMORPHISM are found
        //! @param GRAPH_A should be the smaller of the two graphs for maximum speed
        //! @param GRAPH_B should be the larger of the two graphs for maximum speed
        //! @param MAX_SIZE_OF_ISOMORPHISM the maximum size of the isomorphism, based on whatever criteria are available
        //! @param MIN_SIZE_OF_ISOMORPHISM do not calculate the maximum isomorphism if it would be smaller than this
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @return a vector of maps of graph a indices to graph b indices for the best subisomorphism
        storage::Vector< storage::Map< size_t, size_t> > FindAllIsomorphisms
        (
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B,
          const size_t MAX_SIZE_OF_ISOMORPHISM,
          const size_t MIN_SIZE_OF_ISOMORPHISM,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          const SolutionType original_solution_type( m_SolutionType);
          if( m_SolutionType == e_GreedyUnconnected)
          {
            BCL_MessageVrb
            (
              "Greedy unconnected solution is only supported when finding a single isomorphism. Defaulting to unconnected"
            );
            m_SolutionType = e_Unconnected;
          }
          Setup( GRAPH_A, GRAPH_B, MAX_SIZE_OF_ISOMORPHISM, MIN_SIZE_OF_ISOMORPHISM, MATCHING_VERTICES);

          m_BestIsomorphisms.Reset();
          if( m_IsomorphismSizeUpperBound > 0)
          {
            FindAll(); // find the largest isomorphism
          }

          storage::Vector< storage::Map< size_t, size_t> > isomorphisms( m_BestIsomorphisms.GetSize());

          size_t isomorphism_id( 0);

          // create isomorphism maps for each isomorphism
          for
          (
            storage::List< storage::Vector< size_t> >::const_iterator
              itr_iso( m_BestIsomorphisms.Begin()),
              itr_iso_end( m_BestIsomorphisms.End());
            itr_iso != itr_iso_end;
            ++itr_iso, ++isomorphism_id
          )
          {
            const storage::Vector< size_t> &best_isomorphism( *itr_iso);
            storage::Map< size_t, size_t> &isomorphism( isomorphisms( isomorphism_id));

            for( size_t index( 0); index < m_GraphASize; index++)
            {
              // walk through the vector,  store any vertices that are defined
              if( best_isomorphism( index) != m_UndefinedVertexID)
              {
                isomorphism[ index] = best_isomorphism( index);
              }
            }
          }
          m_SolutionType = original_solution_type;

          return isomorphisms;
        }

      }; // class IntermediateSolution

    //////////
    // data //
    //////////

      //! the two graphs being studied
      util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_GraphA;
      util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_GraphB;

      //! isomorphism between graphs A and B, each pair connects isomorphic vertices
      storage::Vector< storage::Map< t_VertexData, t_VertexData> > m_DataIsomorphisms;

      //! The object used to compare vertices
      const util::ShPtr< util::BinaryFunctionInterface< t_VertexData, t_VertexData, bool> > m_VertexComparer;

      //! The object used to compare edges
      const util::ShPtr< util::BinaryFunctionInterface< t_EdgeData, t_EdgeData, bool> > m_EdgeComparer;

      //! Store the intermediate solution
      //! This is primarily done to enhance the speed of later isomorphisms, since we can avoid reallocating most of the arrays
      IntermediateSolution m_IntermediateSolution;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief copy constructor
      CommonSubgraphIsomorphism( const CommonSubgraphIsomorphism &PARENT) :
        CommonSubgraphIsomorphismBase( PARENT),
        m_GraphA( PARENT.m_GraphA),
        m_GraphB( PARENT.m_GraphB),
        m_DataIsomorphisms( PARENT.m_DataIsomorphisms),
        m_VertexComparer( PARENT.m_VertexComparer),
        m_EdgeComparer( PARENT.m_EdgeComparer),
        m_IntermediateSolution( m_SolutionType, *PARENT.m_EdgeComparer)
      {
      }

      //! @brief constructor from two Graphs and a minimum size for the result, using the default operator==
      //! @param SOLUTION_TYPE whether the isomorphism must be connected
      //! @param VERTEX_COMPARER object that compares vertices
      //! @param EDGE_COMPARER object that compares edges
      CommonSubgraphIsomorphism
      (
        const SolutionType &SOLUTION_TYPE = CommonSubgraphIsomorphismBase::e_Connected,
        const util::ShPtr< util::BinaryFunctionInterface< t_VertexData, t_VertexData, bool> > &VERTEX_COMPARER
          = GetDefaultVertexComparison(),
        const util::ShPtr< util::BinaryFunctionInterface< t_EdgeData, t_EdgeData, bool> > &EDGE_COMPARER
          = GetDefaultEdgeComparison()
      ) :
        CommonSubgraphIsomorphismBase( SOLUTION_TYPE),
        m_GraphA(),
        m_GraphB(),
        m_DataIsomorphisms(),
        m_VertexComparer( VERTEX_COMPARER),
        m_EdgeComparer( EDGE_COMPARER),
        m_IntermediateSolution( m_SolutionType, *EDGE_COMPARER)
      {
      }

      //! @brief Get the default vertex comparison
      static util::ShPtr< util::BinaryFunctionInterface< t_VertexData, t_VertexData, bool> >
      GetDefaultVertexComparison()
      {
        return util::ShPtr< util::BinaryFunctionInterface< t_VertexData, t_VertexData, bool> >
               (
                 new util::BinaryFunctionSTLWrapper< std::equal_to< t_VertexData> >()
               );
      }

      //! @brief get the default edge comparison object
      static util::ShPtr< util::BinaryFunctionInterface< t_EdgeData, t_EdgeData, bool> >
        GetDefaultEdgeComparison()
      {
        return util::ShPtr< util::BinaryFunctionInterface< t_EdgeData, t_EdgeData, bool> >
               (
                 new util::BinaryFunctionSTLWrapper< std::equal_to< t_EdgeData> >()
               );
      }

      //! @brief find the isomorphism between two graphs
      //! @param GRAPH_A, GRAPH_B the graphs to find the subgraph isomorphism for
      //! @param ENFORCE_CONNECTIVITY whether the isomorphism must be connected
      //! @param UPPER_BOUNDS the maximum size to expect the isomorphism to be (often obtained using EstimateUpperBounds)
      //! @param LOWER_BOUNDS do not calculate the maximum isomorphism if it would be smaller than this number
      //! @param MATCHING_VERTICES the initial mapping of vertices in the smaller graph that map to vertices in the larger graph
      //!        If the graphs are equal in size, then MATCHING_VERTICES(i) should be the indices of vertices in graph B
      //!        that can be matched to the i-th vertex in GRAPH_A
      //! @param FIND_ALL true -> find all isomorphisms between the graphs that are at least LOWER_BOUNDS in size
      //!                 false -> only find the largest isomorphism
      void FindIsomorphism
      (
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B,
        const size_t UPPER_BOUNDS,
        const size_t LOWER_BOUNDS = 1,
        const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES =
          storage::Vector< storage::Vector< size_t> >(),
        const bool &FIND_ALL = false
      )
      {
        m_GraphA = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH_A.Clone());
        m_GraphB = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH_B.Clone());

        this->FindIsomorphism( UPPER_BOUNDS, LOWER_BOUNDS, MATCHING_VERTICES, FIND_ALL);
      }

      //! @brief find the isomorphism between the two graphs that are already set
      //! @param UPPER_BOUNDS the maximum size to expect the isomorphism to be (often obtained using EstimateUpperBounds)
      //! @param LOWER_BOUNDS do not calculate the maximum isomorphism if it would be smaller than this number
      //! @param MATCHING_VERTICES the initial mapping of vertices in the smaller graph that map to vertices in the larger graph
      //!        If the graphs are equal in size, then MATCHING_VERTICES(i) should be the indices of vertices in graph B
      //! @param FIND_ALL true -> find all isomorphisms between the graphs that are at least LOWER_BOUNDS in size
      //!                 false -> only find the largest isomorphism
      void FindIsomorphism
      (
        const size_t UPPER_BOUNDS,
        const size_t LOWER_BOUNDS = 1,
        const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES =
          storage::Vector< storage::Vector< size_t> >(),
        const bool &FIND_ALL = false
      )
      {
        if( UPPER_BOUNDS >= LOWER_BOUNDS) // the graphs could be equal
        {
          const ConstGraph< t_VertexData, t_EdgeData>
            &smaller_graph( m_GraphA->GetSize() <= m_GraphB->GetSize() ? *m_GraphA : *m_GraphB);
          const ConstGraph< t_VertexData, t_EdgeData>
            &larger_graph( m_GraphA->GetSize() <= m_GraphB->GetSize() ? *m_GraphB : *m_GraphA);
          // make a vector to hold the matching vertices (in case MATCHING_VERTICES was not given)
          storage::Vector< storage::Vector< size_t> > matching_vertices;

          // get a reference on the matching vertices vector that will be used
          const storage::Vector< storage::Vector< size_t> >
            &matching_vertices_ref( MATCHING_VERTICES.IsEmpty() ? matching_vertices : MATCHING_VERTICES);

          // generate the initial vertex matching matrix, if one was not given
          if( MATCHING_VERTICES.GetSize() == 0)
          {
            matching_vertices.Resize( smaller_graph.GetSize());
            for( size_t index_a( 0), size_a( smaller_graph.GetSize()); index_a < size_a; ++index_a)
            {
              const t_VertexData &vertex_a_datum( smaller_graph.GetVertexData( index_a));

              storage::Vector< size_t> &row( matching_vertices( index_a));
              row.AllocateMemory( larger_graph.GetSize());
              for( size_t index_b( 0), size_b( larger_graph.GetSize()); index_b < size_b; ++index_b)
              {
                if( ( *m_VertexComparer)( larger_graph.GetVertexData( index_b), vertex_a_datum))
                {
                  row.PushBack( index_b);
                }
              }
            }
          }
          else if( MATCHING_VERTICES.GetSize() != smaller_graph.GetSize())
          {
            BCL_Exit
            (
              GetClassIdentifier()
              + " was given a vertex matching matrix that was neither empty nor the size of the smaller graph ("
              + util::Format()( larger_graph.GetSize()) + "), aborting",
              -1
            );
          }

          // find the actual isomorphism
          if( FIND_ALL)
          {
            // find all isomorphisms
            m_Isomorphisms =
              m_IntermediateSolution.FindAllIsomorphisms
              (
                smaller_graph,
                larger_graph,
                UPPER_BOUNDS,
                LOWER_BOUNDS,
                matching_vertices_ref
              );
          }
          else
          {
            // just find the largest isomorphism
            m_Isomorphisms.Reset();
            m_Isomorphisms.PushBack
            (
              m_IntermediateSolution.FindIsomorphism
              (
                smaller_graph,
                larger_graph,
                UPPER_BOUNDS,
                LOWER_BOUNDS,
                matching_vertices_ref
              )
            );
          }

          // swap the graphs if m_GraphA is currently graph b (happens if graph a is larger than graph b)
          if( m_GraphA->GetSize() > m_GraphB->GetSize())
          {
            SwapIsomorphisms();
          }

          // regenerate the data isomorphism
          RegenerateVertexDataIsomorphismMap();
        }
        else
        {
          m_Isomorphisms.Reset();
          m_DataIsomorphisms.Reset();
        }

      }

      //! clone the object
      CommonSubgraphIsomorphism *Clone() const
      {
        return new CommonSubgraphIsomorphism( *this);
      }

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      void SetIsomorphisms( const storage::Vector< storage::Map< size_t, size_t> > &ISOMORPHISMS)
      {
        m_Isomorphisms = ISOMORPHISMS;
        RegenerateVertexDataIsomorphismMap();
      }

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      void SetIsomorphism( const storage::Map< size_t, size_t> &ISOMORPHISM)
      {
        m_Isomorphisms = storage::Vector< storage::Map< size_t, size_t> >( 1, ISOMORPHISM);
        RegenerateVertexDataIsomorphismMap();
      }

      //! @brief Get the isomorphic vertices between graphs A and B
      //! @return the isomorphic vertices
      const storage::Vector< storage::Map< t_VertexData, t_VertexData> > &GetVertexIsomorphisms() const
      {
        return m_DataIsomorphisms;
      }

      //! @brief Get the 1st graph
      //! @return the 1st graph
      const ConstGraph< t_VertexData, t_EdgeData> &GetGraphA() const
      {
        return *m_GraphA;
      }

      //! @brief Get the 2nd graph
      //! @return the 2nd graph
      const ConstGraph< t_VertexData, t_EdgeData> &GetGraphB() const
      {
        return *m_GraphB;
      }

      //! @brief Set m_GraphA
      //! @param GRAPH a graph to set as m_GraphA
      void SetGraphA( const util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > &GRAPH)
      {
        m_GraphA = GRAPH;
      }

      //! @brief Set m_GraphB
      //! @param GRAPH  a graph to set as m_GraphB
      void SetGraphB( const util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > &GRAPH)
      {
        m_GraphB = GRAPH;
      }

      //! @brief Set m_GraphA
      //! @param GRAPH a graph to set as m_GraphA
      void SetGraphA( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_GraphA = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH.Clone(), true);
      }

      //! @brief Set m_GraphB
      //! @param GRAPH  a graph to set as m_GraphB
      void SetGraphB( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_GraphB = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH.Clone(), true);
      }

      //! @brief Set the graphs
      //! @param GRAPH_A one of the graphs
      //! @param GRAPH_B another of the graphs
      void SetGraphs
      (
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B
      )
      {
        m_GraphA = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH_A.Clone());
        m_GraphB = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH_B.Clone());
      }

      //! @brief estimate upper bounds on the size of the isomorphism given the current graphs
      size_t EstimateUpperBounds() const
      {
        return CommonSubgraphIsomorphismBase::EstimateUpperBounds( GetGraphA(), GetGraphB());
      }

      //! @brief get all the common subgraphs of graph a
      //! @return all the common subgraphs of graph a
      storage::Vector< Subgraph< t_VertexData, t_EdgeData> > GetSubgraphIsomorphismsOfGraphA() const
      {
        // create subgraphs from each of the isomorphisms
        storage::Vector< Subgraph< t_VertexData, t_EdgeData> > subgraphs( m_Isomorphisms.GetSize());
        for
        (
          size_t isomorphism_number( 0), number_isomorphisms( m_Isomorphisms.GetSize());
          isomorphism_number < number_isomorphisms;
          ++isomorphism_number
        )
        {
          // GetKeys returns the indices of graph a that are in the isomorphism
          subgraphs( isomorphism_number)
            = Subgraph< t_VertexData, t_EdgeData>( m_GraphA, m_Isomorphisms( isomorphism_number).GetKeysAsVector());
        }
        if( subgraphs.IsEmpty())
        {
          subgraphs.PushBack( Subgraph< t_VertexData, t_EdgeData>( m_GraphA, storage::Vector< size_t>()));
        }
        return subgraphs;
      }

      //! @brief get all the common subgraphs of graph b
      //! @return all the common subgraphs of graph b
      storage::Vector< Subgraph< t_VertexData, t_EdgeData> > GetSubgraphIsomorphismsOfGraphB() const
      {
        // create subgraphs from each of the isomorphisms
        storage::Vector< Subgraph< t_VertexData, t_EdgeData> > subgraphs( m_Isomorphisms.GetSize());
        for
        (
          size_t isomorphism_number( 0), number_isomorphisms( m_Isomorphisms.GetSize());
          isomorphism_number < number_isomorphisms;
          ++isomorphism_number
        )
        {
          // GetMappedValues returns the indices of graph b that are in the isomorphism
          subgraphs( isomorphism_number)
            = Subgraph< t_VertexData, t_EdgeData>( m_GraphB, m_Isomorphisms( isomorphism_number).GetMappedValues());
        }
        if( subgraphs.IsEmpty())
        {
          subgraphs.PushBack( Subgraph< t_VertexData, t_EdgeData>( m_GraphB, storage::Vector< size_t>()));
        }
        return subgraphs;
      }

      //! @brief assignment operator
      //! @param PARENT parent isomorphism
      //! @return a reference to this
      CommonSubgraphIsomorphism &operator =( const CommonSubgraphIsomorphism &PARENT)
      {
        if( this != &PARENT)
        {
          CommonSubgraphIsomorphismBase::operator=( PARENT);
          m_GraphA = PARENT.m_GraphA;
          m_GraphB = PARENT.m_GraphB;
          m_DataIsomorphisms = PARENT.m_DataIsomorphisms;
        }
        return *this;
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
        // read members
        CommonSubgraphIsomorphismBase::Read( ISTREAM);
        m_GraphA = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( new ConstGraph< t_VertexData, t_EdgeData>());
        m_GraphB = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( new ConstGraph< t_VertexData, t_EdgeData>());
        io::Serialize::Read( *m_GraphA, ISTREAM);
        io::Serialize::Read( *m_GraphB, ISTREAM);

        RegenerateVertexDataIsomorphismMap();

        // return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        CommonSubgraphIsomorphismBase::Write( OSTREAM, INDENT) << '\n';
        io::Serialize::Write( *m_GraphA, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( *m_GraphB, OSTREAM, INDENT) << '\n';

        // return
        return OSTREAM;
      }

      //! @brief regenerate m_DataIsomorphisms from m_Isomorphisms and the graphs
      void RegenerateVertexDataIsomorphismMap()
      {
        m_DataIsomorphisms.Resize( m_Isomorphisms.GetSize());
        typename storage::Vector< storage::Map< t_VertexData, t_VertexData> >::iterator
          itr_data_isos( m_DataIsomorphisms.Begin());

        for
        (
          storage::Vector< storage::Map< size_t, size_t> >::const_iterator
            itr_isos( m_Isomorphisms.Begin()), itr_isos_end( m_Isomorphisms.End());
          itr_isos != itr_isos_end;
          ++itr_isos, ++itr_data_isos
        )
        {
          // erase the old isomorphism
          *itr_data_isos = storage::Map< t_VertexData, t_VertexData>();

          // generate the new vertex isomorphism
          for
          (
            typename storage::Map< size_t, size_t>::const_iterator iso( itr_isos->Begin()), iso_end( itr_isos->End());
            iso != iso_end;
            ++iso
          )
          {
            ( *itr_data_isos)[ m_GraphA->GetVertexData( iso->first)] = m_GraphB->GetVertexData( iso->second);
          }
        }
      }

    }; // class CommonSubgraphIsomorphism

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_H_

