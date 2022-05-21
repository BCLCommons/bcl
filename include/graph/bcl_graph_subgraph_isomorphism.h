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

#ifndef BCL_GRAPH_SUBGRAPH_ISOMORPHISM_H_
#define BCL_GRAPH_SUBGRAPH_ISOMORPHISM_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_connectivity.h"
#include "bcl_graph_csi_substructure.h"
#include "bcl_graph_subgraph.h"
#include "bcl_graph_subgraph_isomorphism_base.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SubgraphIsomorphism
    //! @brief Solve the subgraph isomorphism problem on a ConstGraph with arbitrary data types
    //! @details allows user to choose either connected or unconnected solutions
    //!   The algorithm for unconnected csi is given in see Krissinel, Evgeny B., Henrick, Kim,
    //!   "Common subgraph isomorphism detection by backtracking search",
    //!   European Bioinformatics Institute, Genome Campus, Hinxton, Cambridge CB10 ISD, U.K.
    //!   The connected version is a simple refinement of the unconnected version wherein the vertex that is picked
    //!   is required to be connected to the existing isomorphism, if there is one
    //! @see @link example_graph_subgraph_isomorphism.cpp @endlink
    //! @author mendenjl
    //! @date Jun 04, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class SubgraphIsomorphism :
      public SubgraphIsomorphismBase
    {

    private:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class IntermediateSolution
      //! @brief Class used internally to find common substructure(s)
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Jun 04, 2012
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class IntermediateSolution
      {

      protected:

        //! a protected typedef for how the data is stored internally
        //! by putting this in the protected section, it can only be accessed by this class and deriving classes
        typedef storage::Vector< storage::Vector< size_t> > t_MatchingVertexVector;

      //////////
      // data //
      //////////

        size_t m_SubgraphSize;              //!< Size of graph A, cached for convenience
        size_t m_GraphSize;                 //!< Size of graph B, cached for convenience
        size_t m_NumberIsomorphisms;        //!< Actual # of isomorphisms
        size_t m_MaxDisparateIsomorphisms;  //!< Maximum # of disparate isomorphisms ( = BinomialCoeff( graph size, subgraph size)
        const size_t m_UndefinedVertexID;   //!< Undefined size_t for this function

        size_t m_CurrentIsomorphismSize; //!< number of vertices in the current isomorphism

        //! the two graphs being studied
        util::SiPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_Subgraph;
        util::SiPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_Graph;

        //! the largest isomorphism(s) found
        //! m_Subgraph(i) <-> m_Graph( m_BestIsomorphism(i)),
        //! if m_BestIsomorphism(i) != m_UndefinedVertexID
        storage::List< storage::Vector< size_t> > m_BestIsomorphisms;

        //! Sets of disparate sets of vertices that have been found as incidence strings
        storage::Set< std::string> m_DisparateIsomorphisms;

        //! the current isomorphism
        //! m_Subgraph(i) <-> m_Graph( m_CurrentIsomorphism(i)),
        //! if m_CurrentIsomorphism(i) != m_UndefinedVertexID
        storage::Vector< size_t> m_CurrentIsomorphism;

        //! m_CouldAddToIsomorphism(i) == 1 iff vertex i in m_Subgraph is not in the current isomorphism but could be added
        //! m_CouldAddToIsomorphism(i) == 0 iff vertex i cannot be added to the isomorphism for any reason
        storage::Vector< size_t> m_CouldAddToIsomorphism;

        //! vector of the maximum number of MatchingVertexVectors that might be needed
        //! while finding the isomorphism
        //! Storing these vectors means that memory will be allocated only O(m) times, vs. O(e^m)
        //!   if memory were allocated/deallocated within the recursive functions
        //! this is effectively a 3-dimensional vector, initially of size
        //! m_SubgraphSize X m_SubgraphSize X 0.
        storage::Vector< t_MatchingVertexVector> m_MatchingVertexStorage;

        //! Whether the next call to PickVertex should consider reachability
        //! This is generally true except when:
        //! 1. Picking the first vertex
        //! 2. On disconnected graphs, after matching a connected block
        bool m_ConsiderReachability;

        //! a vector containing the number of neighbors of a vertex in m_Subgraph/B that are
        //! currently in the isomorphism (only used/tracked if searching for a connected solution)
        storage::Vector< size_t> m_NeighborsInIsomorphismA;

        //! a vector containing the number of neighbors of a vertex in m_Subgraph/B that are
        //! currently in the isomorphism (only used/tracked if searching for a connected solution)
        storage::Vector< size_t> m_NeighborsInIsomorphismB;

      protected:

        //! @brief create an incidence string, 1 for each vertex that is in GRAPH_B in the isomorphism, - otherwise
        std::string GetCurrentIsomoprhismIncidenceString()
        {
          std::string incidence( m_GraphSize, '-');
          for
          (
            storage::Vector< size_t>::const_iterator
              itr( m_CurrentIsomorphism.Begin()), itr_end( m_CurrentIsomorphism.End());
            itr != itr_end;
            ++itr
          )
          {
            incidence[ *itr] = '1';
          }
          return incidence;
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
          m_CouldAddToIsomorphism( ID_IN_GRAPHA) = 1; // allow ID_IN_GRAPHA to be added to the isomorphism again
          m_CurrentIsomorphism( ID_IN_GRAPHA) = m_UndefinedVertexID; // Map of ID_IN_GRAPHA to the undefined vertex
          m_CurrentIsomorphismSize--; // record the decrement of the current isomorphism size

          const storage::Vector< size_t> &neighbors_a( m_Subgraph->GetNeighborIndices( ID_IN_GRAPHA));

          // decrement count of all neighbors of ID_IN_GRAPHA because they are now connected to by 1 less vertex in A
          for
          (
            size_t target_index( 0), neighbors_a_size( neighbors_a.GetSize());
            target_index != neighbors_a_size;
            ++target_index
          )
          {
            --m_NeighborsInIsomorphismA( neighbors_a( target_index));
          }

          // decrement count of all neighbors of ID_IN_GRAPHB because they are now connected to by 1 less vertex in B
          const storage::Vector< size_t> &neighbors_b( m_Graph->GetNeighborIndices( ID_IN_GRAPHB));
          for
          (
            size_t target_index( 0), neighbors_b_size( neighbors_b.GetSize());
            target_index != neighbors_b_size;
            ++target_index
          )
          {
            --m_NeighborsInIsomorphismB( neighbors_b( target_index));
          }
        }

        //! @brief Add a vertex from graph a and b to the current isomorphism
        //! @param ID_IN_GRAPHA The vertex (index) from graph A to add to the current isomorphism
        //! @param ID_IN_GRAPHB The vertex (index) from graph B to add to the current isomorphism
        void AddToIsomorphism( const size_t &ID_IN_GRAPHA, const size_t &ID_IN_GRAPHB)
        {
          // make sure that the vertex isn't already in the isomorphism.  If it is, something must have
          // went wrong in the algorithm
          m_CouldAddToIsomorphism( ID_IN_GRAPHA) = 0; // forbid ID_IN_GRAPHA from being added to the isomorphism twice
          m_CurrentIsomorphism( ID_IN_GRAPHA) = ID_IN_GRAPHB; // Map of ID_IN_GRAPHA to ID_IN_GRAPHB
          ++m_CurrentIsomorphismSize; // record the increment of the current isomorphism size

          // increment count of all neighbors of ID_IN_GRAPHA because they are now connected to by 1 more vertex in A
          const storage::Vector< size_t> &neighbors_a( m_Subgraph->GetNeighborIndices( ID_IN_GRAPHA));
          for
          (
            size_t target_index( 0), neighbors_a_size( neighbors_a.GetSize());
            target_index != neighbors_a_size;
            ++target_index
          )
          {
            ++m_NeighborsInIsomorphismA( neighbors_a( target_index));
          }

          // increment count of all neighbors of ID_IN_GRAPHB because they are now connected to by 1 more vertex in B
          const storage::Vector< size_t> &neighbors_b( m_Graph->GetNeighborIndices( ID_IN_GRAPHB));
          for
          (
            size_t target_index( 0), neighbors_b_size( neighbors_b.GetSize());
            target_index != neighbors_b_size;
            ++target_index
          )
          {
            ++m_NeighborsInIsomorphismB( neighbors_b( target_index));
          }
        }

        //! @brief The recursive part of the algorithm
        void Backtrack()
        {
          // pick a vertex in graph a and see what it could be isomorphic to in graph b
          const size_t picked_vertex_index( PickVertex());

          // get the potential matches in graph b to the picked_vertex_index from graph A (from the vertex matching)
          const storage::Vector< size_t> &mappable_vertices
          (
            m_MatchingVertexStorage( m_CurrentIsomorphismSize)( picked_vertex_index)
          );

          // check whether adding the next vertex will complete the isomorphism
          if( m_CurrentIsomorphismSize + 1 == m_SubgraphSize)
          {
            // adding any of the mappable vertices will result in the full isomorphism
            // just take the first isomorphism though
            m_CurrentIsomorphism( picked_vertex_index) = mappable_vertices.FirstElement();
            m_BestIsomorphisms.PushBack( m_CurrentIsomorphism);
            ++m_NumberIsomorphisms;
            m_CurrentIsomorphism( picked_vertex_index) = m_UndefinedVertexID;
          }
          else
          {
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
              if( Refine( picked_vertex_index, *itr_mappable_vertices))
              {
                Backtrack();

                // if an isomorphism was found, just return
                if( !m_BestIsomorphisms.IsEmpty())
                {
                  return;
                }
              }
              RemoveFromIsomorphism( picked_vertex_index, *itr_mappable_vertices);
            }
          }
        }

        //! @brief finds all disparate isomorphisms
        void FindDisparate( const size_t &MAX)
        {
          // yep, so pick a vertex in graph a and see what it could be isomorphic to in graph b
          const size_t picked_vertex_index( PickVertex());

          // get the potential matches in graph b to the picked_vertex_index from graph A (from the vertex matching)
          const storage::Vector< size_t> &mappable_vertices
          (
            m_MatchingVertexStorage( m_CurrentIsomorphismSize)( picked_vertex_index)
          );

          // check whether adding the next vertex will complete the isomorphism
          if( m_CurrentIsomorphismSize + 1 == m_SubgraphSize)
          {
            // adding any of the mappable vertices will result in the full isomorphism
            for // add each matching combination
            (
              storage::Vector< size_t>::const_iterator
                itr_mappable_vertices( mappable_vertices.Begin()),
                itr_mappable_vertices_end( mappable_vertices.End());
              itr_mappable_vertices != itr_mappable_vertices_end;
              ++itr_mappable_vertices
            )
            {
              // add picked_vertex_index <-> *itr_mappable_vertices to the isomorphism and see how far that gets us
              m_CurrentIsomorphism( picked_vertex_index) = *itr_mappable_vertices;
              if( m_DisparateIsomorphisms.Insert( GetCurrentIsomoprhismIncidenceString()).second)
              {
                m_BestIsomorphisms.PushBack( m_CurrentIsomorphism);
                ++m_NumberIsomorphisms;
                if( m_NumberIsomorphisms == std::min( m_MaxDisparateIsomorphisms, MAX))
                {
                  return;
                }
              }
            }
            m_CurrentIsomorphism( picked_vertex_index) = m_UndefinedVertexID;
          }
          else
          {
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
              if( Refine( picked_vertex_index, *itr_mappable_vertices))
              {
                FindDisparate( MAX);
                if( m_NumberIsomorphisms == std::min( m_MaxDisparateIsomorphisms, MAX))
                {
                  return;
                }
              }
              RemoveFromIsomorphism( picked_vertex_index, *itr_mappable_vertices);
              // try connecting the picked_vertex to something else if we can
            }
          }
        }

        //! @brief finds all isomorphisms larger than the min isomorphism size (rather than just one as backtrack does)
        void FindAll()
        {
          // yep, so pick a vertex in graph a and see what it could be isomorphic to in graph b
          const size_t picked_vertex_index( PickVertex());

          // get the potential matches in graph b to the picked_vertex_index from graph A (from the vertex matching)
          const storage::Vector< size_t> &mappable_vertices
          (
            m_MatchingVertexStorage( m_CurrentIsomorphismSize)( picked_vertex_index)
          );

          // check whether adding the next vertex will complete the isomorphism
          if( m_CurrentIsomorphismSize + 1 == m_SubgraphSize)
          {
            // adding any of the mappable vertices will result in the full isomorphism
            for // add each matching combination
            (
              storage::Vector< size_t>::const_iterator
                itr_mappable_vertices( mappable_vertices.Begin()),
                itr_mappable_vertices_end( mappable_vertices.End());
              itr_mappable_vertices != itr_mappable_vertices_end;
              ++itr_mappable_vertices
            )
            {
              // add picked_vertex_index <-> *itr_mappable_vertices to the isomorphism and see how far that gets us
              m_CurrentIsomorphism( picked_vertex_index) = *itr_mappable_vertices;
              m_BestIsomorphisms.PushBack( m_CurrentIsomorphism);
            }
            m_NumberIsomorphisms += mappable_vertices.GetSize();
            m_CurrentIsomorphism( picked_vertex_index) = m_UndefinedVertexID;
          }
          else
          {
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
              if( Refine( picked_vertex_index, *itr_mappable_vertices))
              {
                FindAll();
              }
              RemoveFromIsomorphism( picked_vertex_index, *itr_mappable_vertices);
              // try connecting the picked_vertex to something else if we can
            }
          }
        }

        //! @brief Pick a vertex to extend the existing isomorphism with
        //! @return index of the vertex with the smallest non-zero number of mappable vertices
        size_t PickVertex() const
        {
          const t_MatchingVertexVector &vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize));

          size_t best_index( m_UndefinedVertexID);
          if( m_ConsiderReachability)
          {
            for( size_t index( 0), fewest_mappable_vertices( m_UndefinedVertexID); index < m_SubgraphSize; ++index)
            {
              const size_t mappable_vertices( vmv( index).GetSize()); // number of vertices that index could map to in B
              if( m_NeighborsInIsomorphismA( index) && mappable_vertices && mappable_vertices < fewest_mappable_vertices)
              {
                if( mappable_vertices == size_t( 1))
                {
                  // stop if we find a vertex with only 1 potential match.
                  return index;
                }
                // choose the vertex with the least, non-zero number of potential matches in B
                fewest_mappable_vertices = mappable_vertices;
                best_index = index;
              }
            }
          }
          else
          {
            for( size_t index( 0), fewest_mappable_vertices( m_UndefinedVertexID); index < m_SubgraphSize; ++index)
            {
              const size_t mappable_vertices( vmv( index).GetSize()); // number of vertices that index could map to in B
              if( mappable_vertices && mappable_vertices < fewest_mappable_vertices)
              {
                if( mappable_vertices == size_t( 1))
                {
                  // stop if we find a vertex with only 1 potential match.
                  return index;
                }
                // choose the vertex with the least, non-zero number of potential matches in B
                fewest_mappable_vertices = mappable_vertices;
                best_index = index;
              }
            }
          }

          return best_index;
        }

        //! @brief Refine Update the current matching vertex structure with the most newly-added match (from graph a/b)
        //! @param LAST_INDEX_FROM_SUBGRAPH index of the most-recently added vertex from graph a
        //! @param LAST_INDEX_FROM_GRAPH index of the most-recently added vertex from graph b
        //! @return true if all unmapped vertices in the subgraph still have at least one
        bool Refine
        (
          const size_t &LAST_INDEX_FROM_SUBGRAPH,
          const size_t &LAST_INDEX_FROM_GRAPH
        )
        {
          // get a reference to the vertex matching vector that resulted in the current vertex being added to the isomorphism
          // which is the last one on the stack of vertex matching vectors
          t_MatchingVertexVector &old_vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize - 1));

          // get a reference to the current vertex matching vector
          t_MatchingVertexVector &new_vmv( m_MatchingVertexStorage( m_CurrentIsomorphismSize));

          // remove any data from the current vertex matching vector
          for( t_MatchingVertexVector::iterator itr( new_vmv.Begin()), itr_end( new_vmv.End()); itr != itr_end; ++itr)
          {
            itr->Resize( 0);
          }

          // test whether there are still reachable vertices
          bool still_reachable( false);
          bool still_chooseable( false);

          m_ConsiderReachability = true;

          // iterate over all the vertices in graph a
          for
          (
            typename ConstGraph< t_VertexData, t_EdgeData>::AdjacencyMatrixIterator
              itr_a( m_Subgraph->GetAdjacencyMatrixIterator( LAST_INDEX_FROM_SUBGRAPH));
            itr_a.GetTarget() < m_SubgraphSize;
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
              // then add the new vertex in graph b to the VMV.  Ignore LAST_INDEX_FROM_GRAPH,
              // since it would already be added
              const size_t index_in_graph_b( old_row( index_in_vmv));
              if
              (
                m_Graph->GetEdgeData( LAST_INDEX_FROM_GRAPH, index_in_graph_b) == edge_a_data
                && m_NeighborsInIsomorphismA( index_in_graph_a) == m_NeighborsInIsomorphismB( index_in_graph_b)
                && index_in_graph_b != LAST_INDEX_FROM_GRAPH
              )
              {
                // could still add old_row( index_in_vmv) to the isomorphism
                new_row.PushBack( index_in_graph_b);
              }
            }

            if( new_row.IsEmpty())
            {
              return false;
            }
            else if( m_NeighborsInIsomorphismA( index_in_graph_a))
            {
              still_reachable = true;
            }
            else
            {
              still_chooseable = true;
            }
          }

          if( still_reachable)
          {
            return true;
          }
          else if( still_chooseable)
          {
            m_ConsiderReachability = false;
            return true;
          }
          else
          {
            return false;
          }
        }

        //! @brief setup the intermediate solution for a new search
        //! @param SUBGRAPH should be the smaller of the two graphs for maximum speed
        //! @param GRAPH should be the larger of the two graphs for maximum speed
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        void Setup
        (
          const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          m_Subgraph = SUBGRAPH;
          m_Graph = GRAPH;
          BCL_Assert
          (
            SUBGRAPH.GetUnconnectedEdgeData() == GRAPH.GetUnconnectedEdgeData(),
            "CSI requires the graphs to have the same unconnected edge data value"
          );
          m_SubgraphSize = m_Subgraph->GetSize();
          m_GraphSize = m_Graph->GetSize();
          m_CurrentIsomorphismSize = 0;
          m_BestIsomorphisms.Reset();
          m_CurrentIsomorphism.Resize( m_SubgraphSize);
          m_CouldAddToIsomorphism.Resize( m_SubgraphSize);
          m_CurrentIsomorphism.SetAllElements( m_UndefinedVertexID);
          m_CouldAddToIsomorphism.SetAllElements( 1);
          m_DisparateIsomorphisms.Reset();
          m_MaxDisparateIsomorphisms = math::BinomialCoefficient( m_GraphSize, m_SubgraphSize);
          m_ConsiderReachability = false;

          m_MatchingVertexStorage.Resize( m_SubgraphSize); // allocate room for the vertex matching
          m_MatchingVertexStorage( 0) = MATCHING_VERTICES;
          for( size_t index( 1); index < m_SubgraphSize; ++index)
          {
            m_MatchingVertexStorage( index).Resize( m_SubgraphSize);
          }
          m_NeighborsInIsomorphismA.Resize( m_SubgraphSize);
          m_NeighborsInIsomorphismB.Resize( m_GraphSize);
          m_NeighborsInIsomorphismA.SetAllElements( 0);
          m_NeighborsInIsomorphismB.SetAllElements( 0);
        }

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor ( optional in this case)
        IntermediateSolution() :
          m_SubgraphSize( 0),
          m_GraphSize( 0),
          m_UndefinedVertexID( util::GetUndefinedSize_t()),
          m_CurrentIsomorphismSize( 0),
          m_Subgraph(),
          m_Graph()
        {
        }

        //! @brief Runs the common substructure algorithm
        //! @param SUBGRAPH should be the smaller of the two graphs for maximum speed
        //! @param GRAPH should be the larger of the two graphs for maximum speed
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @return a map of graph a indices to graph b indices for the best subisomorphism
        storage::Vector< storage::Vector< size_t> > Find
        (
          const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          Setup( SUBGRAPH, GRAPH, MATCHING_VERTICES);
          Backtrack(); // find one
          return storage::Vector< storage::Vector< size_t> >( m_BestIsomorphisms.Begin(), m_BestIsomorphisms.End());
        }

        //! @brief Runs the substructure algorithm to find all occurrences of subgraph within graph
        //! @param SUBGRAPH should be the smaller of the two graphs for maximum speed
        //! @param GRAPH should be the larger of the two graphs for maximum speed
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @param DISPARATE whether to return only disparate vertex isomorphisms
        //! @return a map of graph a indices to graph b indices for the best subisomorphism
        storage::Vector< storage::Vector< size_t> > FindAll
        (
          const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES
        )
        {
          Setup( SUBGRAPH, GRAPH, MATCHING_VERTICES);
          this->FindAll();
          return storage::Vector< storage::Vector< size_t> >( m_BestIsomorphisms.Begin(), m_BestIsomorphisms.End());
        }

        //! @brief Runs the substructure algorithm to find all occurrences of subgraph within graph
        //! @param SUBGRAPH should be the smaller of the two graphs for maximum speed
        //! @param GRAPH should be the larger of the two graphs for maximum speed
        //! @param MATCHING_VERTICES vertex indices from graph a that compare equal to vertex indices from graph b
        //! @param DISPARATE whether to return only disparate vertex isomorphisms
        //! @return a map of graph a indices to graph b indices for the best subisomorphism
        storage::Vector< storage::Vector< size_t> > FindDisparate
        (
          const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
          const ConstGraph< t_VertexData, t_EdgeData> &GRAPH,
          const storage::Vector< storage::Vector< size_t> > &MATCHING_VERTICES,
          const size_t &MAX
        )
        {
          Setup( SUBGRAPH, GRAPH, MATCHING_VERTICES);
          this->FindDisparate( MAX);
          return storage::Vector< storage::Vector< size_t> >( m_BestIsomorphisms.Begin(), m_BestIsomorphisms.End());
        }

      }; // class IntermediateSolution

    //////////
    // data //
    //////////

      //! the two graphs being studied
      util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_Subgraph;
      util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_Graph;

      //! isomorphism between graphs A and B, each pair connects isomorphic vertices
      storage::Vector< storage::Vector< t_VertexData> > m_DataIsomorphisms;

      //! Store the intermediate solution
      //! This is primarily done to enhance the speed of later isomorphisms, since we can avoid reallocating most of the arrays
      IntermediateSolution m_IntermediateSolution;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      //! @brief copy constructor
      SubgraphIsomorphism() :
        SubgraphIsomorphismBase(),
        m_Subgraph(),
        m_Graph(),
        m_DataIsomorphisms(),
        m_IntermediateSolution()
      {
      }

      //! @brief copy constructor
      SubgraphIsomorphism( const SubgraphIsomorphism &PARENT) :
        SubgraphIsomorphismBase( PARENT),
        m_Subgraph( PARENT.m_Subgraph),
        m_Graph( PARENT.m_Graph),
        m_DataIsomorphisms( PARENT.m_DataIsomorphisms),
        m_IntermediateSolution()
      {
      }

      //! clone the object
      SubgraphIsomorphism *Clone() const
      {
        return new SubgraphIsomorphism( *this);
      }

    /////////////////
    // Data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      void SetIsomorphisms( const storage::Vector< storage::Vector< size_t> > &ISOMORPHISMS)
      {
        m_Isomorphisms = ISOMORPHISMS;
        RegenerateVertexDataIsomorphismMap();
      }

      //! @brief Get the isomorphic vertices between graphs A and B
      //! @return the isomorphic vertices
      const storage::Vector< storage::Vector< t_VertexData> > &GetVertexIsomorphisms() const
      {
        return m_DataIsomorphisms;
      }

      //! @brief Get the 1st graph
      //! @return the 1st graph
      const ConstGraph< t_VertexData, t_EdgeData> &GetSubgraph() const
      {
        return *m_Subgraph;
      }

      //! @brief Get the 2nd graph
      //! @return the 2nd graph
      const ConstGraph< t_VertexData, t_EdgeData> &GetGraph() const
      {
        return *m_Graph;
      }

      //! @brief Set m_Subgraph
      //! @param GRAPH a graph to set as m_Subgraph
      void SetSubgraph( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_Subgraph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH.Clone(), true);
      }

      //! @brief Set m_Graph
      //! @param GRAPH  a graph to set as m_Graph
      void SetGraph( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_Graph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( GRAPH.Clone(), true);
      }

      //! @brief Set m_Subgraph; but assume another object owns it
      //! @param GRAPH a graph to set as m_Subgraph
      void SetSubgraphExternalOwnership( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_Subgraph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( &GRAPH, false);
      }

      //! @brief Set m_Graph
      //! @param GRAPH  a graph to set as m_Graph
      void SetGraphExternalOwnership( const ConstGraph< t_VertexData, t_EdgeData> &GRAPH)
      {
        m_Graph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( &GRAPH, false);
      }

      //! @brief Set the graphs
      //! @param SUBGRAPH the subgraph
      //! @param GRAPH the other graph
      void SetGraphs
      (
        const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH
      )
      {
        SetSubgraph( SUBGRAPH);
        SetGraph( GRAPH);
      }

      //! @brief get all the isomorphisms from the subgraph to the graph
      //! @return all the isomorphisms from the subgraph to the graph
      storage::Vector< Subgraph< t_VertexData, t_EdgeData> > GetSubgraphIsomorphisms() const
      {
        // create subgraphs from each of the isomorphisms
        storage::Vector< Subgraph< t_VertexData, t_EdgeData> > subgraphs;
        subgraphs.AllocateMemory( m_Isomorphisms.GetSize());
        for
        (
          storage::Vector< storage::Vector< size_t> >::const_iterator
            itr( m_Isomorphisms.Begin()), itr_end( m_Isomorphisms.End());
          itr != itr_end;
          ++itr
        )
        {
          subgraphs.PushBack( Subgraph< t_VertexData, t_EdgeData>( m_Graph, *itr));
        }
        return subgraphs;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief find the isomorphism between two graphs
      //! @param SUBGRAPH, GRAPH the graphs to find the subgraph isomorphism for
      void FindIsomorphism
      (
        const ConstGraph< t_VertexData, t_EdgeData> &SUBGRAPH,
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH
      )
      {
        SetSubgraph( SUBGRAPH);
        SetGraph( GRAPH);
        this->FindIsomorphism();
      }

      //! @brief find the isomorphism between the two graphs that are already set
      //! @return true iff an isomorphism is found
      bool FindIsomorphism()
      {
        return InternalFindIsomorphism( false, false);
      }

      //! @brief find all possible isomorphisms between the two graphs that are already set
      //! @return true iff an isomorphism is found
      bool FindAllIsomorphisms()
      {
        return InternalFindIsomorphism( true, false);
      }

      //! @brief find all disparate isomorphisms (disparate isomorphisms have unique vertex sets) between the two graphs that are already set
      //! @return true iff an isomorphism is found
      //! @note this should be used whenver symmetry can be ignored but otherwise all isomorphisms are desired
      bool FindDisparateIsomorphisms()
      {
        return InternalFindIsomorphism( true, true);
      }

      //! @brief assignment operator
      //! @param PARENT parent isomorphism
      //! @return a reference to this
      SubgraphIsomorphism &operator =( const SubgraphIsomorphism &PARENT)
      {
        if( this != &PARENT)
        {
          SubgraphIsomorphismBase::operator=( PARENT);
          m_Graph = PARENT.m_Graph;
          m_Subgraph = PARENT.m_Subgraph;
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
        SubgraphIsomorphismBase::Read( ISTREAM);
        bool subg_isdefined( false), graph_is_defined( false);
        io::Serialize::Read( subg_isdefined, ISTREAM);
        io::Serialize::Read( graph_is_defined, ISTREAM);
        if( subg_isdefined)
        {
          m_Subgraph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( new ConstGraph< t_VertexData, t_EdgeData>());
          io::Serialize::Read( *m_Subgraph, ISTREAM);
        }
        if( graph_is_defined)
        {
          m_Graph = util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> >( new ConstGraph< t_VertexData, t_EdgeData>());
          io::Serialize::Read( *m_Graph, ISTREAM);
        }

        if( graph_is_defined && subg_isdefined)
        {
          RegenerateVertexDataIsomorphismMap();
        }
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
        SubgraphIsomorphismBase::Write( OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Subgraph.IsDefined() && m_Subgraph.IsOwner(), OSTREAM, INDENT) << '\t'
                       << ( m_Graph.IsDefined() && m_Graph.IsOwner()) << '\n';
        if( m_Subgraph.IsDefined() && m_Subgraph.IsOwner())
        {
          io::Serialize::Write( *m_Subgraph, OSTREAM, INDENT) << '\n';
        }
        if( m_Graph.IsDefined() && m_Graph.IsOwner())
        {
          io::Serialize::Write( *m_Graph, OSTREAM, INDENT) << '\n';
        }

        // return
        return OSTREAM;
      }

      //! @brief find the isomorphism between the two graphs that are already set
      //! @param FIND_ALL true -> find all isomorphisms between the graphs that are at least LOWER_BOUNDS in size
      //!                 false -> only find the largest isomorphism
      //! @param FIND_DISPARATE whether to only find disparate vertex sets, e.g. ignore symmetry
      //! @return true iff an isomorphism is found
      bool InternalFindIsomorphism( const bool &FIND_ALL, const bool &FIND_DISPARATE)
      {
        BCL_Assert( !FIND_DISPARATE || FIND_ALL, "Find all must be given with find disparate");
        // reset any existing isomorphisms
        m_Isomorphisms.Reset();
        m_DataIsomorphisms.Reset();

        // return if there is no way that subgraph is in m_Graph
        if
        (
          m_Subgraph->GetSize() > m_Graph->GetSize()
          || m_Subgraph->NumEdges() > m_Graph->NumEdges()
          || !CSISubstructure::IsContainedIn( m_Subgraph->GetVertexTypeCountMap(), m_Graph->GetVertexTypeCountMap())
          || !CSISubstructure::IsContainedIn( m_Subgraph->GetEdgeTypeCountMap(), m_Graph->GetEdgeTypeCountMap())
        )
        {
          return false;
        }

        // make a vector to hold the matching vertices (in case MATCHING_VERTICES was not given)
        storage::Vector< storage::Vector< size_t> > matching_vertices
        (
          CSISubstructure::GetVertexMatchingMatrixForSubgraph( *m_Subgraph, *m_Graph)
        );

        // if there are no matching vertices remaining, the vertex isomorphism will not work
        if( matching_vertices.IsEmpty() || matching_vertices.LastElement().IsEmpty())
        {
          if( m_Subgraph->GetSize() == size_t( 0))
          {
            m_Isomorphisms.PushBack();
            m_DataIsomorphisms.PushBack();
            return true;
          }
          return false;
        }
        static util::Stopwatch s_iso( "Isomorphism search", util::Time( 1, 0), util::Message::e_Verbose, true, false);
        s_iso.Start();
        // find the isomorphism(s)
        if( FIND_ALL)
        {
          if( FIND_DISPARATE)
          {
            storage::Set< std::string> incidence_strings;
            storage::Vector< storage::Vector< size_t> > new_iso
            (
              m_IntermediateSolution.Find( *m_Subgraph, *m_Graph, matching_vertices)
            );
            bool can_continue( true);
            while( new_iso.GetSize() && !new_iso.FirstElement().IsEmpty() && can_continue)
            {
              std::string incidence( m_Graph->GetVertices().GetSize(), '0');

              size_t index( 0);
              for
              (
                auto itrf( new_iso.FirstElement().Begin()), itrfend( new_iso.FirstElement().End());
                itrf != itrfend;
                ++itrf, ++index
              )
              {
                matching_vertices( index).RemoveElements( matching_vertices( index).Find( *itrf), 1);
                incidence[ *itrf] = '1';
                if( matching_vertices( index).IsEmpty())
                {
                  can_continue = false;
                }
              }
              if( incidence_strings.Insert( incidence).second)
              {
                m_Isomorphisms.PushBack( new_iso.FirstElement());
              }
              if( can_continue)
              {
                new_iso = m_IntermediateSolution.Find( *m_Subgraph, *m_Graph, matching_vertices);
              }
            }
          }
          else
          {
            m_Isomorphisms = m_IntermediateSolution.FindAll( *m_Subgraph, *m_Graph, matching_vertices);
          }
        }
        else
        {
          m_Isomorphisms = m_IntermediateSolution.Find( *m_Subgraph, *m_Graph, matching_vertices);
        }
        s_iso.Stop();
        // regenerate the data isomorphism
        RegenerateVertexDataIsomorphismMap();

        return m_Isomorphisms.GetSize();
      }

      //! @brief regenerate m_DataIsomorphisms from m_Isomorphisms and the graphs
      void RegenerateVertexDataIsomorphismMap()
      {
        m_DataIsomorphisms.Resize( m_Isomorphisms.GetSize());
        typename storage::Vector< storage::Vector< t_VertexData> >::iterator itr_data_isos( m_DataIsomorphisms.Begin());

        for
        (
          storage::Vector< storage::Vector< size_t> >::const_iterator
            itr_isos( m_Isomorphisms.Begin()), itr_isos_end( m_Isomorphisms.End());
          itr_isos != itr_isos_end;
          ++itr_isos, ++itr_data_isos
        )
        {
          itr_data_isos->Reset();
          itr_data_isos->AllocateMemory( itr_isos->GetSize());
          // generate the new vertex isomorphism
          for
          (
            storage::Vector< size_t>::const_iterator iso( itr_isos->Begin()), iso_end( itr_isos->End());
            iso != iso_end;
            ++iso
          )
          {
            itr_data_isos->PushBack( m_Graph->GetVertexData( *iso));
          }
        }
      }

    }; // class SubgraphIsomorphism

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_SUBGRAPH_ISOMORPHISM_H_

