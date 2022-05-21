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

#ifndef BCL_GRAPH_CONST_GRAPH_H_
#define BCL_GRAPH_CONST_GRAPH_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_own_ptr.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <numeric>

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstGraph
    //! @brief a high-performance graph class that has a number of vertices defined at construction.
    //! @details ConstGraph provides
    //!   1. O(1) time lookups for whether two vertices are connected and the type of edge that connects them
    //!   2. O(1) iteration through connections for a given vertex
    //!   3. O(1) time checks for whether a particular vertex index is in the graph
    //!   4. O(N log N) time checks for the distance from any particular vertex to all other vertices in a random graph
    //!      (O(N^2) worst-case performance on distance-from-one-to-all)
    //! @note By default, ConstGraph uses O(#edges+#vertices) memory
    //!
    //! @param t_VertexData - the type of data associated with each vertex. Must be comparable via ==
    //!        t_VertexData should not contain any information about its edges
    //!
    //! @param t_EdgeData   - the type of data associated with each edge.  Must be comparable via ==
    //!        To reduce memory usage, t_EdgeData should not contain any information about the vertices it connects to
    //!
    //! @see @link example_graph_const_graph.cpp @endlink
    //! @author mendenjl
    //! @date March 3, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class ConstGraph :
      public util::ObjectInterface
    {
    public:

      // Typedefs for how the graph stores its data
      // This makes it trivial to test different data structures and shortens what would otherwise be 50+ character
      // long type names, which are generally clunky and hard to read
      typedef storage::Vector< t_EdgeData> t_EdgeDataOfVertex;               //!< Edge data for a given vertex
      typedef storage::Vector< t_EdgeDataOfVertex> t_EdgeDataContainer;      //!< Container of Vertex Edge Data containers
      typedef storage::Vector< size_t> t_EdgeTargetsOfVertex;                //!< Neighbor IDs for a given vertex
      typedef storage::Vector< t_EdgeTargetsOfVertex> t_EdgeTargetContainer; //!< Container of Vertex Neighbor containers
      typedef t_EdgeData t_EdgeDatum;

    private:

    //////////
    // data //
    //////////

      size_t                         m_Size;                 //!< Size of the graph; cached
      storage::Vector< t_VertexData> m_Vertices;             //!< vertex data
      t_EdgeDataContainer            m_EdgeDataVectors;      //!< data for connected edges indexed by vertex
      t_EdgeTargetContainer          m_EdgeTargetVectors;    //!< connected edge targets indexed by vertex
      t_EdgeData                     m_UnconnectedEdgeValue; //!< value of m_EdgeDataMatrix(A,B) if A and B are unconnected
      size_t                         m_NumEdges;             //!< number of edges in the graph
      bool                           m_Directed;             //!< false if functions should maintain an undirected invariant

      //! edge data indexed by vertex from, vertex to
      //! Stored only on user request for this matrix, so has to be mutable
      mutable util::OwnPtr< linal::Matrix< t_EdgeData> > m_EdgeDataMatrix;

      //! all edges. This will be created and managed only if GetEdges is called
      mutable util::OwnPtr< storage::Vector< t_EdgeData> > m_EdgeDataVector;
      mutable util::OwnPtr< storage::Map< t_VertexData, size_t> > m_VertexDataCountMap;
      mutable util::OwnPtr< storage::Map< t_EdgeData, size_t> > m_EdgeDataCountMap;

      //! Mutex for creation of mutable data elements above
      mutable sched::Mutex m_Mutex;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, creates an empty graph
      ConstGraph() :
        m_Size( 0),
        m_NumEdges( 0),
        m_Directed( false)
      {
      }

      //! @brief construct from vertices and edges
      //! @param VERTEX_DATA a vector of t_VertexData containing the comparator for vertices
      //! @param EDGE_DATA a matrix of t_EdgeData, whose values indicate the type of edge between each pair of vertices
      //! @param UNCONNECTED_EDGE_VALUE value used to denote an unconnected edge in the matrix
      //! @param DIRECTED whether the graph is to be directed
      //! @param MAKE_UNDIRECTED boolean to whether to make the edge data matrix undirected even if the given matrix is not symmetric
      ConstGraph
      (
        const storage::Vector< t_VertexData> &VERTEX_DATA,
        const storage::Vector< UndirectedEdge< t_EdgeData> > &EDGE_DATA,
        const t_EdgeData &UNCONNECTED_EDGE_VALUE = t_EdgeData()
      ) :
        m_Size( VERTEX_DATA.GetSize()),
        m_Vertices( VERTEX_DATA),
        m_EdgeDataVectors( m_Size),
        m_EdgeTargetVectors( m_Size),
        m_UnconnectedEdgeValue( UNCONNECTED_EDGE_VALUE),
        m_NumEdges( 2 * EDGE_DATA.GetSize()),
        m_Directed( false)
      {
        if( EDGE_DATA.IsEmpty())
        {
          return;
        }

        // test whether edge data is sorted
        bool is_edge_data_sorted( true);
        for( size_t edge_number( 1), n_edges( EDGE_DATA.GetSize()); edge_number < n_edges; ++edge_number)
        {
          if( EDGE_DATA( edge_number) < EDGE_DATA( edge_number - 1))
          {
            is_edge_data_sorted = false;
            break;
          }
        }
        storage::Vector< UndirectedEdge< t_EdgeData> > edge_data_sort;
        if( !is_edge_data_sorted)
        {
          edge_data_sort = EDGE_DATA;
          edge_data_sort.Sort( std::less< UndirectedEdge< t_EdgeData> >());
        }

        const storage::Vector< UndirectedEdge< t_EdgeData> >
          &sorted_edges( is_edge_data_sorted ? EDGE_DATA : edge_data_sort);

        // construct m_EdgeTargetVectors and m_EdgeDataVectors from m_EdgeDataMatrix
        for
        (
          typename storage::Vector< UndirectedEdge< t_EdgeData> >::const_iterator
            itr( sorted_edges.Begin()), itr_end( sorted_edges.End());
          itr != itr_end;
          ++itr
        )
        {
          BCL_Assert( itr->GetIndexHigh() < m_Size, "Graph constructor called with edges to non-existant vertices!");
          m_EdgeTargetVectors( itr->GetIndexLow()).PushBack( itr->GetIndexHigh());
          m_EdgeDataVectors( itr->GetIndexLow()).PushBack( itr->GetEdgeData());
          m_EdgeTargetVectors( itr->GetIndexHigh()).PushBack( itr->GetIndexLow());
          m_EdgeDataVectors( itr->GetIndexHigh()).PushBack( itr->GetEdgeData());
        }

        // cache the adjacency matrix automatically if it consumes no more than 2x the amount of memory as the
        // adjacency list.  This ensures that very sparse graphs (e.g. chemical graphs) do not take up excessive space
        const size_t bytes_used_edges
        (
          ( sizeof( t_EdgeData) + sizeof( size_t)) * m_NumEdges + 2 * sizeof( storage::Vector< size_t>) * m_Size
        );
        const size_t bytes_used_adjacency_matrix
        (
          sizeof( t_EdgeData) * m_Size * m_Size + sizeof( linal::Matrix< t_EdgeData>)
        );
        if( bytes_used_edges * 2 > bytes_used_adjacency_matrix)
        {
          CacheAdjacencyMatrix();
        }
      }

      //! @brief construct from vertices and edges
      //! @param VERTEX_DATA a vector of t_VertexData containing the comparator for vertices
      //! @param EDGE_DATA a matrix of t_EdgeData, whose values indicate the type of edge between each pair of vertices
      //! @param UNCONNECTED_EDGE_VALUE value used to denote an unconnected edge in the matrix
      //! @param DIRECTED whether the graph is to be directed
      //! @param MAKE_UNDIRECTED boolean to whether to make the edge data matrix undirected even if the given matrix is not symmetric
      ConstGraph
      (
        const storage::Vector< t_VertexData> &VERTEX_DATA,
        const linal::Matrix< t_EdgeData> &EDGE_DATA,
        const t_EdgeData &UNCONNECTED_EDGE_VALUE = t_EdgeData(),
        const bool DIRECTED = false,
        const bool MAKE_UNDIRECTED = true
      ) :
        m_Size( VERTEX_DATA.GetSize()),
        m_Vertices( VERTEX_DATA),
        m_EdgeDataVectors( m_Size),
        m_EdgeTargetVectors( m_Size),
        m_UnconnectedEdgeValue( UNCONNECTED_EDGE_VALUE),
        m_NumEdges( 0)
      {
        // ensure that EDGE_DATA is a square matrix and that all the sizes are consistent
        BCL_Assert
        (
          m_Size == EDGE_DATA.GetNumberRows() && m_Size == EDGE_DATA.GetNumberCols(),
          "ConstGraph Constructor called with different size vertex/edge vectors or matrices"
        );

        // construct m_EdgeTargetVectors and m_EdgeDataVectors from m_EdgeDataMatrix
        for( size_t i( 0); i < m_Size; i++)
        {
          for( size_t j( 0); j < m_Size; j++)
          {
            if( EDGE_DATA( i, j) != m_UnconnectedEdgeValue)
            {
              m_EdgeTargetVectors( i).PushBack( j);
              m_EdgeDataVectors( i).PushBack( EDGE_DATA( i, j));
              ++m_NumEdges;
            }
            else if( MAKE_UNDIRECTED && EDGE_DATA( j, i) != m_UnconnectedEdgeValue)
            {
              m_EdgeTargetVectors( i).PushBack( j);
              m_EdgeDataVectors( i).PushBack( EDGE_DATA( j, i));
              ++m_NumEdges;
            }
          }
        }

        // Determine whether the graph is directed or not
        m_Directed = ( DIRECTED || !m_Size || !EDGE_DATA.IsSymmetric()) && !MAKE_UNDIRECTED;

        // check that if the graph was supposed to be undirected, it really was
        BCL_Assert
        (
          !( !DIRECTED && m_Size && m_Directed),
          " Error: Given edge data matrix was not symmetric but it should have been: "
          + util::Format()( EDGE_DATA)
        );

        // cache the adjacency matrix automatically if it consumes no more than 2x the amount of memory as the
        // adjacency list.  This ensures that very sparse graphs (e.g. chemical graphs) do not take up excessive space
        const size_t bytes_used_edges
        (
          ( sizeof( t_EdgeData) + sizeof( size_t)) * m_NumEdges + 2 * sizeof( storage::Vector< size_t>) * m_Size
        );
        const size_t bytes_used_adjacency_matrix
        (
          sizeof( t_EdgeData) * m_Size * m_Size + sizeof( linal::Matrix< t_EdgeData>)
        );
        if( bytes_used_edges * 2 > bytes_used_adjacency_matrix)
        {
          CacheAdjacencyMatrix();
        }
      }

      //! clone the object
      ConstGraph *Clone() const
      {
        return new ConstGraph( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @return true if the graph is directed
      bool IsDirected() const
      {
        return m_Directed;
      }

      //! @return sets the graph as directed
      void SetDirected()
      {
        m_Directed = true;
      }

      //! @return true if the graph is undirected
      bool IsUndirected() const
      {
        return !m_Directed;
      }

      //! @return the number of edges
      size_t NumEdges() const
      {
        return m_NumEdges;
      }

      //! @brief return the edge data matrix
      //! @return the edge data matrix
      const linal::Matrix< t_EdgeData> &GetEdgeDataMatrix() const
      {
        CacheAdjacencyMatrix();
        return *m_EdgeDataMatrix;
      }

      //! @brief GetSize
      //! @return the number of vertices in the matrix
      const size_t &GetSize() const
      {
        return m_Size;
      }

      //! @brief Given a vertex index, return its data
      const t_VertexData &GetVertexData( const size_t &VERTEX) const
      {
        return m_Vertices( VERTEX);
      }

      //! @brief Given a vertex index, return its data
      const storage::Vector< t_VertexData> &GetVertices() const
      {
        return m_Vertices;
      }

      //! @brief Given a vertex index, return its neighbors data
      const t_EdgeDataOfVertex &GetNeighborData( const size_t &VERTEX) const
      {
        return m_EdgeDataVectors( VERTEX);
      }

      //! @brief Given a vertex index, return its neighbors indices
      const t_EdgeTargetsOfVertex &GetNeighborIndices( const size_t &VERTEX) const
      {
        return m_EdgeTargetVectors( VERTEX);
      }

      //! @brief GetEdgeData Get the data corresponding to the edge from index VERTEX_A to VERTEX_B
      //! @param VERTEX_A index of the first vertex
      //! @param VERTEX_B index of the second vertex
      //! @return data corresponding to the edge
      const t_EdgeData &GetEdgeData( const size_t &VERTEX_A, const size_t &VERTEX_B) const
      {
        if( m_EdgeDataMatrix.IsDefined())
        {
          // call the adjacency matrix, if it is already available
          return GetEdgeDataMatrix()( VERTEX_A, VERTEX_B);
        }

        const storage::Vector< size_t> &neighbors_a( m_EdgeTargetVectors( VERTEX_A));
        const size_t index( neighbors_a.Find( VERTEX_B));
        return index < neighbors_a.GetSize() ? m_EdgeDataVectors( VERTEX_A)( index) : m_UnconnectedEdgeValue;
      }

      //! @brief retrieve the value used in the adjacency matrix to represent edges that are not connected
      //! @return the value used in the adjacency matrix to represent edges that are not connected
      const t_EdgeData &GetUnconnectedEdgeData() const
      {
        return m_UnconnectedEdgeValue;
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class AdjacencyMatrixIterator
      //! @brief acts like an iterator through an adjacency matrix, without requiring the actual adjacency matrix
      //! @author mendenjl
      //! @date May 29, 2012
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class BCL_API AdjacencyMatrixIterator
      {

      private:

      //////////
      // data //
      //////////

        size_t            m_Index;           //!< Logical index currently on
        const size_t     *m_ItrTarget;       //!< Pointer to next target
        const size_t     *m_ItrTargetEnd;    //!< End of targets
        const t_EdgeData *m_ItrData;         //!< Iterator on data
        const t_EdgeData &m_UnconnectedData; //!< Unconnected data type

      public:

        //! @brief construct iterator from target and data vectors
        //! @param TARGETS     vector of target indices
        //! @param DATA        vector of data indices
        //! @param UNCONNECTED value to return for non-connections
        AdjacencyMatrixIterator
        (
          const t_EdgeTargetsOfVertex &TARGETS,
          const t_EdgeDataOfVertex &DATA,
          const t_EdgeData &UNCONNECTED
        ) :
          m_Index( 0),
          m_ItrTarget( TARGETS.Begin().base()),
          m_ItrTargetEnd( TARGETS.End().base()),
          m_ItrData( DATA.Begin().base()),
          m_UnconnectedData( UNCONNECTED)
        {
        }

        //! @brief pre increment - sets m_Index to the next size_t of the adjacency matrix
        //! @return reference to this iterator
        AdjacencyMatrixIterator &operator ++()
        {
          ++m_Index;
          if( m_ItrTarget != m_ItrTargetEnd && m_Index > *m_ItrTarget)
          {
            ++m_ItrTarget;
            ++m_ItrData;
          }
          return *this;
        }

        //! @brief get the target indice
        //! @return the current target indice
        const size_t &GetTarget() const
        {
          return m_Index;
        }

        //! @brief get the data from the current target
        //! @return the data from the current target
        const t_EdgeData &GetData() const
        {
          if( m_ItrTarget == m_ItrTargetEnd || m_Index != *m_ItrTarget)
          {
            return m_UnconnectedData;
          }
          return *m_ItrData;
        }

      }; // class AdjacencyMatrixIterator

    ////////////////
    // operations //
    ////////////////

      //! @brief get an adjacency matrix iterator for a given vertex
      //! @param VERTEX the vertex for which to obtain the iterator
      //! @return an iterator that acts like it was on an adjacency matrix iterator
      AdjacencyMatrixIterator GetAdjacencyMatrixIterator( const size_t &VERTEX) const
      {
        return
          AdjacencyMatrixIterator( m_EdgeTargetVectors( VERTEX), m_EdgeDataVectors( VERTEX), m_UnconnectedEdgeValue);
      }

      //! @brief Cache complete adjacency matrix for O(1) calls to GetEdgeData()
      //! results in 10-20% faster largest common substructure search, but at a cost of nearly 10x memory for sparse,
      //! chemical graphs
      void CacheAdjacencyMatrix() const
      {
        m_Mutex.Lock();
        if( m_EdgeDataMatrix.IsDefined())
        {
          m_Mutex.Unlock();
          // adjacency matrix is already cached, return
          return;
        }

        m_EdgeDataMatrix =
          util::OwnPtr< linal::Matrix< t_EdgeData> >
          (
            new linal::Matrix< t_EdgeData>( m_Size, m_Size, m_UnconnectedEdgeValue)
          );

        // get a reference on the matrix
        linal::Matrix< t_EdgeData> &matrix( *m_EdgeDataMatrix);

        // iterate over vertices
        for( size_t vertex_number( 0); vertex_number < m_Size; ++vertex_number)
        {
          // iterate over neighbors and data
          typename t_EdgeDataOfVertex::const_iterator
            itr_data( m_EdgeDataVectors( vertex_number).Begin()),
            itr_data_end( m_EdgeDataVectors( vertex_number).End());

          for
          (
            typename t_EdgeTargetsOfVertex::const_iterator itr_target( m_EdgeTargetVectors( vertex_number).Begin());
            itr_data != itr_data_end;
            ++itr_data, ++itr_target
          )
          {
            // set up the matrix
            matrix( vertex_number, *itr_target) = *itr_data;
          }
        }
        m_Mutex.Unlock();
      }

      //! @brief Uncache complete adjacency matrix to save memory
      void UncacheAdjacencyMatrix() const
      {
        m_Mutex.Lock();
        m_EdgeDataMatrix = util::OwnPtr< linal::Matrix< t_EdgeData> >();
        m_Mutex.Unlock();
      }

      //! @brief checks if two Vertices are connected
      //! @param VERTEX_INDEX_A index of the first Vertex
      //! @param VERTEX_INDEX_B index of the second Vertex
      bool AreConnected( const size_t &VERTEX_INDEX_A, const size_t &VERTEX_INDEX_B) const
      {
        return GetEdgeData( VERTEX_INDEX_A, VERTEX_INDEX_B) != m_UnconnectedEdgeValue;
      }

      //! @brief checks if two Vertices are connected
      //! @param VERTEX_DATA_A first Vertex
      //! @param VERTEX_DATA_B second Vertex
      bool AreVerticesConnected( const t_VertexData &VERTEX_DATA_A, const t_VertexData &VERTEX_DATA_B) const
      {
        // find the index for the first vertex
        size_t vertex_index_a( m_Vertices.Find( VERTEX_DATA_A));

        // find the index for the second vertex
        size_t vertex_index_b( m_Vertices.Find( VERTEX_DATA_B));

        // return false if either passed vertex is not in the vector
        if( vertex_index_a == m_Vertices.GetSize() || vertex_index_b == m_Vertices.GetSize())
        {
          return false;
        }

        // end
        return AreConnected( vertex_index_a, vertex_index_b);
      }

      //! @brief get the subgraph induced by a set of indices
      //! @param SUBGRAPH_INDICES the indices to get the subgraph of
      //! @return the resulting subgraph
      ConstGraph< t_VertexData, t_EdgeData> GetSubgraph( const storage::Vector< size_t> &SUBGRAPH_INDICES) const
      {
        return GetSubgraph( SUBGRAPH_INDICES, util::SiPtr< storage::Vector< size_t> >());
      }

      //! @brief get the subgraph induced by a set of indices
      //! @param SUBGRAPH_INDICES the indices to get the subgraph of
      //! @param STORE_MAPPING stores the mapping of the old indices to the new indices
      //! @return the resulting subgraph
      ConstGraph< t_VertexData, t_EdgeData> GetSubgraph
      (
        const storage::Vector< size_t> &SUBGRAPH_INDICES,
        util::SiPtr< storage::Vector< size_t> > STORE_MAPPING
      ) const
      {
        const size_t new_size( SUBGRAPH_INDICES.GetSize());

        // make sure that all SUBGRAPH_INDICES does not contain any index twice
        BCL_Assert
        (
          new_size == storage::Set< size_t>( SUBGRAPH_INDICES.Begin(), SUBGRAPH_INDICES.End()).GetSize(),
          "Indices given to GetSubgraph must be unique!"
        );

        // create a mapping from old index to new index
        storage::Vector< size_t> old_to_new_index( m_Size, util::GetUndefined< size_t>());
        for( size_t index( 0), size( SUBGRAPH_INDICES.GetSize()); index < size; ++index)
        {
          old_to_new_index( SUBGRAPH_INDICES( index)) = index;
        }

        // make a new vector to store the data
        storage::Vector< t_VertexData> new_vertex_data( new_size);

        // setup the new vertex data
        {
          typename storage::Vector< t_VertexData>::iterator itr_new_data( new_vertex_data.Begin());
          for
          (
            storage::Vector< size_t>::const_iterator itr( SUBGRAPH_INDICES.Begin()), itr_end( SUBGRAPH_INDICES.End());
            itr != itr_end;
            ++itr, ++itr_new_data
          )
          {
            // store the old vertex datum
            *itr_new_data = m_Vertices( *itr);
          }
        }

        if( m_Directed || m_EdgeDataMatrix.IsDefined())
        {
          // make a matrix to store the new edge data
          linal::Matrix< t_EdgeData> new_edge_data( new_size, new_size, m_UnconnectedEdgeValue);

          // store the row we're current at in the new_edge_data matrix
          size_t row( 0);

          // now that old_index_to_new is setup, walk through the edges of each vertex in the subgraph
          for
          (
            storage::Vector< size_t>::const_iterator itr( SUBGRAPH_INDICES.Begin()), itr_end( SUBGRAPH_INDICES.End());
            itr != itr_end;
            ++itr, ++row
          )
          {
            // iterate over vertices
            // iterate over neighbors and data
            typename t_EdgeDataOfVertex::const_iterator
              itr_data( m_EdgeDataVectors( *itr).Begin()),
              itr_data_end( m_EdgeDataVectors( *itr).End());

            for
            (
              typename t_EdgeTargetsOfVertex::const_iterator itr_target( m_EdgeTargetVectors( *itr).Begin());
              itr_data != itr_data_end;
              ++itr_data, ++itr_target
            )
            {
              if( util::IsDefined( old_to_new_index( *itr_target)))
              {
                // set up the matrix
                new_edge_data( row, old_to_new_index( *itr_target)) = *itr_data;
              }
            }
          }

          if( STORE_MAPPING.IsDefined())
          {
            *STORE_MAPPING = old_to_new_index;
          }

          return ConstGraph< t_VertexData, t_EdgeData>
                 (
                   new_vertex_data,
                   new_edge_data,
                   m_UnconnectedEdgeValue,
                   m_Directed,
                   false
                 );
        }
        // make a matrix to store the new edge data
        storage::Vector< UndirectedEdge< t_EdgeData> > new_edge_data;
        const size_t estimated_edges( m_NumEdges * math::Sqr( double( new_vertex_data.GetSize()) / double( m_Size)));
        new_edge_data.AllocateMemory( estimated_edges);

        // store the row we're current at in the new_edge_data matrix
        size_t row( 0);

        // now that old_index_to_new is setup, walk through the edges of each vertex in the subgraph
        for
        (
          storage::Vector< size_t>::const_iterator itr( SUBGRAPH_INDICES.Begin()), itr_end( SUBGRAPH_INDICES.End());
          itr != itr_end;
          ++itr, ++row
        )
        {
          // iterate over vertices
          // iterate over neighbors and data
          typename t_EdgeDataOfVertex::const_iterator
            itr_data( m_EdgeDataVectors( *itr).Begin()),
            itr_data_end( m_EdgeDataVectors( *itr).End());

          for
          (
            typename t_EdgeTargetsOfVertex::const_iterator itr_target( m_EdgeTargetVectors( *itr).Begin());
            itr_data != itr_data_end;
            ++itr_data, ++itr_target
          )
          {
            const size_t new_index( old_to_new_index( *itr_target));
            if( util::IsDefined( new_index) && row <= new_index)
            {
              // set up the matrix
              new_edge_data.PushBack
              (
                UndirectedEdge< t_EdgeData>( row, old_to_new_index( *itr_target), *itr_data)
              );
            }
          }
        }

        if( STORE_MAPPING.IsDefined())
        {
          *STORE_MAPPING = old_to_new_index;
        }

        return ConstGraph< t_VertexData, t_EdgeData>
               (
                 new_vertex_data,
                 new_edge_data,
                 m_UnconnectedEdgeValue
               );
      }

      //! @brief output simple connectivity information for each vertex
      std::string GetBasicConnectivity() const
      {
        std::string return_value;
        // iterate through the vertices
        for( size_t index_a( 0); index_a < m_Size; index_a++)
        {
          return_value += util::Format()( index_a) + " <-> ";
          for
          (
            t_EdgeTargetsOfVertex::const_iterator
              itr_edges( m_EdgeTargetVectors( index_a).Begin()),
              itr_edges_end( m_EdgeTargetVectors( index_a).End());
            itr_edges != itr_edges_end;
            ++itr_edges
          )
          {
            return_value += util::Format()( *itr_edges) + ", ";
          }

          return_value += '\n';
        }

        return return_value;
      }

      //! @brief output simple connectivity information for each vertex
      void ReadBasicConnectivity( std::istream &STREAM)
      {
        storage::Vector< storage::Vector< std::string> > lines( util::SplittedStringLineListFromIStream( STREAM, " <->,"));
        const size_t size( lines.GetSize());
        m_Size = size;
        m_Vertices = storage::Vector< t_VertexData>( size);
        m_EdgeDataVectors = t_EdgeDataContainer( m_Size);
        m_EdgeTargetVectors = t_EdgeTargetContainer( m_Size);
        m_UnconnectedEdgeValue = util::GetUndefined< t_EdgeDatum>();
        m_NumEdges = 0;
        m_Directed = true;

        // iterate through the vertices
        for( size_t index_a( 0); index_a < m_Size; index_a++)
        {
          for
          (
            storage::Vector< std::string>::const_iterator
              itr( lines( index_a).Begin() + 1), itr_end( lines( index_a).End());
            itr != itr_end;
            ++itr
          )
          {
            AddEdge( index_a, util::ConvertStringToNumericalValue< size_t>( *itr), t_EdgeData());
          }
        }
      }

      //! @brief ChangeEdge is like the operator( ROW, COLUMN) of a matrix. Matrix A( FROM, TO) = EDGE_DATA means
      //! @brief to ConstGraph A.ChangeEdge(FROM, TO, EDGE_DATA), while maintaining the undirected invariant if needed
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void ChangeEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        if( EDGE_DATA == m_UnconnectedEdgeValue) // remove any existing edge
        {
          RemoveEdge( FROM, TO);
        }
        else if( GetEdgeData( FROM, TO) == m_UnconnectedEdgeValue) // No existing edge, add a new one
        {
          InternalAddEdge( FROM, TO, EDGE_DATA);
        }
        else // Edge exists and EDGE_DATA is not unconnected.  Edit the current edge.
        {
          InternalEditEdge( FROM, TO, EDGE_DATA);
        }
      }

      //! @brief remove any edge between FROM and TO
      //! @param FROM the index of the from vertex
      void RemoveAllEdges( const size_t &FROM)
      {
        t_EdgeTargetsOfVertex neighbors( GetNeighborIndices( FROM));
        for
        (
          t_EdgeTargetsOfVertex::const_iterator
            itr( neighbors.Begin()), itr_end( neighbors.End());
          itr != itr_end;
          ++itr
        )
        {
          RemoveEdge( FROM, *itr);
        }
      }

      //! @brief remove any edge between FROM and TO
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      void RemoveEdge( const size_t &FROM, const size_t &TO)
      {
        if( GetEdgeData( FROM, TO) != m_UnconnectedEdgeValue)
        {
          const size_t index_in_from( m_EdgeTargetVectors( FROM).Find( TO));
          m_EdgeTargetVectors( FROM).RemoveElements( index_in_from, 1);
          m_EdgeDataVectors( FROM).RemoveElements( index_in_from, 1);
          m_NumEdges--;
          if( m_EdgeDataMatrix.IsDefined())
          {
            SetEdgeDataMatrix( FROM, TO, m_UnconnectedEdgeValue);
          }
          if( IsUndirected())
          {
            if( m_EdgeDataMatrix.IsDefined())
            {
              SetEdgeDataMatrix( TO, FROM, m_UnconnectedEdgeValue);
            }

            const size_t index_in_to( m_EdgeTargetVectors( TO).Find( FROM));
            m_EdgeTargetVectors( TO).RemoveElements( index_in_to, 1);
            m_EdgeDataVectors( TO).RemoveElements( index_in_to, 1);
            m_NumEdges--;
          }
          m_EdgeDataVector = util::OwnPtr< storage::Vector< t_EdgeData> >();
          m_EdgeDataCountMap = util::OwnPtr< storage::Map< t_EdgeData, size_t> >();
        }
      }

      //! @brief add an edge between FROM and TO if one does not already exist.  Change the edge if it already exists
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void AddEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        if( EDGE_DATA != m_UnconnectedEdgeValue && GetEdgeData( FROM, TO) == m_UnconnectedEdgeValue)
        {
          InternalAddEdge( FROM, TO, EDGE_DATA);
          m_EdgeDataVector = util::OwnPtr< storage::Vector< t_EdgeData> >();
          m_EdgeDataCountMap = util::OwnPtr< storage::Map< t_EdgeData, size_t> >();
        }
      }

      //! @brief Edit the edge between FROM and TO if it exists and the new edge is not the unconnected edge value
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void EditEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        if( EDGE_DATA != m_UnconnectedEdgeValue && GetEdgeData( FROM, TO) != m_UnconnectedEdgeValue)
        {
          InternalEditEdge( FROM, TO, EDGE_DATA);
          m_EdgeDataVector = util::OwnPtr< storage::Vector< t_EdgeData> >();
          m_EdgeDataCountMap = util::OwnPtr< storage::Map< t_EdgeData, size_t> >();
        }
      }

      //! @brief Edit the vertex data
      //! @param VERTEX the index of the vertex
      //! @param VERTEX_DATA the data to use
      void EditVertex( const size_t &VERTEX, const t_VertexData &VERTEX_DATA)
      {
        m_Vertices( VERTEX) = VERTEX_DATA;
        m_Mutex.Lock();
        m_VertexDataCountMap = util::OwnPtr< storage::Map< t_VertexData, size_t> >();
        m_Mutex.Unlock();
      }

      //! @brief get all the edges in the graph
      //! @return a vector containing copies of all the edges
      const storage::Vector< t_EdgeData> &GetEdges() const
      {
        m_Mutex.Lock();
        if( !m_EdgeDataVector.IsDefined())
        {
          m_EdgeDataVector = util::OwnPtr< storage::Vector< t_EdgeData> >( new storage::Vector< t_EdgeData>());
          m_EdgeDataVector->AllocateMemory( NumEdges());
          // iterate through the vertices
          for( size_t index_a( 0), size( GetSize()); index_a < size; index_a++)
          {
            m_EdgeDataVector->Append( m_EdgeDataVectors( index_a));
          }
        }
        m_Mutex.Unlock();

        return *m_EdgeDataVector;
      }

      const storage::Map< t_EdgeData, size_t> &GetEdgeTypeCountMap() const
      {
        m_Mutex.Lock();
        if( !m_EdgeDataCountMap.IsDefined())
        {
          m_EdgeDataCountMap = util::OwnPtr< storage::Map< t_EdgeData, size_t> >( new storage::Map< t_EdgeData, size_t>());

          // iterate through the vertices
          for( size_t index_a( 0), size( GetSize()); index_a < size; ++index_a)
          {
            for
            (
              t_EdgeTargetsOfVertex::const_iterator
                itr_edges( m_EdgeDataVectors( index_a).Begin()),
                itr_edges_end( m_EdgeDataVectors( index_a).End());
              itr_edges != itr_edges_end;
              ++itr_edges
            )
            {
              ( *m_EdgeDataCountMap)[ *itr_edges]++;
            }
          }
        }
        m_Mutex.Unlock();

        return *m_EdgeDataCountMap;
      }

      const storage::Map< t_VertexData, size_t> &GetVertexTypeCountMap() const
      {
        m_Mutex.Lock();
        if( !m_VertexDataCountMap.IsDefined())
        {
          m_VertexDataCountMap = util::OwnPtr< storage::Map< t_VertexData, size_t> >( new storage::Map< t_VertexData, size_t>());

          // iterate through the vertices
          for( size_t index_a( 0), size( GetSize()); index_a < size; ++index_a)
          {
            ( *m_VertexDataCountMap)[ m_Vertices( index_a)]++;
          }
        }
        m_Mutex.Unlock();
        return *m_VertexDataCountMap;
      }

    ///////////////
    // operators //
    ///////////////

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
        io::Serialize::Read( m_Vertices, ISTREAM);
        m_Size = m_Vertices.GetSize();

        io::Serialize::Read( m_EdgeDataVectors, ISTREAM);
        io::Serialize::Read( m_EdgeTargetVectors, ISTREAM);
        io::Serialize::Read( m_UnconnectedEdgeValue, ISTREAM);
        io::Serialize::Read( m_Directed, ISTREAM);
        UncacheAdjacencyMatrix();

        // count the edges in the matrix that are not equal to the undirected edge value
        m_NumEdges = 0;
        for( size_t i( 0), size( m_EdgeDataVectors.GetSize()); i < size; ++i)
        {
          m_NumEdges += m_EdgeDataVectors( i).GetSize();
        }
        //return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write
        io::Serialize::Write( m_Vertices, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_EdgeDataVectors, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_EdgeTargetVectors, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_UnconnectedEdgeValue, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Directed, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the edge between FROM and TO to have data type EDGE_DATA, if the matrix exists
      //! @note  protected because this function assumes no edge between FROM and TO exists
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void SetEdgeDataMatrix( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        if( m_EdgeDataMatrix.IsDefined())
        {
          ( *m_EdgeDataMatrix)( FROM, TO) = EDGE_DATA;
        }
      }

      //! @brief add the edge between FROM and TO, assuming that one exists
      //! @note  protected because this function assumes no edge between FROM and TO exists
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void InternalEditEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        m_EdgeDataVectors( FROM)( m_EdgeTargetVectors( FROM).Find( TO)) = EDGE_DATA;
        SetEdgeDataMatrix( FROM, TO, EDGE_DATA);

        if( IsUndirected())
        {
          m_EdgeDataVectors( TO)( m_EdgeTargetVectors( TO).Find( FROM)) = EDGE_DATA;
          SetEdgeDataMatrix( TO, FROM, EDGE_DATA);
        }
      }

      //! @brief add an edge between FROM and TO
      //! @note  protected because this function assumes no edge between FROM and TO exists
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void InternalAddEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        InternalAddDirectedEdge( FROM, TO, EDGE_DATA);
        if( IsUndirected())
        {
          InternalAddDirectedEdge( TO, FROM, EDGE_DATA);
        }
      }

      //! @brief add a directed edge between FROM and TO
      //! @note  protected because this function assumes no edge between FROM and TO exists
      //! @param FROM the index of the from vertex
      //! @param TO the index of the to vertex
      //! @param EDGE_DATA the new datum to use for the edge
      void InternalAddDirectedEdge( const size_t &FROM, const size_t &TO, const t_EdgeData &EDGE_DATA)
      {
        storage::Vector< size_t> &targets( m_EdgeTargetVectors( FROM));
        size_t insert_pos( 0), size( targets.GetSize());
        for( ; insert_pos < size && targets( insert_pos) < TO; ++insert_pos)
        {
        }
        targets.InsertElement( targets[ insert_pos], TO);
        m_EdgeDataVectors( FROM).InsertElement( m_EdgeDataVectors( FROM)[ insert_pos], EDGE_DATA);
        m_NumEdges++;
        SetEdgeDataMatrix( FROM, TO, EDGE_DATA);
      }

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief MakeRandomUndirectedGraph
      //! @param GRAPH_SIZE number of vertices in the graph
      //! @param EDGE_COLOR_RATIOS a vector of doubles indicating the ratios of each color, e.g. {1.0, 2.0, 4.0}
      //!        should give 4x as edges colored 2 as edges colored 0
      //! @param VERTEX_COLOR_RATIOS the ratio of colors of for vertices
      //! @param DIAGONAL_COLOR the color for edges to/from a vertex to itself (usually 0)
      //! @return a graph with the given edge, vertex, and diagonal edge color ratios
      static ConstGraph MakeRandomUndirectedGraph
      (
        const size_t GRAPH_SIZE,
        const storage::Vector< double> &EDGE_COLOR_RATIOS,
        const storage::Vector< double> &VERTEX_COLOR_RATIOS,
        const size_t &DIAGONAL_COLOR = 0
      )
      {
        const size_t size( GRAPH_SIZE);
        if( size == 0)
        {
          return ConstGraph();
        }

        const size_t total_elements( size * size);
        const size_t total_off_diagonal_elements( total_elements - size);
        const size_t undirected_off_diagonal_elements( total_off_diagonal_elements / 2);

        linal::Matrix< size_t>    edge_colors( size, size, DIAGONAL_COLOR);
        storage::Vector< size_t> vertex_colors( size);
        if( size > 1)
        {
          storage::Vector< size_t> edge_color_count( EDGE_COLOR_RATIOS.GetSize());
          // convert EDGE_COLOR_RATIOS to the expected number of edges
          const double edge_renormalization_factor
          (
            double( undirected_off_diagonal_elements)
            /
            std::max( std::accumulate( EDGE_COLOR_RATIOS.Begin(), EDGE_COLOR_RATIOS.End(), 0.0), std::numeric_limits< double>::epsilon())
          );
          for( size_t i( 0); i < edge_color_count.GetSize(); i++)
          {
            // calculate how many edges to make in color i, err on the side of too many
            // because we will eventually shuffle the array and ignore the excess elements, if there
            // are any.
            edge_color_count( i) = std::ceil( edge_renormalization_factor * EDGE_COLOR_RATIOS( i) + 0.5);
          }

          storage::Vector< size_t> off_diag_edge_colors;
          for( size_t i( 0); i < edge_color_count.GetSize(); i++)
          {
            if( edge_color_count( i))
            {
              off_diag_edge_colors.Append( storage::Vector< size_t>( edge_color_count( i), i));
            }
          }

          off_diag_edge_colors.Shuffle();

          // walk through the matrix, setting each off-diagonal element to the next edge color
          for
          (
            size_t matrix_i( 0), edge_color_i = 0;
            edge_color_i < undirected_off_diagonal_elements;
            matrix_i++
          )
          {
            for( size_t matrix_j( matrix_i + 1); matrix_j < size; matrix_j++, edge_color_i++)
            {
              edge_colors[ matrix_i][ matrix_j] = off_diag_edge_colors( edge_color_i);
              edge_colors[ matrix_j][ matrix_i] = off_diag_edge_colors( edge_color_i);
            }
          }
        }

        {
          // similarly, choose the vertex colors
          storage::Vector< size_t> vertex_color_count( VERTEX_COLOR_RATIOS.GetSize());
          const double vertex_renormalization_factor
          (
            double( size)
            /
            std::max
            (
              std::accumulate( VERTEX_COLOR_RATIOS.Begin(), VERTEX_COLOR_RATIOS.End(), 0.0),
              std::numeric_limits< double>::epsilon()
            )
          );
          for( size_t i( 0); i < VERTEX_COLOR_RATIOS.GetSize(); i++)
          {
            // calculate how many edges to make in color i, err on the side of too many
            // because we will eventually shuffle the array and ignore the excess elements, if there
            // are any.
            vertex_color_count( i) = std::ceil( vertex_renormalization_factor * VERTEX_COLOR_RATIOS( i) + 0.5);
          }

          vertex_colors.Resize( 0);
          for( size_t i( 0); i < vertex_color_count.GetSize(); i++)
          {
            if( vertex_color_count( i))
            {
              vertex_colors.Append( storage::Vector< size_t>( vertex_color_count( i), i));
            }
          }

          vertex_colors.Shuffle();
          vertex_colors.Resize( size);
        }
        return ConstGraph( vertex_colors, edge_colors, DIAGONAL_COLOR);
      }

      //! @brief randomly reorder the graph (used for testing isomorphism classes)
      //! @return the actual reordering performed
      storage::Vector< size_t> Shuffle()
      {
        // Given a symmetric matrix, swap Rows and Columns together into a random order
        // Thus making a random isomorph of the matrix.
        // Returns a vector (new_indices) of the new order of indices
        // such that
        // Matrix M1 = some symmetric matrix
        // Matrix M2 = M1
        // storage::Vector< size_t> new_indices( M1.ShuffleRowColumns())
        // now M1(A,B) == M2(new_indices(A),new_indices(B)) for all A & B < size

        // make a new vector to store the indices in random order
        storage::Vector< size_t> new_indices( m_Size);

        // fill the vector with 0,1,2...m_NumberOfRows-1
        for( size_t index = 0; index < m_Size; ++index)
        {
          new_indices( index) = index;
        }

        if( m_Size <= size_t( 1)) // need at least two elements to shuffle
        {
          return new_indices;
        }

        CacheAdjacencyMatrix();

        // the last element would just get swapped with itself, so stop just before we get to the last element
        for( size_t i( 0), total_size_minus_one( m_Size - 1); i < total_size_minus_one; ++i)
        {
          // pick a new index from here to the end of the array.
          // No need to consider swapping with previous elements, since this index (i)
          // already may have been swapped with one of them.
          const size_t target_index( random::GetGlobalRandom().Random< size_t>( i, total_size_minus_one));
          if( target_index != i) // avoid redundant and potentially dangerous swaps
          {
            std::swap( new_indices( i), new_indices( target_index));
            m_EdgeDataMatrix->SwapRows( target_index, i);
            m_EdgeDataMatrix->SwapCols( target_index, i);
            std::swap( m_Vertices( i), m_Vertices( target_index));
          }
        }

        *this =
          ConstGraph< t_VertexData, t_EdgeData>
          (
            m_Vertices,
            *m_EdgeDataMatrix,
            m_UnconnectedEdgeValue,
            m_Directed,
            false
          );

        return new_indices;
      }

    }; // template class Graph

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_CONST_GRAPH_H_

