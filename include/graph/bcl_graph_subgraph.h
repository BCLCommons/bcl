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

#ifndef BCL_GRAPH_SUBGRAPH_H_
#define BCL_GRAPH_SUBGRAPH_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_const_graph.h"
#include "bcl_graph_undirected_edge.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Subgraph
    //! @brief provides special methods for handling graphs that are subsets of other graphs
    //! @param t_VertexData - the type of data associated with each vertex. Must be comparable via ==
    //! @param t_EdgeData   - the type of data associated with each edge.  Must be comparable via ==
    //!
    //! @see @link example_graph_subgraph.cpp @endlink
    //! @author mendenjl
    //! @date Jan 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class Subgraph :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > m_Graph;     //!< the complete graph
      storage::Vector< size_t>                              m_VertexIds; //!< the vertex indices of the complete graph that are part of the subgraph

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Subgraph()
      {
      }

      //! @brief construct from a graph and selected vertices
      //! @param GRAPH the graph to consider
      //! @param VERTEX_INDICES the vertex indices of the graph that are part of the subgraph
      Subgraph
      (
        const util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > &GRAPH,
        const storage::Vector< size_t> &VERTEX_IDS
      ) :
        m_Graph( GRAPH),
        m_VertexIds( VERTEX_IDS)
      {
        // check that all vertex ids are < the size of the graph
        const size_t graph_size( GetSizeOfParentGraph());

        if( m_VertexIds.IsEmpty())
        {
          return;
        }

        // check that indices are valid
        const size_t max_index( *std::max_element( m_VertexIds.Begin(), m_VertexIds.End()));
        BCL_Assert
        (
          max_index < graph_size,
          "vertex index " + util::Format()( max_index) + " from graph with only " + util::Format()( graph_size) +
          " indices"
        );
      }

      //! @brief construct from a graph and selected vertices
      //! @param GRAPH the graph to consider
      //! @param VERTEX_INDICES the vertex indices of the graph that are part of the subgraph
      Subgraph
      (
        const util::OwnPtr< ConstGraph< t_VertexData, t_EdgeData> > &GRAPH,
        const storage::Vector< size_t> &VERTEX_IDS
      ) :
        m_Graph( GRAPH),
        m_VertexIds( VERTEX_IDS)
      {
        // check that all vertex ids are < the size of the graph
        const size_t graph_size( GetSizeOfParentGraph());

        if( m_VertexIds.IsEmpty())
        {
          return;
        }

        // check that indices are valid
        const size_t max_index( *std::max_element( m_VertexIds.Begin(), m_VertexIds.End()));
        BCL_Assert
        (
          max_index < graph_size,
          "vertex index " + util::Format()( max_index) + " from graph with only " + util::Format()( graph_size) +
          " indices"
        );
      }

      //! clone the object
      Subgraph *Clone() const
      {
        return new Subgraph( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetSize returns the number of vertices in the subgraph
      //! @return the number of vertices in the subgraph
      size_t GetSize() const
      {
        return m_VertexIds.GetSize();
      }

      //! @brief Get a pointer to the internally held graph
      const util::OwnPtr< const ConstGraph< t_VertexData, t_EdgeData> > &GetParentGraphPtr() const
      {
        return m_Graph;
      }

      //! @brief Given the indices of the vertex involved in the subgraph
      const storage::Vector< size_t> &GetVertexIndices() const
      {
        return m_VertexIds;
      }

      //! @brief Set the indices of the vertex involved in the subgraph
      void SetVertexIndices( const storage::Vector< size_t> &VERTICES)
      {
        m_VertexIds = VERTICES;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a new graph from the vertex ids
      //! @return the a new graph generated from this subgraph
      ConstGraph< t_VertexData, t_EdgeData> ToGraph() const
      {
        return m_Graph.IsDefined() ? m_Graph->GetSubgraph( m_VertexIds) : ConstGraph< t_VertexData, t_EdgeData>();
      }

      //! @brief construct a new graph from the vertex ids
      //! @return the a new graph generated from this subgraph
      operator ConstGraph< t_VertexData, t_EdgeData>() const
      {
        return ToGraph();
      }

      //! @brief construct a subgraph containing the vertices that are not in this subgraph, and the edges between them
      //! @return a subgraph containing the vertices that are not in this subgraph, and the edges between them
      Subgraph< t_VertexData, t_EdgeData> GetComplement() const
      {
        storage::Vector< size_t> complementary_ids;
        complementary_ids.AllocateMemory( GetSizeOfParentGraph() - m_VertexIds.GetSize());

        const storage::Set< size_t> current_id_set( m_VertexIds.Begin(), m_VertexIds.End());

        // store all the vertices that are not in this subgraph
        for( size_t i( 0), size( GetSizeOfParentGraph()); i < size; ++i)
        {
          if( !current_id_set.Contains( i))
          {
            complementary_ids.PushBack( i);
          }
        }

        // get that subgraph of a
        return Subgraph< t_VertexData, t_EdgeData>( m_Graph, complementary_ids);
      }

      //! @brief construct list of undirected edges in the subgraph
      //! @return the list of undirected edges in the subgraph
      storage::Vector< UndirectedEdge< t_EdgeData> > GetUndirectedEdges() const
      {
        BCL_Assert( !m_Graph->IsDirected(), "Cannot get an undirected edge from a directed graph");

        // make a vector that will hold 1 if a vertex from graph was in the subgraph
        storage::Vector< size_t> in_subgraph( GetSizeOfParentGraph(), size_t( 0));

        for
        (
          storage::Vector< size_t>::const_iterator itr( m_VertexIds.Begin()), itr_end( m_VertexIds.End());
          itr != itr_end;
          ++itr
        )
        {
          in_subgraph( *itr) = 1;
        }

        // make a list to hold the edges in the subgraph
        storage::List< UndirectedEdge< t_EdgeData> > subgraph_edges;

        for
        (
          storage::Vector< size_t>::const_iterator itr_vertex( m_VertexIds.Begin()), itr_vertex_end( m_VertexIds.End());
          itr_vertex != itr_vertex_end;
          ++itr_vertex
        )
        {
          typename storage::Vector< t_EdgeData>::const_iterator itr_data( m_Graph->GetNeighborData( *itr_vertex).Begin());
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( m_Graph->GetNeighborIndices( *itr_vertex).Begin()),
              itr_neighbors_end( m_Graph->GetNeighborIndices( *itr_vertex).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors, ++itr_data
          )
          {
            // if both vertices are in the subgraph
            if( in_subgraph( *itr_neighbors) == 1 && *itr_vertex < *itr_neighbors)
            {
              // store the edge
              subgraph_edges.PushBack( UndirectedEdge< t_EdgeData>( *itr_vertex, *itr_neighbors, *itr_data));
            }
          }
        }

        return storage::Vector< UndirectedEdge< t_EdgeData> >( subgraph_edges.Begin(), subgraph_edges.End());
      }

      //! @brief construct list of edges (by indices) in the subgraph
      //! @return the list of edges (by indices) in the subgraph
      storage::List< storage::Pair< size_t, size_t> > GetEdgeIndices() const
      {
        // make a vector that will hold 1 if a vertex from graph was in the subgraph
        storage::Vector< size_t> in_subgraph( GetSizeOfParentGraph(), size_t( 0));

        for
        (
          storage::Vector< size_t>::const_iterator itr( m_VertexIds.Begin()), itr_end( m_VertexIds.End());
          itr != itr_end;
          ++itr
        )
        {
          in_subgraph( *itr) = 1;
        }

        // make a list to hold the edges in the subgraph
        storage::List< storage::Pair< size_t, size_t> > subgraph_edges;

        for
        (
          storage::Vector< size_t>::const_iterator itr_vertex( m_VertexIds.Begin()), itr_vertex_end( m_VertexIds.End());
          itr_vertex != itr_vertex_end;
          ++itr_vertex
        )
        {
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( m_Graph->GetNeighborIndices( *itr_vertex).Begin()),
              itr_neighbors_end( m_Graph->GetNeighborIndices( *itr_vertex).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            // if both vertices are in the subgraph
            if( in_subgraph( *itr_neighbors) == 1)
            {
              // store the edge
              subgraph_edges.PushBack( storage::Pair< size_t, size_t>( *itr_vertex, *itr_neighbors));
            }
          }
        }

        return subgraph_edges;
      }

      //! @brief construct a vector containing the ids of vertices whose edges all stay within the subgraph
      //! @return a vector containing the ids of vertices whose edges all stay within the subgraph
      storage::Vector< size_t> GetIdsOfInteriorVertices() const
      {
        // make a vector that holds 1 if a vertex from graph was in the subgraph
        storage::Vector< size_t> in_subgraph( GetSizeOfParentGraph(), size_t( 0));

        for
        (
          storage::Vector< size_t>::const_iterator itr( m_VertexIds.Begin()), itr_end( m_VertexIds.End());
          itr != itr_end;
          ++itr
        )
        {
          in_subgraph( *itr) = 1;
        }

        // make a list to hold the ids of all vertices whose edges remain within the subgraph
        storage::Vector< size_t> enclosed_vertices;
        enclosed_vertices.AllocateMemory( m_VertexIds.GetSize());
        for
        (
          storage::Vector< size_t>::const_iterator itr_vertex( m_VertexIds.Begin()), itr_vertex_end( m_VertexIds.End());
          itr_vertex != itr_vertex_end;
          ++itr_vertex
        )
        {
          bool is_fully_enclosed( true);
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( m_Graph->GetNeighborIndices( *itr_vertex).Begin()),
              itr_neighbors_end( m_Graph->GetNeighborIndices( *itr_vertex).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            // if both vertices are in the subgraph
            if( in_subgraph( *itr_neighbors) != 1)
            {
              // then the vertex is not fully enclosed, so we needn't look any farther
              is_fully_enclosed = false;
              break;
            }
          }
          if( is_fully_enclosed)
          {
            enclosed_vertices.PushBack( *itr_vertex);
          }
        }

        return enclosed_vertices;
      }

      //! @brief construct a list of edges (by indices) that go between the subgraph and the supergraph
      //! @return the list of edges (by indices) that go between the subgraph and the supergraph, with vertices in the subgraph
      //!         listed first
      storage::List< storage::Pair< size_t, size_t> > GetAdjacentEdgeIndices() const
      {
        const size_t graph_size( GetSizeOfParentGraph());

        // make a vector that will hold 1 iff a vertex from the parent graph is in the subgraph
        storage::Vector< size_t> in_subgraph( graph_size, size_t( 0));

        for
        (
          storage::Vector< size_t>::const_iterator itr( m_VertexIds.Begin()), itr_end( m_VertexIds.End());
          itr != itr_end;
          ++itr
        )
        {
          in_subgraph( *itr) = 1;
        }

        // make a list to hold the edges that have only one vertex in the subgraph
        storage::List< storage::Pair< size_t, size_t> > adjacent_edges;

        for( size_t vertex_number( 0); vertex_number < graph_size; ++vertex_number)
        {
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( m_Graph->GetNeighborIndices( vertex_number).Begin()),
              itr_neighbors_end( m_Graph->GetNeighborIndices( vertex_number).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            // if one vertex was in the subgraph but the other was not
            if( in_subgraph( vertex_number) != in_subgraph( *itr_neighbors))
            {
              if( in_subgraph( vertex_number))
              {
                // store the edge
                adjacent_edges.PushBack( storage::Pair< size_t, size_t>( vertex_number, *itr_neighbors));
              }
            }
          }
        }

        return adjacent_edges;
      }

      //! @brief construct a list of edges (by indices) that go between the subgraph and the supergraph
      //! @return a vector of vectors. Outer index is the index into the subgraph, inner index is index of the original graph
      storage::Vector< storage::Vector< size_t> > GetOrderedAdjacentEdgeIndices() const
      {
        const size_t graph_size( GetSizeOfParentGraph());

        // make a vector that will hold 1 iff a vertex from the parent graph is in the subgraph
        storage::Vector< size_t> in_subgraph( graph_size, size_t( 0));

        for
        (
          storage::Vector< size_t>::const_iterator itr( m_VertexIds.Begin()), itr_end( m_VertexIds.End());
          itr != itr_end;
          ++itr
        )
        {
          in_subgraph( *itr) = 1;
        }
        storage::Vector< storage::Vector< size_t> > adjacent_edges( m_VertexIds.GetSize());

        size_t i( 0);
        for
        (
          storage::Vector< size_t>::const_iterator itr_vertex( m_VertexIds.Begin()), itr_vertex_end( m_VertexIds.End());
          itr_vertex != itr_vertex_end;
          ++itr_vertex, ++i
        )
        {
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_neighbors( m_Graph->GetNeighborIndices( *itr_vertex).Begin()),
              itr_neighbors_end( m_Graph->GetNeighborIndices( *itr_vertex).End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            // if one vertex was in the subgraph but the other was not
            if( !in_subgraph( *itr_neighbors))
            {
              // store the edge
              adjacent_edges( i).PushBack( *itr_neighbors);
            }
          }
        }
        return adjacent_edges;
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
        io::Serialize::Read( m_Graph    , ISTREAM);
        io::Serialize::Read( m_VertexIds, ISTREAM);

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
        io::Serialize::Write( m_Graph    , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_VertexIds, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the size of the parent graph, if defined
      //! @return the size of the parent graph, if defined
      size_t GetSizeOfParentGraph() const
      {
        return m_Graph.IsDefined() ? m_Graph->GetSize() : 0;
      }

    }; // template class Graph

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_SUBGRAPH_H_

