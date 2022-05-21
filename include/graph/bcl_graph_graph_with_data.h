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

#ifndef BCL_GRAPH_GRAPH_WITH_DATA_H_
#define BCL_GRAPH_GRAPH_WITH_DATA_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_vertex_with_data.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GraphWithData
    //! @brief A graph with a unique data on each vertex (t_VertexDataType) along with edges that can also store data
    //!
    //! @tparam t_VertexDataType data stored in each vertex, needs to be unique to for each vertex
    //! @tparam t_EdgeDataType   data stored in each edge
    //!
    //! @details NOTE: before this class is used, the respective s_Instance member variables of Vertex and Edge need to
    //! be initialized for each datatype that you will be using the class template Graph with. This is demonstrated in
    //! the example file for this class. If this is not done, outputting the Graph via BCLMessage will produce the
    //! following error:
    //!
    //! <span style="color:red">
    //! =crt=bcl::util=> there is no such enum with the name |bcl::graph::Vertex|
    //! <br>=crt=bcl::util=> entry with name: "bcl::graph::Vertex" does not exist => reading to ShPtr will be
    //! impossible
    //! </span>
    //!
    //! @see @link example_graph_graph_with_data.cpp @endlink
    //! @author selicd, karakam, weinerbe
    //! @date Dec 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class GraphWithData :
      public util::ObjectInterface
    {

    public:

      // Typedefs to reduce external dependence on data structures
      typedef VertexWithData< t_VertexDataType, t_EdgeDataType> VertexType;
      typedef typename VertexType::EdgeType                     EdgeType;
      typedef util::ShPtrVector< VertexType>                    VertexContainerType;
      typedef typename VertexType::EdgeContainerType            EdgeContainerType;

    private:

    //////////
    // data //
    //////////

      //! list of all vertices in the graph
      VertexContainerType m_Vertices;

      //! true if the graph directed, otherwise false
      bool m_Directed;

    public:

      typedef storage::Map< util::SiPtr< const VertexType>, size_t> VertexToIndexMap;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, creates an empty graph
      GraphWithData( const bool DIRECTED = false) :
        m_Vertices(),
        m_Directed( DIRECTED)
      {
      }

      //! @brief construct from a list of edges and directed
      GraphWithData
      (
        const VertexContainerType &VERTICES,
        const bool DIRECTED = false
      ) :
        m_Vertices( VERTICES),
        m_Directed( DIRECTED)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Graph
      GraphWithData< t_VertexDataType, t_EdgeDataType> *Clone() const
      {
        return new GraphWithData< t_VertexDataType, t_EdgeDataType>( *this);
      }

      //! @brief creates a new GraphWithData by hard-copying the vertices
      //! @return a new GraphWithData
      GraphWithData< t_VertexDataType, t_EdgeDataType> *HardCopy() const
      {
        return new GraphWithData< t_VertexDataType, t_EdgeDataType>( m_Vertices.HardCopy(), m_Directed);
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

      //! @brief return a const reference to the shared pointer list containing the vertices
      //! @return a const reference to the shared pointer list containing the vertices
      const VertexContainerType &GetVertices() const
      {
        return m_Vertices;
      }

      //! @brief return a non-const reference to the shared pointer list containing the vertices
      //! @return a non-const reference to the shared pointer list containing the vertices
      VertexContainerType &GetVertices()
      {
        return m_Vertices;
      }

      //! @brief return number of vertices
      //! @return number of vertices
      size_t GetNumberVertices() const
      {
        return m_Vertices.GetSize();
      }

      //! @brief return the number of vertices that have the requested number of edges
      //! @param NUMBER_OF_EDGES number of edges that requested vertices will contain
      //! @return number vertices that have the requested number of edges
      size_t GetNumberVerticesWithNEdges( const size_t NUMBER_OF_EDGES) const
      {
        // initialize counter
        size_t counter( 0);

        // iterate through the vertices
        for
        (
          typename VertexContainerType::const_iterator vertex_itr( m_Vertices.Begin()), vertex_itr_end( m_Vertices.End());
          vertex_itr != vertex_itr_end; ++vertex_itr
        )
        {
          // if the number of edges matches the requested size
          if( ( *vertex_itr)->GetDegree() == NUMBER_OF_EDGES)
          {
            ++counter;
          }
        }

        // end
        return counter;
      }

      //! @brief return all vertices that have the requested number of edges
      //! @param NUMBER_OF_EDGES number of edges that requested vertices will contain
      //! @return SiPtrVector of vertices that have the requested number of edges
      util::SiPtrVector< const VertexType> GetVerticesWithNEdges( const size_t NUMBER_OF_EDGES) const
      {
        // initialize the SiPtrVector
        util::SiPtrVector< const VertexType> vertices;

        // iterate through the vertices
        for
        (
          typename VertexContainerType::const_iterator vertex_itr( m_Vertices.Begin()), vertex_itr_end( m_Vertices.End());
          vertex_itr != vertex_itr_end; ++vertex_itr
        )
        {
          // if the number of edges matches the requested size
          if( ( *vertex_itr)->GetDegree() == NUMBER_OF_EDGES)
          {
            // add to the vector
            vertices.PushBack( **vertex_itr);
          }
        }

        // end
        return vertices;
      }

      //! @brief return all vertex data in a SiPtrVector
      //! @return SiPtrVector that contains all vertices data
      util::SiPtrVector< const t_VertexDataType> GetAllVertexData() const
      {
        util::SiPtrVector< const t_VertexDataType> vertex_data_vector; // initialize 

        // iterate over vertices
        for
        (
          typename util::ShPtrVector< VertexType>::const_iterator
            vertex_itr( m_Vertices.Begin()), vertex_itr_end( m_Vertices.End());
          vertex_itr != vertex_itr_end; ++vertex_itr
        )
        {
          vertex_data_vector.PushBack( util::ToSiPtr( ( *vertex_itr)->GetData())); // pushback into vertex_data_vector
        }

        return vertex_data_vector;
      }
      
      //! @brief return all vertex pairs between edges exist
      //! @return Vector of VectorND which contains the verteces that are connected
      storage::Vector< storage::VectorND< 2, VertexType> > GetEdgeVertices() const
      {
        storage::Vector< storage::VectorND< 2, VertexType> > edge_verteces;
        
        for
        ( 
          typename VertexContainerType::const_iterator itr_a( GetVertices().Begin()), itr_end( GetVertices().End());
          itr_a != itr_end; ++itr_a
        )
        {
          for( typename VertexContainerType::const_iterator itr_b( itr_a + 1); itr_b != itr_end; ++itr_b)
          {
            const VertexType &vertex_a( **itr_a), &vertex_b( **itr_b);
            const EdgeType &edge( vertex_a.FindEdge( vertex_b));
                
            if( edge.IsDefined())
            {
              edge_verteces.PushBack( storage::VectorND< 2, VertexType>( vertex_a, vertex_b));
            }
          }
        }
        
        return edge_verteces;
      }

      //! @brief return whether the graph is directed
      //! @return whether the graph is directed
      bool IsDirected() const
      {
        return m_Directed;
      }

      //! @brief checks that each vertex has at least one edge. If the graph has isolated subgraphs,
      //!        this will still return true
      //! @return whether each vertex has at least one edge
      bool IsConnected() const
      {
        // initialize bool to check whether any vertices are unconnected
        bool is_connected( true);

        // if the graph has more than one vertex
        if( m_Vertices.GetSize() > 1)
        {
          // iterate through the graph vertices
          for
          (
             typename VertexContainerType::const_iterator vertex_itr( m_Vertices.Begin()),
               vertex_itr_end( m_Vertices.End());
            vertex_itr != vertex_itr_end; ++vertex_itr
          )
          {
            // if the number of edges is zero, set the bool to false and break
            if( ( *vertex_itr)->GetDegree() == 0)
            {
              is_connected = false;
              break;
            }
          }
        }

        // end
        return is_connected;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns a shared pointer to a vertex containing a specified instance of data if it exists
      //! @param VERTEX_DATA vertex data to search for
      //! @return ShPtr to the VertexWithData
      const util::ShPtr< VertexType> &FindVertex( const t_VertexDataType &VERTEX_DATA) const
      {
        // a NULL pointer to return if nothing is found
        static util::ShPtr< VertexType> sp_undefined_vertex;

        // iterate through the vertex list in search of the designated value
        for
        (
          typename VertexContainerType::const_iterator itr( m_Vertices.Begin()), itr_end( m_Vertices.End());
          itr != itr_end;
          ++itr
        )
        {
          // if the current vertex contains the data being sought, return it
          if( ( *itr)->GetData() == VERTEX_DATA)
          {
            return *itr;
          }
        }

        // if not found return undefined vertex
        return sp_undefined_vertex;
      }

      //! @brief return if the given vertex is already stored
      //! @param VERTEX Vertex of interest
      //! @return if the given vertex is already stored
      bool HasVertex( const VertexType &VERTEX) const
      {
        return FindVertex( VERTEX.GetData()).IsDefined();
      }

      //! @brief check if a Vertex with the specified value exits
      //! @param VERTEX_DATA value of the Vertex
      //! @return true if vertex was found, false otherwise
      bool HasVertex( const t_VertexDataType &VERTEX_DATA) const
      {
        return FindVertex( VERTEX_DATA).IsDefined();
      }

      //! @brief add the specified vertex
      //! @param SP_VERTEX the vertex to add
      //! @return whether insertion was successful
      bool AddVertex( const util::ShPtr< VertexType> &SP_VERTEX)
      {
        // if the given vertex does not exist
        if( !HasVertex( *SP_VERTEX))
        {
          // add the vertex to the vertex list
          m_Vertices.PushBack( SP_VERTEX);

          // return success
          return true;
        }

        // return failure
        return false;
      }

      //! @brief add the specified vertex
      //! @param VERTEX_DATA the vertex to add
      //! @return //! @return whether insertion was successful
      bool AddVertex( const t_VertexDataType &VERTEX_DATA)
      {
        // if the vertex already exists, return false but do nothing else
        if( HasVertex( VERTEX_DATA))
        {
          return false;
        }

        // if the vertex list is empty, also assign the vertex to m_Root
        util::ShPtr< VertexType> new_vertex( new VertexType( VERTEX_DATA));

        // add the vertex to the vertex list
        m_Vertices.PushBack( new_vertex);

        return true;
      }

      //! @brief get the edge connecting two vertices
      //! @param VERTEX_A first Vertex
      //! @param VERTEX_B the second Vertex
      const EdgeType &FindEdge
      (
        const VertexType &VERTEX_A,
        const VertexType &VERTEX_B
      ) const
      {
        // find the corresponding edge
        return VERTEX_A.FindEdge( VERTEX_B);
      }

      //! @brief return the number of Edges
      size_t GetNumberEdges() const
      {
        // initialize edge number counter
        size_t nr_edges( 0);

        // iterate through the Vertices
        for
        (
          typename VertexContainerType::const_iterator itr( m_Vertices.Begin()), itr_end( m_Vertices.End());
          itr != itr_end; ++itr
        )
        {
          // update the counter
          nr_edges += ( *itr)->GetEdges().GetSize();
        }

        // end
        return nr_edges;
      }

      //! @brief Add an edge
      //! @param VERTEX_DATA_A start vertex
      //! @param VERTEX_DATA_B end vertex
      //! @param EDGE_DATA edge data
      //! @return true if succeeded, false if failed
      bool AddEdge
      (
        const t_VertexDataType &VERTEX_DATA_A,
        const t_VertexDataType &VERTEX_DATA_B,
        const t_EdgeDataType &EDGE_DATA
      )
      {
        // get shared pointers to the two vertices because it's inefficient to call FindVertex a million times
        util::ShPtr< VertexType> vertex_a( FindVertex( VERTEX_DATA_A));
        util::ShPtr< VertexType> vertex_b( FindVertex( VERTEX_DATA_B));

        // check if vertices exist or are already Connected
        if( !vertex_a.IsDefined() || !vertex_b.IsDefined() || vertex_a->FindEdge( VERTEX_DATA_B).IsDefined())
        {
          return false;
        }

        // check if the edge already exists, provided that the vertex contains edges in the first place
        // add the edge
        vertex_a->AddEdge( vertex_b, EDGE_DATA);

        // if the graph is an undirected graph, also insert an edge in the other direction
        if( !m_Directed)
        {
          vertex_b->AddEdge( vertex_a, EDGE_DATA);
        }

        return true;
      }

      //! @brief delete an Edge between two Vertices
      //! @param VERTEX_DATA_A value of the origin Vertex
      //! @param VERTEX_DATA_B value of the destination Vertex
      //! @return whether deletion succeeded
      bool DeleteEdge( const t_VertexDataType &VERTEX_DATA_A, const t_VertexDataType &VERTEX_DATA_B)
      {
        // Only try to delete the edge if the vertices are connected
        if( AreVerticesConnected( VERTEX_DATA_A, VERTEX_DATA_B))
        {
          // find the first vertex
          util::ShPtr< VertexType> sp_vertex_a( FindVertex( VERTEX_DATA_A));
          sp_vertex_a->DeleteEdge( VERTEX_DATA_B);

          // if the Graph is undirected, delete the other shared_ptr of the edge
          if( !m_Directed)
          {
            util::ShPtr< VertexType> sp_vertex_b( FindVertex( VERTEX_DATA_B));
            sp_vertex_b->DeleteEdge( VERTEX_DATA_A);
          }
          // deletions are complete so return true
          return true;
        }
        // if the vertices were not connected to start with, then return false
        return false;
      }

      //! @brief removes a Vertex from the Graph, also removing all Edges to and from it
      //! @param VERTEX a Vertex needing deletion
      void DeleteVertex( const util::ShPtr< VertexType> &VERTEX)
      {
        for
        (
          typename VertexContainerType::iterator
            itr( m_Vertices.Begin()),
            itr_end( m_Vertices.End());
          itr != itr_end;
          ++itr
        )
        {
          ( *itr)->DeleteEdge( VERTEX->GetData());
        }
        m_Vertices.Remove( m_Vertices.Begin() + m_Vertices.GetIndex( VERTEX));
      } // DeleteVertex

      //! @brief checks if two Vertices are connected
      //! @param VERTEX_DATA_A value of the first Vertex
      //! @param VERTEX_DATA_B value of the second Vertex
      bool AreVerticesConnected( const t_VertexDataType &VERTEX_DATA_A, const t_VertexDataType &VERTEX_DATA_B) const
      {
        // try to locate both vertices with the given data
        const util::ShPtr< VertexType> vertex_a( FindVertex( VERTEX_DATA_A));
        const util::ShPtr< VertexType> vertex_b( FindVertex( VERTEX_DATA_B));

        // return if they exist and they are actually connected
        return vertex_a.IsDefined() && vertex_b.IsDefined() && vertex_a->FindEdge( VERTEX_DATA_B).IsDefined();
      }

      //! @brief deletes all vertices that have no edges
      void DeleteUnconnectedVertices()
      {
        // initialize new vertex vector
        VertexContainerType vertices;

        // iterate over the vertices
        for
        (
          typename VertexContainerType::const_iterator itr( m_Vertices.Begin()), itr_end( m_Vertices.End());
          itr != itr_end; ++itr
        )
        {
          // if the vertex has edges
          if( ( *itr)->GetDegree() > 0)
          {
            // add it
            vertices.PushBack( *itr);
          }
        }

        // update the vertices
        m_Vertices = vertices;
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
        // read in whether the graph is directed
        io::Serialize::Read( m_Directed, ISTREAM);

        // read in the basic vertex data
        io::Serialize::Read( m_Vertices, ISTREAM);

        //return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write whether the graph is directed
        io::Serialize::Write( m_Directed, OSTREAM, INDENT);

        // write the basic vertex data
        io::Serialize::Write( m_Vertices, OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class Graph

    // instantiate s_Instance
    template< typename t_VertexDataType, typename t_EdgeDataType>
    const util::SiPtr< const util::ObjectInterface> GraphWithData< t_VertexDataType, t_EdgeDataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new GraphWithData< t_VertexDataType, t_EdgeDataType>())
    );
  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_GRAPH_WITH_DATA_H_

