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

#ifndef BCL_GRAPH_VERTEX_WITH_DATA_H_
#define BCL_GRAPH_VERTEX_WITH_DATA_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_edge_with_data.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VertexWithData
    //! @brief is a template class for representing vertices in a graph that also store an associated data
    //! Each Vertex can be connected via Edges to any number of other Vertices
    //!
    //! @see @link example_graph_vertex_with_data.cpp @endlink
    //! @author selicd, karakam, weinerbe
    //! @date 12/04/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class VertexWithData :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

    public:

    //////////////
    // typedefs //
    //////////////

      // so the user does not have to know what type of edge we are using
      typedef EdgeWithData< t_VertexDataType, t_EdgeDataType> EdgeType;

      // so the user does not have to know what type of container these edges are in
      typedef storage::List< EdgeType> EdgeContainerType;

    private:

    //////////
    // data //
    //////////

      //! vertex data
      t_VertexDataType m_Data;

      //! list of connected edges
      EdgeContainerType m_Edges;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, needed to define s_Instance
      VertexWithData() :
        m_Data(),
        m_Edges()
      {
      }

      //! @brief constructor that takes an instance of Vertex data
      //! @param DATA the data contained in the vertex
      VertexWithData( const t_VertexDataType &DATA) :
        m_Data( DATA),
        m_Edges()
      {
      }

      //! @brief Clone function
      //! @return pointer to new Vertex
      VertexWithData< t_VertexDataType, t_EdgeDataType> *Clone() const
      {
        return new VertexWithData< t_VertexDataType, t_EdgeDataType>( *this);
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

      //! @brief return reference to the data contained in the vertex
      //! @return reference to the data contained in the vertex
      const t_VertexDataType &GetData() const
      {
        return m_Data;
      }

      //! @brief return const-reference to the list of edges
      //! @return const-reference to the list of edges
      const EdgeContainerType &GetEdges() const
      {
        return m_Edges;
      }

      //! @brief return non-const reference to the list of edges
      //! @return non-const reference to the list of edges
      EdgeContainerType &GetEdges()
      {
        return m_Edges;
      }

      //! @brief return reference to the degree of the vertex, the number of edges
      //! @return reference to the degree of the vertex, the number of edges
      size_t GetDegree() const
      {
        return m_Edges.GetSize();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a new edge to the vertex
      //! @param VERTEX connecting vertex
      //! @param EDGE_DATA data for the edge
      void AddEdge
      (
        util::ShPtr< VertexWithData< t_VertexDataType, t_EdgeDataType> > &VERTEX,
        const t_EdgeDataType &EDGE_DATA
      )
      {
        // create the new edge and add it to the list of edges
        m_Edges.PushBack( EdgeType( VERTEX, EDGE_DATA));
      }

      //! @brief delete the edge with the specified data
      //! @param VERTEX_DATA data of the vertex to be deleted
      void DeleteEdge( const t_VertexDataType &VERTEX_DATA)
      {
        // iterate over edges
        for
        (
          typename EdgeContainerType::iterator iter( m_Edges.Begin()), iter_end( m_Edges.End());
          iter != iter_end; ++iter
        )
        {
          // if the edge pointing to the vertex with the specified data is found
          if( iter->GetTarget()->GetData() == VERTEX_DATA)
          {
            // delete this edge and return
            m_Edges.Remove( iter);
            return;
          }
        }
      }

      //! @brief delete an edge that is directed towards the specified vertex
      //! @param VERTEX connecting vertex of interest
      void DeleteEdge( const util::ShPtr< VertexWithData< t_VertexDataType, t_EdgeDataType> > &VERTEX)
      {
        // iterate over edges
        for
        (
          typename EdgeContainerType::iterator iter( m_Edges.Begin()), iter_end( m_Edges.End());
          iter != iter_end;
        )
        {
          // if this edge points to the specified vertex
          if( iter->GetTarget() == VERTEX)
          {
            // remove it and assign the iterator
            iter = m_Edges.Remove( iter);
          }
          else
          {
            // otherwise just increment
            ++iter;
          }
        }
      }

      //! @brief find the edge that points to the given VERTEX
      //! @param VERTEX Vertex of interest
      //! @return the edge that points to the given VERTEX
      const EdgeType &FindEdge( const VertexWithData< t_VertexDataType, t_EdgeDataType> &VERTEX) const
      {
        // iterate through the Edges of the Vertex
        for
        (
          typename EdgeContainerType::const_iterator iter( m_Edges.Begin()), iter_end( m_Edges.End());
          iter != iter_end; ++iter
        )
        {
          // if the requested data is found
          if( iter->GetTarget().GetPointer() == &VERTEX)
          {
            // return this edge
            return *iter;
          }
        }

        // if no such edge exists return invalid Edge
        return EdgeType::GetUndefined();
      }

      //! @brief return ShPtr to the edge from this vertex that points to the the vertex with the given data
      //! @param VERTEX_DATA datum in the Vertex to be found
      //! @return requested edge
      const EdgeType &FindEdge( const t_VertexDataType &VERTEX_DATA) const
      {
        // iterate through the Edges of the Vertex
        for
        (
          typename EdgeContainerType::const_iterator iter( m_Edges.Begin()), iter_end( m_Edges.End());
          iter != iter_end; ++iter
        )
        {
          // if the requested data is found
          if( iter->GetTarget()->GetData() == VERTEX_DATA)
          {
            // return this edge
            return *iter;
          }
        }

        // if no such edge exists return invalid Edge
        return EdgeType::GetUndefined();
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      //! @note  edges are NOT read with this functions, otherwise there would be an infinite
      //! @note  loop where each edge prints out each vertex, which prints out each edge, etc.
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Data, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @note  edges are NOT printed with this functions, otherwise there is an infinite
      //! @note  loop where each edge prints out each vertex, which prints out each edge, etc.
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class VertexWithData

    // instantiate s_Instance of Vertex
    template< typename t_VertexDataType, typename t_EdgeDataType>
    const util::SiPtr< const util::ObjectInterface> VertexWithData< t_VertexDataType, t_EdgeDataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VertexWithData< t_VertexDataType, t_EdgeDataType>())
    );

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_VERTEX_WITH_DATA_H_ 

