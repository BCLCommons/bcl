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

#ifndef BCL_GRAPH_EDGE_WITH_DATA_H_
#define BCL_GRAPH_EDGE_WITH_DATA_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_vertex_with_data.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EdgeWithData
    //! @brief is a template class that represents an edge in a graph that can contain a value ( color)
    //! @details Each Edge is connected to a single Vertex and can be colored.
    //!
    //! @tparam t_VertexDataType type of the data associated with the target vertex
    //! @tparam t_EdgeDataType   type of the data associated with this edge
    //!
    //! @see @link example_graph_edge_with_data.cpp @endlink
    //! @author selicd, karakam, weinerbe
    //! @date 12/04/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class EdgeWithData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! pointer to the vertex the edge is directed toward
      util::SiPtr< const VertexWithData< t_VertexDataType, t_EdgeDataType> > m_Target;

      //! the data associated with the edge
      t_EdgeDataType m_Data;

    public:

      //! @brief return reference to an undefined EdgeWithData
      //! @return reference to an undefined EdgeWithData
      static const EdgeWithData< t_VertexDataType, t_EdgeDataType> &GetUndefined()
      {
        // initialize static instance
        static const EdgeWithData< t_VertexDataType, t_EdgeDataType> s_undefined_edge;

        // end
        return s_undefined_edge;
      }

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, needed to define s_Instance
      EdgeWithData() :
        m_Target(),
        m_Data()
      {
      }

      //! @brief constructor that takes a connected Vertex and an optional color
      //! @param SP_VERTEX the connected Vertex
      //! @param EDGE_DATA Data associated with the edge
      EdgeWithData
      (
        const util::ShPtr< VertexWithData< t_VertexDataType, t_EdgeDataType> > &SP_VERTEX,
        const t_EdgeDataType &EDGE_DATA
      ) :
        m_Target( SP_VERTEX),
        m_Data( EDGE_DATA)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Edge
      EdgeWithData< t_VertexDataType, t_EdgeDataType> *Clone() const
      {
        return new EdgeWithData< t_VertexDataType, t_EdgeDataType>( *this);
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

      //! @brief return a reference to edge data
      //! @return reference to edge data
      const t_EdgeDataType &GetData() const
      {
        return m_Data;
      }

      //! @brief return const reference to the vertex the edge is directed toward
      //! @return const reference to the vertex the edge is directed toward
      const util::SiPtr< const VertexWithData< t_VertexDataType, t_EdgeDataType> > &GetTarget() const
      {
        return m_Target;
      }

//      //! @brief return non-const reference to the vertex the edge is directed toward
//      //! @return non-const reference to the vertex the edge is directed toward
//      util::SiPtr< VertexWithData< t_VertexDataType, t_EdgeDataType> > &GetTarget()
//      {
//        return m_Target;
//      }

      //! @brief returns whether this edge is defined
      //! @return whether this edge is defined
      bool IsDefined() const
      {
        return m_Target.IsDefined();
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
        // an edge is meaningless without its vertices, but edges cannot "own" vertices, otherwise
        // there will be memory leaks when we delete the graph with data.  so read/write are unusable here
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // an edge is meaningless without its vertices, but edges cannot "own" vertices, otherwise
        // there will be memory leaks when we delete the graph with data.  so read/write are unused here
        return OSTREAM;
      }

    }; // template class EdgeWithData

    // instantiate s_Instance
    template< typename t_VertexDataType, typename t_EdgeDataType>
    const util::SiPtr< const util::ObjectInterface> EdgeWithData< t_VertexDataType, t_EdgeDataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new EdgeWithData< t_VertexDataType, t_EdgeDataType>())
    );

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_EDGE_WITH_DATA_H_ 
