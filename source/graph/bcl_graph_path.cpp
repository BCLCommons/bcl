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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "graph/bcl_graph_path.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

    //! @brief default constructor, creates an empty path
    Path::Path() :
      m_Undirected( false)
    {
    }

    //! @brief construct from vertices and edges
    //! @param GRAPH_SIZE the size of the graph from which the path is a part
    //! @param VERTICES a vector containing the vertices in the order they were visited
    //! @param UNDIRECTED whether the path is undirected
    //! @brief initialize a path with a certain graph size
    Path::Path
    (
      const size_t &GRAPH_SIZE,
      const storage::Vector< size_t> &VERTICES,
      const bool &UNDIRECTED
    ) :
      m_Vertices( VERTICES),
      m_VertexBitMap( GRAPH_SIZE, '0'),
      m_Undirected( UNDIRECTED)
    {
      // if undirected, then reverse the path if necessary
      if( m_Undirected && m_Vertices.GetSize() > size_t( 1) && m_Vertices.FirstElement() > m_Vertices.LastElement())
      {
        std::reverse( m_Vertices.Begin(), m_Vertices.End());
      }

      // set up the bit map with '1' for every vertex that has been visited
      for
      (
        storage::Vector< size_t>::const_iterator itr_path( VERTICES.Begin()), itr_path_end( VERTICES.End());
        itr_path != itr_path_end;
        ++itr_path
      )
      {
        BCL_Assert
        (
          *itr_path < GRAPH_SIZE,
          " path contained vertex with index (" + util::Format()( *itr_path) +
          ")  greater than graph size (" + util::Format()( GRAPH_SIZE) + " )"
        );

        m_VertexBitMap[ *itr_path] = '1';
      }
    }

    //! @brief create a new path by connecting two existing paths and a mask for which vertices to fuse
    //! @param PATH_A, PATH_B the two paths to combine
    //! @param OVERLAP the path along which to combine the vertices
    //! e.g Path([1,2,3], [0,1,2]) -> 0,1,2,1,2,3
    //!     Path([1,2,3], [0,1,2], [1,2]) -> 0,1,2,3
    //!     Path([1,2,3], [0,1,2], [2]) -> 0,1,1,2,3
    Path::Path( const Path &PATH_A, const Path &PATH_B, const Path &OVERLAP) :
      m_VertexBitMap( PATH_A.m_VertexBitMap),
      m_Undirected( PATH_A.m_Undirected)
    {
      BCL_Assert
      (
        ( PATH_A.m_Undirected && PATH_B.m_Undirected)
        || ( !PATH_A.m_Undirected && !PATH_B.m_Undirected),
        "Cannot combined directed paths with undirected paths"
      );

      // add all the vertices in PATH_B to the vertex bit map as well
      for( size_t path_b_index( 0), path_b_size( PATH_B.GetSize()); path_b_index < path_b_size; ++path_b_index)
      {
        m_VertexBitMap[ PATH_B.m_Vertices( path_b_index)] = '1';
      }

      if( !OVERLAP.GetSize())
      {
        m_Vertices = PATH_A.GetVertices();
        m_Vertices.InsertElements( m_Vertices.End(), PATH_B.GetVertices());
        return;
      }

      // make enough room for the combined path in m_Vertices
      m_Vertices.Resize( PATH_B.GetSize() + PATH_A.GetSize() - OVERLAP.GetSize());

      // get an iterator to the first element
      storage::Vector< size_t>::iterator new_path_itr( m_Vertices.Begin());

      if( PATH_A.StartsWith( OVERLAP))
      {
        // copy in reverse
        new_path_itr = std::copy( PATH_A.ReverseBegin(), PATH_A.ReverseEnd(), new_path_itr);
      }
      else
      {
        // copy forward
        new_path_itr = std::copy( PATH_A.Begin(), PATH_A.End(), new_path_itr);
      }

      if( PATH_B.StartsWith( OVERLAP))
      {
        // copy forward
        std::copy( PATH_B.Begin() + OVERLAP.GetSize(), PATH_B.End(), new_path_itr);
      }
      else
      {
        // copy in reverse
        std::copy( PATH_B.ReverseBegin() + OVERLAP.GetSize(), PATH_B.ReverseEnd(), new_path_itr);
      }
    }

    //! clone the object
    Path *Path::Clone() const
    {
      return new Path( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @return class name of the object behind a pointer or the current object
    const std::string &Path::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @return true if the path is directed
    bool Path::IsDirected() const
    {
      return !m_Undirected;
    }

    //! @return true if the path is undirected
    bool Path::IsUndirected() const
    {
      return m_Undirected;
    }

    //! @brief return the edge ids
    //! @return the edge ids
    const storage::Vector< size_t> &Path::GetVertices() const
    {
      return m_Vertices;
    }

    //! @brief GetSize
    //! @return the number of vertices in the path
    size_t Path::GetSize() const
    {
      return m_Vertices.GetSize();
    }

    //! @brief get the first element of the path
    //! @return the first element of the path
    size_t Path::FirstElement() const
    {
      return m_Vertices.FirstElement();
    }

    //! @brief GetSize
    //! @return the number of vertices in the path
    size_t Path::LastElement() const
    {
      return m_Vertices.LastElement();
    }

    //! @brief get an iterator to the beginning of the path
    //! @return an iterator to the beginning of the path
    Path::const_iterator Path::Begin() const
    {
      return m_Vertices.Begin();
    }

    //! @brief get an iterator to the end of the path
    //! @return an iterator to the end of the path
    Path::const_iterator Path::End() const
    {
      return m_Vertices.End();
    }

    //! @brief get an iterator to the reverse begining of the path
    //! @return an iterator to the reverse begining of the path
    Path::const_reverse_iterator Path::ReverseBegin() const
    {
      return m_Vertices.ReverseBegin();
    }

    //! @brief get an iterator to the reverse end of the path
    //! @return an iterator to the reverse end of the path
    Path::const_reverse_iterator Path::ReverseEnd() const
    {
      return m_Vertices.ReverseEnd();
    }

    //! @brief Test whether the path contains a given vertex
    //! @return true if the path contains a particular vertex
    bool Path::Contains( const size_t &VERTEX) const
    {
      return VERTEX < m_VertexBitMap.size() && m_VertexBitMap[ VERTEX] == '1';
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine if the path ends at VERTEX
    //! @return true if the path ends on VERTEX
    bool Path::EndsAt( const size_t &VERTEX) const
    {
      return GetSize() > 0 && ( VERTEX == m_Vertices.FirstElement() || VERTEX == m_Vertices.LastElement());
    }

    //! @brief count the number of times vertices visited more than once in the path
    //! @return the number of times vertices visited more than once in the path
    size_t Path::CountVertexRepetitions() const
    {
      const size_t number_unique_vertices( std::count( m_VertexBitMap.begin(), m_VertexBitMap.end(), '1'));

      return m_Vertices.GetSize() - number_unique_vertices;
    }

    //! @brief determine if the paths overlap anywhere other than at their ends
    //! @return true if the paths overlap anywhere other than at their ends
    bool Path::Crosses( const Path &PATH) const
    {
      BCL_Assert
      (
        m_VertexBitMap.size() == PATH.m_VertexBitMap.size(),
        "Can't test for crossing of paths from different graphs"
      );

      for( size_t vertex( 1), number_vertices( m_Vertices.GetSize() - 1); vertex < number_vertices; ++vertex)
      {
        if( PATH.m_VertexBitMap[ m_Vertices( vertex)] == '1')
        {
          return true;
        }
      }

      return false;
    }

    //! @brief true if one path connects to another, meaning, the paths contain a common starting or ending vertex
    bool Path::Connects( const Path &PATH) const
    {
      if( IsDirected())
      {
        if
        (
          LastElement() == PATH.FirstElement()
          || FirstElement() == PATH.LastElement()
        )
        {
          return true;
        }
      }
      return EndsAt( PATH.FirstElement()) || EndsAt( PATH.LastElement());
    }

    //! @brief test whether the end of this path is the same as PATH
    //! @param PATH the path to search for
    //! @return whether the path starts with PATH
    //! @note This functions treats both paths as directed
    bool Path::StartsWith( const Path &PATH) const
    {
      if( PATH.GetSize() > GetSize())
      {
        return false;
      }
      return std::equal( PATH.Begin(), PATH.End(), m_Vertices.Begin())
             ||
             ( PATH.IsUndirected() && std::equal( PATH.ReverseBegin(), PATH.ReverseEnd(), m_Vertices.Begin()));
    }

    //! @brief search for a path in this path
    //! @param PATH the path to search for, must be directed
    //! @return whether the path starts with PATH
    //! @note This functions treats both paths as directed
    bool Path::EndsWith( const Path &PATH) const
    {
      if( PATH.GetSize() > GetSize())
      {
        return false;
      }
      return std::equal( PATH.ReverseBegin(), PATH.ReverseEnd(), m_Vertices.ReverseBegin())
             ||
             ( PATH.IsUndirected() && std::equal( PATH.Begin(), PATH.End(), m_Vertices.ReverseBegin()));
    }

    //! @brief determine whether this path covers another path (that is, visits the same vertices)
    //! @return true if this path visits all the vertices in PATH
    bool Path::Covers( const Path &PATH) const
    {
      for( size_t vertex( 0), number_vertices( PATH.GetSize()); vertex < number_vertices; ++vertex)
      {
        if( !Contains( PATH.m_Vertices( vertex)))
        {
          return false;
        }
      }

      return true;
    }

    //! @brief determine whether two paths are equivalent tours
    //! Equivalent tours begin and end at the same point, are the same length, and visit the same set of vertices
    bool Path::EquivalentTour( const Path &PATH) const
    {
      return GetSize() == PATH.GetSize()
             && FirstElement() == PATH.FirstElement()
             && LastElement() == PATH.LastElement()
             && Covers( PATH);
    }

    //! @return true if two paths are identical
    bool Path::Identical( const Path &PATH) const
    {
      return GetSize() == PATH.GetSize() && StartsWith( PATH);
    }

    //! @brief Count the number of overlapping vertices
    size_t Path::CountOverlappingVertices( const Path &PATH) const
    {
      if( GetSize() > PATH.GetSize())
      {
        return PATH.CountOverlappingVertices( *this);
      }
      size_t overlap( 0);
      for( size_t i( 0), sz( m_Vertices.GetSize()); i < sz; ++i)
      {
        if( PATH.Contains( m_Vertices( i)))
        {
          ++overlap;
        }
      }
      return overlap;
    }

    //! @brief reverse the path
    void Path::Reverse()
    {
      if( !m_Undirected) // only reverse directed paths
      {
        std::reverse( m_Vertices.Begin(), m_Vertices.End());
      }
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Path::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Vertices, ISTREAM);
      io::Serialize::Read( m_VertexBitMap, ISTREAM);
      io::Serialize::Read( m_Undirected, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Path::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Vertices, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_VertexBitMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Undirected, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the size of parent graph
    //! @return return size of parent graph
    size_t Path::GetGraphSize() const
    {
      return m_VertexBitMap.size();
    }

    //! @brief reorder vertices
    //! @param VERTICES the new ordering of the vertices; must not change which vertices are visited
    void Path::ReorderVertices( const storage::Vector< size_t> &PATH)
    {
      m_Vertices = PATH;
    }

  } // namespace graph
} // namespace bcl
