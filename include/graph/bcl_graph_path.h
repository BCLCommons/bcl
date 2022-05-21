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

#ifndef BCL_GRAPH_PATH_H_
#define BCL_GRAPH_PATH_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Path
    //! @brief a high-performance path class.
    //! @details Paths have a defined order in which they visit vertices and can be undirected
    //! Path provides
    //!   1. O(1) time lookups for whether a vertex is in the graph
    //!   2. O(1) iteration through connections for a given vertex
    //!   3. Capacity to handle directed as well as undirected paths (paths which should always be ordered w/ the 1st vertex <= the last)
    //!
    //! @see @link example_graph_path.cpp @endlink
    //! @author mendenjl
    //! @date September 2, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Path :
      public virtual util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      storage::Vector< size_t> m_Vertices;       //!< the order of the vertices
      std::string              m_VertexBitMap;   //!< m_VertexBitMap[ x] == '1' iff the path contains x, '0'
      bool                     m_Undirected;     //!< whether the path can be reversed so as to maintain the invariant that the first vertex in the path is always < the last vertex

    public:

      typedef storage::Vector< size_t>::const_iterator         const_iterator;
      typedef storage::Vector< size_t>::const_reverse_iterator const_reverse_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, creates an empty path
      Path();

      //! @brief construct from vertices and edges
      //! @param GRAPH_SIZE the size of the graph from which the path is a part
      //! @param VERTICES a vector containing the vertices in the order they were visited
      //! @param UNDIRECTED whether the path is undirected
      //! @brief initialize a path with a certain graph size
      Path
      (
        const size_t &GRAPH_SIZE,
        const storage::Vector< size_t> &VERTICES,
        const bool &UNDIRECTED
      );

      //! @brief create a new path by connecting two existing paths and a mask for which vertices to fuse
      //! @param PATH_A, PATH_B the two paths to combine
      //! @param OVERLAP the path along which to combine the vertices
      //! e.g Path([1,2,3], [0,1,2]) -> 0,1,2,1,2,3
      //!     Path([1,2,3], [0,1,2], [1,2]) -> 0,1,2,3
      //!     Path([1,2,3], [0,1,2], [2]) -> 0,1,1,2,3
      Path( const Path &PATH_A, const Path &PATH_B, const Path &OVERLAP);

      //! clone the object
      Path *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const;

      //! @return true if the path is directed
      bool IsDirected() const;

      //! @return true if the path is undirected
      bool IsUndirected() const;

      //! @brief return the edge ids
      //! @return the edge ids
      const storage::Vector< size_t> &GetVertices() const;

      //! @brief GetSize
      //! @return the number of vertices in the path
      size_t GetSize() const;

      //! @brief get the first element of the path
      //! @return the first element of the path
      size_t FirstElement() const;

      //! @brief GetSize
      //! @return the number of vertices in the path
      size_t LastElement() const;

      //! @brief get an iterator to the beginning of the path
      //! @return an iterator to the beginning of the path
      const_iterator Begin() const;

      //! @brief get an iterator to the end of the path
      //! @return an iterator to the end of the path
      const_iterator End() const;

      //! @brief get an iterator to the reverse begining of the path
      //! @return an iterator to the reverse begining of the path
      const_reverse_iterator ReverseBegin() const;

      //! @brief get an iterator to the reverse end of the path
      //! @return an iterator to the reverse end of the path
      const_reverse_iterator ReverseEnd() const;

      //! @brief Test whether the path contains a given vertex
      //! @return true if the path contains a particular vertex
      bool Contains( const size_t &VERTEX) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief determine if the path ends at VERTEX
      //! @return true if the path ends on VERTEX
      bool EndsAt( const size_t &VERTEX) const;

      //! @brief count the number of times vertices visited more than once in the path
      //! @return the number of times vertices visited more than once in the path
      size_t CountVertexRepetitions() const;

      //! @brstorage::Vector< size_t>::const_iterator itr_aief determine if the paths overlap anywhere other than at their ends
      //! @return true if the paths overlap anywhere other than at their ends
      bool Crosses( const Path &PATH) const;

      //! @brief true if one path connects to another, meaning, the paths contain a common starting or ending vertex
      bool Connects( const Path &PATH) const;

      //! @brief determine whether this path covers another path (that is, visits the same vertices)
      //! @return true if this path visits all the vertices in PATH
      bool Covers( const Path &PATH) const;

      //! @brief determine whether two paths are equivalent tours
      //! Equivalent tours begin and end at the same point, are the same length, and visit the same set of vertices
      bool EquivalentTour( const Path &PATH) const;

      //! @return true if two paths are identical
      bool Identical( const Path &PATH) const;

      //! @brief reverse the path
      void Reverse();

      //! @brief test whether the start of this path is the same as PATH
      //! @param PATH the path to search for
      //! @return whether the path starts with PATH
      //! @note This functions treats both paths as directed
      bool StartsWith( const Path &PATH) const;

      //! @brief search for a path in this path
      //! @param PATH the path to search for, must be directed
      //! @return whether the path starts with PATH
      //! @note This functions treats both paths as directed
      bool EndsWith( const Path &PATH) const;

      //! @brief Count the number of overlapping vertices
      size_t CountOverlappingVertices( const Path &PATH) const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the size of parent graph
      //! @return return size of parent graph
      size_t GetGraphSize() const;

      //! @brief reorder vertices
      //! @param VERTICES the new ordering of the vertices; must not change which vertices are visited
      void ReorderVertices( const storage::Vector< size_t> &PATH);

    }; // class Path

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_PATH_H_

