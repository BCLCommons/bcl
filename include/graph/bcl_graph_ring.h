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

#ifndef BCL_GRAPH_RING_H_
#define BCL_GRAPH_RING_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_path.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Ring
    //! @brief a high-performance Ring class.
    //! @details Rings have a defined order in which they visit vertices and are directed
    //! Ring provides
    //!   1. O(1) time lookups for whether a vertex is in the ring
    //!   2. Provision for fusing two rings along an edge ( this is helpful for determining aromaticity of fused rings
    //!      even if component rings are not aromatic e.g. azulene)
    //!
    //! @see @link example_graph_ring.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 12, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Ring :
      protected Path,
      public virtual util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

    public:

      typedef Path::const_iterator         const_iterator;
      typedef Path::const_reverse_iterator const_reverse_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, creates an empty Ring
      Ring();

      //! @brief construct from vertices and edges
      //! @param GRAPH_SIZE the size of the graph from which the Ring is a part
      //! @param VERTICES a vector containing the vertices in the order they were visited
      Ring( const size_t &GRAPH_SIZE, const storage::Vector< size_t> &VERTICES);

      //! clone the object
      Ring *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class CircularIterator
      //! @brief can iterate through the ring (left is positive direction), cannot change ring elements
      //! @author kothiwsk, mendenjl
      //! @date Dec 12, 2011
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class BCL_API CircularIterator
      {

      private:

      //////////
      // data //
      //////////

        const size_t *m_Current; //!< pointer to the ring element currently iterated at
        const size_t *m_Begin;   //!< pointer to the start of the ring element
        const size_t *m_End;     //!< pointer to the end of the ring; enables wrap around
        bool m_IsForward;        //!< Is this iterator going forward in the ring

      public:

        //! @brief default constructor
        CircularIterator();

        //! @brief construct iterator from ring info, index, and direction
        //! @param RING ring the iterator should operate over
        //! @param INDEX iterator points to vertex indexed in vertex vector of ring
        //! @param IS_FORWARD determines the direction in which iterator moves
        CircularIterator( const storage::Vector< size_t> &RING, const size_t &INDEX, const bool &IS_FORWARD = true);

        //! @brief construct from ring info and iterator
        //! @param RING ring the iterator should operate over
        //! @param ITR iterator points to desired vertex
        CircularIterator( const storage::Vector< size_t> &RING, const storage::Vector< size_t>::const_iterator &ITR);

        //! @brief checks if iterator is defined
        //! @return true if iterator is defined
        bool IsDefined() const;

        //! @brief pre increment - sets m_Current to the next size_t of the ring
        //! @return reference to this iterator
        CircularIterator &operator ++();

        //! @brief pre decrement - sets m_Current to the previous size_t of the ring
        //! @return reference to this iterator
        CircularIterator &operator --();

        //! @brief inequality comparison
        //! @param ITERATOR right hand side iterator
        //! @return true, if m_Current point to different size_ts
        bool operator !=( const CircularIterator &ITERATOR) const;

        //! @brief equality comparison
        //! @param ITERATOR right hand side iterator
        //! @return true, if m_Current point to the same size_ts
        bool operator ==( const CircularIterator &ITERATOR) const;

        //! @brief equality comparison
        //! @param ITERATOR right hand side iterator
        //! @return true, if m_Current point to the same size_ts
        bool operator ==( const const_iterator &ITERATOR) const;

        //! @brief operator *
        //! @return reference to ring element Iterator is pointing to
        const size_t &operator *() const;

      }; // class CircularIterator

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const;

      //! @brief return the ordered vertices in the ring
      //! @return the ordered vertices in the ring
      const storage::Vector< size_t> &GetVertices() const;

      //! @brief GetSize
      //! @return the number of vertices in the ring
      size_t GetSize() const;

      //! @brief get the first element of the ring
      //! @return the first element of the ring
      size_t FirstElement() const;

      //! @brief get the last element of the ring
      //! @return the last element of the ring
      size_t LastElement() const;

      //! @brief get an iterator to the beginning of the ring
      //! @return an iterator to the beginning of the ring
      const_iterator Begin() const;

      //! @brief get an iterator to the end of the ring
      //! @return an iterator to the end of the ring
      const_iterator End() const;

      //! @brief get an iterator to the reverse beginning of the ring
      //! @return an iterator to the reverse beginning of the ring
      const_reverse_iterator ReverseBegin() const;

      //! @brief get an iterator to the reverse end of the ring
      //! @return an iterator to the reverse end of the ring
      const_reverse_iterator ReverseEnd() const;

      //! @brief Test whether the Ring contains a given vertex
      //! @return true if the Ring contains a particular vertex
      bool Contains( const size_t &VERTEX) const;

      //! @brief get the vertices that this Ring has in common with another ring
      //! @return the path that this ring has in common with another ring
      Path GetOverlap( const Ring &RING) const;

    ////////////////
    // operations //
    ////////////////

      //! @return true if two Rings are identical
      bool Identical( const Ring &Ring) const;

      //! @brief remove path from this ring
      //! @param PATH remove a path from Ring
      //! @return new path obtained from this ring after path is removed
      Path Remove( const Path &PATH) const;

      //! @brief find whether the a path is part of this ring
      //! @param PATH the path to search for
      //! @return a circular-iterator that points to vertex that is part of path in this ring
      CircularIterator Find( const Path &PATH) const;

      //! @brief fuse two rings
      //! @param RING_A first ring that needs to be fused to the second
      //! @param RING_B second ring that needs to be fused to the fused
      //! @return fused ring obtained by fusing RING_A with RING_B
      static Ring FuseRings( const Ring &RING_A, const Ring &RING_B);

      //! @brief Get the nominal size of a fusion of RING_A with RING_B
      //! @param RING_A first ring
      //! @param RING_B second ring
      //! @return The nominal size of a fusion of RING_A with RING_B. The actual size will be smaller (0) if the rings do
      //!         not share u unique, single, contiguous overlapping path
      static size_t GetNominalFusionSize( const Ring &RING_A, const Ring &RING_B);

    ///////////////
    // operators //
    ///////////////

      //! @brief test whether one ring is less than another, needed for ordering rings in a set or map
      //! @return bool, whether one ring's vertices are less than another
      bool operator <( const Ring &RING) const;

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

    }; // class Ring

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_RING_H_

