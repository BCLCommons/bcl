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
#include "graph/bcl_graph_ring.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

    //! @brief default constructor, creates an empty Ring
    Ring::Ring() :
      Path()
    {
    }

    //! @brief construct from vertices and edges
    //! @param GRAPH_SIZE the size of the graph from which the Ring is a part
    //! @param VERTICES a vector containing the vertices in the order they were visited
    //! @brief initialize a Ring with a certain graph size
    Ring::Ring
    (
      const size_t &GRAPH_SIZE,
      const storage::Vector< size_t> &VERTICES
    ) :
      Path( GRAPH_SIZE, VERTICES, false)
    {
      if( VERTICES.GetSize() <= 1)
      {
        return;
      }

      BCL_Assert( !Path::CountVertexRepetitions(), "Should be no repeated vertices for a ring!");

      // canonicalize the path; the first vertex in VERTICES should be the one with the smallest index
      // the second vertex should be the smaller of the two indices links to vertices
      const size_t min_index( math::Statistics::MinimumIndex( VERTICES.Begin(), VERTICES.End()));

      if( min_index != 0 || VERTICES( 1) > VERTICES.LastElement())
      {
        // construct a new vector of vertices starting from just after min index till the last element of this ring's indices
        storage::Vector< size_t> reordered( Path::GetVertices(), min_index + 1);

        // append vertices from beginning to the desired index
        reordered.Append( storage::Vector< size_t>( Path::GetVertices(), 0, min_index));

        // test whether the second vertex is greater than the last; if so, reverse
        if( reordered.FirstElement() > reordered.LastElement())
        {
          std::reverse( reordered.Begin(), reordered.End());
        }

        // last, insert the min element again at the beginning of the vector
        reordered.InsertElements( 0, VERTICES( min_index), 1);

        // reorder the vertices accordingly
        Path::ReorderVertices( reordered);
      }
    }

    //! clone the object
    Ring *Ring::Clone() const
    {
      return new Ring( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief default constructor
    Ring::CircularIterator::CircularIterator() :
      m_Current( NULL),
      m_Begin( NULL),
      m_End( NULL),
      m_IsForward( true)
    {
    }

    //! @brief construct iterator from ring info, index, and direction
    //! @param RING ring the iterator should operate over
    //! @param INDEX iterator points to vertex indexed in vertex vector of ring
    //! @param IS_FORWARD determines the direction in which iterator moves
    Ring::CircularIterator::CircularIterator
    (
      const storage::Vector< size_t> &RING,
      const size_t &INDEX,
      const bool &IS_FORWARD
    ) :
      m_Current( RING.Begin().base() + INDEX),
      m_Begin( IS_FORWARD ? RING.Begin().base() : RING.End().base() - 1),
      m_End( IS_FORWARD ? RING.End().base() : RING.Begin().base() - 1),
      m_IsForward( IS_FORWARD)
    {
    }

    //! @brief construct from ring info and iterator
    //! @param RING ring the iterator should operate over
    //! @param ITR iterator points to desired vertex
    Ring::CircularIterator::CircularIterator
    (
      const storage::Vector< size_t> &RING,
      const storage::Vector< size_t>::const_iterator &ITR
    ) :
      m_Current( ITR.base()),
      m_Begin( RING.Begin().base()),
      m_End( RING.End().base()),
      m_IsForward( true)
    {
    }

    //! @brief checks if iterator is defined
    //! @return true if iterator is defined
    bool Ring::CircularIterator::IsDefined() const
    {
      return m_Current != m_End;
    }

    //! @brief pre increment - sets m_Current to the left ring element of the current ring
    //! @return Iterator pointing to the next left element
    Ring::CircularIterator &Ring::CircularIterator::operator ++()
    {
      m_IsForward ? ++m_Current : --m_Current;

      if( m_Current == m_End)
      {
        m_Current = m_Begin;
      }

      return *this;
    }

    //! @brief pre decrement - sets m_Current to the previous size_t of the ring
    //! @return reference to this iterator
    Ring::CircularIterator &Ring::CircularIterator::operator --()
    {
      if( m_Current == m_Begin)
      {
        m_Current = m_End;
      }

      m_IsForward ? --m_Current : ++m_Current;

      return *this;

    }

    //! @brief inequality comparison
    //! @param ITERATOR right hand side iterator
    //! @return true, if m_Current point to different size_t
    bool Ring::CircularIterator::operator !=( const Ring::CircularIterator &ITERATOR) const
    {
      return m_Current != ITERATOR.m_Current;
    }

    //! @brief equality comparison
    //! @param ITERATOR right hand side iterator
    //! @return true, if m_Current point to the same size_t
    bool Ring::CircularIterator::operator ==( const Ring::CircularIterator &ITERATOR) const
    {
      return m_Current == ITERATOR.m_Current;
    }

    //! @brief equality comparison
    //! @param ITERATOR right hand side iterator
    //! @return true, if m_Current point to the same size_t
    bool Ring::CircularIterator::operator ==( const const_iterator &ITERATOR) const
    {
      return m_Current == ITERATOR.base();
    }

    //! @brief operator *
    //! @return reference to index Iterator is pointing to
    const size_t &Ring::CircularIterator::operator *() const
    {
      return *m_Current;
    }

    //! @return class name of the object behind a pointer or the current object
    const std::string &Ring::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the ordered vertices in the ring
    //! @return the ordered vertices in the ring
    const storage::Vector< size_t> &Ring::GetVertices() const
    {
      return Path::GetVertices();
    }

    //! @brief GetSize
    //! @return the number of vertices in the ring
    size_t Ring::GetSize() const
    {
      return Path::GetSize();
    }

    //! @brief get the first element of the ring
    //! @return the first element of the ring
    size_t Ring::FirstElement() const
    {
      return Path::FirstElement();
    }

    //! @brief get the last element of the ring
    //! @return the last element of the ring
    size_t Ring::LastElement() const
    {
      return Path::LastElement();
    }

    //! @brief get an iterator to the beginning of the ring
    //! @return an iterator to the beginning of the ring
    Ring::const_iterator Ring::Begin() const
    {
      return Path::Begin();
    }

    //! @brief get an iterator to the end of the ring
    //! @return an iterator to the end of the ring
    Ring::const_iterator Ring::End() const
    {
      return Path::End();
    }

    //! @brief get an iterator to the reverse beginning of the ring
    //! @return an iterator to the reverse beginning of the ring
    Ring::const_reverse_iterator Ring::ReverseBegin() const
    {
      return Path::ReverseBegin();
    }

    //! @brief get an iterator to the reverse end of the ring
    //! @return an iterator to the reverse end of the ring
    Ring::const_reverse_iterator Ring::ReverseEnd() const
    {
      return Path::ReverseEnd();
    }

    //! @brief Test whether the Ring contains a given vertex
    //! @return true if the Ring contains a particular vertex
    bool Ring::Contains( const size_t &VERTEX) const
    {
      return Path::Contains( VERTEX);
    }

    //! @brief get the vertices that this Ring has in common with another ring
    //! @return the path that this ring has in common with another ring
    Path Ring::GetOverlap( const Ring &RING) const
    {
      // set the first_match to point to iterator pointing to end of vertex vector
      const_iterator first_match( End());

      // if the given ring does not contain the first element of this ring, search which vertex of this ring is common with given ring
      if( !RING.Contains( FirstElement()))
      {
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          if( RING.Contains( *itr))
          {
            // set first_match to the vertex that is common to both the rings and break after finding first match
            first_match = itr;
            break;
          }
        }

        storage::Vector< size_t> overlapping_indices;

        // if no common vertex is found, then first_match points to end of vertex vector
        if( first_match == End())
        {
          // if no common vertex is found return an empty path
          return Path( GetGraphSize(), overlapping_indices, false);
        }

        // if a common vertex is found then start iterating over vertices and store those that are common to both rings
        const size_t match_position( RING.GetVertices().Find( *first_match));
        CircularIterator itr_ring( RING.GetVertices(), match_position, true);
        for( const_iterator itr_end( End()); first_match != itr_end && *first_match == *itr_ring; ++first_match, ++itr_ring)
        {
          // pushback common vertices into overlapping_indices
          overlapping_indices.PushBack( *first_match);
        }
        if( overlapping_indices.GetSize() == 1)
        {
          ----itr_ring;
          for( const_iterator itr_end( End()); first_match != itr_end && *first_match == *itr_ring; ++first_match, --itr_ring)
          {
            // pushback common vertices into overlapping_indices
            overlapping_indices.PushBack( *first_match);
          }
        }

        // return path containing indices that is common to both rings
        return Path( GetGraphSize(), overlapping_indices, false);
      }

      // if first element of this ring is contained in the given ring then find the first vertex that is unique to one of the rings
      const_iterator first_nonmatch( End());
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( !RING.Contains( *itr))
        {
          first_nonmatch = itr;
          break;
        }
      }

      // if all vertices of RING are in this ring, then check if the order is correct
      if( first_nonmatch == End())
      {
        if( Identical( RING)) // correct order, return the whole path
        {
          return Path( GetGraphSize(), GetVertices(), false);
        }
        first_nonmatch = Begin();
      }

      storage::Vector< size_t> overlapping_indices;

      // if a an uncommon index is found then get an circular iterator pointing to first unique vertex found after search
      CircularIterator itr( GetVertices(), first_nonmatch);
      // search for next common index
      while( !RING.Contains( *itr))
      {
        ++itr;
      }

      // if a common vertex is found then start iterating over vertices and store those that are common to both rings
      const size_t match_position( RING.GetVertices().Find( *itr));
      CircularIterator itr_ring( RING.GetVertices(), match_position, true);
      for( size_t index( 0), size_this( GetSize()); index < size_this && *itr == *itr_ring; ++itr, ++itr_ring, ++index)
      {
        // pushback common vertices into overlapping_indices
        overlapping_indices.PushBack( *itr);
      }
      if( overlapping_indices.GetSize() == 1)
      {
        ----itr_ring;
        for( size_t index( 1), size_this( GetSize()); index < size_this && *itr == *itr_ring; ++itr, --itr_ring, ++index)
        {
          // pushback common vertices into overlapping_indices
          overlapping_indices.PushBack( *itr);
        }
      }

      // return a path containing the vertices
      return Path( GetGraphSize(), overlapping_indices, false);
    }

  ////////////////
  // operations //
  ////////////////

    //! @return true if two Rings are identical
    bool Ring::Identical( const Ring &RING) const
    {
      return
        RING.GetSize() == GetSize()
        && RING.GetGraphSize() == GetGraphSize()
        && RING.GetVertices() == GetVertices();
    }

    //! @brief remove path from this ring
    //! @param PATH remove a path from Ring
    //! @return new Path obtained from Ring after Path is removed
    Path Ring::Remove( const Path &PATH) const
    {
      storage::Vector< size_t> new_path;

      // get the reverse iterator that points to vertex in the ring that is first vertex in this path
      CircularIterator iterator( Find( PATH));

      // if vertex is null, that means path was not found in ring and so return the path of this ring
      if( !iterator.IsDefined())
      {
        return *this;
      }

      // get a copy of iterator
      CircularIterator itr_begin( iterator);

      // skip vertices that are present in the given path
      for( size_t index( 0), path_size( PATH.GetSize()); index < path_size; ++iterator, ++index)
      {
      }

      // once vertices that are found in the given path are skipped, store remaining vertices in a vector
      while( iterator != itr_begin)
      {
        new_path.PushBack( *iterator);
        ++iterator;
      }

      // return new path constructed from new vertex vector
      return Path( GetGraphSize(), new_path, false);
    }

    //! @brief find whether the a path is part of this ring
    //! @param PATH the path to search for
    //! @return a circular-iterator that points to vertex that is part of path in this ring
    Ring::CircularIterator Ring::Find( const Path &PATH) const
    {
      // if given path is of size 0 then return a null iterator
      if( PATH.GetSize() == 0 || PATH.GetSize() > GetSize())
      {
        return CircularIterator();
      }

      // find the index of first element of given path in this ring
      const size_t find_position( GetVertices().Find( PATH.FirstElement()));

      // get circular iterator pointing to the vertex of this ring that find_position specifies
      CircularIterator itr( GetVertices(), find_position, true);

      // check if itr is defined
      if( !itr.IsDefined())
      {
        return itr;
      }

      // get a circulariterator for the given path which points to first vertex
      CircularIterator itr_path_for( PATH.GetVertices(), 0, true);

      // itr points to first vertex of the given path. increment itr until we reach end of the given path or until we find an vertex in given path that
      // does not exist in this ring
      for( size_t index( 0), path_size( PATH.GetSize()); index < path_size && *itr == *itr_path_for; ++itr_path_for, ++itr, ++index)
      {
      }
      // if the path is found return the forward iterator
      if( itr_path_for == PATH.Begin())
      {
        return CircularIterator( GetVertices(), find_position, true);
      }

      // get a forward iterator starting from first vertex in the given path
      itr_path_for = CircularIterator( PATH.GetVertices(), 0);

      // iterate in the other direction from before on this ring
      itr = CircularIterator( GetVertices(), find_position, false);

      // search in the other direction
      for( size_t index( 0), path_size( PATH.GetSize()); index < path_size && *itr == *itr_path_for; ++itr_path_for, ++itr, ++index)
      {
      }

      // if path is found return the reverse iterator
      if( itr_path_for == PATH.Begin())
      {
        return CircularIterator( GetVertices(), find_position, false);
      }

      // if the complete path not found in the ring then return a null iterator
      return CircularIterator();
    }

    //! @brief fuse two rings
    //! @param RING_A first ring that needs to be fused to the second
    //! @param RING_B second ring that needs to be fused to the fused
    //! @return fused ring obtained by fusing RING_A with RING_B
    Ring Ring::FuseRings( const Ring &RING_A, const Ring &RING_B)
    {
      // check that both the rings have same graph size
      BCL_Assert
      (
        RING_A.GetGraphSize() == RING_B.GetGraphSize(),
        "Can't fuse Rings from different graphs"
      );

      // get overlapping path between two rings
      Path overlap_path( RING_A.GetOverlap( RING_B));

      // if overlap is less than two vertices, then rings cannot be fused
      if( overlap_path.GetSize() < size_t( 2))
      {
        return Ring();
      }

      // for a bridge system retain only the base vertices
      Path bridge_a
      (
        RING_A.GetGraphSize(),
        storage::Vector< size_t>( overlap_path.GetVertices(), 1, overlap_path.GetSize() - 1),
        false
      );

      // remove overlap path from both ring A and B; the remaining elements of a and b are now aligned in the ring along the vector given by overlap path
      Path fuseable_component_a( RING_A.Remove( overlap_path));
      Path fuseable_component_b( RING_B.Remove( overlap_path));

      // because a and b are now aligned along the overlapping path, one ring is going clockwise, the other counterclockwise, so they cannot be combined without
      // reversing one of the paths.  e.g.
      //       5
      // 0   1   6
      // 2   3   7
      //   4
      // the two rings are
      // 0 1 3 4 2
      // 1 5 6 7 3
      // after removal of the overlapping path ( 1 -> 3 ), the rings are then
      // 4 2 0 (which is clockwise)
      // 7 6 5 (which is counterclockwise)
      // The combined path is of course 0 2 4 3  7 6 5 1
      // so we have to reverse one of the paths to be able to combine them

      // reverse one path
      fuseable_component_a.Reverse();

      // add the endpoints of the overlapping path back onto the copy_b path; note that they too have to be reversed in order
      storage::Vector< size_t> ring_b_fusable_component_indices( 1, overlap_path.LastElement());
      ring_b_fusable_component_indices.Append( fuseable_component_b.GetVertices());
      ring_b_fusable_component_indices.PushBack( overlap_path.FirstElement());
      fuseable_component_b = Path( RING_A.GetGraphSize(), ring_b_fusable_component_indices, false);

      // create a path that fuses the rings together
      Path fused_rings( fuseable_component_a, fuseable_component_b, Path());

      // check whether the fusable components are disjoint (no common vertices; usual case) or
      // non-disjoint, which only happens if combining, e.g. a ring made of three basic rings in a row with the interior ring
      // in this case, there is no valid single fused ring, so return an empty ring
      if( fused_rings.CountVertexRepetitions() != size_t( 0))
      {
        return Ring();
      }

      // return the new ring
      return Ring( RING_A.GetGraphSize(), fused_rings.GetVertices());
    }

    //! @brief Get the nominal size of a fusion of RING_A with RING_B
    //! @param RING_A first ring
    //! @param RING_B second ring
    //! @return The nominal size of a fusion of RING_A with RING_B. The actual size will be smaller (0) if the rings do
    //!         not share u unique, single, contiguous overlapping path
    size_t Ring::GetNominalFusionSize( const Ring &RING_A, const Ring &RING_B)
    {
      size_t overlap( RING_A.CountOverlappingVertices( RING_B));
      if( overlap < size_t( 2))
      {
        // no fusion is possible as no more than a single vertex is shared
        return 0;
      }
      return RING_A.GetSize() + RING_B.GetSize() - 2 * ( overlap - 1);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief test whether one ring is less than another, needed for ordering rings in a set or map
    //! @return bool, whether one ring's vertices are less than another
    bool Ring::operator <( const Ring &RING) const
    {
      // smaller rings < larger rings
      if( GetSize() < RING.GetSize())
      {
        return true;
      }
      else if( GetSize() > RING.GetSize())
      {
        return false;
      }

      // rings are equal sized; compare vertices
      for( size_t i( 0), size( GetSize()); i < size; ++i)
      {
        const size_t this_vertex( GetVertices()( i));
        const size_t other_vertex( RING.GetVertices()( i));
        if( this_vertex < other_vertex)
        {
          return true;
        }
        else if( this_vertex > other_vertex)
        {
          return false;
        }
      }

      // rings are identical
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Ring::Read( std::istream &ISTREAM)
    {
      // read members
      Path::Read( ISTREAM);
      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Ring::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      Path::Write( OSTREAM, INDENT);
      //return
      return OSTREAM;
    }

  } // namespace graph
} // namespace bcl
