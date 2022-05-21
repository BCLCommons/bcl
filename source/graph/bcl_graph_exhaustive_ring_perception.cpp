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
#include "graph/bcl_graph_exhaustive_ring_perception.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_graph_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ExhaustiveRingPerception::s_Instance
    (
      GetObjectInstances().AddInstance( new ExhaustiveRingPerception())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExhaustiveRingPerception::ExhaustiveRingPerception() :
      m_Paths(),
      m_Rings()
    {
    }

    //! @brief constructor from a graph
    //! @brief GRAPH graph to search rings in
    ExhaustiveRingPerception::ExhaustiveRingPerception
    (
      const GraphWithData< size_t, size_t> &GRAPH
    ) :
      m_GraphSize( GRAPH.GetVertices().GetSize())
    {
      // initialize m_Paths with all edges
      Initialize( GRAPH);

      // remove all vertices, update paths, collect rings
      while( !m_Paths.IsEmpty())
      {
        // in the paper it's suggested to remove vertices with smaller amount of edges first
        // such a heuristic should go here instead of taking the 'next' vertex
        Remove( m_Paths.FirstElement().FirstElement());
      }
    }

    //! @brief constructor from a graph
    //! @brief GRAPH graph to search rings in
    ExhaustiveRingPerception::ExhaustiveRingPerception
    (
      const ConstGraph< size_t, size_t> &GRAPH
    ) :
      m_GraphSize( GRAPH.GetSize())
    {
      // initialize m_Paths with all edges
      Initialize( GRAPH);

      // remove all vertices, update paths, collect rings
      while( !m_Paths.IsEmpty())
      {
        // in the paper it's suggested to remove vertices with smaller amount of edges first
        // such a heuristic should go here instead of taking the 'next' vertex
        Remove( m_Paths.FirstElement().FirstElement());
      }
    }

    //! @brief initialize m_Paths with all edges
    void ExhaustiveRingPerception::Initialize( const GraphWithData< size_t, size_t> &GRAPH)
    {
      // initialize m_Paths
      for
      (
        GraphWithData< size_t, size_t>::VertexContainerType::const_iterator
          itr_graph_vertex( GRAPH.GetVertices().Begin()),
          itr_graph_vertex_end( GRAPH.GetVertices().End());
        itr_graph_vertex != itr_graph_vertex_end;
        ++itr_graph_vertex
      )
      {
        for
        (
          GraphWithData< size_t, size_t>::EdgeContainerType::const_iterator
            itr_graph_edge( ( *itr_graph_vertex)->GetEdges().Begin()),
            itr_graph_edge_end( ( *itr_graph_vertex)->GetEdges().End());
          itr_graph_edge != itr_graph_edge_end;
          ++itr_graph_edge
        )
        {
          storage::Vector< size_t> initial_path( 2);
          initial_path( 0) = ( *itr_graph_vertex)->GetData();
          initial_path( 1) = itr_graph_edge->GetTarget()->GetData();

          if( GRAPH.IsDirected() || initial_path( 0) <= initial_path( 1))
          {
            m_Paths.PushBack( initial_path);
          }
        }
      }
    }

    //! @brief initialize m_Paths with all edges
    void ExhaustiveRingPerception::Initialize
    (
      const ConstGraph< size_t, size_t> &GRAPH
    )
    {
      const size_t size( GRAPH.GetSize());

      // determine the vertex girths; e.g. the shortest cycle containing a particular vertex
      storage::Vector< size_t> shortest_cycles( GRAPH.GetSize());
      for( size_t vertex_number = 0; vertex_number < size; ++vertex_number)
      {
        shortest_cycles( vertex_number) = Connectivity::LengthOfSmallestCycleWithVertex( GRAPH, vertex_number);
      }

      const bool is_directed( !GRAPH.IsUndirected()); // determine whether the graph is undirected
      // initialize m_Paths
      for
      (
        size_t source_vertex_number = 0, graph_size = GRAPH.GetSize();
        source_vertex_number != graph_size;
        ++source_vertex_number
      )
      {
        if( shortest_cycles( source_vertex_number) > graph_size) // this vertex is not part of any cycle
        {
          continue;
        }
        const ConstGraph< size_t, size_t>::t_EdgeTargetsOfVertex &neighborhood
        (
          GRAPH.GetNeighborIndices( source_vertex_number)
        );
        for
        (
          ConstGraph< size_t, size_t>::t_EdgeTargetsOfVertex::const_iterator
            itr_edge_target( neighborhood.Begin()),
            itr_edge_target_end( neighborhood.End());
          itr_edge_target != itr_edge_target_end;
          ++itr_edge_target
        )
        {
          if( shortest_cycles( *itr_edge_target) > graph_size) // the edge target vertex is not part of any cycle
          {
            continue;
          }
          if( is_directed || source_vertex_number <= *itr_edge_target)
          {
            storage::Vector< size_t> initial_path( 2);
            initial_path( 0) = source_vertex_number;
            initial_path( 1) = *itr_edge_target;
            m_Paths.PushBack( initial_path);
          }
        }
      }

    }

    //! @brief remove a vertex, update paths, collect rings
    //! @param VERTEX vertex to be removed
    void ExhaustiveRingPerception::Remove( const size_t VERTEX)
    {
      BCL_MessageDbg
      (
        "Now " + util::Format()( m_Paths.GetSize()) + " paths and " + util::Format()( m_Rings.GetSize()) + " rings"
      );

      // storage for the updated paths
      storage::List< storage::Vector< size_t> > paths_with_vertex;

      // transfer any paths that end in VERTEX into paths_with_vertex
      for
      (
        storage::List< storage::Vector< size_t> >::iterator
          itr_paths( m_Paths.Begin()),
          itr_paths_end( m_Paths.End());
        itr_paths != itr_paths_end;
        // increment performed in loop
      )
      {
        if( VERTEX == itr_paths->FirstElement() || VERTEX == itr_paths->LastElement())
        {
          storage::List< storage::Vector< size_t> >::iterator old_itr_paths( itr_paths);
          ++itr_paths;
          paths_with_vertex.InternalData().splice( paths_with_vertex.End(), m_Paths.InternalData(), old_itr_paths);
        }
        else
        {
          ++itr_paths;
        }
      }

      // make an incidence vector for PATH_A
      // vertex_is_in_path[ x] -> tells whether vertex #x is in the path
      std::vector< bool> vertex_is_in_path( m_GraphSize);

      // update paths
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr_paths( paths_with_vertex.Begin()),
          itr_paths_end( paths_with_vertex.End());
        itr_paths != itr_paths_end;
        ++itr_paths
      )
      {
        const size_t path_size( itr_paths->GetSize() - 2);
        SetupAdjacencyVector( vertex_is_in_path, *itr_paths);

        storage::List< storage::Vector< size_t> >::const_iterator itr_paths_match( itr_paths);
        ++itr_paths_match;

        while( itr_paths_match != itr_paths_end)
        {
          if
          (
            path_size + itr_paths_match->GetSize() <= m_GraphSize // eliminate obvious overlaps
            && !Overlap( vertex_is_in_path, *itr_paths_match)
          )
          {
            // combine the paths and add them to m_Paths
            m_Paths.InsertElement( m_Paths.Begin(), CombinePaths( VERTEX, *itr_paths, *itr_paths_match));

            storage::List< storage::Vector< size_t> >::iterator itr_begin( m_Paths.Begin());

            // check whether the new ring was actually a path
            if( itr_begin->FirstElement() == itr_begin->LastElement())
            {
              // yep, so splice it into rings instead (this removes it from m_Paths)
              m_Rings.InternalData().splice( m_Rings.End(), m_Paths.InternalData(), itr_begin);
            }
          }

          ++itr_paths_match;
        }
      }
    } // Remove

    //! @brief Combine two paths at their common vertex
    //! @param PATH_A first path to be combined
    //! @param PATH_B second path to be combined
    //! @return combined path
    storage::Vector< size_t> ExhaustiveRingPerception::CombinePaths
    (
      const size_t COMMON_VERTEX,
      const storage::Vector< size_t> &PATH_A,
      const storage::Vector< size_t> &PATH_B
    ) const
    {
      storage::Vector< size_t> new_path( PATH_B.GetSize() + PATH_A.GetSize() - 1);

      storage::Vector< size_t>::iterator new_path_itr( new_path.Begin());

      // copy path A such that the common vertex is last
      if( COMMON_VERTEX == PATH_A.LastElement())
      { // copy forward
        new_path_itr = std::copy( PATH_A.Begin(), --PATH_A.End(), new_path_itr);
      }
      else // COMMON_VERTEX == PATH_A.FirstElement()
      { // copy in reverse
        new_path_itr = std::copy( PATH_A.ReverseBegin(), --PATH_A.ReverseEnd(), new_path_itr);
      }

      // copy path B such that the common vertex is first
      if( COMMON_VERTEX == PATH_B.FirstElement()) // follow with path B, omitting the common element
      { //copy forward
        std::copy( PATH_B.Begin(), PATH_B.End(), new_path_itr);
      }
      else
      { // copy in reverse
        std::copy( PATH_B.ReverseBegin(), PATH_B.ReverseEnd(), new_path_itr);
      }

      return new_path;
    } // CombinePaths

    //! @brief Check if two paths overlap
    //! @param VERTEX_IS_IN_PATH_A adjacency vector for path a
    //! @param PATH_B second path to be combined
    //! @return if paths overlap
    bool ExhaustiveRingPerception::Overlap
    (
      const std::vector< bool> &VERTEX_IS_IN_PATH_A,
      const storage::Vector< size_t> &PATH_B
    ) const
    {
      // check all internal elements in PATH_A for occurrence in PATH_B
      for
      (
        storage::Vector< size_t>::const_iterator
          itr_path( PATH_B.Begin()),
          itr_path_end( PATH_B.End());
        itr_path != itr_path_end;
        ++itr_path
      )
      {
        if( VERTEX_IS_IN_PATH_A[ *itr_path])
        {
          return true;
        }
      }

      return false;
    }

    //! @brief set the adjacency vector with a given path
    //! @param VERTEX_IS_IN_PATH the vector to setup such that VERTEX_IS_IN_PATH[ x] == true iff x is in PATH
    //! @param PATH the path to use in setting up the adjacency vector
    void ExhaustiveRingPerception::SetupAdjacencyVector
    (
      std::vector< bool> &VERTEX_IS_IN_PATH,
      const storage::Vector< size_t> &PATH
    ) const
    {
      VERTEX_IS_IN_PATH.assign( m_GraphSize, false);
      for
      (
        storage::Vector< size_t>::const_iterator itr_path( ++PATH.Begin()), itr_path_end( --PATH.End());
        itr_path != itr_path_end;
        ++itr_path
      )
      {
        VERTEX_IS_IN_PATH[ *itr_path] = true;
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExhaustiveRingPerception::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Paths, ISTREAM);
      io::Serialize::Read( m_Rings, ISTREAM);
      io::Serialize::Read( m_GraphSize, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ExhaustiveRingPerception::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Paths, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Rings, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GraphSize, OSTREAM, INDENT);

      return OSTREAM;
    }
  } // namespace graph
} // namespace bcl
