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
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CommonSubgraphIsomorphismBase::s_Instance
    (
      GetObjectInstances().AddInstance( new CommonSubgraphIsomorphismBase())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a solution type (connected by default)
    //! @param SOLUTION_TYPE whether the isomorphism must be connected
    CommonSubgraphIsomorphismBase::CommonSubgraphIsomorphismBase( const SolutionType &SOLUTION_TYPE) :
      m_SolutionType( SOLUTION_TYPE)
    {
    }

    //! clone the object
    CommonSubgraphIsomorphismBase *CommonSubgraphIsomorphismBase::Clone() const
    {
      return new CommonSubgraphIsomorphismBase( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @return class name of the object behind a pointer or the current object
    const std::string &CommonSubgraphIsomorphismBase::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the isomorphisms between graphs A and B
    //! @return the isomorphisms
    const storage::Vector< storage::Map< size_t, size_t> > &CommonSubgraphIsomorphismBase::GetIsomorphisms() const
    {
      return m_Isomorphisms;
    }

    //! @brief Get the largest common subgraph isomorphism between graphs A and B
    //! @return the largest common subgraph isomorphism
    const storage::Map< size_t, size_t> &CommonSubgraphIsomorphismBase::GetIsomorphism() const
    {
      static storage::Map< size_t, size_t> empty;
      return m_Isomorphisms.IsEmpty() ? empty : m_Isomorphisms.FirstElement();
    }

    //! @brief Get whether this isomorphism had to be connected
    //! @return bool, true indicates this isomorphism had to be connected
    CommonSubgraphIsomorphismBase::SolutionType CommonSubgraphIsomorphismBase::GetSolutionType() const
    {
      return m_SolutionType;
    }

    //! @brief set the isomorphisms
    //! @param ISOMORPHISMS the new isomorphisms to use
    void CommonSubgraphIsomorphismBase::SetIsomorphisms
    (
      const storage::Vector< storage::Map< size_t, size_t> > &ISOMORPHISMS
    )
    {
      m_Isomorphisms = ISOMORPHISMS;
    }

    //! @brief set the isomorphism, assuming there is only one
    //! @param ISOMORPHISM the new isomorphism to use
    void CommonSubgraphIsomorphismBase::SetIsomorphism( const storage::Map< size_t, size_t> &ISOMORPHISM)
    {
      m_Isomorphisms.Reset();
      m_Isomorphisms.PushBack( ISOMORPHISM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CommonSubgraphIsomorphismBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Isomorphisms, ISTREAM);

      size_t enum_value;
      ISTREAM >> enum_value;
      m_SolutionType = SolutionType( enum_value);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &CommonSubgraphIsomorphismBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Isomorphisms, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolutionType, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief swaps the isomorphisms
    //! thus makes indices in graph B be the keys to the isomorphism which map to indices in graph a
    void CommonSubgraphIsomorphismBase::SwapIsomorphisms()
    {
      for
      (
        storage::Vector< storage::Map< size_t, size_t> >::iterator
          itr_isos( m_Isomorphisms.Begin()), itr_isos_end( m_Isomorphisms.End());
        itr_isos != itr_isos_end;
        ++itr_isos
      )
      {
        // swap the isomorphism; used if graph_a was originally the larger graph

        // create a  new isomorphism, which will hold the swapped values
        storage::Map< size_t, size_t> new_isomorphism;
        for
        (
          storage::Map< size_t, size_t>::const_iterator iso( itr_isos->Begin()), iso_end( itr_isos->End());
          iso != iso_end;
          ++iso
        )
        {
          new_isomorphism[ iso->second] = iso->first;
        }

        // swap the isomorphism maps
        itr_isos->Swap( new_isomorphism);
      }
    }

  } // namespace graph
} // namespace bcl
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
#include "graph/bcl_graph.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace graph
} // namespace bcl
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
#include "graph/bcl_graph_edge_cover_ring_perception.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EdgeCoverRingPerception::s_Instance
    (
      GetObjectInstances().AddInstance( new EdgeCoverRingPerception())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EdgeCoverRingPerception::EdgeCoverRingPerception() :
      m_Rings()
    {
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EdgeCoverRingPerception::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Rings, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &EdgeCoverRingPerception::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Rings, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace graph
} // namespace bcl
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
#include "graph/bcl_graph_subgraph_isomorphism_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SubgraphIsomorphismBase::s_Instance
    (
      GetObjectInstances().AddInstance( new SubgraphIsomorphismBase())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! clone the object
    SubgraphIsomorphismBase *SubgraphIsomorphismBase::Clone() const
    {
      return new SubgraphIsomorphismBase( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @return class name of the object behind a pointer or the current object
    const std::string &SubgraphIsomorphismBase::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the isomorphisms between the sub and super graphs
    //! @return the isomorphisms
    const storage::Vector< storage::Vector< size_t> > &SubgraphIsomorphismBase::GetIsomorphisms() const
    {
      return m_Isomorphisms;
    }

    //! @brief Get the first subgraph isomorphism
    //! @return the first subgraph isomorphism
    const storage::Vector< size_t> &SubgraphIsomorphismBase::GetIsomorphism() const
    {
      static storage::Vector< size_t> empty;
      return m_Isomorphisms.IsEmpty() ? empty : m_Isomorphisms.FirstElement();
    }

    //! @brief Get the first subgraph isomorphism
    //! @return the first subgraph isomorphism
    const storage::Vector< size_t> SubgraphIsomorphismBase::GetInverseIsomorphism() const
    {
      static storage::Vector< size_t> empty;
      if( m_Isomorphisms.IsEmpty())
      {
        return empty;
      }

      const storage::Vector< size_t> &isomorphism( m_Isomorphisms.FirstElement());
      storage::Vector< size_t> inverse( isomorphism.GetSize());

      size_t count( 0);
      for
      (
        storage::Vector< size_t>::const_iterator itr( isomorphism.Begin()), itr_end( isomorphism.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        inverse( *itr) = count;
      }
      return inverse;
    }

    //! @brief set the isomorphisms
    //! @param ISOMORPHISMS the new isomorphisms to use
    void SubgraphIsomorphismBase::SetIsomorphisms
    (
      const storage::Vector< storage::Vector< size_t> > &ISOMORPHISMS
    )
    {
      m_Isomorphisms = ISOMORPHISMS;
    }

    //! @brief prune isomorphisms to those that are unique
    //! @param VERTICES_DESIRED_TO_DIFFER vertex indices of subgraph that should differ to consider an isomorphism unique
    //!        from one found already
    void SubgraphIsomorphismBase::PruneIsomorphismsToThoseUniqueInField( const storage::Vector< size_t> &VERTICES_THAT_SHOULD_DIFFER)
    {
      storage::Set< storage::Vector< size_t> > iso_uniq;
      auto itr_place( m_Isomorphisms.Begin()), itr_next( m_Isomorphisms.Begin());
      for( auto itr_end( m_Isomorphisms.End()); itr_next != itr_end; ++itr_next)
      {
        storage::Vector< size_t> iso_key( *itr_next);
        iso_key.Reorder( VERTICES_THAT_SHOULD_DIFFER);
        if( iso_uniq.Insert( iso_key).second)
        {
          if( itr_place != itr_next)
          {
            *itr_place = *itr_next;
          }
          ++itr_place;
        }
      }
      const size_t n_size( std::distance( m_Isomorphisms.Begin(), itr_place));
      m_Isomorphisms.Resize( n_size);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SubgraphIsomorphismBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Isomorphisms, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SubgraphIsomorphismBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Isomorphisms, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace graph
} // namespace bcl
