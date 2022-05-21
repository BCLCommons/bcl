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

#ifndef BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_BASE_H_
#define BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_BASE_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_const_graph.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CommonSubgraphIsomorphismBase
    //! @brief a common base class for all CommonSubgraphIsomorphism that provides data and access methods
    //!
    //! @see @link example_graph_common_subgraph_isomorphism_base.cpp @endlink
    //! @author mendenjl
    //! @date 10/29/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CommonSubgraphIsomorphismBase :
      public util::ObjectInterface
    {

    ///////////
    // Enums //
    ///////////

    public:

      // the type of solution desired
      enum SolutionType
      {
        e_Unconnected, //!< Unconnected - solution may contain components that are not connected
        e_Connected,   //!< Connected - solution will contain only the largest connected component
        e_GreedyUnconnected //!< Uses largest connected components iteratively to produce the unconnected isomorphism
                            //! This option is usually several orders of magnitude faster than Unconnected but may
                            //! occasionally miss a slightly larger isomorphism
      };

    protected:

    //////////
    // data //
    //////////

      //! all isomorphisms between graphs A and B, each pair connects isomorphic vertices from graph a to graph b
      storage::Vector< storage::Map< size_t, size_t> > m_Isomorphisms;

      //! true if the solution must be connected (i.e. all vertices in the isomorphism mutually reachable)
      SolutionType m_SolutionType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a solution type (connected by default)
      //! @param SOLUTION_TYPE whether the isomorphism must be connected
      CommonSubgraphIsomorphismBase( const SolutionType &SOLUTION_TYPE = e_Connected);

      //! clone the object
      CommonSubgraphIsomorphismBase *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const;

      //! @brief Get the isomorphisms between graphs A and B
      //! @return the isomorphisms
      const storage::Vector< storage::Map< size_t, size_t> > &GetIsomorphisms() const;

      //! @brief Get the largest common subgraph isomorphism between graphs A and B
      //! @return the largest common subgraph isomorphism
      const storage::Map< size_t, size_t> &GetIsomorphism() const;

      //! @brief Get whether this isomorphism had to be connected
      //! @return bool, true indicates this isomorphism had to be connected
      SolutionType GetSolutionType() const;

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      virtual void SetIsomorphisms( const storage::Vector< storage::Map< size_t, size_t> > &ISOMORPHISMS);

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      virtual void SetIsomorphism( const storage::Map< size_t, size_t> &ISOMORPHISM);

      //! @brief estimates the maximum size of the isomorphism by considering vertex types
      //! @note  edge colors are ignored
      //! @param GRAPH_A First graph
      //! @param GRAPH_B Second graph
      //! @return an estimate of the maximum size of the isomorphism
      //! @note EstimateUpperBounds assumes that vertex data can be sorted, which requires that
      //!       A<<B && B<C implies A<C.  If this is not the case for t_VertexData, then this function cannot
      //!       be used to estimate the maximum size of the isomorphism
      template< typename t_VertexData, typename t_EdgeData>
      static size_t EstimateUpperBounds
      (
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_A,
        const ConstGraph< t_VertexData, t_EdgeData> &GRAPH_B
      )
      {
        // get the vertex counts for graph a and b
        const storage::Map< t_VertexData, size_t> &vertex_data_a( GRAPH_A.GetVertexTypeCountMap());
        const storage::Map< t_VertexData, size_t> &vertex_data_b( GRAPH_B.GetVertexTypeCountMap());

        // make a variable to store how many atoms are in common in both graphs
        size_t sum_common_vertices( 0);

        // find out how many vertices the graphs have in common by walking through the map and
        // comparing counts
        // e.g. if graph a has 5 of a particular vertex, and graph b has 4, then the graphs can have no more than 4
        // of that vertex type in common.
        auto itr_b( vertex_data_b.Begin()), itr_b_end( vertex_data_b.End());
        for
        (
          typename storage::Map< t_VertexData, size_t>::const_iterator
            itr( vertex_data_a.Begin()), itr_end( vertex_data_a.End());
          itr != itr_end;
          ++itr
        )
        {
          while( itr_b != itr_b_end && itr_b->first < itr->first)
          {
            ++itr_b;
          }
          if( itr_b->first == itr->first)
          {
            sum_common_vertices += std::min( itr_b->second, itr->second);
          }
        }

        return sum_common_vertices;
      } // EstimateUpperBounds

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

      //! @brief swaps the isomorphisms
      //! thus makes indices in graph B be the keys to the isomorphism which map to indices in graph a
      void SwapIsomorphisms();

    };

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_COMMON_SUBGRAPH_ISOMORPHISM_BASE_H_

