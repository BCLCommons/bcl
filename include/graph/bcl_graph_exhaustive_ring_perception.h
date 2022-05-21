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

#ifndef BCL_GRAPH_EXHAUSTIVE_RING_PERCEPTION_H_
#define BCL_GRAPH_EXHAUSTIVE_RING_PERCEPTION_H_

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
    //! @class ExhaustiveRingPerception
    //! @brief determines all possible rings in a graph.
    //!
    //! @see Hanser, Th., Jauffret, Ph., Kaufmann, G., "A New Algorithm for Exhaustive Ring Perception in
    //! a Molecular Graph", Laboratoire de Modeles Informatiques Appliques a la Synthese, URA 405 du CNRS, Universite
    //! Louis Pasteur, 67000 Strasbourg, France, Received May 15, 1996
    //! @see @link example_graph_exhaustive_ring_perception.cpp @endlink
    //! @author mueller, mendenjl
    //! @date 01/11/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ExhaustiveRingPerception :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! paths between vertices in the graph
      storage::List< storage::Vector< size_t> > m_Paths;

      //! rings in the graph
      storage::List< storage::Vector< size_t> > m_Rings;

      //! size of the graph
      size_t m_GraphSize;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExhaustiveRingPerception();

      //! @brief constructor from a graph (either GraphWithData or ConstGraph)
      //! @param GRAPH graph for exhaustive ring detection
      //! @note prefer using ConstGraph here, it is often many times faster than GraphWithData
      ExhaustiveRingPerception( const GraphWithData< size_t, size_t> &GRAPH);
      ExhaustiveRingPerception( const ConstGraph< size_t, size_t> &GRAPH);

      //! clone the object
      ExhaustiveRingPerception *Clone() const
      {
        return new ExhaustiveRingPerception( *this);
      }

      //! @return class name of the object behind a pointer or the current object
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief Get the paths between vertices in the graph
      //! @return list of paths
      const storage::List< storage::Vector< size_t> > &GetPaths() const
      {
        return m_Paths;
      }

      //! @brief Get all rings in the graph
      //! @return list of rings (all closed paths)
      const storage::List< storage::Vector< size_t> > &GetRings() const
      {
        return m_Rings;
      }

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    /////////////
    // methods //
    /////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief initialize m_Paths with all edges from GRAPH (which is either a GraphWithData or a ConstGraph)
      void Initialize( const GraphWithData< size_t, size_t> &GRAPH);
      void Initialize( const ConstGraph< size_t, size_t> &GRAPH);

      //! @brief Remove vertices, store new paths and rings
      //! @param VERTEX vertex to be removed
      void Remove( const size_t VERTEX);

      //! @brief Combine two paths at their common vertex
      //! @param COMMON_VERTEX vertex to join paths at
      //! @param PATH_A first path to be combined
      //! @param PATH_B second path to be combined
      //! @return combined path
      storage::Vector< size_t> CombinePaths
      (
        const size_t COMMON_VERTEX,
        const storage::Vector< size_t> &PATH_A,
        const storage::Vector< size_t> &PATH_B
      ) const;

      //! @brief Check if two paths overlap
      //! @param VERTEX_IS_IN_PATH_A adjacency vector for path a
      //! @param PATH_B second path to be combined
      //! @return if paths overlap
      bool Overlap
      (
        const std::vector< bool> &VERTEX_IS_IN_PATH_A,
        const storage::Vector< size_t> &PATH_B
      ) const;

      //! @brief set the adjacency vector with a given path
      //! @param VERTEX_IS_IN_PATH the vector to setup such that VERTEX_IS_IN_PATH[ x] == true iff x is in PATH
      //! @param PATH the path to use in setting up the adjacency vector
      void SetupAdjacencyVector
      (
        std::vector< bool> &VERTEX_IS_IN_PATH,
        const storage::Vector< size_t> &PATH
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;
    }; // class ExhaustiveRingPerception

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_EXHAUSTIVE_RING_PERCEPTION_H_

