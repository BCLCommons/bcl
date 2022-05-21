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
#include "chemistry/bcl_chemistry_fragment_split_ecfp_fragments.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitECFPFragments::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitECFPFragments)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, sets default steps to 4
    FragmentSplitECFPFragments::FragmentSplitECFPFragments( const size_t &STEPS) :
        m_Steps( STEPS)
    {
    }

    //! @brief copy constructor
    //! @return a pointer to a copy of this class
    FragmentSplitECFPFragments *FragmentSplitECFPFragments::Clone() const
    {
      return new FragmentSplitECFPFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentSplitECFPFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitECFPFragments::GetAlias() const
    {
      static const std::string s_name( "ECFPFragments");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitECFPFragments::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the minimum size of fragments
    //! @return the minimum size of fragments
    const size_t FragmentSplitECFPFragments::GetMinSize() const
    {
      return 0;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitECFPFragments::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {

      //! This algorithm starts at each atom of a molecule and makes fragments that build out one bond at a time from
      //! that atom.  In other words, it begins at an atom, considers it a "fragment".  It then gets the neighbors of that
      //! atom (including the original atom) and considers those a "fragment".  It then steps out one more bond etc...
      //! until the molecule runs out of atoms or the desired number of steps has been achieved

      size_t vertices( MOLECULE_GRAPH.GetSize());
      util::OwnPtr< ConformationGraphConverter::t_AtomGraph> graph_ptr( &MOLECULE_GRAPH, false);

      storage::Set< storage::Vector< size_t> > index_set;
      storage::List< storage::Vector< size_t> > indices;

      size_t number_steps( std::min( m_Steps, vertices));

      // Begin fragmenting at each atom
      for( size_t vertex( 0); vertex < vertices; ++vertex)
      {
        // used vertices so that we don't repeat
        storage::Vector< size_t> vertex_used( vertices, size_t( 0));

        // fragment that will be grown out from this vertex
        storage::Vector< size_t> fragment;
        fragment.AllocateMemory( vertices);

        fragment.PushBack( vertex);
        vertex_used( vertex) = 1;

        for( size_t step( 0); step <= number_steps; ++step)
        {

          // Sort the vector by index, so that adding to the index_set will find identical vectors
          fragment.Sort( std::less< size_t>());

          // Insert it into the fragment set.  If we have seen it the fragment before any expansion that is done past this point
          // will also have been seen (by design of the algorithm), so break the loop here
          if( !index_set.Insert( fragment).second)
          {
            break;
          }

          // Add the fragment to the list of usable fragments
          indices.PushBack( fragment);

          // If we have maxed out our vertices then stop
          if( fragment.GetSize() == vertices)
          {
            break;
          }

          // Subgraph for getting the next layer of vertices
          graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> frag_subgraph
          (
            graph_ptr,
            fragment
          );

          // Add the next layer of vertices to the fragment
          storage::List< storage::Pair< size_t, size_t> > neighbors( frag_subgraph.GetAdjacentEdgeIndices());
          for
          (
            storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_neighbor( neighbors.Begin()),
              itr_neighbor_end( neighbors.End());
            itr_neighbor != itr_neighbor_end;
            ++itr_neighbor
          )
          {
            if( !vertex_used( itr_neighbor->Second()))
            {
              fragment.PushBack( itr_neighbor->Second());
              vertex_used( itr_neighbor->Second()) = 1;
            }
          }
        }
      }

      return indices;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitECFPFragments::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "splits molecules into fragments similar to those used for extended connectivity fingerprints "
        "(see http://pubs.acs.org/doi/abs/10.1021/ci100050t)"
      );

      parameters.AddInitializer
      (
        "steps",
        "how many bonds to go out from each atom",
        io::Serialization::GetAgent( &m_Steps),
        "4"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
