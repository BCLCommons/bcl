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
#include "chemistry/bcl_chemistry_fragment_split_linear_fragments.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitLinearFragments::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitLinearFragments)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, sets default steps to 4
    FragmentSplitLinearFragments::FragmentSplitLinearFragments( const size_t &STEPS) :
        m_Steps( STEPS)
    {
    }

    //! @brief copy constructor
    //! @return a pointer to a copy of this class
    FragmentSplitLinearFragments *FragmentSplitLinearFragments::Clone() const
    {
      return new FragmentSplitLinearFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentSplitLinearFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitLinearFragments::GetAlias() const
    {
      static const std::string s_name( "LinearFragments");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitLinearFragments::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the minimum size of fragments
    //! @return the minimum size of fragments
    const size_t FragmentSplitLinearFragments::GetMinSize() const
    {
      return 0;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief a recursive algorithm for getting the fragments of each molecule
    //! @param MOLECULE_GRAPH the atom graph of the molecule of interest
    //! @param CURRENT_CHAIN a set of indices that represent the current linear fragment up to the current vertex
    //! @param CURRENT_INDEX the index that should be added to CURRENT_CHAIN, if possible
    //! @param FRAGMENTS a set of vectors that represent all unique linear fragments that have been seen
    //! @param MAX_STEPS the maximum number of bonds to go out; default  = 7 akin to openbabel
    //! @param IGNORE_H whether to ignore hydrogens or not
    void FragmentSplitLinearFragments::GetFragments
    (
      const ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH,
      const storage::Set< size_t> &CURRENT_CHAIN,
      const size_t &CURRENT_INDEX,
      storage::Set< storage::Vector< size_t> > &FRAGMENTS,
      const size_t &MAX_STEPS,
      const bool &IGNORE_H
    ) const
    {
      // if the next atom is a hydrogen and it should be ignored, don't do anything
      if( IGNORE_H && MOLECULE_GRAPH.GetVertexData( CURRENT_INDEX)->GetElementType() == GetElementTypes().e_Hydrogen)
      {
        return;
      } 

      // copy the old chain to a new set, then add CURRENT_INDEX to it
      storage::Set< size_t> new_chain;
      new_chain.InsertElements( CURRENT_CHAIN.Begin(), CURRENT_CHAIN.End());
      new_chain.Insert( CURRENT_INDEX);

      // convert the set of indices into a vector and add it to the list of fragments
      storage::Vector< size_t> fragment( new_chain.Begin(), new_chain.End());

      // if the fragment has been seen, or we have gone out too many steps then stop here
      if( !FRAGMENTS.Insert( fragment).second || new_chain.GetSize() > MAX_STEPS)
      {
        return;
      }

      // get neighbor indices of the current index, and expand the fragment to each
      storage::Vector< size_t> neighbor_indices( MOLECULE_GRAPH.GetNeighborIndices( CURRENT_INDEX));
      for( size_t n( 0), end_n( neighbor_indices.GetSize()); n < end_n; ++n)
      {
        if( !CURRENT_CHAIN.Contains( neighbor_indices( n)))
        {
          GetFragments( MOLECULE_GRAPH, new_chain, neighbor_indices( n), FRAGMENTS, MAX_STEPS);
        }
      }
    }

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitLinearFragments::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {

      // the number of vertices
      size_t n_vertices( MOLECULE_GRAPH.GetSize());

      // a set of vectors that represent fragments of the molecule
      storage::Set< storage::Vector< size_t> > index_set;

      // Begin fragment generation at each atom
      for( size_t vertex( 0); vertex < n_vertices; ++vertex)
      {
        // start with an empty chain, and grow out from each vertex
        storage::Set< size_t> new_chain;
        GetFragments( MOLECULE_GRAPH, new_chain, vertex, index_set, m_Steps); 
      }

      // convert the set into a list so that it can be returned
      storage::List< storage::Vector< size_t> > indices( index_set.Begin(), index_set.End());

      return indices;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitLinearFragments::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Splits a molecule into a linear fragments (i.e. fragments that progress 1 atom out at a time without branching)."
        "  This implementation is similar to the openbabel FP2 fingerprint"
      );

      parameters.AddInitializer
      (
        "steps",
        "how many bonds to step out from each starting atom",
        io::Serialization::GetAgent( &m_Steps),
        "7"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
