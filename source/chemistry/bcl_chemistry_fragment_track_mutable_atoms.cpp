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
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "find/bcl_find_collector_interface.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "random/bcl_random_uniform_distribution.h"
#include "storage/bcl_storage_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // intiialize static member data
    FragmentComplete FragmentTrackMutableAtoms::s_MutableFragment = FragmentComplete();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentTrackMutableAtoms::FragmentTrackMutableAtoms()
    {
    }

    //! @brief clone constructor
    FragmentTrackMutableAtoms *FragmentTrackMutableAtoms::Clone() const
    {
      return new FragmentTrackMutableAtoms( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentTrackMutableAtoms::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief a function that returns the non-mutable base fragment provided a mutable fragment/atoms
    //! @param FRAGMENT the small molecule of interest
    //! @param MUTABLE_FRAGMENT sub-fragment of small molecule that can be mutated
    //! @param MUTABLE_ATOMS the atoms within the sub-fragment that can be mutated
    //! @param COMPLEMENT get the complement indices instead
    //! @return a fragment that can be modified by mutate interface classes
    FragmentComplete FragmentTrackMutableAtoms::GetBaseFragment
    (
      const FragmentComplete &FRAGMENT,
      const FragmentComplete &MUTATABLE_FRAGMENT,
      const storage::Vector< size_t> &MUTABLE_ATOMS,
      const bool COMPLEMENT,
      const ConformationGraphConverter &GRAPH_MAKER,
      const graph::CommonSubgraphIsomorphismBase::SolutionType &GRAPH_SOLUTION_TYPE
    )
    {
        if( MUTATABLE_FRAGMENT.GetSize() && MUTABLE_ATOMS.GetSize())
        {
          BCL_MessageStd
          (
            "You cannot simultaneously specify a mutable substructure and mutable atom indices. "
            "Please choose one, or specify neither and allow all atoms to be mutable"
          );
          return FragmentComplete();
        }

      // vertices to keep from original FRAGMENT
      storage::Set< size_t> unique_a_vertices;

      // vertices of the mutable fragment
      graph::Subgraph< size_t, size_t> complement_fragment_subgraph;

      // if we pass a mutable fragment, use to find base fragment
      if( MUTATABLE_FRAGMENT.GetSize())
      {
        // set fragment to our data member
        s_MutableFragment = MUTATABLE_FRAGMENT;
      }
      else if( MUTABLE_ATOMS.GetSize())
      {
        storage::Vector< size_t> mutable_atoms( MUTABLE_ATOMS);
        // make a fragment from the provided indices
        AtomVector< AtomComplete> mutable_frag_atoms( FRAGMENT.GetAtomVector());

        storage::Vector< size_t> atom_indices( storage::CreateIndexVector( FRAGMENT.GetSize()));
        size_t n_removed( 0);
        for( size_t i( 0); i < mutable_atoms.GetSize(); ++i, ++n_removed)
        {
          mutable_atoms( i) = mutable_atoms( i) - n_removed;
          atom_indices.RemoveElements( mutable_atoms( i), size_t( 1));
        }
        mutable_frag_atoms.Reorder( atom_indices);

        s_MutableFragment = FragmentComplete( mutable_frag_atoms, "");
        return s_MutableFragment;
      }
      else
      {
        // return empty molecule
        return FragmentComplete();
      }

      // make a graph of our input molecule
      graph::ConstGraph< size_t, size_t> fragment_graph( GRAPH_MAKER( FRAGMENT)), mutable_fragment_graph( GRAPH_MAKER( s_MutableFragment));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > fragment_graph_ptr( &fragment_graph, false); // for common subgraph isomorphism
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mutable_fragment_graph_ptr( &mutable_fragment_graph, false); // for common subgraph isomorphism

      // find the largest common subgraph with the mutable_fragment
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( GRAPH_SOLUTION_TYPE);
      isomorphism.SetGraphA( fragment_graph_ptr);
      isomorphism.SetGraphB( mutable_fragment_graph_ptr);
      isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

      // get the remainder of the input molecule after subtracting the common atoms from the mutable_fragment
      graph::Subgraph< size_t, size_t> fragment_subgraph
      (
        isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement()
      );

      // get the complement atom indices instead
      if( COMPLEMENT)
      {
        fragment_subgraph = fragment_subgraph.GetComplement();
      }

      // convert subgraph into small molecule
      unique_a_vertices.InsertElements( fragment_subgraph.GetVertexIndices().Begin(), fragment_subgraph.GetVertexIndices().End());
      AtomVector< AtomComplete> new_atom_vec( FRAGMENT.GetAtomVector());
      new_atom_vec.Reorder( storage::Vector< size_t>( unique_a_vertices.Begin(), unique_a_vertices.End()));
      FragmentComplete base_mol( new_atom_vec, "");
      return base_mol;
    }

    //! @brief a function that randomly picks a fragment to be the mutable region and returns the non-mutable base
    //! @param FRAGMENT the small molecule of interest
    //! @param SPLITTER the splitter to use to make the fragmnts
    //! @return a fragment that can be modified by mutate interface classes
    FragmentComplete FragmentTrackMutableAtoms::GetRandomBaseFragment
    (
      const FragmentComplete &FRAGMENT,
      const util::Implementation< FragmentSplitInterface> &SPLITTER
    )
    {
      // split fragment into ensemble components
      FragmentEnsemble split_frags;
      if( SPLITTER.IsDefined())
      {
        BCL_MessageStd( "Splitting molecule with " + util::Format()( SPLITTER->GetAlias()) + " splitter!");
        split_frags = ( *SPLITTER)( FRAGMENT);
      }
      else
      {
        BCL_MessageStd( "Splitter not defined; defaulting to Ring splitter");
        FragmentSplitRings ring_splitter( true, size_t( 3));
        split_frags = ring_splitter( FRAGMENT);
      }
      split_frags.SaturateWithH();

      // choose a random split fragment to be mutable
      split_frags.Shuffle();
      return GetBaseFragment( FRAGMENT, split_frags.GetMolecules().FirstElement(), storage::Vector< size_t>());
    }

    //! @brief a function that randomly picks a fragment from an allowed fragment/atoms
    //! to be the mutable region and returns the non-mutable base
    //! @param FRAGMENT the small molecule of interest
    //! @param SPLITTER the splitter to use to make the fragmnts
    //! @param MUTABLE_FRAGMENT sub-fragment of small molecule that can be mutated
    //! @param MUTABLE_ATOMS the atoms within the sub-fragment that can be mutated
    //! @return a fragment that can be modified by mutate interface classes
    FragmentComplete FragmentTrackMutableAtoms::GetRandomRestrictedBaseFragment
    (
      const FragmentComplete &FRAGMENT,
      const FragmentSplitInterface &SPLITTER,
      const FragmentComplete &MUTATABLE_FRAGMENT,
      const storage::Vector< size_t> &MUTABLE_ATOMS
    )
    {
      // TODO
      return FragmentComplete();
    }

    //! @brief a function that randomly picks an atom from a list of allowed element types for transformation
    //! @param FRAGMENT the small molecule of interest
    //! @param MUTABLE_ATOMS the atoms indicating the fragment that can be mutated
    //! @return a pointer to the selected atom
    storage::Vector< size_t> FragmentTrackMutableAtoms::SetMutableElements
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< ElementType> &MUTABLE_ELEMENTS
    )
    {
      // bail early
      if( !MUTABLE_ELEMENTS.GetSize())
      {
        return storage::Vector< size_t>();
      }

      // only consider elements that are both allowed and in our fragment
      storage::Vector< size_t> mutable_atoms;
      for
      (
          auto atom_itr( FRAGMENT.GetAtomVector().Begin()), atom_itr_end( FRAGMENT.GetAtomVector().End());
          atom_itr != atom_itr_end;
          ++atom_itr
      )
      {
        for( size_t e( 0), e_sz( MUTABLE_ELEMENTS.GetSize()); e < e_sz; ++e)
        {
          if( atom_itr->GetElementType() == MUTABLE_ELEMENTS( e))
          {
            // if this atom is an allowed element, include it in mutable atoms
            mutable_atoms.PushBack( FRAGMENT.GetAtomVector().GetAtomIndex( *atom_itr));
            break;
          }
        }
      }
      return mutable_atoms;
    }

    //! @brief a function that randomly picks an atom from a molecular substructure for transformation
    //! @param FRAGMENT the small molecule of interest
    //! @param REF_FRAGMENT the starting molecule
    //! @return a pointer to the selected atom
    storage::Vector< size_t> FragmentTrackMutableAtoms::SetMutableFragments
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const bool COMPLEMENT,
      const ConformationGraphConverter &GRAPH_MAKER,
      const graph::CommonSubgraphIsomorphismBase::SolutionType &GRAPH_SOLUTION_TYPE
    )
    {
      // bail early
      if( !MUTABLE_FRAGMENTS.GetSize())
      {
        return storage::Vector< size_t>();
      }

      // initialize container for unique atom indices across all fragments
      storage::Set< size_t> mutable_atoms_indices_unique;

      // create graph objects for fragment of interest
      graph::ConstGraph< size_t, size_t> fragment_graph( GRAPH_MAKER( FRAGMENT));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > fragment_graph_ptr( &fragment_graph, false);

      // prepare fragment isomorphism search with mutable fragments
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( GRAPH_SOLUTION_TYPE);
      isomorphism.SetGraphA( fragment_graph_ptr);

      // loop over all fragments and identify mutable atom indices from mutable fragments
      for
      (
          auto frag_itr( MUTABLE_FRAGMENTS.Begin()), frag_itr_end( MUTABLE_FRAGMENTS.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // make a graph of our input molecule
        graph::ConstGraph< size_t, size_t> mutable_fragment_graph( GRAPH_MAKER( *frag_itr));
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > mutable_fragment_graph_ptr( &mutable_fragment_graph, false);

        // find the largest common subgraph with the starting_fragment
        isomorphism.SetGraphB( mutable_fragment_graph_ptr);
        isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

        // get the remainder of the input molecule after subtracting the common atoms from the starting_fragment
        graph::Subgraph< size_t, size_t> fragment_subgraph
        (
          isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement()
        );

        // get the complement atom indices instead
        if( COMPLEMENT)
        {
          fragment_subgraph = fragment_subgraph.GetComplement();
        }

        // get the indices of the mutable atoms
        storage::Vector< size_t> mutable_atoms_indices( fragment_subgraph.GetVertexIndices());

        // collect the unique atoms from current fragment
        mutable_atoms_indices_unique.InsertElements( mutable_atoms_indices.Begin(), mutable_atoms_indices.End());
      }

      // return all unique mutable atom indices as a vector
      return storage::Vector< size_t>( mutable_atoms_indices_unique.Begin(), mutable_atoms_indices_unique.End());
    }

    //! @brief a function that randomly picks an atom for transformation
    //! @param FRAGMENT the small molecule of interest
    //! @param MUTABLE_ATOMS the atoms indicating the fragment that can be mutated
    //! @return a pointer to the selected atom
    util::SiPtr< const AtomConformationalInterface> FragmentTrackMutableAtoms::GetAtomFromMutable
    (
      const FragmentComplete &FRAGMENT,
      const bool ALL_ATOMS,
      const storage::Vector< size_t> &MUTABLE_ATOMS,
      const storage::Vector< ElementType> &MUTABLE_ELEMENTS,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &FIXED_ATOMS,
      const storage::Vector< ElementType> &FIXED_ELEMENTS,
      const FragmentEnsemble &FIXED_FRAGMENTS,
      const bool COMPLEMENT_MUTABLE_FRAGMENTS,
      const bool COMPLEMENT_FIXED_FRAGMENTS,
      const ConformationGraphConverter &MUTABLE_GRAPH_MAKER,
      const ConformationGraphConverter &FIXED_GRAPH_MAKER
    )
    {
      // ultimately we will grab atom from input fragment
      iterate::Generic< const AtomConformationalInterface> atom_vec( FRAGMENT.GetAtomsIterator());
      storage::Set< size_t> mutable_atoms_unique;

      // if we are setting all atoms to mutable then we can skip the other add steps
      if( ALL_ATOMS)
      {
        for( size_t i( 0), sz( FRAGMENT.GetSize()); i < sz; ++i)
        {
          mutable_atoms_unique.InsertElements( i);
        }
      }
      // otherwise add mutable atoms on a case-by-case basis
      else
      {
        // add individually specified mutable atoms
        mutable_atoms_unique.InsertElements( MUTABLE_ATOMS.Begin(), MUTABLE_ATOMS.End());
        // in addition to any mutable atoms passed as arguments, add atom indices of allowed mutable elements
        storage::Vector< size_t> mutable_elements( SetMutableElements( FRAGMENT, MUTABLE_ELEMENTS));
        if( mutable_elements.GetSize())
        {
          mutable_atoms_unique.InsertElements( mutable_elements.Begin(), mutable_elements.End());
        }
        // add atom indices of allowed mutable fragments
        storage::Vector< size_t> mutable_fragments
        (
          SetMutableFragments
          (
            FRAGMENT,
            MUTABLE_FRAGMENTS,
            COMPLEMENT_MUTABLE_FRAGMENTS,
            MUTABLE_GRAPH_MAKER
          )
        );
        if( mutable_fragments.GetSize())
        {
          mutable_atoms_unique.InsertElements( mutable_fragments.Begin(), mutable_fragments.End());
        }
      }

      // change the status of some mutable atoms to fixed
      if( mutable_atoms_unique.GetSize())
      {
        // remove individually specified fixed atoms from mutable atom list
        mutable_atoms_unique.EraseKeys( FIXED_ATOMS.Begin(), FIXED_ATOMS.End());
        // remove atom indices of fixed elements
        storage::Vector< size_t> fixed_elements( SetMutableElements( FRAGMENT, FIXED_ELEMENTS));
        if( fixed_elements.GetSize())
        {
          mutable_atoms_unique.EraseKeys( fixed_elements.Begin(), fixed_elements.End());
        }

        // remove atom indices of fixed fragments
        storage::Vector< size_t> fixed_fragments
        (
          SetMutableFragments
          (
            FRAGMENT,
            FIXED_FRAGMENTS,
            COMPLEMENT_FIXED_FRAGMENTS,
            FIXED_GRAPH_MAKER
          )
        );
        if( fixed_fragments.GetSize())
        {
          mutable_atoms_unique.EraseKeys( fixed_fragments.Begin(), fixed_fragments.End());
        }
      }

      // return a randomly selected atom from the available mutable atom indices
      storage::Vector< size_t> mutable_atoms( mutable_atoms_unique.Begin(), mutable_atoms_unique.End());
      mutable_atoms.Shuffle();
      atom_vec.GotoPosition( *mutable_atoms.Begin());
      return *atom_vec;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentTrackMutableAtoms::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentTrackMutableAtoms::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
