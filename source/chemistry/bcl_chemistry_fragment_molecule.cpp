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
#include "chemistry/bcl_chemistry_fragment_molecule.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from list of FragmentConstitutionShared
    //! @param CONFORMATIONS shptr to conformation statistics
    FragmentMolecule::FragmentMolecule
    (
      util::ShPtr< ConfigurationSet> &CONFIGURATION, const size_t MAX_ROT
    )
    : m_ConfigurationSet( CONFIGURATION),
      m_MaxRot( MAX_ROT)
    {
    }

    //! @brief construct from list of FragmentConstitutionShared
    //! @param CONSTITUTION_SET shptr to constitution set
    FragmentMolecule::FragmentMolecule
    (
      util::ShPtr< ConstitutionSet> &CONSTITUTION_SET, const size_t MAX_ROT
    ) :
      m_ConstitutionSet( CONSTITUTION_SET),
      m_MaxRot( MAX_ROT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FragmentMolecule
    FragmentMolecule *FragmentMolecule::Clone() const
    {
      return new FragmentMolecule( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentMolecule::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentMolecule::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentMolecule::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    namespace
    {
      //! @brief remove bonds that are simple terminal bonds that hang off a ring. The positions of such atoms is trivial to
      //!        reconstruct
      //! @param MOL molecule to remove the unininteresting terminal bonds from rings from
      //! @return molecule without any terminal bonds off of rings
      FragmentComplete RemoveTerminalBondsOfRings( const FragmentComplete &FRAGMENT)
      {
        typename ConformationGraphConverter::t_AtomGraph atom_graph
        (
          ConformationGraphConverter::CreateGraphWithAtoms( FRAGMENT)
        );
        storage::Vector< size_t> indices_to_keep;
        indices_to_keep.AllocateMemory( FRAGMENT.GetSize());
        size_t i( 0);
        for( auto itr( FRAGMENT.GetAtomsIterator()); itr.NotAtEnd(); ++itr, ++i)
        {
          if
          (
            itr->GetBonds().GetSize() != size_t( 1)
          )
          {
            indices_to_keep.PushBack( i);
          }
          else if
          (
            itr->GetBonds().Begin()->GetTargetAtom().CountNonValenceBondsWithProperty
            (
              ConfigurationalBondTypeData::e_IsInRing,
              size_t( 1)
            ) >= size_t( 2)
          )
          {
            if
            (
              itr->GetBonds().Begin()->GetTargetAtom().CountNonValenceBondsWithProperty
              (
                ConfigurationalBondTypeData::e_IsAromatic,
                size_t( 1)
              ) < size_t( 2)
              && itr->GetAtomType() == GetAtomTypes().H_S
            )
            {
              // bonds off of non-aromatic rings. Only keep the hydrogens for exact substitution pattern matching, which
              // can influence the most common rotamer of the ring
              indices_to_keep.PushBack( i);
            }
          }
          else if
          (
            itr->GetAtomType() != GetAtomTypes().H_S
            || itr->GetBonds().Begin()->GetTargetAtom().GetBonds().GetSize()
               > itr->GetBonds().Begin()->GetTargetAtom().GetNumberCovalentlyBoundHydrogens() + size_t( 1)
          )
          {
            // skip terminal hydrogens to reduce computational complexity
            indices_to_keep.PushBack( i);
          }
        }
        return indices_to_keep.GetSize() == FRAGMENT.GetSize()
               ? FRAGMENT
               : FragmentComplete( ConformationGraphConverter::CreateAtomsFromGraph( atom_graph.GetSubgraph( indices_to_keep)), FRAGMENT.GetName());
      }
    }

    //! @brief fragment a molecule/fragment that is passed
    //! @param MOLECULE molecule that needs to be fragmented
    void FragmentMolecule::operator()( const FragmentComplete &MOLECULE)
    {
      FragmentComplete reduced( RemoveTerminalBondsOfRings( MOLECULE));
      if( reduced.GetSize() != MOLECULE.GetSize())
      {
        return operator()( reduced);
      }
      m_ConsiderConfigurations = m_ConfigurationSet.IsDefined();

      // skip molecules with less than 3 atoms
      if( MOLECULE.GetNumberAtoms() < 3)
      {
        return;
      }

      // create a graph of the molecule
      const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> molecule_graph
      (
        ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE)
      );

      // make a sub-fragment of the molecule which is passed in
      SubFragment molecule( MOLECULE);

      // make a storage list for storing sub-fragments
      storage::List< SubFragment> stored_fragments;
      stored_fragments.PushBack( molecule);

      // make a storage list for storing unique fragment which will be obtained from the molecule which is passed in
      storage::List< SubFragment> unique_fragments;
      unique_fragments.PushBack( MOLECULE);

      // unique fragments will be determined by sets of atoms of the MOLECULE from which fragment is constructed
      // make a storage set to store set of atom indices from which unique fragments are constructed
      storage::Set< storage::Set< size_t> > node_iso
      (
        storage::Set< size_t>
        (
          molecule.GetThisToNode().Begin(),
          molecule.GetThisToNode().End()
        )
      );

      double rotatable_bond_count( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( MOLECULE)( 0));
      if( rotatable_bond_count < m_MaxRot)
      {
        if( m_ConsiderConfigurations)
        {
          m_ConfigurationSet->Insert( FragmentConfigurationShared( MOLECULE)).second;
        }
        if( m_ConstitutionSet.IsDefined())
        {
          m_ConstitutionSet->Insert( FragmentConstitutionShared( MOLECULE)).second;
        }
      }

      // call algorithm which fragments the given MOLECULE
      Initialize
      (
        molecule_graph,
        stored_fragments,
        unique_fragments,
        node_iso
      );

    }

    //! @brief the function that calls the fragmentation algorithm
    //! @param FRAGMENTS list containing sub fragments
    //! @param UNIQUE_FRAGMENTS list containing unique FragmentComplete
    //! @param NODE_ISO a set containing atom indices of original molecule being fragmented, from which unique FragmentComplete are constructed
    //! @param return true if
    bool FragmentMolecule::Initialize
    (
      const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
      storage::List< SubFragment> &FRAGMENTS,
      storage::List< SubFragment> &UNIQUE_FRAGMENTS,
      storage::Set< storage::Set< size_t> > &NODE_ISO
    )
    {
      util::Stopwatch iso( "Fragment", util::Time( 1, 0), util::Message::e_Verbose, false, true);
      while( !FRAGMENTS.IsEmpty())
      {
        double stored_time( iso.GetTotalTime().GetSecondsFractional());
        if( stored_time > double( 1500))
        {
          BCL_MessageStd( "Skipping molecule after processing for : " + util::Format()( stored_time) + " seconds");
          break;
        }
        SubFragment first_subfragment( FRAGMENTS.FirstElement());
        const FragmentComplete molecule
        (
          ConformationGraphConverter::CreateAtomsFromGraph( MOLECULE_GRAPH.GetSubgraph( first_subfragment.GetThisToNode()), false),
          std::string()
        );
        // get list of breakable bonds of the subfragment
        storage::List< storage::Vector< sdf::BondInfo> > breakable_bonds( GetBreakableConnectivityPairs( molecule));
        // check for fragments with less than 3 atoms
        // call the function that does the fragmentation
        Fragmentation
        (
          MOLECULE_GRAPH,
          molecule,
          first_subfragment,
          breakable_bonds,
          FRAGMENTS,
          UNIQUE_FRAGMENTS,
          NODE_ISO
        );

        FRAGMENTS.PopFront();
      }
      return true;
    }

    namespace
    {
      void RemoveVectorDirect( storage::Vector< size_t> &TO_REMOVE_FROM, const storage::Vector< size_t> &REMOVE)
      {
        for( auto itr( REMOVE.Begin()), itr_end( REMOVE.End()); itr != itr_end; ++itr)
        {
          size_t find_pos( TO_REMOVE_FROM.Find( *itr));
          if( find_pos < TO_REMOVE_FROM.GetSize())
          {
            TO_REMOVE_FROM.RemoveElements( find_pos, size_t( 1));
          }
        }
      }
    }

    //! @brief fragment a molecule/fragment that is passed.
    //! @param MOLECULE the sub fragment which will be fragmented
    //! @param EDEGES_TO_BREAK a list of edges which are not in a ring and can be broken
    //! @param FRAGMENTS list containing sub fragments
    //! @param UNIQUE_FRAGMENTS list containing unique FragmentComplete
    //! @param NODE_ISO a set containing atom indices of original molecule being fragmented, from which unique FragmentComplete are constructed
    void FragmentMolecule::Fragmentation
    (
      const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
      const FragmentComplete &MOLECULE,
      const SubFragment &SUB_FRAGMENT,
      storage::List< storage::Vector< sdf::BondInfo> > &EDGES_TO_BREAK,
      storage::List< SubFragment> &FRAGMENTS,
      storage::List< SubFragment> &UNIQUE_FRAGMENTS,
      storage::Set< storage::Set< size_t> > &NODE_ISO
    )
    {
      // get graph if fragment exists
      if( !FRAGMENTS.IsEmpty())
      {
        graph::ConstGraph< size_t, size_t> fragment_graph( ConformationGraphConverter()( MOLECULE));
        storage::Vector< storage::Vector< size_t> > neighboring_h( MOLECULE.GetSize());
        for( size_t atom_i( 0), atom_e( MOLECULE.GetSize()); atom_i < atom_e; ++atom_i)
        {
          for
          (
            auto itr( MOLECULE.GetAtomVector()( atom_i).GetBonds().Begin()),
                 itr_end( MOLECULE.GetAtomVector()( atom_i).GetBonds().End());
            itr != itr_end;
            ++itr
          )
          {
            if( itr->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
            {
              neighboring_h( atom_i).PushBack( MOLECULE.GetAtomIndex( itr->GetTargetAtom()));
            }
          }
        }

        // make a container to store all the fragments that will be obtained by fragmentation of all bonds of the
        // SUB_FRAGMENT. This list is pruned to remove fragments that are completely within some other fragment(in the same list).
        storage::List< SubFragment> list_fragments;

        // go through each edge that can be broken to get fragments
        for
        (
          storage::List< storage::Vector< sdf::BondInfo> >::const_iterator
            itr( EDGES_TO_BREAK.Begin()), itr_end( EDGES_TO_BREAK.End());
          itr != itr_end;
          ++itr
        )
        {
          // if removing non ring bond
          if( itr->GetSize() == 1)
          {
            // get the indices of atoms in the bond
            const size_t first_atom_index( itr->FirstElement().GetAtomIndexLow());
            const size_t second_atom_index( itr->FirstElement().GetAtomIndexHigh());

            // remove the edge by passing indices of atoms
            fragment_graph.RemoveEdge( first_atom_index, second_atom_index);

            if( !fragment_graph.GetNeighborIndices( first_atom_index).IsEmpty())
            {
              // get fragment A which contains the first atom after edge is broken
              storage::Vector< size_t> vertices_reachable_first( *graph::Connectivity::GetVerticesReachableFrom( fragment_graph, first_atom_index));
              RemoveVectorDirect( vertices_reachable_first, neighboring_h( first_atom_index));

              // make a sub fragment of the fragment A
              SubFragment component_a( SUB_FRAGMENT, vertices_reachable_first);

              // insert the subfragment into list_fragment
              FragmentElimination( MOLECULE_GRAPH, list_fragments, component_a);
            }

            if( !fragment_graph.GetNeighborIndices( second_atom_index).IsEmpty())
            {
              // get fragment B which contains the second atom after edge is broken
              storage::Vector< size_t> vertices_reachable_second( *graph::Connectivity::GetVerticesReachableFrom( fragment_graph, second_atom_index));
              RemoveVectorDirect( vertices_reachable_second, neighboring_h( second_atom_index));

              // make a sub fragment of the fragment B
              SubFragment component_b( SUB_FRAGMENT, vertices_reachable_second);

              // insert the subfragment into list_fragment
              FragmentElimination( MOLECULE_GRAPH, list_fragments, component_b);
            }

            // add the edge which was removed back into the molecule
            fragment_graph.AddEdge( first_atom_index, second_atom_index, itr->FirstElement().GetConfigurationalBondType().GetIndex());
          }
          // removing non ring bonds to get substituents of rings
          else if( itr->GetSize() > 1)
          {
            // get the ring atom to which a non ring bond is attached. Remove the ring edges of the atom
            size_t first_low( itr->FirstElement().GetAtomIndexLow());
            size_t common_atom;
            if( first_low == itr->LastElement().GetAtomIndexLow() || first_low == itr->LastElement().GetAtomIndexHigh())
            {
              common_atom = first_low;
            }
            else
            {
              common_atom = itr->FirstElement().GetAtomIndexHigh();
            }

            // remove the ring edges and get the connected vertices to connected atom
            for
            (
              storage::Vector< sdf::BondInfo>::const_iterator
                itr_bonds( itr->Begin()), itr_bonds_end( itr->End());
              itr_bonds != itr_bonds_end;
              ++itr_bonds
            )
            {
              fragment_graph.RemoveEdge( itr_bonds->GetAtomIndexLow(), itr_bonds->GetAtomIndexHigh());
            }

            // get fragment A which contains the first atom after edge is broken
            storage::Vector< size_t> vertices_reachable_first( *graph::Connectivity::GetVerticesReachableFrom( fragment_graph, common_atom));
            RemoveVectorDirect( vertices_reachable_first, neighboring_h( common_atom));

            // make a sub fragment of the fragment A
            SubFragment component_a( SUB_FRAGMENT, vertices_reachable_first);

            // insert the subfragment into list_fragment
            FragmentElimination( MOLECULE_GRAPH, list_fragments, component_a);

            // add the edges back into the molecule
            for
            (
              storage::Vector< sdf::BondInfo>::const_iterator
              itr_bonds_add( itr->Begin()), itr_bonds_add_end( itr->End());
              itr_bonds_add != itr_bonds_add_end;
              ++itr_bonds_add
            )
            {
              fragment_graph.AddEdge
              (
                itr_bonds_add->GetAtomIndexLow(),
                itr_bonds_add->GetAtomIndexHigh(),
                itr_bonds_add->GetConfigurationalBondType().GetIndex()
              );
            }
          }
        }

        // insert fragments from the list_fragments to the list of unique fragments
        UniqueSetInsert( list_fragments, UNIQUE_FRAGMENTS, NODE_ISO);

        // append fragments obtained from SUB_FRAGMENT to the list of sub fragments
        FRAGMENTS.Append( list_fragments);
      }
    }

    namespace
    {
      //! @brief simple function to hash atom/bond type as a single size_t; Because this function only will deal with
      //!        terminal atom types, chirality is safely ignored
      size_t HashAtomAndTerminalChainBondType
      (
        const BondConformational &BOND,
        const bool &CONSIDER_CHIRALITY
      )
      {
        const size_t n_atom_types( GetAtomTypes().GetEnumCount());
        if( !CONSIDER_CHIRALITY)
        {
          return
            n_atom_types * BOND.GetBondType()->WithoutIsometry().GetIndex() + BOND.GetTargetAtom().GetAtomType().GetIndex();
        }
        return
            s_NumberChiralities * ( n_atom_types * BOND.GetBondType().GetIndex() + BOND.GetTargetAtom().GetAtomType().GetIndex())
            + size_t( BOND.GetTargetAtom().GetChirality());
      }
    }

    //! @brief get a list of bonds that need to be broken i.e. all non ring bonds and a pair of ring bonds( usually)
    //! @brief that need to broken to get bonds sticking out of a ring.
    //! @param MOLECULE the molecule to be fragmented
    //! @return list of bonds that need be broken
    storage::List< storage::Vector< sdf::BondInfo> > FragmentMolecule::GetBreakableConnectivityPairs
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // get the bond information
      storage::Vector< sdf::BondInfo> con_pairs( MOLECULE.GetBondInfo());

      // create a list to store bonds that can be broken
      storage::List< sdf::BondInfo> con_tuple_tobreak_in_chains;

      // create a vector to store bond info for each atom of a molecule, separating ring and chain bonds
      storage::Vector< storage::Vector< sdf::BondInfo> > chain_bond_info_vector( MOLECULE.GetNumberAtoms());
      storage::Vector< storage::Vector< sdf::BondInfo> > ring_bond_info_vector( MOLECULE.GetNumberAtoms());

      // pushback bonds that are not in a ring
      for( size_t i( 0), n_bonds( con_pairs.GetSize()); i < n_bonds; ++i)
      {
        const sdf::BondInfo &con_pair( con_pairs( i));
        // never break hydrogens
        if( MOLECULE.GetAtomVector()( con_pair.GetAtomIndexHigh()).GetAtomType() == GetAtomTypes().H_S)
        {
          continue;
        }
        if( MOLECULE.GetAtomVector()( con_pair.GetAtomIndexLow()).GetAtomType() == GetAtomTypes().H_S)
        {
          continue;
        }
        if( !con_pair.GetConfigurationalBondType()->IsBondInRing())
        {
          con_tuple_tobreak_in_chains.PushBack( con_pair);
          chain_bond_info_vector( con_pair.GetAtomIndexLow()).PushBack( con_pair);
          chain_bond_info_vector( con_pair.GetAtomIndexHigh()).PushBack( con_pair);
        }
        else
        {
          ring_bond_info_vector( con_pair.GetAtomIndexLow()).PushBack( con_pair);
          ring_bond_info_vector( con_pair.GetAtomIndexHigh()).PushBack( con_pair);
        }
      }

      // remove all but one equivalent bonds to atoms that are terminal
      storage::Set< sdf::BondInfo> bonds_to_ignore;
      size_t symmetries_product( 1);
      for( size_t i( 0), n_atoms( MOLECULE.GetNumberAtoms()); i < n_atoms; ++i)
      {
        // always skip terminal atoms or atoms that are completely in-rings, or atoms with only one chain substituent
        if( chain_bond_info_vector( i).GetSize() <= size_t( 1))
        {
          continue;
        }

        // create a map from bondtype, atomtype to count

        // indices of atoms with each bond type hash
        storage::Map< size_t, storage::Vector< sdf::BondInfo> > atom_bond_types_by_hash;
        size_t bond_id( 0);
        size_t n_inserted( 0);
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bond( MOLECULE.GetAtomVector()( i).GetBonds().Begin()),
            itr_bond_end( MOLECULE.GetAtomVector()( i).GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond, ++bond_id
        )
        {
          const ConfigurationalBondType bond_type( itr_bond->GetBondType());
          if( bond_type->IsBondInRing())
          {
            continue;
          }
          const size_t target_atom_index( MOLECULE.GetAtomIndex( itr_bond->GetTargetAtom()));
          if
          (
            chain_bond_info_vector( target_atom_index).GetSize() != size_t( 1)
            || !ring_bond_info_vector( target_atom_index).IsEmpty()
          )
          {
            continue;
          }
          ++n_inserted;
          // create hash for the atom / bond type
          size_t bond_hash( HashAtomAndTerminalChainBondType( *itr_bond, m_ConsiderConfigurations));
          atom_bond_types_by_hash[ bond_hash].PushBack
          (
            sdf::BondInfo( i, target_atom_index, bond_type)
          );
        }

        // if all the hashes were unique, continue
        if( n_inserted <= size_t( 1) || atom_bond_types_by_hash.GetSize() == n_inserted)
        {
          continue;
        }

        // find vectors in the map that have size > 1
        for
        (
          storage::Map< size_t, storage::Vector< sdf::BondInfo> >::const_iterator
            itr_map( atom_bond_types_by_hash.Begin()), itr_map_end( atom_bond_types_by_hash.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          if( itr_map->second.GetSize() == size_t( 1))
          {
            continue;
          }
          // only retain the first element of the vector, remove others from the list of connectivity pairs
          bonds_to_ignore.InsertElements( itr_map->second.Begin() + 1, itr_map->second.End());
          symmetries_product *= itr_map->second.GetSize();
        }
      }
      storage::List< storage::Vector< sdf::BondInfo> > con_tuple_tobreak;
      // iterate through list of bond info, for each bond not in bonds to ignore, add it to its own vector in
      // con_tuple to break
      for
      (
        storage::List< sdf::BondInfo>::const_iterator
          itr_bond( con_tuple_tobreak_in_chains.Begin()), itr_bond_end( con_tuple_tobreak_in_chains.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        if( !bonds_to_ignore.Contains( *itr_bond))
        {
          con_tuple_tobreak.PushBack( storage::Vector< sdf::BondInfo>( size_t( 1), *itr_bond));
        }
      }

      // now iterate through bonds of each atom and get ring bonds or non ring bonds an atom has
      for
      (
        storage::Vector< storage::Vector< sdf::BondInfo> >::const_iterator
          itr_ring( ring_bond_info_vector.Begin()),
          itr_ring_end( ring_bond_info_vector.End()),
          itr_chain( chain_bond_info_vector.Begin());
        itr_ring != itr_ring_end;
        ++itr_ring, ++itr_chain
      )
      {
        // for atoms with both ring and chain bonds, allow removal of the ring
        if( !itr_chain->IsEmpty() && !itr_ring->IsEmpty())
        {
          // if the central atom is ring, remove ring bonds
          con_tuple_tobreak.PushBack( *itr_ring);
        }
      }

      return con_tuple_tobreak;
    }

    //! @brief inserts fragments from list of sub fragments into list of unique_fragments and update fragment vertices
    //! that have been seen
    //! @param FRAGMENTS list of sub fragments that need to be fragmented
    //! @param UNIQUE_FRAGMENTS list to store unique fragment in FragmentComplete form
    //! @param NODE_ISO set to store set of vertices(isomorphic to fragment) which have been seen in fragments
    void FragmentMolecule::UniqueSetInsert
    (
      storage::List< SubFragment> &FRAGMENTS,
      storage::List< SubFragment> &UNIQUE_FRAGMENTS,
      storage::Set< storage::Set< size_t> > &NODE_ISO
    )
    {
      // insert sub fragments into a list of unique_fragments if it hasn't been seen earlier otherwise remove it from
      // the list of sub fragments
      for
      (
        storage::List< SubFragment>::iterator itr( FRAGMENTS.Begin()), itr_end( FRAGMENTS.End());
        itr != itr_end;
      )
      {
        if( !NODE_ISO.Contains( itr->GetThisToNodeSet()))
        {
          UNIQUE_FRAGMENTS.PushBack( *itr);
          NODE_ISO.Insert( itr->GetThisToNodeSet());
          ++itr;
        }
        else
        {
          itr = FRAGMENTS.Remove( itr);
        }
      }
      return;
    }

    //! @brief function performs checks to see if subfragment passed in contained in existing sub fragments
    //! @param FRAGMENT_LIST list containing sub fragments which are not part of any other fragment in the list
    //! @param FRAGMENT sub fragment which need to be checked for being contained in any other fragments in the fragment list
    void FragmentMolecule::FragmentElimination
    (
      const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
      storage::List< SubFragment> &FRAGMENT_LIST,
      SubFragment &FRAGMENT
    )
    {
      // if fragment has less than 2 atoms then no need to add it to the list
      if( FRAGMENT.GetThisToNode().GetSize() < 4)
      {
        return;
      }
      // get the isomorphism of the fragment to its parent molecule for checking if the same set of molecules have been
      // seen in earlier fragments
      const storage::Set< size_t> &new_set( FRAGMENT.GetThisToParentSet());

      // set a bool to true to keep track of whether to insert the fragment into the list or not
      bool insert( true);

      // check whether vertex of the new sub fragment is contained in existing sub fragments
      for
      (
        storage::List< SubFragment>::iterator itr( FRAGMENT_LIST.Begin()), itr_end( FRAGMENT_LIST.End());
        itr != itr_end;
      )
      {
        // get the vertex isomorphism of sub fragment to the parent fragment
        const storage::Set< size_t> &this_set( itr->GetThisToParentSet());

        // check for set equality, if equal that mean fragment has already been seen
        if( this_set == new_set)
        {
          insert = false;
          break;
        }

        // if sizes of sets are equal that means none of the fragments is contained in another
        else if( this_set.GetSize() == new_set.GetSize())
        {
          ++itr;
          continue;
        }

        // if one fragment is smaller than other, then check if any one contains the other
        // create a vector to store vertex intersection of two fragments being compared
        std::vector< size_t> v( std::min( this_set.GetSize(), new_set.GetSize()));

        std::vector< size_t>::iterator it;

        // find set intersection
        it = std::set_intersection( this_set.Begin(), this_set.End(), new_set.Begin(), new_set.End(), v.begin());

        // get the intersection size
        size_t intersection_size( it - v.begin());

        // if intersection size equal to the size of any fragment being compared, that means the fragment is completely
        // contained in the other
        if( intersection_size == this_set.GetSize())
        {
          itr = FRAGMENT_LIST.Remove( itr);
          continue;
        }
        else if( intersection_size == new_set.GetSize())
        {
          insert = false;
          break;
        }
        ++itr;
      }

      if( !insert)
      {
        return;
      }

      // cache the index to the fragment before inserting it into conformation set
      FragmentComplete molecule
      (
        ConformationGraphConverter::CreateAtomsFromGraph( MOLECULE_GRAPH.GetSubgraph( FRAGMENT.GetThisToNode()), false),
        std::string()
      );

      // skip molecules < four atoms; which could not have a dihedral
      if( molecule.GetNumberAtoms() < 4)
      {
        return;
      }

      PriorityDihedralAngles pda;
      auto dih_edges( pda.GetDihedralEdges( molecule));
      size_t n_dih( 0);
      for( auto itr_di( dih_edges.Begin()), itr_dih_end( dih_edges.End()); itr_di != itr_dih_end; ++itr_di)
      {
        if( !itr_di->GetEdgeData()->IsBondInRing())
        {
          ++n_dih;
        }
      }

      double rotatable_bond_count( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( molecule)( 0));

      if( n_dih <= size_t( 4) && rotatable_bond_count < m_MaxRot)
      {
        bool is_new_scaffold( false);
        if( m_ConsiderConfigurations)
        {
          is_new_scaffold = m_ConfigurationSet->Insert( FragmentConfigurationShared( molecule)).second;
        }
        if( m_ConstitutionSet.IsDefined())
        {
          is_new_scaffold = m_ConstitutionSet->Insert( FragmentConstitutionShared( molecule)).second;
        }
        // if insert is set to true that means the fragment is unique and is not completely contained in any other fragment
        // in the list
        if( is_new_scaffold)
        {
          FRAGMENT_LIST.PushBack( FRAGMENT);
        }
      }
      else
      {
        FRAGMENT_LIST.PushBack( FRAGMENT);
      }

    }

  } // namespace chemistry
} // namespace bcl
