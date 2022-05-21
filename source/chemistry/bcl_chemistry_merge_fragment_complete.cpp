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
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_valence_handler.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MergeFragmentComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new MergeFragmentComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new MergeFragmentComplete
    MergeFragmentComplete *MergeFragmentComplete::Clone() const
    {
      return new MergeFragmentComplete( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MergeFragmentComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief merge two fragment given common vertices for the fragments
    //! @param MOLECULE_A small Molecule to be merged with MOLECULE_B
    //! @param MOLECULE_B small Molecule to be merged with MOLECULE_A
    //! @param COMMON_INDICES atoms that are common between MOLECULE_A (keys) and MOLECULE_B(mapped values)
    //! @param APPENDED_ATOMS container to store atoms of molecule B that have been appended to molecule A
    //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
    storage::Pair< bool, FragmentComplete> MergeFragmentComplete::MergeFragments
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Map< size_t, size_t> &COMMON_INDICES,
      storage::Vector< size_t> &APPENDED_ATOMS
    )
    {

      // check for empty molecule
      if( MOLECULE_A.GetNumberAtoms() == size_t( 0))
      {
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_B);
      }

      // check for empty molecule
      if( MOLECULE_B.GetNumberAtoms() == size_t( 0))
      {
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
      }

      // if common indices are empty then return atom a
      if( COMMON_INDICES.IsEmpty())
      {
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
      }

      // create molecule which will be assembled
      FragmentComplete assemble_molecule( MOLECULE_A);

      // get indices of atoms that are common to molecules A and B
      storage::Vector< size_t> values_b( COMMON_INDICES.GetMappedValues());

      // get atoms that are in b but not in a
      values_b.Sort( std::less< size_t>());
      values_b.PushBack( MOLECULE_B.GetNumberAtoms());
      storage::Vector< size_t>::const_iterator itr_b( values_b.Begin());

      storage::Vector< size_t> indices_in_b_not_a;
      indices_in_b_not_a.AllocateMemory( MOLECULE_B.GetNumberAtoms() - COMMON_INDICES.GetSize());
      for( size_t i( 0), size_b( MOLECULE_B.GetNumberAtoms()); i < size_b; ++i)
      {
        if( i == *itr_b)
        {
          ++itr_b;
        }
        else
        {
          indices_in_b_not_a.PushBack( i);
        }
      }

      // get graph for molecule b
      graph::ConstGraph< size_t, size_t> graph_b
      (
        ConformationGraphConverter
        (
          ConformationGraphConverter::e_AtomType,
          ConfigurationalBondTypeData::e_BondOrderOrAromatic
        )( MOLECULE_B)
      );

      // get atom graph of molecule b
      graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> atom_graph_b
      (
        ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE_B)
      );

      // iterate through common indices and get neighbors that are in common between molecule A and B or neighbors that
      // are not common
      for
      (
        storage::Map< size_t, size_t>::const_iterator
          itr_common_b( COMMON_INDICES.Begin()), itr_common_b_end( COMMON_INDICES.End());
        itr_common_b != itr_common_b_end;
        ++itr_common_b
      )
      {
        // get neighbors of common atom B
        const storage::Vector< size_t> &neighbors_common_b( graph_b.GetNeighborIndices( itr_common_b->second));

        // get neighbors of atoms that exist only in molecule b and exist only in molecule a
        storage::Vector< size_t> neighbors_in_only_b;
        storage::Vector< size_t> neighbors_in_only_a;
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_neighbor( neighbors_common_b.Begin()), itr_neighbor_end( neighbors_common_b.End());
          itr_neighbor != itr_neighbor_end;
          ++itr_neighbor
        )
        {
          // get atoms that are in b but not in a
          if( indices_in_b_not_a.Find( *itr_neighbor) < indices_in_b_not_a.GetSize())
          {
            neighbors_in_only_b.PushBack( *itr_neighbor);
          }
          else
          {
            neighbors_in_only_a.PushBack( *itr_neighbor);
          }
        }
        // if no neighbors exist that are contained only in molecule b
        if( neighbors_in_only_b.IsEmpty())
        {
          continue;
        }

        // get part of graph_b that is unique to molecule b connected to atom
        storage::Vector< size_t> edges_removed;
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_neighbor_a( neighbors_in_only_a.Begin()), itr_neighbor_a_end( neighbors_in_only_a.End());
          itr_neighbor_a != itr_neighbor_a_end;
          ++itr_neighbor_a
        )
        {
          edges_removed.PushBack( graph_b.GetEdgeData( itr_common_b->second, *itr_neighbor_a));
          graph_b.RemoveEdge( itr_common_b->second, *itr_neighbor_a);
        }

        storage::Vector< size_t> connected_vertices_b
        (
          *graph::Connectivity::GetVerticesReachableFrom( graph_b, itr_common_b->second)
        );

        // add the edges back to the graph
        storage::Vector< size_t>::const_iterator itr_edges_removed( edges_removed.Begin());
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_neighbor_a( neighbors_in_only_a.Begin()), itr_neighbor_a_end( neighbors_in_only_a.End());
          itr_neighbor_a != itr_neighbor_a_end;
          ++itr_neighbor_a, ++itr_edges_removed
        )
        {
          graph_b.AddEdge( itr_common_b->second, *itr_neighbor_a, *itr_edges_removed);
        }

        // get subgraph of fragment in the molecule of interest
        graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> subgraph
        (
          util::OwnPtr< graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> >( &atom_graph_b, false),
          connected_vertices_b
        );

        // get the graph and molecule from the atom subgraphs
        graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> graph_part( subgraph.ToGraph());

        // get the atoms itr_common_b->second and first neighbor from subgraph
        storage::Vector< util::SiPtr< const AtomConformationalInterface> > vertices_vector( graph_part.GetVertices());
        storage::Vector< size_t> neighbor_atom_b;
        size_t common_atom( util::GetUndefinedSize_t());
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_neigbors_in_only_b( neighbors_in_only_b.Begin()), itr_neigbors_in_only_b_end( neighbors_in_only_b.End());
          itr_neigbors_in_only_b != itr_neigbors_in_only_b_end;
          ++itr_neigbors_in_only_b
        )
        {
          size_t count( 0);
          for
          (
              storage::Vector< util::SiPtr< const AtomConformationalInterface> >::const_iterator
              itr_atoms( vertices_vector.Begin()), itr_atoms_end( vertices_vector.End());
              itr_atoms != itr_atoms_end;
              ++itr_atoms, ++count
          )
          {
            if( &**itr_atoms == &*atom_graph_b.GetVertexData( itr_common_b->second))
            {
              common_atom = count;
            }
            if( &**itr_atoms == &*atom_graph_b.GetVertexData( *itr_neigbors_in_only_b))
            {
              neighbor_atom_b.PushBack( count);
            }
          }
        }

        // get subgraph of fragment in the molecule of interest
        graph::Subgraph< size_t, size_t> subgraph_simple
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &graph_b, false), connected_vertices_b
        );

        graph::ConstGraph< size_t, size_t> subgraph_simple_part( subgraph_simple.ToGraph());
        storage::Vector< size_t>::const_iterator itr_real_neigh( neighbors_in_only_b.Begin());
        bool perpendicular_ring( false);
        for
        (
          storage::Vector< size_t>::const_iterator itr_neigh( neighbor_atom_b.Begin()), itr_neigh_end( neighbor_atom_b.End());
          itr_neigh != itr_neigh_end;
          ++itr_neigh, ++itr_real_neigh
        )
        {
          size_t cur_edge_data( subgraph_simple_part.GetEdgeData( *itr_neigh, common_atom));
          subgraph_simple_part.RemoveEdge( *itr_neigh, common_atom);
          storage::Vector< size_t> cur_conn_vertices( *graph::Connectivity::GetVerticesReachableFrom( subgraph_simple_part, *itr_neigh));
          subgraph_simple_part.AddEdge( *itr_neigh, common_atom, cur_edge_data);

          if( cur_conn_vertices.Find( common_atom) == cur_conn_vertices.GetSize())
          {
            cur_conn_vertices.PushBack( common_atom);
          }

          if( cur_conn_vertices.GetSize() == connected_vertices_b.GetSize())
          {
            storage::Vector< size_t> cur_conn_vertices_orig( cur_conn_vertices.GetSize());
            size_t vert_count( 0);
            for
            (
              storage::Vector< size_t>::const_iterator itr( cur_conn_vertices.Begin()), itr_end( cur_conn_vertices.End());
                itr != itr_end;
              ++itr, ++vert_count
            )
            {
              cur_conn_vertices_orig( vert_count) = connected_vertices_b( *itr);
            }

            graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> subgraph_second
            (
              util::OwnPtr< graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> >( &graph_part, false),
              cur_conn_vertices
            );
            // get the graph and molecule from the atom subgraphs
            graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> graph_part_second( subgraph_second.ToGraph());

            FragmentComplete part_second( ConformationGraphConverter::CreateAtomsFromGraph( graph_part_second), "");

            const AtomVector< AtomComplete> part_b_atomvec( part_second.GetAtomVector());

            // get the atoms itr_common_b->second and first neighbor from subgraph
            storage::Vector< util::SiPtr< const AtomConformationalInterface> > vertices_vector_second( graph_part_second.GetVertices());
            size_t neighbor_atom_b_second( util::GetUndefinedSize_t());
            size_t common_atom_second( util::GetUndefinedSize_t());
            size_t count_second( 0);
            for
            (
              storage::Vector< util::SiPtr< const AtomConformationalInterface> >::const_iterator
                itr_atoms( vertices_vector_second.Begin()), itr_atoms_end( vertices_vector_second.End());
              itr_atoms != itr_atoms_end;
              ++itr_atoms, ++count_second
            )
            {
              if( &**itr_atoms == &*graph_part.GetVertexData( common_atom))
              {
                common_atom_second = count_second;
              }
              if( &**itr_atoms == &*graph_part.GetVertexData( *itr_neigh))
              {
                neighbor_atom_b_second = count_second;
              }
            }
            if( part_second.GetAtomVector()( common_atom_second).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1))
            {
              perpendicular_ring = true;
              // get molecule from subgraph and atoms that need to be overlapped between this molecule-b part and molecule a
              storage::Map< size_t, size_t> overlapping_part;
              overlapping_part.Insert( storage::Pair< size_t, size_t>( itr_common_b->first, common_atom_second));

              // get the result of merger of molecule A and B
              storage::Pair< bool, FragmentComplete> merged_parts
              (
                MergeFragmentParts
                (
                  assemble_molecule,
                  part_second,
                  overlapping_part,
                  part_second.GetAtomVector()( neighbor_atom_b_second),
                  graph_b.GetEdgeData( itr_common_b->second, *itr_real_neigh),
                  perpendicular_ring
                )
              );

              if( merged_parts.First())
              {
                assemble_molecule = merged_parts.Second();
                // update the mapping between molecule which is being built and part that has already been built
                APPENDED_ATOMS.Append( cur_conn_vertices_orig);
              }
              else
              {
                return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
              }
              break;
            }
          }
        }
        if( !perpendicular_ring)
        {
          itr_real_neigh = neighbors_in_only_b.Begin();
          for
          (
            storage::Vector< size_t>::const_iterator itr_neigh( neighbor_atom_b.Begin()), itr_neigh_end( neighbor_atom_b.End());
            itr_neigh != itr_neigh_end;
            ++itr_neigh, ++itr_real_neigh
          )
          {
            size_t cur_edge_data( subgraph_simple_part.GetEdgeData( *itr_neigh, common_atom));
            subgraph_simple_part.RemoveEdge( *itr_neigh, common_atom);
            storage::Vector< size_t> cur_conn_vertices( *graph::Connectivity::GetVerticesReachableFrom( subgraph_simple_part, *itr_neigh));
            subgraph_simple_part.AddEdge( *itr_neigh, common_atom, cur_edge_data);

            if( cur_conn_vertices.Find( common_atom) == cur_conn_vertices.GetSize())
            {
              cur_conn_vertices.PushBack( common_atom);
            }

            storage::Vector< size_t> cur_conn_vertices_orig( cur_conn_vertices.GetSize());
            size_t vert_count( 0);
            for
            (
              storage::Vector< size_t>::const_iterator itr( cur_conn_vertices.Begin()), itr_end( cur_conn_vertices.End());
                itr != itr_end;
              ++itr, ++vert_count
            )
            {
              cur_conn_vertices_orig( vert_count) = connected_vertices_b( *itr);
            }

            graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> subgraph_second
            (
              util::OwnPtr< graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> >( &graph_part, false),
              cur_conn_vertices
            );
            // get the graph and molecule from the atom subgraphs
            graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> graph_part_second( subgraph_second.ToGraph());

            FragmentComplete part_second( ConformationGraphConverter::CreateAtomsFromGraph( graph_part_second), "");

            const AtomVector< AtomComplete> part_b_atomvec( part_second.GetAtomVector());

            // get the atoms itr_common_b->second and first neighbor from subgraph
            storage::Vector< util::SiPtr< const AtomConformationalInterface> > vertices_vector_second( graph_part_second.GetVertices());
            size_t neighbor_atom_b_second( util::GetUndefinedSize_t());
            size_t common_atom_second( util::GetUndefinedSize_t());
            size_t count_second( 0);
            for
            (
              storage::Vector< util::SiPtr< const AtomConformationalInterface> >::const_iterator
                itr_atoms( vertices_vector_second.Begin()), itr_atoms_end( vertices_vector_second.End());
              itr_atoms != itr_atoms_end;
              ++itr_atoms, ++count_second
            )
            {
              if( &**itr_atoms == &*graph_part.GetVertexData( common_atom))
              {
                common_atom_second = count_second;
              }
              if( &**itr_atoms == &*graph_part.GetVertexData( *itr_neigh))
              {
                neighbor_atom_b_second = count_second;
              }
            }

            // get molecule from subgraph and atoms that need to be overlapped between this molecule-b part and molecule a
            storage::Map< size_t, size_t> overlapping_part;
            overlapping_part.Insert( storage::Pair< size_t, size_t>( itr_common_b->first, common_atom_second));

            // get the result of merger of molecule A and B
            storage::Pair< bool, FragmentComplete> merged_parts
            (
              MergeFragmentParts
              (
                assemble_molecule,
                part_second,
                overlapping_part,
                part_second.GetAtomVector()( neighbor_atom_b_second),
                graph_b.GetEdgeData( itr_common_b->second, *itr_real_neigh),
                perpendicular_ring
              )
            );
            if( merged_parts.First())
            {
              assemble_molecule = merged_parts.Second();
              // update the mapping between molecule which is being built and part that has already been built
              APPENDED_ATOMS.Append( cur_conn_vertices_orig);
            }
            else
            {
              return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
            }
          }
        }
      }
      return storage::Pair< bool, FragmentComplete>( true, assemble_molecule);
    }

    //! @brief merge two fragment given common vertices for the fragments
    //! @param MOLECULE_A small Molecule to be merged with MOLECULE_B assembled
    //! @param MOLECULE_B small Molecule to be merged with MOLECULE_A assembled
    //! @param COMMON_INDICES atoms that are common between MOLECULE_A (keys) and MOLECULE_B(mapped values)
    //! @param NEIGHBOR_ATOM atom in molecule B which is bonded to an atom in molecule A (which is not in B)
    //! @param EDGE_DATA edge data between neighbor atom and an atom in A to which it is bonded
    //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
    storage::Pair< bool, FragmentComplete> MergeFragmentComplete::MergeFragmentParts
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Map< size_t, size_t> &COMMON_INDICES,
      const AtomComplete &NEIGHBOR_ATOM,
      const size_t EDGE_DATA,
      const bool PERPENDICULAR_RING
    )
    {
      // get atom vector for atom A
      AtomVector< AtomComplete> atom_vector_a( MOLECULE_A.GetAtomVector());

      // get atom indices that are common in a vector for both the molecules
      storage::Vector< size_t> common_atoms_a;
      storage::Vector< size_t> common_atoms_b;
      for
      (
        storage::Map< size_t, size_t>::const_iterator
          itr( COMMON_INDICES.Begin()), itr_end( COMMON_INDICES.End());
        itr != itr_end;
        ++itr
      )
      {
        common_atoms_a.PushBack( itr->first);
        common_atoms_b.PushBack( itr->second);
      }

      // get atom indices of molecule B that are not to be overlapped with atoms of molecule A
      storage::Vector< size_t> values_b( COMMON_INDICES.GetMappedValues());
      values_b.Sort( std::less< size_t>());
      values_b.PushBack( MOLECULE_B.GetNumberAtoms());
      storage::Vector< size_t>::const_iterator itr_b( values_b.Begin());

      storage::Vector< size_t> indices_in_b_not_a;
      indices_in_b_not_a.AllocateMemory( MOLECULE_B.GetNumberAtoms() - COMMON_INDICES.GetSize());
      for( size_t i( 0), size_b( MOLECULE_B.GetNumberAtoms()); i < size_b; ++i)
      {
        if( i == *itr_b)
        {
          ++itr_b;
        }
        else
        {
          indices_in_b_not_a.PushBack( i);
        }
      }

      // get new coordinates of molecule B after it has been overlapped with molecule A
      storage::Vector< linal::Vector3D> new_coordinates_b
      (
        GetTransformedCoordinates
        (
          MOLECULE_A,
          MOLECULE_B,
          common_atoms_a,
          common_atoms_b,
          NEIGHBOR_ATOM,
          EDGE_DATA,
          PERPENDICULAR_RING
        )
      );

      // if new_coordinates are empty then something went wrong, return molecule A
      if( new_coordinates_b.IsEmpty())
      {
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
      }

      // get atoms of b and set the coordinates of the atoms
      AtomVector< AtomComplete> atom_vector_b( MOLECULE_B.GetAtomVector());

      storage::Vector< linal::Vector3D>::const_iterator itr_coord( new_coordinates_b.Begin());
      for
      (
        AtomVector< AtomComplete>::iterator itr( atom_vector_b.Begin()), itr_end( atom_vector_b.End());
        itr != itr_end;
        ++itr, ++itr_coord
      )
      {
        itr->SetPosition( *itr_coord);
      }

      // create bond info for the merged molecule A and B. Bond info for A part remains the same but for B it will change
      storage::Vector< size_t> vector_b_hash( MOLECULE_B.GetNumberAtoms(), util::GetUndefined< size_t>());

      for
      (
        storage::Map< size_t, size_t>::const_iterator itr( COMMON_INDICES.Begin()), itr_end( COMMON_INDICES.End());
        itr != itr_end;
        ++itr
      )
      {
        vector_b_hash( itr->second) = itr->first;
      }

      size_t new_id_b( MOLECULE_A.GetNumberAtoms());

      for
      (
        storage::Vector< size_t>::iterator itr( vector_b_hash.Begin()), itr_end( vector_b_hash.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !util::IsDefined( *itr))
        {
          *itr = new_id_b++;
        }
      }

      storage::Vector< sdf::BondInfo> bond_info( atom_vector_b.GetBondInfo());

      storage::Vector< sdf::BondInfo> new_bond_info;
      new_bond_info.AllocateMemory( bond_info.GetSize());

      const size_t size_a( MOLECULE_A.GetNumberAtoms());
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator itr( bond_info.Begin()), itr_end( bond_info.End());
        itr != itr_end;
        ++itr
      )
      {
        size_t new_index_low( vector_b_hash( itr->GetAtomIndexLow()));
        size_t new_index_high( vector_b_hash( itr->GetAtomIndexHigh()));
        sdf::BondInfo updated_bond_info( new_index_low, new_index_high, itr->GetConfigurationalBondType());

        if( updated_bond_info.GetAtomIndexLow() >= size_a)
        {
          continue;
        }

        if( updated_bond_info.GetAtomIndexHigh() < size_a)
        {
          continue;
        }

        new_bond_info.PushBack( updated_bond_info);
      }
      atom_vector_b.Reorder( indices_in_b_not_a);
      atom_vector_a.AddAtomsWithConnectivity( atom_vector_b, new_bond_info);

      FragmentComplete new_molecule( atom_vector_a, "");
      return storage::Pair< bool, FragmentComplete>( true, new_molecule);
    }

    //! @brief connect two fragment given connection points for the fragments
    //! @param MOLECULE_A small Molecule to be connected with MOLECULE_B assembled
    //! @param MOLECULE_B small Molecule to be connected with MOLECULE_A assembled
    //! @param BOND_TYPE the bond that is to be used to connect
    //! @param INDICES_TO_CONNECT atoms of MOLECULE A (keys) that need to be connected to atoms of MOLECULE B(mapped values)
    //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
    storage::Pair< bool, FragmentComplete> MergeFragmentComplete::MergeFragments
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const ConfigurationalBondType &BOND_TYPE,
      const storage::Pair< size_t, size_t> &INDICES_TO_CONNECT
    )
    {
      // get transformed coordinates for molecule B so it can be bonded to molecule A
      storage::Vector< linal::Vector3D> new_coordinates_b
      (
        GetTransformedCoordinates
        (
          MOLECULE_A,
          MOLECULE_B,
          BOND_TYPE,
          INDICES_TO_CONNECT
        )
      );

      // if new coordinates is empty then return false and return an empty molecule
      if( new_coordinates_b.IsEmpty())
      {
        BCL_MessageStd( "Molecule B empty coords!");
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
      }

      // get atoms of molecule B and set their coordinates to the new coordinates
      AtomVector< AtomComplete> atom_vector_b( MOLECULE_B.GetAtomVector());

      storage::Vector< linal::Vector3D>::const_iterator itr_coord( new_coordinates_b.Begin());
      for
      (
        AtomVector< AtomComplete>::iterator itr( atom_vector_b.Begin()), itr_end( atom_vector_b.End());
        itr != itr_end;
        ++itr, ++itr_coord
      )
      {
        itr->SetPosition( *itr_coord);
      }

      storage::Vector< sdf::AtomInfo> atom_info_a( MOLECULE_A.GetAtomInfo());
      storage::Vector< sdf::AtomInfo> atom_info_b( atom_vector_b.GetAtomInfo());

      // create atom info for the merged molecules
      atom_info_a.Append( atom_info_b);

      // reassign bonds of molecule b since the atoms (indices) they connect change in the merged molecule
      storage::Vector< size_t> vector_b_hash( MOLECULE_B.GetNumberAtoms(), util::GetUndefined< size_t>());
      size_t new_id_b( MOLECULE_A.GetNumberAtoms());
      for
      (
        storage::Vector< size_t>::iterator itr( vector_b_hash.Begin()), itr_end( vector_b_hash.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !util::IsDefined( *itr))
        {
          *itr = new_id_b++;
        }
      }

      // if number of atoms in A and B are equal to 1 then just connect the atoms
      if( MOLECULE_A.GetNumberAtoms() == 1 && MOLECULE_B.GetNumberAtoms() == 1)
      {
        storage::Vector< sdf::BondInfo> bond_info( 1, sdf::BondInfo( 0, 1, BOND_TYPE));

        return storage::Pair< bool, FragmentComplete>( true, FragmentComplete( AtomVector< AtomComplete>( atom_info_a, bond_info), ""));
      }

      // if number of atoms in B is 1 then appending the single atom to molecule A is trivial
      if( MOLECULE_B.GetNumberAtoms() == 1)
      {
        storage::Vector< sdf::BondInfo> bond_info_a( MOLECULE_A.GetBondInfo());
        bond_info_a.PushBack( sdf::BondInfo( INDICES_TO_CONNECT.First(), MOLECULE_A.GetNumberAtoms(), BOND_TYPE));
        return storage::Pair< bool, FragmentComplete>( true, FragmentComplete( AtomVector< AtomComplete>( atom_info_a, bond_info_a), ""));
      }

      // if number of atoms in b is greater than 1 then when A and B are connected, atoms of B have to be added after
      // atoms of A and numbered accordingly.
      storage::Vector< sdf::BondInfo> bond_info_a;

      if( MOLECULE_A.GetNumberAtoms() == 1)
      {
        bond_info_a.PushBack( sdf::BondInfo( 0, vector_b_hash( INDICES_TO_CONNECT.Second()), BOND_TYPE));
      }
      else
      {
        bond_info_a = MOLECULE_A.GetBondInfo();
        bond_info_a.PushBack( sdf::BondInfo( INDICES_TO_CONNECT.First(), vector_b_hash( INDICES_TO_CONNECT.Second()), BOND_TYPE));
      }

      // update bonds after adding atoms of B to atoms of molecule A
      storage::Vector< sdf::BondInfo> bond_info( atom_vector_b.GetBondInfo());
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator itr( bond_info.Begin()), itr_end( bond_info.End());
        itr != itr_end;
        ++itr
      )
      {
        size_t new_index_low( vector_b_hash( itr->GetAtomIndexLow()));
        size_t new_index_high( vector_b_hash( itr->GetAtomIndexHigh()));
        sdf::BondInfo updated_bond_info( new_index_low, new_index_high, itr->GetConfigurationalBondType());

        bond_info_a.PushBack( updated_bond_info);
      }

      // construct atom_conformation using AtomVector
      AtomVector< AtomComplete> complete_molecule( atom_info_a, bond_info_a);

      // return the assembled molecule
      return storage::Pair< bool, FragmentComplete>( true, FragmentComplete( complete_molecule, ""));
    }

    //! @brief connect two fragment given multiple connection points for the fragments
    //! @details NOTE: this routine sets the coordinates of the merged fragment to ( 0, 0, 0)
    //! @param MOLECULE_A small Molecule to be connected with MOLECULE_B assembled
    //! @param MOLECULE_B small Molecule to be connected with MOLECULE_A assembled
    //! @param BOND_TYPES list of bond types used to connect atoms together
    //! @param INDICES_TO_CONNECT atoms of MOLECULE A (keys) that need to be connected to atoms of MOLECULE B(mapped values)
    //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
    storage::Pair< bool, FragmentComplete> MergeFragmentComplete::MultipointMergeFragments
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::List< ConfigurationalBondType> &BOND_TYPES,
      const storage::List< storage::Pair< size_t, size_t> > &INDICES_TO_CONNECT
    )
    {

      size_t bond_types_list_size( BOND_TYPES.GetSize());
      size_t indices_list_size( INDICES_TO_CONNECT.GetSize());

      // make sure there are the same number of bond types and index entries
      if( bond_types_list_size == 0 || indices_list_size == 0 || bond_types_list_size != indices_list_size)
      {
        BCL_MessageStd( "MergeFragments: bond type and index lists are not the same sizes");
        return storage::Pair< bool, FragmentComplete>( false, MOLECULE_A);
      }

      AtomVector< AtomComplete> mol_atoms( MOLECULE_A.GetAtomVector());

      storage::Vector< sdf::BondInfo> connecting_bonds;
      connecting_bonds.AllocateMemory( bond_types_list_size);

      size_t atom_num_offset( MOLECULE_A.GetNumberAtoms());
      storage::List< ConfigurationalBondType>::const_iterator itr_bond_type( BOND_TYPES.Begin());
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_bond( INDICES_TO_CONNECT.Begin()),
          itr_bond_end( INDICES_TO_CONNECT.End());
        itr_bond != itr_bond_end;
        ++itr_bond, ++itr_bond_type
      )
      {
        connecting_bonds.PushBack
        (
          sdf::BondInfo
          (
            itr_bond->First(),
            itr_bond->Second() + atom_num_offset,
            *itr_bond_type
          )
        );
      }

      mol_atoms.AddAtomsWithConnectivity
      (
        MOLECULE_B.GetAtomVector(),
        connecting_bonds
      );

      // Reset all of the positions to (0, 0, 0) since there is no sensible way to set coordinates here
      linal::Vector3D zero( 0.0, 0.0, 0.0);
      for( size_t atom_no( 0), end_atom_no( mol_atoms.GetSize()); atom_no < end_atom_no; ++atom_no)
      {
        mol_atoms( atom_no).SetPosition( zero);
      }

      return storage::Pair< bool, FragmentComplete>( true, FragmentComplete( mol_atoms, ""));
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
    std::istream &MergeFragmentComplete::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MergeFragmentComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief transformed coordinates of molecule B such that it can be merged with molecule A
    //! @param MOLECULE_A small molecule to be merged with MOLECULE_B assembled
    //! @param MOLECULE_B small molecule to be merged with MOLECULE_A assembled
    //! @param COMMON_ATOMS_A atoms of molecule A that are common with those in molecule B
    //! @param COMMON_ATOMS_B atoms of molecule B that are common with those in molecule A
    //! @param NEIGHBOR_ATOM atom in molecule B which is bonded to an atom in molecule A (which is not in B)
    //! @param EDGE_DATA edge data between neighbor atom and an atom in A to which it is bonded
    //! @return transformed coordinates of molecule B
    storage::Vector< linal::Vector3D> MergeFragmentComplete::GetTransformedCoordinates
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Vector< size_t> COMMON_ATOMS_A,
      const storage::Vector< size_t> COMMON_ATOMS_B,
      const AtomComplete &NEIGHBOR_ATOM,
      const size_t EDGE_DATA,
      const bool PERPENDICULAR_RING
    )
    {
      // get coordinates of molecule A and B
      const util::SiPtrVector< const linal::Vector3D> coordinates_a( MOLECULE_A.GetAtomCoordinates());
      const util::SiPtrVector< const linal::Vector3D> coordinates_b( MOLECULE_B.GetAtomCoordinates());

      storage::Vector< linal::Vector3D> transformed_coordinates_b;

      // if only one atom is common between two molecules that need to be merged
      // get atom which is common for both molecule A and B
      const AtomComplete &atom_of_interest_aa( MOLECULE_A.GetAtomVector()( COMMON_ATOMS_A.FirstElement()));
      const AtomComplete &atom_of_interest_ba( MOLECULE_B.GetAtomVector()( COMMON_ATOMS_B.FirstElement()));

      const double bond_length
      (
        BondLengths::GetBondLength( atom_of_interest_ba.GetAtomType(), EDGE_DATA, NEIGHBOR_ATOM.GetAtomType())
      );

      // determine valence coordinates for both the atoms
      storage::Vector< linal::Vector3D> valence_coord_ab( ValenceHandler::DetermineCoordinates( atom_of_interest_aa));
      storage::Vector< linal::Vector3D> valence_coord_bb( ValenceHandler::DetermineCoordinates( atom_of_interest_ba));

      // declare variables to hold coordinates of
      linal::Vector3D ideal_coord_ab;
      linal::Vector3D ideal_coord_bb;
      if( PERPENDICULAR_RING && atom_of_interest_aa.GetBonds().GetSize() > size_t( 1) && atom_of_interest_ba.GetBonds().GetSize() > size_t( 1))
      {
        linal::Vector3D first_atom( atom_of_interest_ba.GetBonds()( 0).GetTargetAtom().GetPosition());
        linal::Vector3D second_atom( atom_of_interest_ba.GetBonds()( 1).GetTargetAtom().GetPosition());

        util::SiPtrVector< const linal::Vector3D> molecule_b_coords;
        molecule_b_coords.PushBack( atom_of_interest_ba.GetPosition());
        molecule_b_coords.PushBack( first_atom);
        molecule_b_coords.PushBack( second_atom);
        double ring_angle( linal::ProjAngle( atom_of_interest_ba.GetPosition(), first_atom, second_atom));

        // helper coordinates
        linal::Vector3D foot_point_first
        (
          linal::CoordinatesTrigonal
          (
            atom_of_interest_aa.GetPosition(),
            atom_of_interest_aa.GetBonds()( 0).GetTargetAtom().GetPosition(),
            atom_of_interest_aa.GetBonds()( 1).GetTargetAtom().GetPosition(),
            linal::Distance( atom_of_interest_ba.GetPosition(), first_atom) * std::cos( 54.75 / 180 * math::g_Pi)
          )
        );
        // helper coordinates
        linal::Vector3D foot_point_second
        (
          linal::CoordinatesTrigonal
          (
            atom_of_interest_aa.GetPosition(),
            atom_of_interest_aa.GetBonds()( 0).GetTargetAtom().GetPosition(),
            atom_of_interest_aa.GetBonds()( 1).GetTargetAtom().GetPosition(),
            linal::Distance( atom_of_interest_ba.GetPosition(), second_atom) * std::cos( 54.75 / 180 * math::g_Pi)
          )
        );
        linal::Vector3D offset_first
        (
          linal::Distance( atom_of_interest_ba.GetPosition(), first_atom) * std::sin( ring_angle / 2.0) *
          linal::CrossProduct
          (
            atom_of_interest_aa.GetBonds()( 0).GetTargetAtom().GetPosition() - atom_of_interest_aa.GetPosition(),
            atom_of_interest_aa.GetBonds()( 1).GetTargetAtom().GetPosition() - atom_of_interest_aa.GetPosition()
          ).Normalize()
        );
        linal::Vector3D offset_second
        (
          linal::Distance( atom_of_interest_ba.GetPosition(), second_atom) * std::sin( ring_angle / 2.0) *
          linal::CrossProduct
          (
            atom_of_interest_aa.GetBonds()( 0).GetTargetAtom().GetPosition() - atom_of_interest_aa.GetPosition(),
            atom_of_interest_aa.GetBonds()( 1).GetTargetAtom().GetPosition() - atom_of_interest_aa.GetPosition()
          ).Normalize()
        );
        linal::Vector3D first_position( foot_point_first + offset_first);
        linal::Vector3D second_position( foot_point_second - offset_second);

        util::SiPtrVector< const linal::Vector3D> ideal_coords;
        ideal_coords.PushBack( atom_of_interest_aa.GetPosition());
        ideal_coords.PushBack( first_position);
        ideal_coords.PushBack( second_position);

        math::TransformationMatrix3D transform_matrix
        (
          quality::RMSD::SuperimposeCoordinates
          (
            ideal_coords, molecule_b_coords
          )
        );
        transformed_coordinates_b = TransformCoordinates( coordinates_b, transform_matrix);
        return transformed_coordinates_b;
      }
      if( !valence_coord_ab.IsEmpty())
      {
        ideal_coord_ab = valence_coord_ab.FirstElement();
        ideal_coord_ab = ( ( ideal_coord_ab - atom_of_interest_aa.GetPosition()) * bond_length) + atom_of_interest_aa.GetPosition();
      }

      if( !valence_coord_bb.IsEmpty())
      {
        ideal_coord_bb = valence_coord_bb.FirstElement();
        ideal_coord_bb = ( ( ideal_coord_bb - atom_of_interest_ba.GetPosition()) * bond_length) + atom_of_interest_ba.GetPosition();
      }

      // if common atom in molecule A is in ring
      if( atom_of_interest_aa.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1))
      {
        if( !valence_coord_ab.IsEmpty())
        {
          ideal_coord_bb = NEIGHBOR_ATOM.GetPosition();
        }
        else
        {
          return storage::Vector< linal::Vector3D>();
        }
      }
      // if common atom in molecule B is in ring
      else if( atom_of_interest_ba.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1))
      {
        if( !valence_coord_bb.IsEmpty())
        {
          ideal_coord_ab = atom_of_interest_aa.GetBonds().FirstElement().GetTargetAtom().GetPosition();
          double bond_len( linal::Distance( ideal_coord_ab, atom_of_interest_aa.GetPosition()));
          ideal_coord_bb = valence_coord_bb.FirstElement();
          ideal_coord_bb = ( ( ideal_coord_bb - atom_of_interest_ba.GetPosition()) * bond_len) + atom_of_interest_ba.GetPosition();

        }
        else
        {
          return storage::Vector< linal::Vector3D>();
        }
      }
      else
      {
        if( !valence_coord_bb.IsEmpty() && atom_of_interest_ba.GetBonds().GetSize() > 1)
        {
          ideal_coord_ab = atom_of_interest_aa.GetBonds().FirstElement().GetTargetAtom().GetPosition();
          double bond_len( linal::Distance( ideal_coord_ab, atom_of_interest_aa.GetPosition()));
          ideal_coord_bb = valence_coord_bb.FirstElement();
          ideal_coord_bb = ( ( ideal_coord_bb - atom_of_interest_ba.GetPosition()) * bond_len) + atom_of_interest_ba.GetPosition();
        }
        else if( !valence_coord_ab.IsEmpty() && atom_of_interest_ba.GetBonds().GetSize() == 1)
        {
          ideal_coord_bb = NEIGHBOR_ATOM.GetPosition();
        }
        else if( !valence_coord_bb.IsEmpty())
        {
          ideal_coord_ab = atom_of_interest_aa.GetBonds().FirstElement().GetTargetAtom().GetPosition();
          double bond_len( linal::Distance( ideal_coord_ab, atom_of_interest_aa.GetPosition()));
          ideal_coord_bb = valence_coord_bb.FirstElement();
          ideal_coord_bb = ( ( ideal_coord_bb - atom_of_interest_ba.GetPosition()) * bond_len) + atom_of_interest_ba.GetPosition();
        }
        else
        {
          return storage::Vector< linal::Vector3D>();
        }
      }
      math::TransformationMatrix3D transform_matrix
      (
        math::TransformationMatrix3D
        (
          coord::LineSegment3D( atom_of_interest_aa.GetPosition(), ideal_coord_ab),
          coord::LineSegment3D( atom_of_interest_ba.GetPosition(), ideal_coord_bb)
        )
      );
      transformed_coordinates_b = TransformCoordinates( coordinates_b, transform_matrix);

      return transformed_coordinates_b;
    }

    //! @brief transformed coordinates of molecule B such that it can be connected with molecule A
    //! @param MOLECULE_A small molecule to be merged with MOLECULE_B assembled
    //! @param MOLECULE_B small molecule to be merged with MOLECULE_A assembled
    //! @param BOND_TYPE bond type to be used to connect the molecules
    //! @param VERTICES_TO_CONNECT atoms of MOLECULE A (key) that need to be connected to atoms of MOLECULE A (mapped value)
    //! @return transformed coordinates of molecule B
    storage::Vector< linal::Vector3D> MergeFragmentComplete::GetTransformedCoordinates
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const ConfigurationalBondType &BOND_TYPE,
      const storage::Pair< size_t, size_t> &VERTICES_TO_CONNECT
    )
    {
      // get coordinates of molecules A and B
      const util::SiPtrVector< const linal::Vector3D> coordinates_a( MOLECULE_A.GetAtomCoordinates());
      const util::SiPtrVector< const linal::Vector3D> coordinates_b( MOLECULE_B.GetAtomCoordinates());

      // get coordinates that are defined and that are going to connect molecule A and B
      util::SiPtrVector< const linal::Vector3D> coordinates_a_copy( coordinates_a);
      coordinates_a_copy.Reorder( storage::Vector< size_t>( 1, VERTICES_TO_CONNECT.First()));

      // get coordinates that are defined and that are going to connect molecule A and B
      util::SiPtrVector< const linal::Vector3D> coordinates_b_copy( coordinates_b);
      coordinates_b_copy.Reorder( storage::Vector< size_t>( 1, VERTICES_TO_CONNECT.Second()));

      // create a vector to store transformed coordinates
      storage::Vector< linal::Vector3D> transformed_coordinates_b;

      // get atoms that need to be connected
      const AtomComplete &atom_of_interest_aa( MOLECULE_A.GetAtomVector()( VERTICES_TO_CONNECT.First()));
      const AtomComplete &atom_of_interest_bb( MOLECULE_B.GetAtomVector()( VERTICES_TO_CONNECT.Second()));

      double bond_length
      (
        BondLengths::GetBondLength
        (
          atom_of_interest_aa.GetAtomType(),
          BOND_TYPE->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
          atom_of_interest_bb.GetAtomType()
        )
      );

      // get where a new atom will be placed if it has to be added
      storage::Vector< linal::Vector3D> ideal_coords_ab( ValenceHandler::DetermineCoordinates( atom_of_interest_aa));
      storage::Vector< linal::Vector3D> ideal_coords_ba( ValenceHandler::DetermineCoordinates( atom_of_interest_bb));

      // if any of the molecules doesnt have valencies
      if( ideal_coords_ab.IsEmpty() || ideal_coords_ba.IsEmpty())
      {
        return storage::Vector< linal::Vector3D>();
      }

      linal::Vector3D ideal_coord_ab( ideal_coords_ab( random::GetGlobalRandom().Random( ideal_coords_ab.GetSize() - 1)));
      linal::Vector3D ideal_coord_ba( ideal_coords_ba( random::GetGlobalRandom().Random( ideal_coords_ba.GetSize() - 1)));

      ideal_coord_ba = ( ( ideal_coord_ba - atom_of_interest_bb.GetPosition()) * bond_length) + atom_of_interest_bb.GetPosition();
      ideal_coord_ab = ( ( ideal_coord_ab - atom_of_interest_aa.GetPosition()) * bond_length) + atom_of_interest_aa.GetPosition();

      // align the bond-valence coords for the two molecules
      transformed_coordinates_b = TransformCoordinates
          (
            coordinates_b,
            math::TransformationMatrix3D
            (
              coord::LineSegment3D( atom_of_interest_aa.GetPosition(), ideal_coord_ab),
              coord::LineSegment3D( ideal_coord_ba, atom_of_interest_bb.GetPosition())
            )
          );
      return transformed_coordinates_b;
    }

    //! @brief transformed coordinates after application of a given transformation matrix to a given set of coordinates
    //! @param COORDINATES coordinates that need to be transformed
    //! @param TRANSFORMATION_MATRIX transformation matrix that needs to be applied to coordinates of interest
    //! @return transformed coordinates
    storage::Vector< linal::Vector3D> MergeFragmentComplete::TransformCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
    )
    {
      storage::Vector< linal::Vector3D> transformed_coordinates;
      transformed_coordinates.AllocateMemory( COORDINATES.GetSize());

      for
      (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr_transform( COORDINATES.Begin()), itr_transform_end( COORDINATES.End());
          itr_transform != itr_transform_end;
          ++itr_transform
      )
      {
        linal::Vector3D curr_coord( **itr_transform);
        transformed_coordinates.PushBack( curr_coord.Transform( TRANSFORMATION_MATRIX));
      }
      return transformed_coordinates;
    }

    //! @brief transformed coordinates after application of a given transformation matrix to a given set of coordinates
    //! @param COORDINATES coordinates that need to be transformed
    //! @param TRANSFORMATION_MATRIX transformation matrix that needs to be applied to coordinates of interest
    //! @return transformed coordinates
    storage::Vector< linal::Vector3D> MergeFragmentComplete::TransformCoordinates
    (
      const storage::Vector< linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
    )
    {
      storage::Vector< linal::Vector3D> transformed_coordinates;
      transformed_coordinates.AllocateMemory( COORDINATES.GetSize());

      for
      (
          storage::Vector< linal::Vector3D>::const_iterator
          itr_transform( COORDINATES.Begin()), itr_transform_end( COORDINATES.End());
          itr_transform != itr_transform_end;
          ++itr_transform
      )
      {
        linal::Vector3D curr_coord( *itr_transform);
        transformed_coordinates.PushBack( curr_coord.Transform( TRANSFORMATION_MATRIX));
      }
      return transformed_coordinates;
    }

  } // namespace chemistry

} // namespace bcl
