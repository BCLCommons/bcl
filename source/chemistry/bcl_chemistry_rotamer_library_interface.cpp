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
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @return returns default rotamer library to be used
    const std::string &RotamerLibraryInterface::GetDefault()
    {
      static std::string s_default( "cod");
      return s_default;
    }

    namespace
    {
      //! @brief given two sorted vectors, test that one is a subset of the other
      //! @param TEST_SUBSET the potential subset of test-set
      //! @param TEST_SET the potential super-set of test-subset
      //! @return rue if every element in TEST_SUBSET is in TEST_SET under the precondition that both vectors are sorted initially
      bool IsSubsetOfForSortedVectors( const storage::Vector< size_t> &TEST_SUBSET, const storage::Vector< size_t> &TEST_SET)
      {
        if( TEST_SUBSET.GetSize() > TEST_SET.GetSize())
        {
          return false;
        }
        if( TEST_SUBSET.GetSize() == TEST_SET.GetSize())
        {
          return TEST_SUBSET == TEST_SET;
        }
        // handle trivial case of empty subset
        if( TEST_SUBSET.IsEmpty())
        {
          return true;
        }

        // handle case that the last element of the subset is larger
        // This saves the need to check for the end of TEST_SET, below
        if( TEST_SUBSET.LastElement() > TEST_SET.LastElement())
        {
          return false;
        }

        auto itr_set( TEST_SET.Begin());
        for
        (
          auto itr_subset( TEST_SUBSET.Begin()), itr_subset_end( TEST_SUBSET.End());
          itr_subset != itr_subset_end;
          ++itr_subset, ++itr_set
        )
        {
          while( *itr_set < *itr_subset)
          {
            ++itr_set;
          }
          if( *itr_set > *itr_subset)
          {
            return false;
          }
        }
        return true;
      }

      void SortConformationsBySizeThenScaffoldCount( FragmentEnsemble &CONFORMATIONS)
      {
        std::vector< std::pair< std::pair< size_t, size_t>, storage::List< FragmentComplete>::iterator> > ens_list;
        ens_list.reserve( CONFORMATIONS.GetSize());
        for( auto itr( CONFORMATIONS.Begin()), itr_end( CONFORMATIONS.End()); itr != itr_end; ++itr)
        {
          ens_list.push_back
          (
            std::make_pair
            (
              std::pair< size_t, size_t>( itr->GetNumberAtoms(), itr->GetStoredProperties().GetMDLPropertyAsVector( "ScaffoldCount")( 0)),
              itr
            )
          );
        }
        std::stable_sort
        (
          ens_list.begin(),
          ens_list.end(),
          [](
              const std::pair< std::pair< size_t, size_t>, storage::List< FragmentComplete>::iterator> &A,
              const std::pair< std::pair< size_t, size_t>, storage::List< FragmentComplete>::iterator> &B
          )->bool
          {
            return A.first.first < B.first.first
                   || ( A.first.first == B.first.first && A.first.second > B.first.second);
          }
        );
        std::list< FragmentComplete> new_frag_ensemble;
        for( auto itr( ens_list.begin()), itr_end( ens_list.end()); itr != itr_end; ++itr)
        {
          new_frag_ensemble.splice( new_frag_ensemble.end(), CONFORMATIONS.GetMolecules().InternalData(), itr->second);
        }
        std::swap( new_frag_ensemble, CONFORMATIONS.GetMolecules().InternalData());
      }
    }

    //! @brief create the substructure tree of unique constitutions for the  given rotamer library
    //! @param CONFORMATIONS
    void RotamerLibraryInterface::Create( FragmentEnsemble &CONFORMATIONS) const
    {
      const size_t conformations_size( CONFORMATIONS.GetSize());

      // Sort the conformations, first by number of atoms (smallest to largest) then by scaffold count (largest to smallest)
      // This ensures that root nodes are earlier in the list than child nodes, and that the list is generally in order of
      // decreasing frequency among molecules of the same size
      SortConformationsBySizeThenScaffoldCount( CONFORMATIONS);

      // create a constitution set to get unique constitutions
      ConstitutionSet constitution_set;
      const size_t n_conformations( CONFORMATIONS.GetSize());
      BCL_MessageStd( "Finding all constitutions");
      size_t cntr( 0);
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_conformations( CONFORMATIONS.Begin()), itr_conformations_end( CONFORMATIONS.End());
        itr_conformations != itr_conformations_end;
        ++itr_conformations, ++cntr
      )
      {
        FragmentConstitutionShared constitution( *itr_conformations);
        constitution_set.Insert( constitution);
        util::GetLogger().LogStatus
        (
          "Finding all constitutions " + util::Format()( cntr) + " / " + util::Format()( n_conformations)
        );
      }

      // get unique constitutions from the given rotamer library
      const util::ShPtrList< FragmentConstitutionShared> &unique_constitutions( constitution_set.GetConstitutions());

      BCL_MessageStd( "Creating graphs for each constitution");

      // graph maker to create graphs of unique constitutions which will be used for substructure search
      ConstitutionGraphConverter graph_maker;

      // allocate memory to store graphs
      const size_t n_uniq_constitutions( unique_constitutions.GetSize());
      storage::Vector< graph::ConstGraph< size_t, size_t> > search_substructures;
      search_substructures.AllocateMemory( n_uniq_constitutions);
      storage::Vector< size_t> scaffold_counts;
      scaffold_counts.AllocateMemory( n_uniq_constitutions);

      // create and store graphs of unique constitutions in the vector
      for
      (
        util::ShPtrList< FragmentConstitutionShared>::const_iterator
          itr_unique_constitutions( unique_constitutions.Begin()), itr_unique_constitutions_end( unique_constitutions.End());
        itr_unique_constitutions != itr_unique_constitutions_end;
        ++itr_unique_constitutions
      )
      {
        search_substructures.PushBack( graph_maker( **itr_unique_constitutions)); // Substructures to look for
        scaffold_counts.PushBack( ( *itr_unique_constitutions)->GetMDLPropertyAsVector( "ScaffoldCount")( 0));
      }
      util::ShPtrVector< FragmentConstitutionShared> constitution_vector( unique_constitutions.Begin(), unique_constitutions.End());

      BCL_MessageStd( "Detecting all root nodes and children");

      // children to all their root nodes
      storage::Vector< storage::Vector< size_t> > children_to_roots( n_uniq_constitutions);

      // nodes to their immediate, most important children
      storage::Vector< storage::Set< size_t> > child_nodes( n_uniq_constitutions);

      storage::Vector< size_t> root_nodes;
      storage::Vector< size_t> depth( n_uniq_constitutions, size_t( 0));
      storage::Vector< size_t> best_parent_root( n_uniq_constitutions, util::GetUndefined< size_t>());
      storage::Vector< size_t> best_parent( n_uniq_constitutions, util::GetUndefined< size_t>());

      // for each molecule, molecules at or above it in depth
      storage::Vector< storage::Vector< size_t> > required_graph( n_uniq_constitutions);
      // do the sub structure search
      {
        storage::Vector< size_t> needed_root( n_uniq_constitutions, size_t( 0));
        storage::List< size_t> root_nodes_list;
        size_t molecule_index( 0), n_roots( 0);
        for
        (
          storage::Vector< graph::ConstGraph< size_t, size_t> >::iterator
            itr_frags( search_substructures.Begin()), itr_frags_end( search_substructures.End());
          itr_frags != itr_frags_end;
          ++itr_frags, ++molecule_index
        )
        {
          // set the graph ownership to itr_frags
          graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
          isomorphism.SetGraphExternalOwnership( *itr_frags);

          // temporary list to store indices of child fragments
          storage::List< size_t> mol_to_roots;
          for
          (
            storage::List< size_t>::const_iterator itr_root_id( root_nodes_list.Begin()), itr_root_id_end( root_nodes_list.End());
            itr_root_id != itr_root_id_end;
            ++itr_root_id
          )
          {
            // set subgraph ownership to itr_frag_a
            isomorphism.SetSubgraphExternalOwnership( search_substructures( *itr_root_id));

            // if frag_a is substructure of frag, then store the index of frag_a in
            if( isomorphism.FindIsomorphism())
            {
              mol_to_roots.PushBack( *itr_root_id);
            }
          }
          if( mol_to_roots.IsEmpty())
          {
            root_nodes_list.PushBack( molecule_index);
            best_parent_root( molecule_index) = molecule_index;
            ++n_roots;
          }
          else
          {
            children_to_roots( molecule_index).InternalData().assign( mol_to_roots.Begin(), mol_to_roots.End());
            for
            (
              auto itr_m( children_to_roots( molecule_index).Begin()), itr_m_end( children_to_roots( molecule_index).End());
              itr_m != itr_m_end;
              ++itr_m
            )
            {
              needed_root( *itr_m) = size_t( 1);
            }
            // walk back through the list of molecules, find the first one that this molecule could be attached to
            size_t prev_mol_index( molecule_index);
            size_t previous_root( children_to_roots( molecule_index).LastElement());
            while( --prev_mol_index > previous_root)
            {
              if( depth( prev_mol_index) == size_t( 0))
              {
                // root that this node doesn't connect to
                continue;
              }
              if( IsSubsetOfForSortedVectors( children_to_roots( prev_mol_index), children_to_roots( molecule_index)))
              {
                isomorphism.SetSubgraphExternalOwnership( search_substructures( prev_mol_index));
                if( isomorphism.FindIsomorphism())
                {
                  break;
                }
              }
            }
            size_t final_best_prev_mol( prev_mol_index);
            needed_root( best_parent_root( prev_mol_index)) = 0;
            required_graph( molecule_index).PushBack( prev_mol_index);
            while( --prev_mol_index != util::GetUndefinedSize_t())
            {
               if
               (
                 needed_root( best_parent_root( prev_mol_index))
                 && depth( prev_mol_index) <= depth( molecule_index)
                 && IsSubsetOfForSortedVectors( children_to_roots( prev_mol_index), children_to_roots( molecule_index))
               )
               {
                 isomorphism.SetSubgraphExternalOwnership( search_substructures( prev_mol_index));
                 if( isomorphism.FindIsomorphism())
                 {
                   if( scaffold_counts( prev_mol_index) < scaffold_counts( final_best_prev_mol))
                   {
                     final_best_prev_mol = prev_mol_index;
                   }
                   required_graph( molecule_index).PushBack( prev_mol_index);
                   needed_root( best_parent_root( prev_mol_index)) = 0;
                 }
               }
            }
            best_parent( molecule_index) = final_best_prev_mol;
            child_nodes( final_best_prev_mol).Insert( molecule_index);
            depth( molecule_index) = depth( final_best_prev_mol) + 1;
            best_parent_root( molecule_index) = best_parent_root( final_best_prev_mol);
            required_graph( molecule_index).RemoveElements( required_graph( molecule_index).Find( final_best_prev_mol), size_t( 1));
            for
            (
              auto itr_m( children_to_roots( molecule_index).Begin()), itr_m_end( children_to_roots( molecule_index).End());
              itr_m != itr_m_end;
              ++itr_m
            )
            {
              BCL_Assert( !needed_root( *itr_m), "This shouldn't happen");
            }
            // for each root that this molecule has, find the predecessor to this node in search order.
          }
          util::GetLogger().LogStatus
          (
            "Looking for roots "
            + util::Format()( molecule_index)
            + " / "
            + util::Format()( n_uniq_constitutions)
            + " " + util::Format().W( 4)( double( molecule_index) * 100.0 / double( n_uniq_constitutions)) + " % "
            + util::Format()( n_roots) + " roots found"
          );
        }

        root_nodes.InternalData().assign( root_nodes_list.Begin(), root_nodes_list.End());
      }

      // Attach each molecule to the nearest prior molecule in the ensemble (ordered by scaffold counts) that it is part of

      // Record the minimal set of molecules that must have been matched earlier in the search in order for this molecule to
      // be matched. This set of molecules is composed of all matched terminii/leaf nodes on the search across the existing
      // trees

      // Only form the link parent->child for the lowest-scaffold-numbered parent of the molecule

      BCL_MessageStd( "Looking for all molecules within all other molecules");

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // Reorder constitution ensemble such that root nodes are first followed by child nodes.
      // The file will be arranged so that the tree is traversed in breadth-first search manner
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // convert to vector for random access

      util::ShPtrVector< FragmentConstitutionShared> reordered_list;
      reordered_list.AllocateMemory( constitution_vector.GetSize());

      storage::Vector< size_t> constitution_tracker, inv_constitution_tracker( constitution_vector.GetSize());
      constitution_tracker.AllocateMemory( constitution_vector.GetSize());

      storage::List< size_t> current_nodes( root_nodes.Begin(), root_nodes.End());

      std::string have_seen( constitution_vector.GetSize(), '0');
      while( !current_nodes.IsEmpty())
      {
        storage::List< size_t> new_cur_nodes;
        for
        (
            storage::List< size_t>::const_iterator
            itr_cur( current_nodes.Begin()), itr_cur_end( current_nodes.End());
            itr_cur != itr_cur_end;
            ++itr_cur
        )
        {
          if( have_seen[ *itr_cur] == '0')
          {
            have_seen[ *itr_cur] = '1';
            inv_constitution_tracker( *itr_cur) = reordered_list.GetSize();
            reordered_list.PushBack( constitution_vector( *itr_cur));
            constitution_tracker.PushBack( *itr_cur);
            const storage::Set< size_t> &nodes_for_next( child_nodes( *itr_cur));
            new_cur_nodes.InsertElements( new_cur_nodes.End(), nodes_for_next.Begin(), nodes_for_next.End());
          }
        }
        current_nodes = new_cur_nodes;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // Since the constitutions have been reordered, update the child nodes for the tree
      // Also reorder the constitution graphs accordingly
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      child_nodes.Reorder( constitution_tracker);
      best_parent.Reorder( constitution_tracker);
      size_t pos( 0);
      for
      (
        storage::Vector< storage::Set< size_t> >::iterator
          itr_child( child_nodes.Begin()), itr_child_end( child_nodes.End());
        itr_child != itr_child_end;
        ++itr_child, ++pos
      )
      {
        if( util::IsDefined( best_parent( pos)))
        {
          best_parent( pos) = inv_constitution_tracker( best_parent( pos));
        }
        storage::Set< size_t> &this_child_nodes( *itr_child);
        storage::Set< size_t> this_reordered;
        for
        (
          storage::Set< size_t>::const_iterator
            itr( this_child_nodes.Begin()), itr_end( this_child_nodes.End());
          itr != itr_end;
          ++itr
        )
        {
          this_reordered.Insert( inv_constitution_tracker( *itr));
        }
        this_child_nodes = this_reordered;
      }

      search_substructures.Reorder( constitution_tracker);

      BCL_MessageStd( "Mapping conformations to constitutions");

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // Get mapping of conformations to reordered constitutions
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a vector to store indices of child fragments
      storage::Vector< storage::Set< size_t> > exact_indices_vector( constitution_vector.GetSize());

      // do the sub structure search
      size_t conformation_index( 0);
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_conformations( CONFORMATIONS.Begin()), itr_conformations_end( CONFORMATIONS.End());
        itr_conformations != itr_conformations_end;
        ++itr_conformations, ++conformation_index
      )
      {
        // create a graph for this molecule
        graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( FragmentConstitutionShared( *itr_conformations)));

        // set the graph ownership to itr_frags
        graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
        isomorphism.SetSubgraphExternalOwnership( mol_graph);

        size_t constitution_index( 0);
        for
        (
          storage::Vector< graph::ConstGraph< size_t, size_t> >::const_iterator
            itr_frags( search_substructures.Begin()), itr_frags_end( search_substructures.End());
          itr_frags != itr_frags_end;
          ++itr_frags, ++constitution_index
        )
        {
          if( mol_graph.GetSize() != itr_frags->GetSize() || mol_graph.NumEdges() != itr_frags->NumEdges())
          {
            continue;
          }

          // set subgraph ownership to itr_frag_a
          isomorphism.SetGraphExternalOwnership( *itr_frags);
          // if frag_a is substructure of frag, then store the index of frag_a in
          if( isomorphism.FindIsomorphism())
          {
            exact_indices_vector( constitution_index).Insert( conformation_index);
          }
        }
      }

      BCL_MessageStd( "Writing data");

      CreateImpl( root_nodes.GetSize(), child_nodes, reordered_list, exact_indices_vector, CONFORMATIONS, best_parent);
    }

    //! @brief get constitutions that are root in the substructure tree i.e. they are not contained in any other fragments
    graph::ConstGraph< size_t, size_t> RotamerLibraryInterface::CreateConstitutionGraphFromString
    (
      const storage::Vector< std::string> &GRAPH_VERTEX_EDGE_DATA,
      const storage::Vector< std::string> &MOLECULE_INFO
    ) const
    {
      size_t number_atoms( util::ConvertStringToNumericalValue< size_t>( GRAPH_VERTEX_EDGE_DATA( 0)));
      size_t number_bonds( util::ConvertStringToNumericalValue< size_t>( GRAPH_VERTEX_EDGE_DATA( 1)));

      // get the adjacency list from the second line.
      storage::Vector< graph::UndirectedEdge< size_t> > adjacency_list;
      adjacency_list.AllocateMemory( number_bonds);

      // get vertices from the second line
      storage::Vector< size_t> vertices( number_atoms);

      // keep track of which atoms have been seen and inserted in the vertices vector
      storage::Set< size_t> atoms_seen;

      for( size_t entry_index( 0); entry_index < number_bonds; entry_index++)
      {
        // each edge is represented as set of five integers corresponding to first_atom_index, first_atom_data,
        // second_atom_index, second_atom_data, bond_data
        size_t string_position( entry_index * 5);

        size_t first_atom_index( util::ConvertStringToNumericalValue< size_t>( MOLECULE_INFO( string_position + 0)));
        size_t first_atom_type( util::ConvertStringToNumericalValue< size_t>( MOLECULE_INFO( string_position + 1)));
        size_t second_atom_index( util::ConvertStringToNumericalValue< size_t>( MOLECULE_INFO( string_position + 2)));
        size_t second_atom_type( util::ConvertStringToNumericalValue< size_t>( MOLECULE_INFO( string_position + 3)));
        size_t bond_type( util::ConvertStringToNumericalValue< size_t>( MOLECULE_INFO( string_position + 4)));

        // get edge from the current bond
        graph::UndirectedEdge< size_t> new_edge
        (
          first_atom_index,
          second_atom_index,
          bond_type
        );

        adjacency_list.PushBack( new_edge);

        // check if first_atom has already been seen, if not then insert into vertices vector.
        if( !atoms_seen.Contains( first_atom_index))
        {
          vertices( first_atom_index) = first_atom_type;
          atoms_seen.Insert( first_atom_index);
        }

        if( !atoms_seen.Contains( second_atom_index))
        {
          vertices( second_atom_index) = second_atom_type;
          atoms_seen.Insert( second_atom_index);
        }
      }

      // create graph from the vertices and adjacency list
      graph::ConstGraph< size_t, size_t> constitution_graph
      (
        vertices,
        adjacency_list,
        GetConstitutionalBondTypes().e_Undefined
      );

      // end
      return constitution_graph;
    }

  } // namespace chemistry
} // namespace bcl
