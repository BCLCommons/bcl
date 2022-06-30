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
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_ofstream.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_molecule_fragment_recombination.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"
using bcl::io::DirectoryEntry;

// external includes

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MoleculeFragmentRecombination::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new MoleculeFragmentRecombination)
    );

    //! @brief default constructor
    MoleculeFragmentRecombination::MoleculeFragmentRecombination
    (
      const std::string &FILENAME,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON,
      const ConfigurationalBondTypeData::Data &BOND_COMPARISON,
      const bool &MUTUALLY_MATCHING_ATOMS
    ) :
          m_File( FILENAME),
          m_FileWasRead( false),
          m_AtomComparison( ATOM_COMPARISON),
          m_BondComparison( BOND_COMPARISON),
          m_Converter( ATOM_COMPARISON, BOND_COMPARISON, false),
          m_MutuallyMatchingAtoms( MUTUALLY_MATCHING_ATOMS)
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeFragmentRecombination::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &MoleculeFragmentRecombination::GetAlias() const
    {
      static const std::string s_name( "MoleculeFragmentRecombination");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &MoleculeFragmentRecombination::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the minimum size of fragments
    //! @return the minimum size of fragments
    const size_t MoleculeFragmentRecombination::GetMinSize() const
    {
      return 0;
    }

    //! @brief returns an ensemble of fragments of a molecule
    //! @param CONFORMATION molecule of interest
    //! @return an ensemble of common substructures relative to those in a file
    FragmentEnsemble MoleculeFragmentRecombination::operator()( const ConformationInterface &CONFORMATION) const
    {
      // read in molecules
      ReadFile();

      // initialize output ensemble
      FragmentEnsemble final_ensemble;

      // loop over comparator molecules
      for
      (
          FragmentEnsemble::iterator mol_itr( m_Molecules.Begin()), mol_itr_end( m_Molecules.End());
          mol_itr != mol_itr_end;
          ++mol_itr
      )
      {
        // recombine molecules
        FragmentEnsemble ensemble( this->CompareSubstructures( CONFORMATION, *mol_itr, m_MutuallyMatchingAtoms, false));

        // add to output ensemble only if it is not an empty fragment
        for
        (
            auto new_mol_itr( ensemble.Begin()), new_mol_itr_end( ensemble.End());
            new_mol_itr != new_mol_itr_end;
            ++new_mol_itr
        )
        {
          if( new_mol_itr->GetSize())
          {
            final_ensemble.PushBack( *new_mol_itr);
          }
        }
      }
      return final_ensemble;
    }

    //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
    //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
    //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
    //! @param mol_a mol_a molecule to be perturbed
    //! @param mol_b ensemble of molecules against which mol_a molecule will be compared
    //! @return indices for each atom in mol_a molecule and their corresponding score contributions
    FragmentEnsemble MoleculeFragmentRecombination::CompareSubstructures
    (
      const FragmentComplete &MOL_A,
      const FragmentComplete &MOL_B,
      const bool &MUTUALLY_MATCHING_ATOMS,
      const bool &OUTPUT_INTERMEDIATES
    ) const
    {
      // we will return an ensemble of new molecules
      FragmentEnsemble recombined_mols;

      // open a stream to output intermediate structures
      io::OFStream output_intermediate_orig, output_intermediate, output_intermediate_fragment, subgraph_a_sdf, subgraph_b_sdf;
      if( OUTPUT_INTERMEDIATES)
      {
        io::File::MustOpenOFStream( output_intermediate_orig, "FeatureMapper.input.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( subgraph_a_sdf, "Subgraph_A.sdf");
        io::File::MustOpenOFStream( subgraph_b_sdf, "Subgraph_B.sdf");
      }

      // make modifiable mol_a molecule
      FragmentComplete mol_a( MOL_A);
      FragmentComplete mol_b( MOL_B);

      // write intermediates
      if( OUTPUT_INTERMEDIATES)
      {
        mol_a.WriteMDL( output_intermediate_orig);
        mol_b.WriteMDL( output_intermediate_orig);
      }

      // prepare graph marker with ring splitter
      FragmentSplitRings split_rings( true);
      FragmentGraphMarker graph_marker( m_Converter, split_rings);

      // make a graph of mol_a
      graph::ConstGraph< size_t, size_t> mol_a_graph( m_Converter( mol_a));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_a_graph_ptr( &mol_a_graph, false);

      // generate a graph for mol_b
      graph::ConstGraph< size_t, size_t> mol_b_graph( m_Converter( mol_b));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_b_graph_ptr( &mol_b_graph, false);

      // set mol_a and mol_b to be the first and second graph in the subgraph iso, respectively
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso( graph::CommonSubgraphIsomorphismBase::e_Connected);
      common_subgraph_iso.SetGraphA( mol_a_graph_ptr);
      common_subgraph_iso.SetGraphB( mol_b_graph_ptr);

      // get the subgraph isomorphism (common subgraph between the two molecules) of the mol_a molecule
      if( MUTUALLY_MATCHING_ATOMS)
      {
        // identify matched atoms
        storage::Vector< size_t> matched_atoms_mol_a, matched_atoms_mol_b;
        ConformationComparisonPsiField::GetAlignedAtoms( mol_a, mol_b, matched_atoms_mol_a, matched_atoms_mol_b);

        // make a map from the matched atoms
        storage::Map< size_t, size_t> matching_atoms;
        BCL_Assert( matched_atoms_mol_a.GetSize() == matched_atoms_mol_b.GetSize(), "Unequal number of mutually matching atoms!");
        for
        (
            size_t i( 0), i_sz( matched_atoms_mol_a.GetSize()); i < i_sz; ++i
        )
        {
          matching_atoms.Insert( std::make_pair( matched_atoms_mol_a( i), matched_atoms_mol_b( i)));
        }

        // set the isomorphism from the matched atoms map
        common_subgraph_iso.SetIsomorphism( matching_atoms);
      }
      else
      {
        common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), size_t( 1));
      }

      // get subgraphs
      graph::Subgraph< size_t, size_t> subgraph_mol_a
      (
        common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement()
      );
      graph::Subgraph< size_t, size_t> subgraph_mol_b
      (
        common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB().FirstElement()
      );

      // debug
      if( OUTPUT_INTERMEDIATES)
      {
        AtomVector< AtomComplete> subgraph_a_atom_vec( MOL_A.GetAtomVector());
        AtomVector< AtomComplete> subgraph_b_atom_vec( MOL_B.GetAtomVector());
        subgraph_a_atom_vec.Reorder( subgraph_mol_a.GetVertexIndices());
        subgraph_b_atom_vec.Reorder( subgraph_mol_b.GetVertexIndices());
        FragmentComplete subgraph_a_frag( subgraph_a_atom_vec, "subgraph_a");
        FragmentComplete subgraph_b_frag( subgraph_b_atom_vec, "subgraph_b");
        subgraph_a_frag.WriteMDL(subgraph_a_sdf);
        subgraph_b_frag.WriteMDL(subgraph_b_sdf);
      }

      // create graphs lacking all intra-isomorphism edges
      graph::ConstGraph< size_t, size_t> mol_a_graph_depleted( mol_a_graph);
      graph::ConstGraph< size_t, size_t> mol_b_graph_depleted( mol_b_graph);
      for( size_t atom_i( 0), iso_sz( subgraph_mol_b.GetSize()); atom_i < iso_sz; ++atom_i)
      {
        mol_a_graph_depleted.RemoveAllEdges( subgraph_mol_a.GetVertexIndices()( atom_i));
        mol_b_graph_depleted.RemoveAllEdges( subgraph_mol_b.GetVertexIndices()( atom_i));
      }

      // get the largest isomorphism mapping atoms from mol_a to mol_b
      storage::Map< size_t, size_t> largest_iso( common_subgraph_iso.GetIsomorphism());
      //BCL_Debug( largest_iso);
      //BCL_Debug( largest_iso.GetKeysAsVector()( 18));
      //BCL_Debug( largest_iso.GetMappedValues()( 18));

      // get edges connecting subgraph to supergraph
      storage::Vector< storage::Vector< size_t> > subgraph_mol_a_edges( subgraph_mol_a.GetOrderedAdjacentEdgeIndices());
      storage::Vector< storage::Vector< size_t> > subgraph_mol_b_edges( subgraph_mol_b.GetOrderedAdjacentEdgeIndices());

      // iterate over the atoms connecting the subgraph to the supergraph (indexed by subgraph indices)
      for( size_t atom_i( 0), iso_sz( subgraph_mol_a_edges.GetSize()); atom_i < iso_sz; ++atom_i)
      {
        // if this atom from the subgraph does not connect with the supergraph then skip it
        if( subgraph_mol_a_edges( atom_i).IsEmpty())
        {
          continue;
        }
        //BCL_Debug( atom_i);
        // for each atom in each subgraph, grab all its adjacent edges and their connected atoms using
        storage::Vector< size_t> mol_a_attribute_score;
        for
        (
            auto subgraph_mol_a_itr( subgraph_mol_a_edges( atom_i).Begin()),
            subgraph_mol_a_itr_end( subgraph_mol_a_edges( atom_i).End());
            subgraph_mol_a_itr != subgraph_mol_a_itr_end;
            ++subgraph_mol_a_itr
        )
        {
          if( mol_a_attribute_score.Find( *subgraph_mol_a_itr) < mol_a_attribute_score.GetSize())
          {
            // cycle; already grabbed this vertex, continue
            continue;
          }

          // grab the inverted vertices
          util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_remove
          (
            graph::Connectivity::GetVerticesReachableFrom
            (
              mol_a_graph_depleted,
              *subgraph_mol_a_itr
            )
          );

          mol_a_attribute_score.Append( *reachable_vertices_to_remove);
        }
        //BCL_Debug( mol_a_attribute_score);

        // atoms for the non-subgraph fragment
        AtomVector< AtomComplete> mol_a_attribute( mol_a.GetAtomVector());
        mol_a_attribute.Reorder( mol_a_attribute_score);
        FragmentComplete mol_a_attribute_fc( mol_a_attribute, "");

        // get remainder of mol_a without the selected non-subgraph fragment
        const storage::Vector< size_t> mol_a_reachable
        (
          graph::Subgraph< size_t, size_t>( mol_a_graph_ptr, mol_a_attribute_score).GetComplement().GetVertexIndices()
        );
        //BCL_Debug( mol_a_reachable);

        // build fragment for remainder of mol_a
        AtomVector< AtomComplete> mol_a_keep( mol_a.GetAtomVector());
        mol_a_keep.Reorder( mol_a_reachable);
        FragmentComplete mol_a_keep_fc( mol_a_keep, "");

        // repeat the procedure for mol_b
        storage::Vector< size_t> mol_b_attribute_score;
        for
        (
            auto subgraph_mol_b_itr( subgraph_mol_b_edges( atom_i).Begin()),
            subgraph_mol_b_itr_end( subgraph_mol_b_edges( atom_i).End());
            subgraph_mol_b_itr != subgraph_mol_b_itr_end;
            ++subgraph_mol_b_itr
        )
        {
          if( mol_b_attribute_score.Find( *subgraph_mol_b_itr) < mol_b_attribute_score.GetSize())
          {
            // cycle; already grabbed this vertex, continue
            continue;
          }
          // grab the inverted vertices
          util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_remove
          (
            graph::Connectivity::GetVerticesReachableFrom
            (
              mol_b_graph_depleted,
              *subgraph_mol_b_itr
            )
          );
          mol_b_attribute_score.Append( *reachable_vertices_to_remove);
        }
        //BCL_Debug( mol_b_attribute_score);

        // atoms for the non-subgraph fragment
        AtomVector< AtomComplete> mol_b_attribute( mol_b.GetAtomVector());
        mol_b_attribute.Reorder( mol_b_attribute_score);
        FragmentComplete mol_b_attribute_fc( mol_b_attribute, "");

        // obtain remainder of mol_b
        const storage::Vector< size_t> mol_b_reachable
        (
          graph::Subgraph< size_t, size_t>( mol_b_graph_ptr, mol_b_attribute_score).GetComplement().GetVertexIndices()
        );
        //BCL_Debug( mol_b_reachable);

        // make fragment for remainder of mol_b
        AtomVector< AtomComplete> mol_b_keep( mol_b.GetAtomVector());
        mol_b_keep.Reorder( mol_b_reachable);
        FragmentComplete mol_b_keep_fc( mol_b_keep, "");

        // write intermediates
        if( OUTPUT_INTERMEDIATES)
        {
          mol_a_keep_fc.WriteMDL( output_intermediate);
          mol_a_attribute_fc.WriteMDL( output_intermediate_fragment);
          mol_b_keep_fc.WriteMDL( output_intermediate);
          mol_b_attribute_fc.WriteMDL( output_intermediate_fragment);
        }

        // recombine the pieces
        recombined_mols.PushBack
        (
          FragmentComplete
          (
            this->RecombineComponents
            (
              MOL_A,
              MOL_B,
              mol_a_keep_fc,
              mol_b_keep_fc,
              largest_iso,
              mol_a_attribute_fc,
              mol_a_attribute_score,
              mol_b_attribute_fc,
              mol_b_attribute_score
            )
          )
        );

      } // end looping over the non-subgraph components

      // close these output streams
      io::File::CloseClearFStream( output_intermediate_orig);
      io::File::CloseClearFStream( output_intermediate);
      io::File::CloseClearFStream( output_intermediate_fragment);
      io::File::CloseClearFStream( subgraph_a_sdf);
      io::File::CloseClearFStream( subgraph_b_sdf);

      // be happy, be done
      return recombined_mols;
    }

    FragmentComplete MoleculeFragmentRecombination::RecombineComponents
    (
      const FragmentComplete &PARENT_A,
      const FragmentComplete &PARENT_B,
      const FragmentComplete &BASE_MOL_A,
      const FragmentComplete &BASE_MOL_B,
      const storage::Map< size_t, size_t> &ISO,
      const FragmentComplete &FRAG_A,
      const storage::Vector< size_t> &FRAG_A_COMPONENT_INDICES,
      const FragmentComplete &FRAG_B,
      const storage::Vector< size_t> &FRAG_B_COMPONENT_INDICES
    ) const
    {
      // start with BASE_MOL_A, get its atoms and bonds
      AtomVector< AtomComplete> base_mol_a_vec( BASE_MOL_A.GetAtomVector());
      storage::Vector< sdf::AtomInfo> base_mol_a_atominfo( base_mol_a_vec.GetAtomInfo());
      storage::Vector< sdf::BondInfo> base_mol_a_bondinfo( base_mol_a_vec.GetBondInfo());

      // get the FRAG_B atoms and bonds
      AtomVector< AtomComplete> frag_b_vec( FRAG_B.GetAtomVector());
      storage::Vector< sdf::AtomInfo> frag_b_atominfo( frag_b_vec.GetAtomInfo());
      storage::Vector< sdf::BondInfo> frag_b_bondinfo( frag_b_vec.GetBondInfo());

      // add FRAG_B atoms to BASE_MOL_A
      // loop over atoms in FRAG_B
      for( size_t i( 0); i < frag_b_atominfo.GetSize(); ++i)
      {
        // add the current atom in FRAG_B to the atoms in BASE_MOL_A
        base_mol_a_atominfo.PushBack( frag_b_atominfo( i));
      }

      // add the FRAG_B internal connections to BASE_MOL_A
      for( size_t i( 0); i < frag_b_bondinfo.GetSize(); ++i)
      {
        //BCL_Debug( frag_b_bondinfo( i).GetAtomIndexLow());
        //BCL_Debug( frag_b_bondinfo( i).GetAtomIndexHigh());
        //BCL_Debug( frag_b_bondinfo( i).GetAtomIndexLow() + BASE_MOL_A.GetSize());
        //BCL_Debug( frag_b_bondinfo( i).GetAtomIndexHigh() + BASE_MOL_A.GetSize());
        //BCL_Debug( frag_b_bondinfo( i).GetConfigurationalBondType());
        base_mol_a_bondinfo.PushBack
        (
          sdf::BondInfo
          (
            frag_b_bondinfo( i).GetAtomIndexLow() + BASE_MOL_A.GetSize(),
            frag_b_bondinfo( i).GetAtomIndexHigh() + BASE_MOL_A.GetSize(),
            frag_b_bondinfo( i).GetConfigurationalBondType()
          )
        );
      }

      // add the connections between FRAG_B and BASE_MOL_A
      storage::Set< size_t> conf_moveable_indices;
      for( size_t i( 0); i < frag_b_atominfo.GetSize(); ++i)
      {
        // update the index of the new atom in BASE_MOL_A
        // atoms are always added to the end of the atom vector
        size_t new_atom_index( BASE_MOL_A.GetSize() + i);
        conf_moveable_indices.InsertElement( new_atom_index);
        //BCL_Debug( new_atom_index);

        //BCL_Debug( frag_b_atominfo( i).GetAtomType());
        //BCL_Debug( PARENT_B.GetAtomVector()( FRAG_B_COMPONENT_INDICES( i)).GetAtomType());

        // get the original connectivity info for FRAG_B in its parent molecule
        storage::Vector< sdf::BondInfo> frag_b_original_bondinfo;
        storage::Vector< storage::Pair< size_t, ConfigurationalBondType>> atom_bond_pairs;
        for
        (
            auto bond_itr( PARENT_B.GetAtomVector()( FRAG_B_COMPONENT_INDICES( i)).GetBonds().Begin()),
            bond_itr_end( PARENT_B.GetAtomVector()( FRAG_B_COMPONENT_INDICES( i)).GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          // get the atom index of the connected atom in PARENT_B as well as the bondtype
          atom_bond_pairs.PushBack
          (
            std::make_pair
            (
              PARENT_B.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()),
              bond_itr->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic ?
                  GetConfigurationalBondTypes().e_AromaticBond :
                  bond_itr->GetBondType()
            )
          );
        }
        //BCL_Debug( atom_bond_pairs);

        // now we need to find the common subgraph iso atom in the BASE_MOL_A for the parent atom index from PARENT_B
        for
        (
            auto itr( atom_bond_pairs.Begin()), itr_end( atom_bond_pairs.End());
            itr != itr_end;
            ++itr
        )
        {
          // one strategy is to get the MappedValues, which are the PARENT_B indices, find out the index of atom in that vector
          // which corresponds to the position in the ISO, and then access that position/index in the GetKeys vector so that
          // I can get the corresponding atom in PARENT_A. then based on the atom index(ices) that was/were removed, I can
          // compute the adjusted index of the BASE_MOL_A that I need to connect to.
          // the same thing can be achieved with easier logic by just iterating over the map
          for
          (
              auto iso_itr( ISO.Begin()), iso_itr_end( ISO.End());
              iso_itr != iso_itr_end;
              ++iso_itr
          )
          {
            // if the map value is our PARENT_B atom index
            //BCL_Debug( iso_itr->second);
            if( iso_itr->second == itr->First())
            {
              // get the corresponding atom in PARENT_A
              size_t atom_index( iso_itr->first);
              //BCL_Debug( atom_index);

              // adjust the atom index count if atom indices removed from PARENT_A to make BASE_MOL_A were smaller than atom index
              for( size_t frag_a_index( 0); frag_a_index < FRAG_A_COMPONENT_INDICES.GetSize(); ++frag_a_index)
              {
                if( FRAG_A_COMPONENT_INDICES( frag_a_index) < iso_itr->first)
                {
                  atom_index -= size_t( 1);
                }
              }

              // create the bond connecting the BASE_MOL_A atom_index atom to the added FRAG_B new_atom_index
              base_mol_a_bondinfo.PushBack
              (
                sdf::BondInfo
                (
                  atom_index,
                  new_atom_index,
                  itr->Second()
                )
              );
            }
          }
          // TODO:
          // if I had done this the other way around, I could just ISO.Find(KEY)->second to get the value, then do my subtraction
          // maybe I will need to do both anyway if I want the molecules in both sets of coordinates. if not, maybe will just do this
        }
      }
        // make our new molecule
        AtomVector< AtomComplete> new_mol_vec( base_mol_a_atominfo, base_mol_a_bondinfo);
        FragmentComplete new_mol_unfixed( new_mol_vec, "");

        // debug
//        io::OFStream debug_out;
//        io::File::MustOpenOFStream( debug_out, "NEW_unfixed.sdf", std::ios::app);
//        new_mol_unfixed.WriteMDL(debug_out);
//        io::File::CloseClearFStream( debug_out);

        // let's add the adjacent atom indices to help us fix the conformer
        storage::Set< size_t> adjacent_indices;
        storage::Vector< size_t> initial_indices( conf_moveable_indices.Begin(), conf_moveable_indices.End());
        for( size_t a( 0); a < initial_indices.GetSize(); ++a)
        {
          for
          (
              auto bond_itr( new_mol_vec( initial_indices( a)).GetBonds().Begin()),
              bond_itr_end( new_mol_vec( initial_indices( a)).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            adjacent_indices.InsertElement( new_mol_vec.GetAtomIndex( bond_itr->GetTargetAtom()));
          }
        }
        conf_moveable_indices.InsertElements( adjacent_indices);

        // create graph for the new molecule
        ConformationGraphConverter graph_maker
        (
          ConformationGraphConverter::e_AtomType,
          ConfigurationalBondTypeData::e_IsInRing
        );
        ConformationGraphConverter::t_AtomGraph molecule_graph( ConformationGraphConverter::CreateGraphWithAtoms( new_mol_unfixed));

        // get ring components
        FragmentSplitRings splitter( true, size_t( 3));
        storage::List< storage::Vector< size_t> > ring_components( splitter.GetComponentVertices( new_mol_unfixed, molecule_graph));
        storage::Vector< storage::Vector< size_t>> ring_components_vec( ring_components.Begin(), ring_components.End());

        // decide which rings matter
        storage::Vector< size_t> rings_that_matter;
        size_t ring_index( 0);
        for
        (
            auto ring_comp_itr( ring_components.Begin()), ring_comp_itr_end( ring_components.End());
            ring_comp_itr != ring_comp_itr_end;
            ++ring_comp_itr, ++ring_index
        )
        {
          for
          (
              auto atom_itr( ring_comp_itr->Begin()), atom_itr_end( ring_comp_itr->End());
              atom_itr != atom_itr_end;
              ++atom_itr
          )
          {
            // check against all of my perturbed atoms
            for( size_t i( 0); i < initial_indices.GetSize(); ++i)
            {
              if( *atom_itr == initial_indices( i))
              {
                rings_that_matter.PushBack( ring_index);
                break;
              }
            }
          }
        }

        // add the atoms from the important rings to the moveable indices
        for
        (
            auto important_rings_itr( rings_that_matter.Begin()),
            import_rings_itr_end( rings_that_matter.End());
            important_rings_itr != import_rings_itr_end;
            ++important_rings_itr
        )
        {
          for( size_t atom( 0); atom < ring_components_vec( *important_rings_itr).GetSize(); ++atom)
          {
            conf_moveable_indices.InsertElement( ring_components_vec( *important_rings_itr)( atom));
          }
        }

        // clean the new molecule
        FragmentMapConformer cleaner( "None", false, storage::Vector< size_t>( conf_moveable_indices.Begin(), conf_moveable_indices.End()));
        util::ShPtr< FragmentComplete> new_mol( cleaner.Clean( new_mol_vec, BASE_MOL_A, "None", false));
        if( new_mol.IsDefined())
        {
          return *new_mol;
        }

        // debug
//        io::File::MustOpenOFStream( debug_out, "NEW_fixed.sdf", std::ios::app);
//        new_mol->WriteMDL(debug_out);
//        io::File::CloseClearFStream( debug_out);

      return FragmentComplete();
    }

    //! @brief reads in molecules from a given file if it is necessary
    void MoleculeFragmentRecombination::ReadFile() const
    {
      if( !m_FileWasRead)
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_File);
        m_Molecules.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
        m_FileWasRead = true;
      }
    }

    //! @brief Set the members with LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeFragmentRecombination::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // initialize our graph converter
      m_Converter = ConformationGraphConverter( m_AtomComparison, m_BondComparison);

      // if it exists, we need to delete the combined intermediate file since it appends when it writes
      const std::string intermediate_input_filename( "FeatureMapper.input.intermediates.combined.sdf");
      io::DirectoryEntry intermediate_input_entry( intermediate_input_filename);
      if( intermediate_input_entry.DoesExist())
      {
        intermediate_input_entry.Remove();
      }

      const std::string intermediate_filename( "FeatureMapper.intermediates.combined.sdf");
      io::DirectoryEntry intermediate_entry( intermediate_filename);
      if( intermediate_entry.DoesExist())
      {
        intermediate_entry.Remove();
      }

      // same thing for the intermediate fragments
      const std::string intermediate_fragment_filename( "FeatureMapper.intermediate_frags.combined.sdf");
      io::DirectoryEntry intermediate_frag_entry( intermediate_fragment_filename);
      if( intermediate_frag_entry.DoesExist())
      {
        intermediate_frag_entry.Remove();
      }

      // no need to do this for the per-molecule files since those overwrite
      return true;
    }

    io::Serializer MoleculeFragmentRecombination::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Recombine molecules based on maximum common substructure differences");
      member_data.AddInitializer
      (
        "file",
        "file containing molecules to compare for largest common substructures",
        io::Serialization::GetAgentInputFilename( &m_File)
      );
      member_data.AddInitializer
      (
        "atom comparison",
        "atom data that is compared to determine whether atoms are equivalent",
        io::Serialization::GetAgent( &m_AtomComparison),
        "ElementType"
      );
      member_data.AddInitializer
      (
        "bond comparison",
        "bond data that is compared",
        io::Serialization::GetAgent( &m_BondComparison),
        "BondOrderAmideOrAromaticWithRingness"
      );
      return member_data;
    }

  } // namespace chemistry
} // namespace bcl
