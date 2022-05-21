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
#include <chemistry/bcl_chemistry_configuration_set.h>
#include <chemistry/bcl_chemistry_conformation_comparison_psi_field.h>
#include <chemistry/bcl_chemistry_conformation_graph_converter.h>
#include <chemistry/bcl_chemistry_constitution_set.h>
#include <chemistry/bcl_chemistry_fragment_configuration_shared.h>
#include <chemistry/bcl_chemistry_fragment_constitution_shared.h>
#include <chemistry/bcl_chemistry_fragment_graph_marker.h>
#include <chemistry/bcl_chemistry_fragment_split_interface.h>
#include <chemistry/bcl_chemistry_fragment_split_largest_component.h>
#include <chemistry/bcl_chemistry_fragment_split_rings.h>
#include <graph/bcl_graph_connectivity.h>
#include <graph/bcl_graph_subgraph.h>
#include <io/bcl_io_directory_entry.h>
#include <io/bcl_io_ofstream.h>
BCL_StaticInitializationFiascoFinder

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
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

    //! @brief default constructor
    MoleculeFeatureMapper::MoleculeFeatureMapper() :
          m_Descriptor(),
          m_AtomComparison( ConformationGraphConverter::e_ElementType),
          m_BondComparison( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
    {
    }

    //! @brief constructor
    //! @param SCORE_LABEL the score descriptor to use
    MoleculeFeatureMapper::MoleculeFeatureMapper
    (
      const util::ObjectDataLabel &SCORE_LABEL
    ) :
      m_Descriptor( SCORE_LABEL)
    {
    }

  ////////////////
  // operations //
  ////////////////

    namespace
    {
      //! @brief used for storing information about a molecule's feature map
      struct FeatureMapInfo
      {
        util::SiPtr< const FragmentComplete> m_Molecule;
        graph::ConstGraph< size_t, size_t> m_Graph;
        storage::Vector< size_t> m_ScaffoldIso;
        float m_Score;
        std::vector< std::map< size_t, float> > m_Perturbations;
      };

      //! @brief used to run a thread
      struct Worker
      {
        util::SiPtr< const FragmentComplete> m_Molecule; //! the molecule of interest
        std::map< size_t, float> m_LocalMap;             //! a local mapping of atom effects
        descriptor::CheminfoProperty m_Model;            //! the scoring function
        storage::Vector< size_t> m_PerturbAtoms;         //! the atoms to perturb (remove)
        size_t m_StartIndex;                             //! the starting index
        size_t m_Stride;                                 //! how many indices to skip each time
        bool m_IgnoreH;                                  //! whether to ignore hydrogens
        bool m_SplitLargest;                             //! whether to clean up accidental fragments
        FragmentSplitLargestComponent m_Splitter;        //! removes accidental fragments when removing atoms
        FragmentEnsemble m_Ensemble;

        Worker() :
          m_Molecule(),
          m_LocalMap(),
          m_Model(),
          m_PerturbAtoms(),
          m_StartIndex( 0),
          m_Stride( 0),
          m_SplitLargest( false),
          m_Splitter(),
          m_IgnoreH( false)
        {
        }

        void RunThread()
        {
          // If the molecule or stride hasn't been specified then do nothing
          if( !m_Molecule.IsDefined() || !m_Stride)
          {
            return;
          }

          FragmentComplete temp_parent( *m_Molecule);
          temp_parent.SaturateWithH();
          float parent_score( m_Model->SumOverObject( temp_parent)( 0));
          BCL_MessageDbg( "Parent molecule score: " + util::Format()( parent_score));

          const AtomVector< AtomComplete> &atom_vector( temp_parent.GetAtomVector());
          size_t orig_n_atoms( atom_vector.GetSize());

          //for( size_t atom_index( 0); atom_index < orig_n_atoms; ++atom_index)
          for( size_t index( m_StartIndex), max_index( m_PerturbAtoms.GetSize()); index < max_index; index += m_Stride)
          {
            size_t atom_index( m_PerturbAtoms( index));
            if( m_IgnoreH && atom_vector( atom_index).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
            {
              continue;
            }

            storage::Vector< size_t> keep_atoms;
            keep_atoms.AllocateMemory( orig_n_atoms - 1);
            for( size_t i( 0); i < orig_n_atoms; ++i)
            {
              if( i != atom_index)
              {
                keep_atoms.PushBack( i);
              }
            }

            // make new molecule with atom removed
            AtomVector< AtomComplete> trimmed_vector( atom_vector);
            trimmed_vector.Reorder( keep_atoms);
            FragmentComplete new_mol( trimmed_vector, "");

            // clean molecule to remove small fragments caused by atom deletion
            FragmentComplete new_clean_mol;
            if( m_SplitLargest)
            {
              FragmentEnsemble largest_component( m_Splitter( new_mol));
              new_clean_mol = largest_component.GetMolecules().FirstElement();
            }
            else
            {
              new_clean_mol = new_mol;
            }

            // finalize
            new_clean_mol.SaturateWithH();
            m_Ensemble.PushBack( new_clean_mol);
            float score( m_Model->SumOverObject( new_clean_mol)( 0));
            BCL_MessageDbg( "Perturbed molecule score: " + util::Format()( score));
            //            float score_diff( score - parent_score);
            float score_diff( parent_score - score);
            BCL_MessageDbg( "deltaScore: " + util::Format()( score_diff));
            m_LocalMap[ atom_index] = score_diff;
          }
        }
      };
    }

    //! @return a vector containing score changes when translating (first) and removing (second) atoms in each molecule
    std::vector< std::map< size_t, float> > MoleculeFeatureMapper::Perturb
    (
      const FragmentComplete &MOLECULE,
      const bool &IGNORE_H,
      const bool &SPLIT_LARGEST,
      storage::Vector< size_t> ATOM_INDICES
    ) const
    {
      BCL_Assert( m_Descriptor.IsDefined(), "MoleculeFeatureMapper::Perturb called before descriptor was set");
      const size_t &orig_n_atoms( MOLECULE.GetNumberAtoms());

      // if the indices vector is empty, operate on all atoms
      if( ATOM_INDICES.IsEmpty())
      {
        ATOM_INDICES.AllocateMemory( orig_n_atoms);
        for( size_t i( 0); i < orig_n_atoms; ++i)
        {
          ATOM_INDICES.PushBack( i);
        }
      }

      // set up threads
      size_t n_threads( sched::GetNumberCPUs());
      std::vector< Worker> threads( n_threads);
      for( size_t i( 0); i < n_threads; ++i)
      {
        threads[ i].m_Molecule = util::SiPtr< const FragmentComplete>( &MOLECULE);
        threads[ i].m_IgnoreH = IGNORE_H;
        threads[ i].m_PerturbAtoms = ATOM_INDICES;
        threads[ i].m_Model = m_Descriptor;
        threads[ i].m_StartIndex = i;
        threads[ i].m_Stride = n_threads;
        threads[ i].m_SplitLargest = SPLIT_LARGEST;
      }

      util::ShPtrVector< sched::JobInterface> jobs;
      jobs.AllocateMemory( n_threads);
      const size_t group( 1);

      for( size_t p( 0); p < n_threads; ++p)
      {
        jobs.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::ThunkJob< Worker, void>
            (
              group,
              threads[ p],
              &Worker::RunThread,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );

        // submit the job
        sched::GetScheduler().RunJob( jobs.LastElement());
      }

      // wait for jobs to finish
      for( size_t p( 0); p < n_threads; ++p)
      {
        sched::GetScheduler().Join( jobs( p));
      }

      std::vector< std::map< size_t, float> > return_vals( 1);
      for( size_t t( 0); t < n_threads; ++t)
      {
        return_vals[ 0].insert( threads[ t].m_LocalMap.begin(), threads[ t].m_LocalMap.end());
      }
      return return_vals;
    }

    //! @brief perturb a chemical structure by splitting out fragments and measure the score changes it makes to a model
    //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
    //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
    //! @param MOLECULE parent molecule to be perturbed
    //! @param SPLITTER split interface determining how parent molecule will be perturbed
    //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
    //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
    //! @return indices for each atom in parent molecule and their corresponding score contributions
    std::vector< std::map< size_t, float> > MoleculeFeatureMapper::PerturbByFragment
    (
      const FragmentComplete &MOLECULE,
      const util::Implementation< FragmentSplitInterface> &SPLITTER,
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const bool &MUTUALLY_MATCHING_ATOMS,
      const bool &NORMALIZE,
      const bool &IGNORE_H,
      const bool &AVERAGE,
      const size_t &MOL_INDEX,
      const bool &OUTPUT_INTERMEDIATES
    ) const
    {
      // require that our scorer is defined
      BCL_Assert( m_Descriptor.IsDefined(), "MoleculeFeatureMapper::PerturbByFragment called before descriptor was set");

      // delete these files before appending
      io::OFStream output_intermediate, output_intermediate_fragment;
      if( OUTPUT_INTERMEDIATES)
      {
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags.combined.sdf", std::ios::app);
      }
      io::File::CloseClearFStream( output_intermediate);
      io::File::CloseClearFStream( output_intermediate_fragment);

      // make modifiable parent molecule
      FragmentComplete molecule( MOLECULE);

      // score parent molecule
      FragmentComplete parent_temp( molecule);
      parent_temp.SaturateWithH();
      float parent_score( m_Descriptor->SumOverObject( parent_temp)( 0));
      BCL_MessageDbg( "Parent molecule score: " + util::Format()( parent_score));

      // initialize each atom contribution to zero
      const size_t &orig_n_atoms( molecule.GetNumberAtoms());
      std::map< size_t, float> atom_contributions;
      for( size_t i( 0); i < orig_n_atoms; ++i)
      {
        atom_contributions.insert( std::make_pair( i, float( 0.0)));
      }

      // split fragments but return inverted components
      ConformationGraphConverter::t_AtomGraph atom_graph
      (
        ConformationGraphConverter::CreateGraphWithAtoms( molecule)
      );
      ConformationGraphConverter::t_AtomGraph atom_graph_copy( atom_graph);
      auto split_component_vertices( SPLITTER->GetComponentVertices( molecule, atom_graph));
      FragmentEnsemble fragments( SPLITTER->ConvertComponentsIntoEnsemble( molecule, split_component_vertices, atom_graph_copy, true));

      if( !fragments.GetSize())
      {
        std::vector< std::map< size_t, float> > return_vals( 1);
        return_vals[ 0].insert( atom_contributions.begin(), atom_contributions.end());
        return return_vals;
      }

      // at this point, we now just have to do CompareSubstructuresRigorous with the inverted components as the references
      return CompareSubstructuresRigorous
          (
            molecule,
            fragments,
            ATOM_TYPE,
            BOND_TYPE,
            MUTUALLY_MATCHING_ATOMS,
            NORMALIZE,
            AVERAGE,
            util::GetUndefinedSize_t(),
            OUTPUT_INTERMEDIATES
          );
    }

    //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
    //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
    //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
    //! @param MOLECULE parent molecule to be perturbed
    //! @param ENSEMBLE ensemble of molecules against which parent molecule will be compared
    //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
    //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
    //! @param AVERAGE if true, compute average value per atom across all comparisons; otherwise take cumulative sums
    //! @return indices for each atom in parent molecule and their corresponding score contributions
    std::vector< std::map< size_t, float> > MoleculeFeatureMapper::CompareSubstructuresNaive
    (
      const FragmentComplete &MOLECULE,
      const FragmentEnsemble &ENSEMBLE,
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const bool &NORMALIZE,
      const bool &IGNORE_H,
      const bool &AVERAGE,
      const size_t &MOL_INDEX,
      const bool &OUTPUT_INTERMEDIATES
    ) const
    {
      // require that our scorer is defined
      BCL_Assert( m_Descriptor.IsDefined(), "MoleculeFeatureMapper::CompareSubstructures called before descriptor was set");

      // open a stream to output intermediate structures
      io::OFStream output_intermediate, output_intermediate_fragment;

      // if MOL_INDEX is defined, we can save to different files
      if( OUTPUT_INTERMEDIATES && util::IsDefined( MOL_INDEX))
      {
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates." + util::Format()( MOL_INDEX) + ".sdf.gz");
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags." + util::Format()( MOL_INDEX) + ".sdf.gz");
      }
      else if( OUTPUT_INTERMEDIATES)
      {
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags.combined.sdf", std::ios::app);
      }

      // make modifiable parent molecule
      FragmentComplete molecule( MOLECULE);
      FragmentEnsemble ensemble( ENSEMBLE);

      // score parent molecule
      float parent_score( m_Descriptor->SumOverObject( molecule)( 0));
      BCL_MessageDbg( "Parent molecule score: " + util::Format()( parent_score));
      if( OUTPUT_INTERMEDIATES)
      {
        molecule.WriteMDL( output_intermediate);
      }

      // remove hydrogen atoms so that this doesn't become too slow
      if( IGNORE_H)
      {
        molecule.RemoveH();
      }

      // make a graph of our molecule
      ConformationGraphConverter graph_maker( ATOM_TYPE, BOND_TYPE);
      graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( molecule));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

      // set the parent molecule as the first graph
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso;
      common_subgraph_iso.SetGraphA( mol_graph_ptr);

      // initialize each atom contribution to zero
      const size_t &orig_n_atoms( molecule.GetNumberAtoms());
      std::map< size_t, float> atom_contributions;
      for( size_t i( 0); i < orig_n_atoms; ++i)
      {
        atom_contributions.insert( std::make_pair( i, float( 0.0)));
      }

      // iterate over each non-parent molecule in the ensemble
      size_t mol_index( 0);
      for
      (
          auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        if( mol_index == MOL_INDEX)
        {
          continue;
        }

        // score the other molecule for naive approach
        // we will attribute the score difference more or less equally across the atoms that differ between the two molecules
        float score_diff( util::GetUndefined< float>());
        float score( m_Descriptor->SumOverObject( *mol_itr)( 0));
        BCL_MessageDbg( "Molecule " + util::Format()( mol_index) + " score: " + util::Format()( score));
        score_diff = parent_score - score;
        BCL_MessageDbg( "deltaScore: " + util::Format()( score_diff));

        // can remove hydrogen atoms to speedup substructure search
        if( OUTPUT_INTERMEDIATES)
        {
          mol_itr->WriteMDL( output_intermediate);
        }
        if( IGNORE_H)
        {
          mol_itr->RemoveH();
        }

        // generate a graph for the current non-parent molecule from ensemble
        graph::ConstGraph< size_t, size_t> frag_graph( graph_maker( *mol_itr));
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > frag_graph_ptr( &frag_graph, false);

        // set the second graph to be the fragment graph
        common_subgraph_iso.SetGraphB( frag_graph_ptr);

        // get the subgraph isomorphism (common subgraph between the two molecules)
        common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), 1);
        graph::Subgraph< size_t, size_t> subgraph
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement()
        );

        // subgraph of the reference molecule
        graph::Subgraph< size_t, size_t> subgraph_b
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB().FirstElement()
        );

        // make new mols from subgraphs
        storage::Vector< size_t> indices( subgraph.GetVertexIndices()), indices_b( subgraph_b.GetVertexIndices());
        AtomVector< AtomComplete> subgraph_v( molecule.GetAtomVector()), subgraph_v_b( mol_itr->GetAtomVector());
        subgraph_v.Reorder( indices);
        subgraph_v_b.Reorder( indices_b);
        FragmentComplete new_mol( subgraph_v, ""), new_mol_b( subgraph_v_b, "");

        // output intermediate mols
        if( OUTPUT_INTERMEDIATES)
        {
          new_mol.WriteMDL( output_intermediate);
          new_mol_b.WriteMDL( output_intermediate);
        }

        // invert the subgraph isomorphism to obtain the subgraph containing the differing vertices
        graph::Subgraph<size_t, size_t> complement_subgraph( subgraph.GetComplement());
        graph::Subgraph<size_t, size_t> complement_subgraph_b( subgraph_b.GetComplement());
        storage::Vector< size_t> complement_indices( complement_subgraph.GetVertexIndices());
        storage::Vector< size_t> complement_indices_b( complement_subgraph_b.GetVertexIndices());

        // make fragments from complement subgraphs
        AtomVector< AtomComplete> subgraph_fragment_v( molecule.GetAtomVector());
        AtomVector< AtomComplete> subgraph_fragment_v_b( mol_itr->GetAtomVector());
        subgraph_fragment_v.Reorder( complement_indices);
        subgraph_fragment_v_b.Reorder( complement_indices_b);
        FragmentComplete subgraph_fragment( subgraph_fragment_v, "");
        FragmentComplete subgraph_fragment_b( subgraph_fragment_v_b, "");

        // output fragments
        if( OUTPUT_INTERMEDIATES)
        {
          subgraph_fragment.WriteMDL( output_intermediate_fragment);
          subgraph_fragment_b.WriteMDL( output_intermediate_fragment);
        }

        // save the score differences to the subgraph atoms
        for
        (
            auto itr( complement_indices.Begin()),
            itr_end( complement_indices.End());
            itr != itr_end;
            ++itr
        )
        {
          // evenly distribute score difference across all differing atoms between the two molecules
          if( NORMALIZE)
          {
            // exclude hydrogen atoms from this accounting,
            // so only divide score difference by the heavy atom count of the unshared fragment
            if( IGNORE_H)
            {
              size_t n_hydrogen_atoms( subgraph_fragment.GetNumberHydrogens());
              atom_contributions[ *itr] += score_diff / ( subgraph.GetVertexIndices().GetSize() - n_hydrogen_atoms);
            }
            else
            {
              atom_contributions[ *itr] += score_diff / subgraph.GetVertexIndices().GetSize();
            }
          }
          // all differing atoms are given equal accounting
          else
          {
            atom_contributions[ *itr] += score_diff;
          }
        }
      }

      // get the average contribution of each atom across all comparisons to
      // other molecules in the ensemble
      if( AVERAGE)
      {
        for( size_t i( 0); i < orig_n_atoms; ++i)
        {
          util::IsDefined( MOL_INDEX) ?
              atom_contributions[ i] = atom_contributions[ i] / ( ENSEMBLE.GetSize() - 1) :
              atom_contributions[ i] = atom_contributions[ i] / ENSEMBLE.GetSize();
        }
      }

      // obtain return values
      std::vector< std::map< size_t, float> > return_vals( 1);
      return_vals[ 0].insert( atom_contributions.begin(), atom_contributions.end());

      // close these output streams
      io::File::CloseClearFStream( output_intermediate);
      io::File::CloseClearFStream( output_intermediate_fragment);

      // done
      return return_vals;
    }

    //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
    //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
    //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
    //! @param MOLECULE parent molecule to be perturbed
    //! @param ENSEMBLE ensemble of molecules against which parent molecule will be compared
    //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
    //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
    //! @param AVERAGE if true, compute average value per atom across all comparisons; otherwise take cumulative sums
    //! @return indices for each atom in parent molecule and their corresponding score contributions
    std::vector< std::map< size_t, float> > MoleculeFeatureMapper::CompareSubstructuresRigorous
    (
      const FragmentComplete &MOLECULE,
      const FragmentEnsemble &ENSEMBLE,
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const bool &MUTUALLY_MATCHING_ATOMS,
      const bool &NORMALIZE,
      const bool &AVERAGE,
      const size_t &MOL_INDEX,
      const bool &OUTPUT_INTERMEDIATES
    ) const
    {
      // require that our scorer is defined
      BCL_Assert( m_Descriptor.IsDefined(), "MoleculeFeatureMapper::CompareSubstructures called before descriptor was set");

      // open a stream to output intermediate structures
      io::OFStream output_intermediate_orig, output_intermediate, output_intermediate_fragment;

      // if MOL_INDEX is defined, we can save to different files
      if( OUTPUT_INTERMEDIATES && util::IsDefined( MOL_INDEX))
      {
        io::File::MustOpenOFStream( output_intermediate_orig, "FeatureMapper.input.intermediates." + util::Format()( MOL_INDEX) + ".sdf.gz");
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates." + util::Format()( MOL_INDEX) + ".sdf.gz");
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags." + util::Format()( MOL_INDEX) + ".sdf.gz");
      }
      else if( OUTPUT_INTERMEDIATES)
      {
        io::File::MustOpenOFStream( output_intermediate_orig, "FeatureMapper.input.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate, "FeatureMapper.intermediates.combined.sdf", std::ios::app);
        io::File::MustOpenOFStream( output_intermediate_fragment, "FeatureMapper.intermediate_frags.combined.sdf", std::ios::app);
      }

      // make modifiable parent molecule
      FragmentComplete molecule( MOLECULE);
      FragmentEnsemble ensemble( ENSEMBLE);

      // score parent molecule
      float parent_score( util::GetUndefinedDouble());
      FragmentComplete temp_parent( molecule);
      temp_parent.SaturateWithH();
      parent_score = m_Descriptor->SumOverObject( temp_parent)( 0);
      BCL_MessageDbg( "Parent molecule score: " + util::Format()( parent_score));

      // write intermediates
      if( OUTPUT_INTERMEDIATES)
      {
        molecule.WriteMDL( output_intermediate_orig);
      }

      // make a graph of our parent molecule
      ConformationGraphConverter graph_maker( ATOM_TYPE, BOND_TYPE);
      FragmentSplitRings split_rings( true);
      FragmentGraphMarker graph_marker( graph_maker, split_rings);
      graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( molecule));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

      // set the parent molecule as the first graph for the isomorphism
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso( graph::CommonSubgraphIsomorphismBase::e_Connected);
      common_subgraph_iso.SetGraphA( mol_graph_ptr);

      // initialize each parent molecule atom contribution to zero
      const size_t &orig_n_atoms( molecule.GetNumberAtoms());
      std::map< size_t, float> atom_contributions;
      for( size_t i( 0); i < orig_n_atoms; ++i)
      {
        atom_contributions.insert( std::make_pair( i, float( 0.0)));
      }

      // iterate over each non-parent molecule in the ensemble
      size_t mol_index( 0);
      for
      (
          auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        if( mol_index == MOL_INDEX)
        {
          continue;
        }

        // get the reference molecule score
        float reference_score( util::GetUndefinedDouble());
        FragmentComplete temp_reference( *mol_itr);
        temp_reference.SaturateWithH();
        reference_score = m_Descriptor->SumOverObject( temp_reference)( 0);
        BCL_MessageDbg( "Reference molecule score: " + util::Format()( reference_score));

        if( OUTPUT_INTERMEDIATES)
        {
          mol_itr->WriteMDL( output_intermediate_orig);
        }

        // generate a graph for the current non-parent molecule from ensemble
        graph::ConstGraph< size_t, size_t> frag_graph( graph_maker( *mol_itr));
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > frag_graph_ptr( &frag_graph, false);

        // set the second graph to be the reference molecule graph
        common_subgraph_iso.SetGraphB( frag_graph_ptr);

        // get the subgraph isomorphism (common subgraph between the two molecules) of the parent molecule
        if( MUTUALLY_MATCHING_ATOMS)
        {
          // identify matched atoms
          storage::Vector< size_t> matched_atoms_parent, matched_atoms_reference;
          ConformationComparisonPsiField::GetAlignedAtoms( molecule, *mol_itr, matched_atoms_parent, matched_atoms_reference);

          // make a map from the matched atoms
          storage::Map< size_t, size_t> matching_atoms;
          BCL_Assert( matched_atoms_parent.GetSize() == matched_atoms_reference.GetSize(), "Unequal number of mutually matching atoms!");
          for
          (
              size_t i( 0), i_sz( matched_atoms_parent.GetSize()); i < i_sz; ++i
          )
          {
            matching_atoms.Insert( std::make_pair( matched_atoms_parent( i), matched_atoms_reference( i)));
          }

          // set the isomorphism from the matched atoms map
          common_subgraph_iso.SetIsomorphism( matching_atoms);
        }
        else
        {
          common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), size_t( 1));
        }

        // get parent and reference subgraphs
        graph::Subgraph< size_t, size_t> subgraph_parent
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement()
        );
        graph::Subgraph< size_t, size_t> subgraph_reference
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB().FirstElement()
        );

        // create graphs lacking all intra-isomorphism edges
        graph::ConstGraph< size_t, size_t> frag_graph_depleted( frag_graph);
        graph::ConstGraph< size_t, size_t> molecule_graph_depleted( mol_graph);
        for( size_t atom_i( 0), iso_sz( subgraph_reference.GetSize()); atom_i < iso_sz; ++atom_i)
        {
          molecule_graph_depleted.RemoveAllEdges( subgraph_parent.GetVertexIndices()( atom_i));
          frag_graph_depleted.RemoveAllEdges( subgraph_reference.GetVertexIndices()( atom_i));
        }

        // get the largest isomorphism mapping atoms from parent to reference
        storage::Map< size_t, size_t> largest_iso( common_subgraph_iso.GetIsomorphism());

        // get edges connecting subgraph to supergraph
        storage::Vector< storage::Vector< size_t> > subgraph_parent_edges( subgraph_parent.GetOrderedAdjacentEdgeIndices());
        storage::Vector< storage::Vector< size_t> > subgraph_reference_edges( subgraph_reference.GetOrderedAdjacentEdgeIndices());

        // get subgraph vertices
        auto const &graph_a_subgraph_indices( subgraph_parent.GetVertexIndices());
        auto const &graph_b_subgraph_indices( subgraph_reference.GetVertexIndices());

        // iterate over the atoms connecting the subgraph to the supergraph
        for( size_t atom_i( 0), iso_sz( subgraph_parent_edges.GetSize()); atom_i < iso_sz; ++atom_i)
        {
          // if this atom from the subgraph does not connect wtih the supergraph then skip it
          // i.e. this atom does not have a supergraph component not in the subgraph
          if( subgraph_parent_edges( atom_i).IsEmpty())
          {
            continue;
          }
          // for each atom in each subgraph, grab all its adjacent edges and their connected atoms using graph connectivity getverticesreachablefromdirectededge.
          // stick them all into two vectors (one for each molecule) of sets (one for each atom in each graph that is within the isomorphism)
          // of size_t s or some other suitable non redundant hashing data structure
          storage::Vector< size_t> parent_attribute_score;
          for
          (
              auto subgraph_parent_itr( subgraph_parent_edges( atom_i).Begin()),
              subgraph_parent_itr_end( subgraph_parent_edges( atom_i).End());
              subgraph_parent_itr != subgraph_parent_itr_end;
              ++subgraph_parent_itr
          )
          {
            if( parent_attribute_score.Find( *subgraph_parent_itr) < parent_attribute_score.GetSize())
            {
              // cycle; already grabbed this vertex, continue
              continue;
            }

            // grab the inverted vertices so that we can track score differences easily
            util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_remove
            (
              graph::Connectivity::GetVerticesReachableFrom
              (
                molecule_graph_depleted,
                *subgraph_parent_itr
              )
            );

            parent_attribute_score.Append( *reachable_vertices_to_remove);
          }

          const storage::Vector< size_t> parent_reachable
          (
            graph::Subgraph< size_t, size_t>( mol_graph_ptr, parent_attribute_score).GetComplement().GetVertexIndices()
          );

          // parent
          AtomVector< AtomComplete> parent_keep( molecule.GetAtomVector()), parent_attribute( molecule.GetAtomVector());
          parent_keep.Reorder( parent_reachable);
          parent_attribute.Reorder( parent_attribute_score);
          FragmentComplete parent_keep_fc( parent_keep, "");
          FragmentComplete parent_attribute_fc( parent_attribute, "");

          // score modified parent
          float parent_keep_fc_score( util::GetUndefinedDouble());
          FragmentComplete parent_temp_copy( parent_keep_fc);
          parent_temp_copy.SaturateWithH();
          parent_keep_fc_score = m_Descriptor->SumOverObject( parent_temp_copy)( 0);
          BCL_MessageDbg( "Modified parent molecule score: " + util::Format()( parent_keep_fc_score));

          // for each atom in each subgraph, grab all its adjacent edges and their connected atoms using graph connectivity getverticesreachablefromdirectededge.
          // stick them all into two vectors (one for each molecule) of sets (one for each atom in each graph that is within the isomorphism)
          // of size_t s or some other suitable non redundant hashing data structure
          storage::Vector< size_t> reference_attribute_score;
          for
          (
              auto subgraph_reference_itr( subgraph_reference_edges( atom_i).Begin()),
              subgraph_reference_itr_end( subgraph_reference_edges( atom_i).End());
              subgraph_reference_itr != subgraph_reference_itr_end;
              ++subgraph_reference_itr
          )
          {
            if( reference_attribute_score.Find( *subgraph_reference_itr) < reference_attribute_score.GetSize())
            {
              // cycle; already grabbed this vertex, continue
              continue;
            }
            // grab the inverted vertices so that we can track score differences easily
            util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_remove
            (
              graph::Connectivity::GetVerticesReachableFrom
              (
                frag_graph_depleted,
                *subgraph_reference_itr
              )
            );
            reference_attribute_score.Append( *reachable_vertices_to_remove);
          }

          const storage::Vector< size_t> reference_reachable
          (
            graph::Subgraph< size_t, size_t>( frag_graph_ptr, reference_attribute_score).GetComplement().GetVertexIndices()
          );

          // reference
          AtomVector< AtomComplete> reference_keep( mol_itr->GetAtomVector()), reference_attribute( mol_itr->GetAtomVector());
          reference_keep.Reorder( reference_reachable);
          reference_attribute.Reorder( reference_attribute_score);
          FragmentComplete reference_keep_fc( reference_keep, "");
          FragmentComplete reference_attribute_fc( reference_attribute, "");

          // score modified reference
          float reference_keep_fc_score( util::GetUndefinedDouble());
          FragmentComplete reference_temp_copy( reference_keep_fc);
          reference_temp_copy.SaturateWithH();
          reference_keep_fc_score = m_Descriptor->SumOverObject( reference_temp_copy)( 0);
          BCL_MessageDbg( "Modified reference molecule score: " + util::Format()( reference_keep_fc_score));

          // compute deltaScore
          float score_diff_parent( parent_keep_fc_score - parent_score);
          float score_diff_reference( reference_keep_fc_score - reference_score);
          float score_diff( score_diff_reference - score_diff_parent);
          BCL_MessageDbg( "deltaScore: " + util::Format()( score_diff));

          // write intermediates
          if( OUTPUT_INTERMEDIATES)
          {
            parent_keep_fc.WriteMDL( output_intermediate);
            parent_attribute_fc.WriteMDL( output_intermediate_fragment);
            reference_keep_fc.WriteMDL( output_intermediate);
            reference_attribute_fc.WriteMDL( output_intermediate_fragment);
          }

          // save scores
          for
          (
              auto score_itr( parent_attribute_score.Begin()), score_itr_end( parent_attribute_score.End());
              score_itr != score_itr_end;
              ++score_itr
          )
          {
            atom_contributions[ *score_itr] += score_diff;
          }
        }
      }

      // get the average contribution of each atom across all comparisons to
      // other molecules in the ensemble
      if( AVERAGE)
      {
        for( size_t i( 0); i < orig_n_atoms; ++i)
        {
          util::IsDefined( MOL_INDEX) ?
              atom_contributions[ i] = atom_contributions[ i] / ( ENSEMBLE.GetSize() - 1) :
              atom_contributions[ i] = atom_contributions[ i] / ENSEMBLE.GetSize();
        }
      }

      // done
      std::vector< std::map< size_t, float> > return_vals( 1);
      return_vals[ 0].insert( atom_contributions.begin(), atom_contributions.end());

      // close these output streams
      io::File::CloseClearFStream( output_intermediate_orig);
      io::File::CloseClearFStream( output_intermediate);
      io::File::CloseClearFStream( output_intermediate_fragment);

      // be happy, be done
      return return_vals;
    }

    storage::Vector< storage::Triplet< size_t, size_t, float> > MoleculeFeatureMapper::GetFeatureRMSDOnScaffold
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &MOLS,
      const std::vector< std::vector< std::map< size_t, float> > > &MAPS,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO
    ) const
    {
      if( MAPS.size() != MOLS.GetSize())
      {
        return storage::Vector< storage::Triplet< size_t, size_t, float> >();
      }
      ConformationGraphConverter graph_maker( ATOM_COMPARISON_TYPE, BOND_TYPE_INFO);

      storage::Vector< graph::ConstGraph< size_t, size_t> > graphs( MOLS.GetSize());
      storage::Vector< storage::Vector< size_t> > isomorphisms( graphs.GetSize());

      graph::ConstGraph< size_t, size_t> scaff_graph( graph_maker( SCAFFOLD));

      graph::SubgraphIsomorphism< size_t, size_t> iso_search;
      iso_search.SetSubgraphExternalOwnership( scaff_graph);

      size_t mol_no( 0);
      for
      (
          FragmentEnsemble::const_iterator itr_mol( MOLS.Begin()), itr_mol_end( MOLS.End());
          itr_mol != itr_mol_end;
          ++itr_mol, ++mol_no
      )
      {
        graphs( mol_no) = graph_maker( *itr_mol);
        iso_search.SetGraphExternalOwnership( graphs( mol_no));
        if( iso_search.FindIsomorphism())
        {
          isomorphisms( mol_no) = iso_search.GetIsomorphism();
        }
      }

      storage::Vector< storage::Triplet< size_t, size_t, float> > result_vector;
      result_vector.AllocateMemory( graphs.GetSize() * ( graphs.GetSize() - 1));
      for( size_t i( 0); i < graphs.GetSize(); ++i)
      {
        for( size_t j( i + 1); j < graphs.GetSize(); ++j)
        {
          if( isomorphisms( i).IsEmpty() || isomorphisms( j).IsEmpty())
          {
            continue;
          }

          float comparison
          (
            CompareCommonStructs
            (
              isomorphisms( i),
              MAPS[ i][ 0],
              isomorphisms( j),
              MAPS[ j][ 0]
            )
          );

          if( util::IsDefined( comparison))
          {
            result_vector.PushBack( storage::Triplet< size_t, size_t, float>( i, j, comparison));
          }

        }
      }
      return result_vector;
    }

    float MoleculeFeatureMapper::CompareCommonStructs
    (
      const storage::Vector< size_t> &SCAFF_ISO_MOL_1,
      const std::map< size_t, float> &PERTURBS_MOL_1,
      const storage::Vector< size_t> &SCAFF_ISO_MOL_2,
      const std::map< size_t, float> &PERTURBS_MOL_2
    ) const
    {
      float result( 0);
      if( SCAFF_ISO_MOL_1.GetSize() != SCAFF_ISO_MOL_2.GetSize())
      {
        BCL_MessageStd( "Molecules had different sized scaffold isomorphisms, cannot compare");
        return util::GetUndefined< float>();
      }

      size_t iso_size( SCAFF_ISO_MOL_1.GetSize());

      for( size_t a( 0); a < iso_size; ++a)
      {
        size_t a_1( SCAFF_ISO_MOL_1( a)), a_2( SCAFF_ISO_MOL_2( a));
        std::map< size_t, float>::const_iterator itr_1( PERTURBS_MOL_1.find( a_1));
        std::map< size_t, float>::const_iterator itr_2( PERTURBS_MOL_2.find( a_2));
        if( itr_1 == PERTURBS_MOL_1.end() || itr_2 == PERTURBS_MOL_2.end())
        {
          BCL_MessageStd( "Couldn't find a scaffold atom in the map");
          return util::GetUndefined< float>();
        }
        float diff( itr_1->second - itr_2->second);
        result += std::sqrt( diff * diff);
      }
      return result;
    }

    //! @brief calculates average atom scores over a common scaffold
    std::vector< std::map< size_t, float> > MoleculeFeatureMapper::AverageFeatures
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &MOLECULES,
      const std::vector< std::vector< std::map< size_t, float> > > &MAPS,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO
    ) const
    {

      ConformationGraphConverter graph_maker( ATOM_COMPARISON_TYPE, BOND_TYPE_INFO);

      size_t n_mols( MOLECULES.GetSize());
      bool use_maps( MAPS.size() == n_mols);

      std::vector< FeatureMapInfo> mol_infos( n_mols);
      graph::ConstGraph< size_t, size_t> scaff_graph( graph_maker( SCAFFOLD));
      graph::SubgraphIsomorphism< size_t, size_t> find_iso;
      find_iso.SetSubgraphExternalOwnership( scaff_graph);

      size_t mol_no( 0);
      math::RunningAverage< linal::Vector< float> > avg;
      for
      (
          FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin()), itr_mol_end( MOLECULES.End());
          itr_mol != itr_mol_end;
          ++itr_mol, ++mol_no
      )
      {
        mol_infos[ mol_no].m_Graph = graph_maker( *itr_mol);
        find_iso.SetGraphExternalOwnership( mol_infos[ mol_no].m_Graph);

        if( !find_iso.FindIsomorphism())
        {
          BCL_MessageStd( "Molecule " + util::Format()( mol_no) + " did not contain the scaffold");
          continue;
        }

        mol_infos[ mol_no].m_Molecule = util::SiPtr< const FragmentComplete>( &( *itr_mol));
        mol_infos[ mol_no].m_ScaffoldIso = find_iso.GetIsomorphism();
        mol_infos[ mol_no].m_Score = m_Descriptor->SumOverObject( *itr_mol)( 0);
        mol_infos[ mol_no].m_Perturbations = use_maps ? MAPS[ mol_no] : Perturb( *itr_mol, false, false, mol_infos[ mol_no].m_ScaffoldIso);

        if( !mol_infos[ mol_no].m_Perturbations.empty())
        {
          size_t n( mol_infos[ mol_no].m_ScaffoldIso.GetSize());
          linal::Vector< float> perturbations( n, float( 0));
          for( size_t i( 0); i < n; ++i)
          {
            size_t atom_no( mol_infos[ mol_no].m_ScaffoldIso( i));
            perturbations( i) = mol_infos[ mol_no].m_Perturbations[ 0][ atom_no];
          }
          avg += perturbations;
        }
      }

      linal::Vector< float> mean( avg.GetAverage());
      std::vector< std::map< size_t, float> > return_vals( 1);
      for( size_t i( 0), n( mean.GetSize()); i < n; ++i)
      {
        return_vals[ 0][ i] = mean( i);
      }
      return return_vals;
    }

    //! @brief extracts fragments from a molecule; keeps atoms whose scores are outside of a
    //! given confidence interval
    //! @param MOLECULE the molecule to extract from
    //! @param ATOM_SCORES scores of each atom
    //! @param CONFIDENCE_INTERVAL any scores inside this range are considered statistically insignificant
    //! @return an ensemble of "significant" fragments
    FragmentEnsemble MoleculeFeatureMapper::ExtractFragments
    (
      const FragmentComplete &MOLECULE,
      const std::map< size_t, float> &ATOM_SCORES,
      const math::Range< float> &CONFIDENCE_INTERVAL
    ) const
    {

      FragmentEnsemble fragments;

      if( !MOLECULE.GetNumberAtoms() || ATOM_SCORES.empty())
      {
        return fragments;
      }

      ConformationGraphConverter converter;
      ConformationGraphConverter::t_AtomGraph mol_graph( converter.CreateGraphWithAtoms( MOLECULE));

      storage::Vector< size_t> keep;
      for
      (
          std::map< size_t, float>::const_iterator itr_atom( ATOM_SCORES.begin()), itr_atom_end( ATOM_SCORES.end());
          itr_atom != itr_atom_end;
          ++itr_atom
      )
      {
        if( !CONFIDENCE_INTERVAL.IsWithin( itr_atom->second))
        {
          keep.PushBack( itr_atom->first);
        }
      }

      // Do nothing if nothing was significant
      if( keep.IsEmpty())
      {
        return fragments;
      }

      ConformationGraphConverter::t_AtomGraph subgraph( mol_graph.GetSubgraph( keep));
      storage::List< storage::Vector< size_t> > components
      (
        graph::Connectivity::GetComponents( subgraph)
      );

      for
      (
          storage::List< storage::Vector< size_t> >::const_iterator itr_comp( components.Begin()), itr_comp_end( components.End());
          itr_comp != itr_comp_end;
          ++itr_comp
      )
      {
        ConformationGraphConverter::t_AtomGraph this_frag( subgraph.GetSubgraph( *itr_comp));
        fragments.PushBack
        (
          FragmentComplete
          (
            ConformationGraphConverter::CreateAtomsFromGraph( this_frag),
            ""
          )
        );
      }

      return fragments;
    }

    //! @brief Set the members with LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeFeatureMapper::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
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

    io::Serializer MoleculeFeatureMapper::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "determines the features of a molecule that are significant to QSAR models");
      member_data.AddInitializer
      (
        "model",
        "the QSAR model to use",
        io::Serialization::GetAgent( &m_Descriptor)
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
