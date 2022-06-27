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
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"

// includes, sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_conformation_set.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_perturb_molecule_pose.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_score_function_generic.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_metropolis.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_distribution_interface.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_implementation.h"
// external includes - sorted alphabetically

#undef AddAtom
#undef RemoveAtom
#undef ATOMS

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // initialize static
    sched::Mutex &FragmentStochasticPoseOptimizer::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentStochasticPoseOptimizer::FragmentStochasticPoseOptimizer() :
        m_PropertyScorer( descriptor::CheminfoProperty()),
        m_BFactors( storage::Vector< float>()),
        m_BindingPocket( FragmentComplete()),
        m_VoxelGridPocket( VoxelGridAtom( 4.0)),
        m_MDL( std::string()),
        m_BindingPocketFilename( std::string()),
        m_Iterations( size_t( 100)),
        m_Attempts( size_t( 1)),
        m_RefinementIterations( 100),
        m_Temperature( float( 1.0)),
        m_Opti( opti::e_SmallerIsBetter),
        m_VDWClashCutoff( float( 5.0))
    {
    }

//    //! @brief pose resolver constructor
//    //! @param PROPERTY_SCORER enables pose-dependent optimization of score
//    //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
//    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
//    FragmentStochasticPoseOptimizer::FragmentStochasticPoseOptimizer
//    (
//      const descriptor::CheminfoProperty &PROPERTY_SCORER,
//      const storage::Vector< float> &BFACTORS,
//      const FragmentComplete &BINDING_POCKET,
//      const size_t &ITERATIONS
//    ) :
//        m_PropertyScorer( PROPERTY_SCORER),
//        m_BFactors( BFACTORS),
//        m_BindingPocket( BINDING_POCKET),
//        m_Iterations( ITERATIONS)
//    {
//      // if the b-factors are empty then just assign a vector size according to each protein residue and a value of 1
//      if( !m_BFactors.GetSize())
//      {
//        m_BFactors = storage::Vector< float>( m_BindingPocket.GetSize(), float(2.0));
//      }
//
//      // now make sure that the b factors are actually the correct length (i.e. the user did not pass the wrong length)
//      BCL_Assert(m_BFactors.GetSize() == m_BindingPocket.GetSize(), "B-factor size is not equal to number of atoms in the provided protein pocket");
//    }

    //! @brief constructor with binding pocket to be created from MDL property
    //! @param PROPERTY_SCORER enables pose-dependent optimization of score
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    //! @param MDL the SDF file MDL property specifying the binding pocket filename
    //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
    FragmentStochasticPoseOptimizer::FragmentStochasticPoseOptimizer
    (
      const descriptor::CheminfoProperty &PROPERTY_SCORER,
      const storage::Vector< float> &BFACTORS,
      const std::string &MDL,
      const std::string &BINDING_POCKET_FILENAME,
      const size_t &ITERATIONS,
      const size_t &ATTEMPTS,
      const size_t &REFINEMENT_ITERATIONS,
      const float &TEMPERATURE,
      const opti::Tracker< FragmentComplete, double> &OPTI,
      const float &VDW_CLASH_CUTOFF
    ) :
        m_PropertyScorer( PROPERTY_SCORER),
        m_BFactors( BFACTORS),
        m_BindingPocket( FragmentComplete()),
        m_VoxelGridPocket( VoxelGridAtom( 4.0)),
        m_MDL( MDL),
        m_BindingPocketFilename( BINDING_POCKET_FILENAME),
        m_Iterations( ITERATIONS),
        m_Attempts( ATTEMPTS),
        m_RefinementIterations( REFINEMENT_ITERATIONS),
        m_Temperature( TEMPERATURE),
        m_Opti( OPTI),
        m_VDWClashCutoff( VDW_CLASH_CUTOFF)
    {
      // read in the protein binding pocket
      if
      (
          m_BindingPocketFilename.size() &&
          !m_BindingPocket.GetSize()
      )
      {
        // Read in PDB
        BCL_MessageVrb( "Caching pocket: " + util::Format()( m_BindingPocketFilename));
        const pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
        assemble::ProteinModel protein_model( factory.ProteinModelFromPDBFilename( m_BindingPocketFilename, pdb::Factory::GetSSETypeMinSizes( 0, 0, 0)));

        //instantiate AASideChainFactory no hydrogens include backbone atoms for sidechain placement superimposition
        biol::AASideChainFactory side_chain_factory( false, true);

        // add side chains to model and set it to model
        protein_model = *side_chain_factory.ProteinModelWithSideChains( protein_model);
        AAFragmentComplete aa_fragment( protein_model.GetAminoAcids(), true);
        aa_fragment.SaturateWithH();
        aa_fragment.RemoveAtomsUndefinedPos();

        // get final model
        m_BindingPocket = aa_fragment;

        // Setup the protein binding pocket voxel grid for clash detection
        m_VoxelGridPocket.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( m_BindingPocket.GetAtomsIterator().Begin(), m_BindingPocket.GetAtomsIterator().End()));
      }

      // if the b-factors are empty then just assign a vector size according to each protein residue and a value of 1
      if( !m_BFactors.GetSize())
      {
        m_BFactors = storage::Vector< float>( m_BindingPocket.GetSize(), float( 2.0));
      }

      // now make sure that the b factors are actually the correct length (i.e. the user did not pass the wrong length)
      BCL_Assert( m_BFactors.GetSize() == m_BindingPocket.GetSize(), "B-factor size is not equal to number of atoms in the provided protein pocket");
    }

    //! @brief clone constructor
    FragmentStochasticPoseOptimizer *FragmentStochasticPoseOptimizer::Clone() const
    {
      return new FragmentStochasticPoseOptimizer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentStochasticPoseOptimizer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief resolve 3D ligand clashes with a protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param BINDING_POCKET the protein binding pocket
    util::ShPtr< FragmentEnsemble> FragmentStochasticPoseOptimizer::ResolveClashes
    (
      const FragmentEnsemble &ENSEMBLE
    )
    {

      // Build molecule voxel grid
      VoxelGridAtom voxel_grid_mol( 4.0);
      voxel_grid_mol.Clear();

      // loop over molecules
      storage::Vector< size_t> mols_to_remove;
      size_t mol_index( 0);
      for
      (
          auto mol_itr( ENSEMBLE.Begin()), mol_itr_end( ENSEMBLE.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        // Setup molecule voxel grid
        voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( mol_itr->GetAtomsIterator().Begin(), mol_itr->GetAtomsIterator().End()));

        // determine molecule atom neighboring atoms in the protein
        auto neighbors( voxel_grid_mol.GetNeighborsIn( m_VoxelGridPocket, 4.0));

        storage::Vector< linal::Vector< double>> radii( neighbors.GetSize(), linal::Vector< double>( size_t( 2)));
        size_t neighbors_index( 0);

        // iterate over protein atom neighbors to each molecule atom
        for( auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors, ++neighbors_index)
        {
          //Get distance between neighbor atoms
          const double distance( itr_neighbors->Third());

          //Get b-factor
          double bfactor( m_BFactors( m_BindingPocket.GetAtomIndex( *itr_neighbors->Second())));

          //Compute csd_vdw radius of molecule atom
          const double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
          radii( neighbors_index)( 0) = mol_csd_vdw_radius;

          //Compute csd_vdw radius of pocket atom
          const double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
          radii( neighbors_index)( 1) = pocket_csd_vdw_radius;

          //Compute the sum of the vdw radii of the molecule atom and the pocket atom
          const double min_radius( mol_csd_vdw_radius + pocket_csd_vdw_radius);

          // If distance is less than min_radius then there is a collision
          bool raw_clash;
          if( distance <= ( min_radius - bfactor)) // b-factor is just a default tolerance of 1.0 right now
          {
            raw_clash = true;
          }
          else
          {
            raw_clash = false;
          }

          // stupid check ; make better
          // perhaps something with AtomClashScore, but would need a new constructor to work with atom pairs
          if( raw_clash)
          {
            mols_to_remove.PushBack( mol_index);
            break;
          }
        } // end atom neighbors iterator
      } // end conformer iterator

      // remove molecules with clashes from ensemble
      storage::List< FragmentComplete> good_confs;
      BCL_MessageStd( "Remove clashing confs");
      size_t conf_index( 0);
      for
      (
        auto conf_itr( ENSEMBLE.Begin()), conf_itr_end( ENSEMBLE.End());
        conf_itr != conf_itr_end;
        ++conf_itr, ++conf_index
      )
      {
        bool good_conf( true);
        for
        (
            size_t bad_conf_indices( 0);
            bad_conf_indices < mols_to_remove.GetSize();
            ++bad_conf_indices
        )
        {
          if( conf_index == mols_to_remove( bad_conf_indices))
          {
            good_conf = false;
            break;
          }
        }
        if( good_conf)
        {
          good_confs.PushBack( *conf_itr);
        }
      }
      if( good_confs.GetSize())
      {
        return util::ShPtr< FragmentEnsemble>( new FragmentEnsemble( good_confs));
      }
      else
      {
        return util::ShPtr< FragmentEnsemble>( new FragmentEnsemble());
      }
    }

    //! @brief optimize pose-dependent score of ligand with protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param BINDING_POCKET the protein binding pocket
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::OptimizePoseConf
    (
      const FragmentEnsemble &ENSEMBLE
    )
    {
      // copy ensemble
      FragmentEnsemble ensemble( ENSEMBLE);

      // make sure have an ensemble, or else bail
      if( ensemble.GetSize())
      {
        // iterate over conformational ensemble and score confs
        storage::Vector< double> scores( ensemble.GetSize(), util::GetUndefinedDouble());
        size_t score_index( 0);
        for
        (
            auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
            conf_itr != conf_itr_end;
            ++conf_itr, ++score_index
        )
        {
          // score
          GetMutex().Lock();
          scores( score_index) = m_PropertyScorer->SumOverObject( *conf_itr)( 0);
          GetMutex().Unlock();
          conf_itr->StoreProperty( "FLD_Pose_Score", linal::Vector< double>( 1, -1 * scores( score_index)));
        }
        ensemble.Sort( "FLD_Pose_Score"); // sorts in ascending order

        // iterate over all conformers
        for( auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End()); mol_itr != mol_itr_end; ++mol_itr)
        {
          // return if it does not clash
          if( !ProteinLigandClash( *mol_itr))
          {
            return util::ShPtr< FragmentComplete>( new FragmentComplete( *mol_itr));
          }
        }
        // if they all clash then return empty
        return util::ShPtr< FragmentComplete>( new FragmentComplete());
      }
      else
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete());
      }
    }

    //! @brief optimize pose-dependent score of ligand with protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param PROPERTY_SCORER the score cheminfo property
    //! @param BINDING_POCKET the protein binding pocket
    //! @param BFACTORS the clash tolerance at each protein atom
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::GetBestScoringNonClashingPose
    (
      const FragmentEnsemble &ENSEMBLE
    )
    {
      // copy ensemble
      FragmentEnsemble ensemble( ENSEMBLE);

      // iterate over conformational ensemble
      storage::Vector< double> scores( ensemble.GetSize(), util::GetUndefinedDouble());
      size_t score_index( 0);
      for
      (
          auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
          conf_itr != conf_itr_end;
          ++conf_itr, ++score_index
      )
      {
        // score
        scores( score_index) = m_PropertyScorer->SumOverObject( *conf_itr)( 0);

        // keep track of last best molecule
        FragmentComplete last_best( *conf_itr);

        // minor perturb rotation/translation
        for( size_t iteration( 0); iteration < m_Iterations; ++iteration)
        {
          auto this_one( PerturbMoleculePoseFiveDegreesMax( *conf_itr));
          double temp_score( m_PropertyScorer->SumOverObject( *this_one)( 0));
          if( temp_score < scores( score_index))
          {
            temp_score = scores( score_index);
            last_best = *this_one;
          }
        }
        // save final score to molecule for sorting and reporting
        conf_itr->StoreProperty( "FLD_Pose_Score", linal::Vector< double>( 1, -1 * scores( score_index)));

        // for testing purposes, we will allow sampling across the conformers in the perturber,
        // so we will break here. The alternative would be to do enumerative perturbation on each conformer
        break;
      }
      ensemble.Sort( "FLD_Pose_Score"); // sorts in ascending order

      // Build molecule voxel grid
      VoxelGridAtom voxel_grid_mol( 4.0);
      voxel_grid_mol.Clear();

      // loop over molecules
      size_t mol_index( 0);
      for
      (
          auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        // Setup molecule voxel grid
        voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( mol_itr->GetAtomsIterator().Begin(), mol_itr->GetAtomsIterator().End()));

        // determine molecule atom neighboring atoms in the protein
        auto neighbors( voxel_grid_mol.GetNeighborsIn( m_VoxelGridPocket, 4.0));

        storage::Vector< linal::Vector< double>> radii( neighbors.GetSize(), linal::Vector< double>( size_t( 2)));
        size_t neighbors_index( 0);

        // iterate over protein atom neighbors to each molecule atom
        bool mol_clash( false);
        for( auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors, ++neighbors_index)
        {
          //Get distance between neighbor atoms
          const double distance( itr_neighbors->Third());

          //Get b-factor
          double bfactor( m_BFactors( m_BindingPocket.GetAtomIndex( *itr_neighbors->Second())));

          //Compute csd_vdw radius of molecule atom
          const double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
          radii( neighbors_index)( 0) = mol_csd_vdw_radius;

          //Compute csd_vdw radius of pocket atom
          const double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
          radii( neighbors_index)( 1) = pocket_csd_vdw_radius;

          //Compute the sum of the vdw radii of the molecule atom and the pocket atom
          const double min_radius( mol_csd_vdw_radius + pocket_csd_vdw_radius);

          // If distance is less than min_radius then there is a collision
          bool raw_clash;
          if( distance <= ( min_radius - bfactor)) // b-factor is just a default tolerance of 1.0 right now
          {
            raw_clash = true;
          }
          else
          {
            raw_clash = false;
          }

          // perhaps something with AtomClashScore, but would need a new constructor to work with atom pairs
          if( raw_clash)
          {
            mol_clash = true;
            break;
          }
        } // end atom neighbors iterator

        // if the molecule does not have clashes, then we want it
        if( !mol_clash)
        {
          break;
        }
      } // end conformer iterator

      // get the best molecule
      FragmentComplete best_mol;
      size_t conf_index( 0);
      for
      (
        auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
        conf_itr != conf_itr_end;
        ++conf_itr, ++conf_index
      )
      {
        if( conf_index != mol_index)
        {
          continue;
        }
        best_mol = *conf_itr;
      }

      // make sure we actually found a molecule;
      if( best_mol.GetSize())
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete( best_mol));
      }
      // if not, return empty
      else
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete());
      }
    }

    //! @brief optimize pose-dependent score of ligand with protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param PROPERTY_SCORER the score cheminfo property
    //! @param BINDING_POCKET the protein binding pocket
    //! @param BFACTORS the clash tolerance at each protein atom
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::StochasticPoseOptimization
    (
      const FragmentEnsemble &ENSEMBLE
    )
    {
      // copy ensemble
      FragmentEnsemble ensemble( ENSEMBLE);
      util::ShPtr< FragmentComplete> best_pose( new FragmentComplete( ensemble.GetMolecules().FirstElement()));

      // score
      GetMutex().Lock();
      double best_score( m_PropertyScorer->SumOverObject( *best_pose)( 0));
      GetMutex().Unlock();

      // minor perturb rotation/translation or conformer swap
      for( size_t iteration( 0); iteration < m_Iterations; ++iteration)
      {
        util::ShPtr< FragmentComplete> pose( PerturbMoleculePoseFiveDegreesMax( *best_pose, ensemble));
        pose->ResetCache(); //< ultimately it would be nice to only clear the cached protein-ligand interaction score
        pose->GetStoredPropertiesNonConst().SetMDLProperty( m_MDL, m_BindingPocketFilename);
        GetMutex().Lock();
        double current_score( m_PropertyScorer->SumOverObject( *pose)( 0));
        GetMutex().Unlock();
        pose->StoreProperty( "FLD_Pose_Score", linal::Vector< double>( 1, -1 * current_score));
        if( current_score > best_score)
        {
          best_score = current_score;
          best_pose = pose;
        }
        else
        {
          pose = best_pose;
        }
      }

      // save final score to molecule for sorting and reporting
      best_pose->StoreProperty( "FLD_Pose_Score", linal::Vector< double>( 1, -1 * best_score));
      bool mol_clash( ProteinLigandClash( *best_pose));

      // make sure we actually found a molecule;
      if( best_pose->GetSize() && !mol_clash)
      {
        return best_pose;
      }
      // if not, return empty (i.e. if our best molecule clashes, we just give up)
      else
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete());
      }
    }

    //! @brief MCM sampling of ligand inside of protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param LOCAL disables initial large perturbation
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::MCMDocker
    (
      const FragmentEnsemble &ENSEMBLE,
      const bool &LOCAL
    )
    {
      // choose random conformer from the ensemble
      FragmentEnsemble::const_iterator mol_itr( ENSEMBLE.Begin());
      size_t pos( random::GetGlobalRandom().Random< double>( ENSEMBLE.GetSize()));
      std::advance( mol_itr, pos);

      // copy molecule
      FragmentComplete molecule( *mol_itr);

      // create the temperature control
      util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureDefault( m_Temperature));

      // create the metropolis
      mc::Metropolis< double> metropolis( sp_temperature, true);

      // create the termination criterion
      opti::CriterionCombine< FragmentComplete, double> criterion_combine;

      // create a scoring object
      util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > scorer
      (
        new ScoreFunctionGeneric( m_PropertyScorer)
      );

      // set up our mutater object
      util::ShPtr< math::MutateDecisionNode< FragmentComplete> > mutater
      (
        new math::MutateDecisionNode< FragmentComplete>()
      );

      // if standard dock, make an initial large perturbation
      if( !LOCAL)
      {
        math::RotationMatrix3D rot;
        rot.SetRand( 2.0 * math::g_Pi);
        linal::Vector3D mol_center_coords( molecule.GetCenter());
        molecule.Translate( -mol_center_coords);
        molecule.Rotate( rot);
        molecule.Translate( mol_center_coords);
        const double mag( random::GetGlobalRandom().Random( 3.0));
        linal::Vector3D trans( mag, 0.0, 0.0);
        molecule.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));

        // add the perturb mutate with medium sized moves
        mutater->AddMutate( PerturbMoleculePose( ENSEMBLE), 1.0);

        // insert termination criteria that depends on the total number of MC iterations
        opti::CriterionNumberIterations< FragmentComplete, double> maximum_number_iterations( m_Iterations);
        criterion_combine.InsertCriteria( maximum_number_iterations);

        // insert termination criteria that depends on the total number of unimproved MC iterations
        opti::CriterionUnimproved< FragmentComplete, double> maximum_number_unimproved( m_Iterations * 0.20);
        criterion_combine.InsertCriteria( maximum_number_unimproved);
      }
      // if local dock, do not perform initial large perturbation
      else
      {
        // add the perturb mutate with small sized moves (1 degree rot, 0.1 A trans)
        mutater->AddMutate( PerturbMoleculePose( ENSEMBLE, 0.0174533, 0.1), 1.0);

        // insert termination criteria that depends on the total number of MC iterations
        opti::CriterionNumberIterations< FragmentComplete, double> maximum_number_iterations( m_RefinementIterations);
        criterion_combine.InsertCriteria( maximum_number_iterations);

        // insert termination criteria that depends on the total number of unimproved MC iterations
        opti::CriterionUnimproved< FragmentComplete, double> maximum_number_unimproved( m_RefinementIterations * 0.33333);
        criterion_combine.InsertCriteria( maximum_number_unimproved);
      }

      // Set up MC method
      mc::Approximator< FragmentComplete, double> approximator
      (
        *scorer,
        *mutater,
        metropolis,
        criterion_combine,
        molecule,
        m_Opti
      );

      // run the approximator
      approximator.Approximate();
      FragmentComplete best_mol( approximator.GetTracker().GetBest()->First());
      linal::Vector< double> best_score( 1, approximator.GetTracker().GetBest()->Second());
      best_mol.StoreProperty( "FLD_Score", best_score);

      // debug options
//      io::OFStream debug_mdl_out;
//      io::File::MustOpenOFStream( debug_mdl_out, "MCM_Poses.sdf", std::ios::app);
//      best_mol.WriteMDL( debug_mdl_out);
//      io::File::CloseClearFStream( debug_mdl_out);

      // return to RunMCMDocker
      return util::ShPtr< FragmentComplete>( new FragmentComplete( best_mol));
    }

    //! @brief MCM sampling of ligand inside of protein pocket
    //! @param ENSEMBLE the ensemble of molecule conformers
    //! @param LOCAL disables initial large perturbation
    //! @param REFINE enables local refinement of the best poses
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::RunMCMDocker
    (
      const FragmentEnsemble &ENSEMBLE,
      const bool &LOCAL,
      const bool &REFINE
    )
    {
      // initialize container for MCM poses
      FragmentEnsemble poses;

      // run the MCM docker this many times
      BCL_MessageStd( "Running " + util::Format()( m_Attempts) + " cyles of Monte Carlo - Metropolis!");
      for( size_t attempt( 0); attempt < m_Attempts; ++attempt)
      {
        util::ShPtr< FragmentComplete> mcm_result
        (
          MCMDocker( ENSEMBLE, LOCAL)
        );
        // add the new pose to the set
        BCL_MessageStd( "Cycle " + util::Format()( attempt + 1) + " score: " + util::Format()( mcm_result->GetMDLProperty( "FLD_Score")));
        poses.PushBack( *mcm_result);
      }

      // prepare final poses
      FragmentEnsemble refined_poses;

      // perform local refinement
      if( REFINE)
      {
        //make a local sample confs object
        static RotamerLibraryFile rotamer_library;
        static SampleConformations sample_confs
        (
          rotamer_library,
          "",
          0.0,
          100,
          1000,
          false
        );
        sample_confs.SetSamplingPreferences( false, false, true, false);

        // perform a refinement MCM optimization run
        size_t pose_index( 0);
        for
        (
            auto pose_itr( poses.GetMolecules().Begin()), pose_itr_end( poses.GetMolecules().End());
            pose_itr != pose_itr_end;
            ++pose_itr, ++pose_index
        )
        {
          float score( m_PropertyScorer->SumOverObject( *pose_itr)( 0));
          BCL_MessageStd( "Pre-refinement score: " + util::Format()( score));
          if
          (
              ( m_Opti.GetImprovementType() == opti::e_LargerIsBetter && score >= 1.5) ||
              ( m_Opti.GetImprovementType() == opti::e_SmallerIsBetter && score <= -1.5)
          )
          {
            FragmentEnsemble local_confs( sample_confs( *pose_itr).First());
            util::ShPtr< FragmentComplete> refinement_result
            (
              MCMDocker( local_confs, true)
            );
            // add the new pose to the set
            BCL_MessageStd( "Refinement Cycle " + util::Format()( pose_index + 1) + " score: " + util::Format()( refinement_result->GetMDLProperty( "FLD_Score")));
            refined_poses.PushBack( *refinement_result);
          }
          else
          {
            BCL_MessageStd( "Refinement Cycle " + util::Format()( pose_index + 1) + " skipped due to poor initial score!");
          }
        }
      }
      // output initial docking results without refinement
      else
      {
        refined_poses = poses;
      }

      // choose the best pose
      refined_poses.Sort( "FLD_Score");
      FragmentComplete best_pose;
      if( m_Opti.GetImprovementType() == opti::e_LargerIsBetter)
      {
        best_pose = refined_poses.GetMolecules().LastElement();
      }
      else
      {
        best_pose = refined_poses.GetMolecules().FirstElement();
      }

      // return the best pose
      BCL_MessageStd( "Best score: " + util::Format()( best_pose.GetMDLProperty( "FLD_Score")));
      return util::ShPtr< FragmentComplete>( new FragmentComplete( best_pose));
    }

    //! @brief Detects atomic clashes between a protein and a ligand
    bool FragmentStochasticPoseOptimizer::ProteinLigandClash
    (
      const FragmentComplete &MOLECULE
    )
    {
      // Build Molecule voxel grids
      VoxelGridAtom voxel_grid_mol( 4.0); // consider making this an argument and just clearing before each SetObjects()
      voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( MOLECULE.GetAtomsIterator().Begin(), MOLECULE.GetAtomsIterator().End()));

      // determine molecule atom neighboring atoms in the protein
      auto neighbors( voxel_grid_mol.GetNeighborsIn( m_VoxelGridPocket, 4.0));
      storage::Vector< linal::Vector< double>> radii( neighbors.GetSize(), linal::Vector< double>( size_t( 2)));
      size_t neighbors_index( 0);

      // iterate over protein atom neighbors to each molecule atom
      bool mol_clash( false);
      for
      (
          auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End());
          itr_neighbors != itr_neighbors_end;
          ++itr_neighbors, ++neighbors_index
      )
      {
        //Get distance between neighbor atoms
        const double distance( itr_neighbors->Third());

        //Get b-factor
        double bfactor( m_BFactors( m_BindingPocket.GetAtomIndex( *itr_neighbors->Second())));

        //Compute csd_vdw radius of molecule atom
        const double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
        radii( neighbors_index)( 0) = mol_csd_vdw_radius;

        //Compute csd_vdw radius of pocket atom
        const double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
        radii( neighbors_index)( 1) = pocket_csd_vdw_radius;

        //Compute the sum of the vdw radii of the molecule atom and the pocket atom
        const double min_radius( mol_csd_vdw_radius + pocket_csd_vdw_radius);

        // If distance is less than min_radius then there is a collision
        bool raw_clash;
        if( distance <= ( min_radius - bfactor)) // b-factor is just a default tolerance of 1.0 right now
        {
          raw_clash = true;
        }
        else
        {
          raw_clash = false;
        }

        // perhaps something with AtomClashScore, but would need a new constructor to work with atom pairs
        if( raw_clash)
        {
          mol_clash = true;
          break;
        }
      } // end atom neighbors iterator
      return mol_clash;
    }

    //! @brief Check for a bad ligand conformer VDW score
    //! @param MOLECULE the molecule conformer of interest
    //! @return returns true if the conformer passes (is good), false otherwise
    bool FragmentStochasticPoseOptimizer::GoodLigandVDWScore
    (
      const FragmentComplete &MOLECULE
    )
    {
      // filter molecules with strained 3D conformers
      GetMutex().Lock();
      bool good_conf
      (
        descriptor::GetCheminfoProperties().calc_MoleculeVdwScore->SumOverObject( MOLECULE)( 0) < m_VDWClashCutoff ?
            true :
            false
      );
      GetMutex().Unlock();
      return good_conf;
    }

  ////////////////////////
  // mutation functions //
  ////////////////////////

    //! @brief performs a minor perturbation of the molecule via rotation or translation, then score
    //! @param MOLECULE the molecule whose pose is to be perturbed
    //! @param CONF_ENSEMBLE optionslly a conformer ensemble of the molecule to be perturbed
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::PerturbMoleculePoseFiveDegreesMax
    (
      const FragmentComplete &MOLECULE,
      const FragmentEnsemble &CONF_ENSEMBLE
    )
    {
      // copy the molecule
      FragmentComplete molecule( MOLECULE);

      // rotate molecule randomly with magnitude less than or equal to approximately 5 degree
      double prob( random::GetGlobalRandom().Random< double>( CONF_ENSEMBLE.IsEmpty() ? 0.80 : 1.0));
      if( prob < 0.80)
      {
        // perform a random rotation <= 5 degrees
        math::RotationMatrix3D rot;
        rot.SetRand( 0.087); // in radians
        linal::Vector3D mol_center_coords( molecule.GetCenter());
        molecule.Translate( -mol_center_coords);
        molecule.Rotate( rot);
        molecule.Translate( mol_center_coords);

        // perform a random translation between 0.0 and 0.1 angstroms
        const double mag( random::GetGlobalRandom().Random( 0.1));
        linal::Vector3D trans( mag, 0.0, 0.0);
        molecule.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));
      }

      // swap one conformer for another
      else
      {
        // choose a random conformer
        FragmentEnsemble::const_iterator mol_itr( CONF_ENSEMBLE.Begin());
        size_t pos( random::GetGlobalRandom().Random< double>( CONF_ENSEMBLE.GetSize()));
        std::advance( mol_itr, pos);

        // superimpose with original conf        ormer
        util::SiPtrVector< const linal::Vector3D> scaffold_coords( molecule.GetAtomCoordinates());
        util::SiPtrVector< const linal::Vector3D> molecule_coords( mol_itr->GetAtomCoordinates());
        math::TransformationMatrix3D transform( quality::RMSD::SuperimposeCoordinates( scaffold_coords, molecule_coords));
        molecule = *mol_itr;
        molecule.Transform( transform);
      }

      // return mutate result
      return util::ShPtr< FragmentComplete>( new FragmentComplete( molecule));
    }

    //! @brief takes an ensemble and returns a random conformer
    util::ShPtr< FragmentComplete> FragmentStochasticPoseOptimizer::ConformerSwap
    (
      const FragmentEnsemble &ENSEMBLE
    )
    {
      // choose random conformer
      FragmentEnsemble::const_iterator mol_itr( ENSEMBLE.Begin());
      size_t pos( random::GetGlobalRandom().Random< double>( ENSEMBLE.GetSize()));
      std::advance( mol_itr, pos);

      // return mutate result
      return util::ShPtr< FragmentComplete>( new FragmentComplete( *mol_itr));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return the MDL property associated with this object
    //! @return the MDL string
    const std::string &FragmentStochasticPoseOptimizer::GetMDL() const
    {
      return m_MDL;
    }

    //! @brief return the pocket filename associated with this object
    //! @return the PDB filename
    const std::string &FragmentStochasticPoseOptimizer::GetPocketFilename() const
    {
      return m_BindingPocketFilename;
    }

    //! @brief return the protein pocket PDB associated with this object
    //! @return the protein pocket object
    const FragmentComplete &FragmentStochasticPoseOptimizer::GetPocket() const
    {
      return m_BindingPocket;
    }

    //! brief return the bfactors associated with this object
    //! return the bfactors for each atom in the pocket
    const storage::Vector< float> &FragmentStochasticPoseOptimizer::GetBFactors() const
    {
      return m_BFactors;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentStochasticPoseOptimizer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentStochasticPoseOptimizer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
