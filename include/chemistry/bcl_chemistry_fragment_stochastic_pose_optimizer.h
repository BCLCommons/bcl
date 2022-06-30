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

#ifndef BCL_CHEMISTRY_FRAGMENT_STOCHASTIC_POSE_OPTIMIZER_H_
#define BCL_CHEMISTRY_FRAGMENT_STOCHASTIC_POSE_OPTIMIZER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "opti/bcl_opti_tracker.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentStochasticPoseOptimizer
    //! @brief Performs random perturbations (rotation, translation, conformation change) to a small molecule
    //! ligand in a protein binding pocket and scores the interaction with a pose-dependent deep neural network score
    //! function. If a perturbation improves the score then the move is accepted, otherwise it is rejected.
    //!
    //! RunMCMDOcker trackers the best pose from each MCMDocker MCM run with a set
    //! MCMDocker does the following:
    //! 1. Choose a random conformer (from the global ensemble)
    //! 2. BigTransform (max translate 3.0 angstroms, max rotate 360 degrees)
    //! 3. MCM trials (transform with rot <= 5.0, trans <= 0.1, local confswap)
    //! 4. Return molecule and score to RunMCMDocker
    //! RunMCMDocker saves each outcome as a sorted ConformationSet
    //!
    //! @see @link example_chemistry_fragment_stochastic_pose_optimizer.cpp @endlink
    //! @author brownbp1
    //! @date May 18, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentStochasticPoseOptimizer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! scoring property to be used during the Clean phase
      descriptor::CheminfoProperty m_PropertyScorer;

      //! per-residue flexibility (lower numbers less rigid, higher numbers more rigid)
      storage::Vector< float> m_BFactors;

      //! protein or protein binding pocket for clash resolution and pose-dependent scoring
      FragmentComplete m_BindingPocket;

      //! protein pocket voxel grid
      VoxelGridAtom m_VoxelGridPocket;

      //! MDL property label specifying path to protein binding pocket
      std::string m_MDL;

      //! filename for protein or protein binding pocket
      std::string m_BindingPocketFilename;

      //! optimization iterations
      size_t m_Iterations;

      //! optimization attempts
      size_t m_Attempts;

      //! refinement iterations
      size_t m_RefinementIterations;

      //! metropolis temperature
      float m_Temperature;

      //! improvement type
      opti::Tracker< FragmentComplete, double> m_Opti;

      //! VDW threshold above which conformers will be excluded
      float m_VDWClashCutoff;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentStochasticPoseOptimizer();

//      //! @brief constructor with binding pocket object
//      //! @param PROPERTY_SCORER enables pose-dependent optimization of score
//      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
//      //! @param BINDING_POCKET the protein binding pocket object
//      FragmentStochasticPoseOptimizer
//      (
//        const descriptor::CheminfoProperty &PROPERTY_SCORER,
//        const storage::Vector< float> &BFACTORS,
//        const FragmentComplete &BINDING_POCKET,
//        const size_t &ITERATIONS
//      );

      //! @brief constructor with binding pocket to be created from MDL property
      //! @param PROPERTY_SCORER enables pose-dependent optimization of score
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      //! @param MDL the SDF file MDL property specifying the binding pocket filename
      //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
      FragmentStochasticPoseOptimizer
      (
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const storage::Vector< float> &BFACTORS,
        const std::string &MDL,
        const std::string &BINDING_POCKET_FILENAME,
        const size_t &ITERATIONS,
        const size_t &ATTEMPTS,
        const size_t &REFINEMENT_ITERATIONS,
        const float &TEMPERATURE = 1.0,
        const opti::Tracker< FragmentComplete, double> &OPTI = opti::e_SmallerIsBetter,
        const float &VDW_CLASH_CUTOFF = 5.0
      );

      //! @brief clone constructor
      FragmentStochasticPoseOptimizer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief chooses the first conformer in an ensemble not clashed with the protein
      //! @param ENSEMBLE the ensemble of molecule conformers
      util::ShPtr< FragmentEnsemble> ResolveClashes
      (
        const FragmentEnsemble &ENSEMBLE
      );

      //! @brief optimize pose-dependent score of ligand with protein pocket
      //! @param ENSEMBLE the ensemble of molecule conformers
      util::ShPtr< FragmentComplete> OptimizePoseConf
      (
        const FragmentEnsemble &ENSEMBLE
      );

      //! @brief optimize pose-dependent score of ligand with protein pocket
      //! @param ENSEMBLE the ensemble of molecule conformers
      util::ShPtr< FragmentComplete> GetBestScoringNonClashingPose
      (
        const FragmentEnsemble &ENSEMBLE
      );

      //! @brief optimize pose-dependent score of ligand with protein pocket
      //! @param ENSEMBLE the ensemble of molecule conformers
      util::ShPtr< FragmentComplete> StochasticPoseOptimization
      (
        const FragmentEnsemble &ENSEMBLE
      );

      //! @brief MCM sampling of ligand inside of protein pocket
      //! @param ENSEMBLE the ensemble of molecule conformers
      //! @param LOCAL disables initial large perturbation
      util::ShPtr< FragmentComplete> MCMDocker
      (
        const FragmentEnsemble &ENSEMBLE,
        const bool &LOCAL = false
      );

      //! @brief MCM sampling of ligand inside of protein pocket
      //! @param ENSEMBLE the ensemble of molecule conformers
      //! @param LOCAL disables initial large perturbation
      //! @param REFINE enables local refinement of the best poses
      util::ShPtr< FragmentComplete> RunMCMDocker
      (
        const FragmentEnsemble &ENSEMBLE,
        const bool &LOCAL = false,
        const bool &REFINE = true
      );

      //! @brief Checks for clashes between protein and small molecule
      bool ProteinLigandClash
      (
        const FragmentComplete &MOLECULE
      );

      //! @brief Check for a bad ligand conformer VDW score
      //! @param MOLECULE the molecule conformer of interest
      //! @return returns true if the conformer passes (is good), false otherwise
      bool GoodLigandVDWScore
      (
        const FragmentComplete &MOLECULE
      );

    ////////////////////////
    // mutation functions //
    ////////////////////////

      //! @brief performs a minor perturbation of the molecule via rotation or translation, then score
      //! @param MOLECULE the molecule whose pose is to be perturbed
      //! @param CONF_ENSEMBLE optionally a conformer ensemble of the molecule to be perturbed
      //! @return returns the perturbed molecule
      util::ShPtr< FragmentComplete> PerturbMoleculePoseFiveDegreesMax
      (
        const FragmentComplete &MOLECULE,
        const FragmentEnsemble &CONF_ENSEMBLE = FragmentEnsemble()
      );

      //! @brief takes an ensemble and returns a random conformer
      util::ShPtr< FragmentComplete> ConformerSwap
      (
        const FragmentEnsemble &ENSEMBLE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

      //! @brief return the MDL property associated with this object
      //! @return the MDL string
      const std::string &GetMDL() const;

      //! @brief return the pocket filename associated with this object
      //! @return the PDB filename
      const std::string &GetPocketFilename() const;

      //! @brief return the protein pocket PDB associated with this object
      //! @return the protein pocket object
      const FragmentComplete &GetPocket() const;

      //! brief return the bfactors associated with this object
      //! return the bfactors for each atom in the pocket
      const storage::Vector< float> &GetBFactors() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentStochasticPoseOptimizer

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_STOCHASTIC_POSE_OPTIMIZER_H_
