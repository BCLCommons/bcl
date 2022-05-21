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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_MCM_H_

#define BCL_CHEMISTRY_FRAGMENT_MUTATE_MCM_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_constitution_set.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "bcl_chemistry_search_fragment_library_from_tree.h"
#include "find/bcl_find_pick_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "mc/bcl_mc_approximator.h"
#include "opti/bcl_opti_tracker.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentMutateMCM
    //! @brief Internal MCM engine  to allow selective optimization of a molecular substructure/fragment prior
    //! to returning it to the main Monte Carlo - Metropolis engine in the FocusedLibraryDesign application
    //!
    //! @see @link example_chemistry_fragment_mutate_mcm.cpp @endlink
    //! @author brownbp1
    //! @date April 30, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateMCM :
      public FragmentMutateInterface
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! MCM optimization goal
      opti::Tracker< FragmentComplete, double> m_OptiGoal;

      //! fragment split
      util::Implementation< FragmentSplitInterface> m_Splitter;

      //! rotamer library information
      util::ShPtr< SearchFragmentLibraryFromTree> m_RotamerLibrarySearcher;

      //! pool of fragments to be picked from
      util::ShPtr< FragmentEnsemble> m_FragmentPool;

      //! flags for the internal Monte Carlo-Metropolis
      size_t m_MaxSequentialMutates;
      float m_RingSwapProb;
      float m_CyclizeProb;
      float m_AlchemyProb;
      float m_RemoveAtomProb;
      float m_RemoveBondProb;
      float m_AddMedChemProb;
      float m_FluorinateProb;
      float m_HalogenateProb;
      float m_ExtendWithLinkerProb;

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
      FragmentMutateMCM();

      //! @brief constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      FragmentMutateMCM
      (
        const opti::Tracker< FragmentComplete, double> &OPTI_GOAL,
        const util::Implementation< FragmentSplitInterface> &SPLITTER,
        const util::ShPtr< SearchFragmentLibraryFromTree> &TREE_SEARCH,
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &CORINA_CONFS
      );

      //! @brief constructor with MCM options
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      FragmentMutateMCM
      (
        const opti::Tracker< FragmentComplete, double> &OPTI_GOAL,
        const util::Implementation< FragmentSplitInterface> &SPLITTER,
        const util::ShPtr< SearchFragmentLibraryFromTree> &TREE_SEARCH,
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &CORINA_CONFS,
        const size_t &MAX_SEQUENTIAL_MUTATES,
        const float &RING_SWAP_PROB,
        const float &CYCLIZE_PROB,
        const float &ALCHEMY_PROB,
        const float &REMOVE_ATOM_PROB,
        const float &REMOVE_BOND_PROB,
        const float &ADD_MEDCHEM_PROB,
        const float &FLUORINATE_PROB,
        const float &HALOGENATE_PROB,
        const float &EXTEND_WITH_LINKER_PROB
      );

      //! @brief local mutate pose-sensitive constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
      //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMutateMCM
      (
        const opti::Tracker< FragmentComplete, double> &OPTI_GOAL,
        const util::Implementation< FragmentSplitInterface> &SPLITTER,
        const util::ShPtr< SearchFragmentLibraryFromTree> &TREE_SEARCH,
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA_CONFS
      );

      //! @brief local mutate pose-sensitive constructor with MCM options
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
      //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMutateMCM
      (
        const opti::Tracker< FragmentComplete, double> &OPTI_GOAL,
        const util::Implementation< FragmentSplitInterface> &SPLITTER,
        const util::ShPtr< SearchFragmentLibraryFromTree> &TREE_SEARCH,
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA_CONFS,
        const size_t &MAX_SEQUENTIAL_MUTATES,
        const float &RING_SWAP_PROB,
        const float &CYCLIZE_PROB,
        const float &ALCHEMY_PROB,
        const float &REMOVE_ATOM_PROB,
        const float &REMOVE_BOND_PROB,
        const float &ADD_MEDCHEM_PROB,
        const float &FLUORINATE_PROB,
        const float &HALOGENATE_PROB,
        const float &EXTEND_WITH_LINKER_PROB
      );

      //! @brief clone constructor
      FragmentMutateMCM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an SmallMolecule and returning a grown SmallMolecule
      //! @param FRAGMENT small molecule of interest
      //! @return Constitution after the mutate
      math::MutateResult< FragmentComplete> operator()( const FragmentComplete &FRAGMENT) const;

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief set up mutate object
      util::ShPtr< math::MutateInterface< FragmentComplete> > SetupMutate
      (
        const FragmentComplete &START_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &RESOLVE_CLASHES,
        const bool &CORINA_CONFS,
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const size_t &MAX_SEQUENTIAL_MUTATES,
        const float &RING_SWAP_PROB,
        const float &CYCLIZE_PROB,
        const float &ALCHEMY_PROB,
        const float &REMOVE_ATOM_PROB,
        const float &REMOVE_BOND_PROB,
        const float &ADD_MEDCHEM_PROB,
        const float &FLUORINATE_PROB,
        const float &HALOGENATE_PROB,
        const float &EXTEND_WITH_LINKER_PROB
      ) const;

      //! @brief set up score object
      util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > SetupScore
      (
        const descriptor::CheminfoProperty &PROPERTY
      ) const;

      //! @brief set up an MCM approximator object
      //! @return returns an approximator
      mc::Approximator< FragmentComplete, double> SetupMCM
      (
        const opti::Tracker< FragmentComplete, double> &OPTI_GOAL,
        const FragmentComplete &START_MOL,
        const size_t &ITERATIONS,
        const size_t &MAX_UNIMPROVED,
        const size_t &MAX_SKIPPED,
        const float &TEMP_ACCEPT_START,
        const float &TEMP_ACCEPT_END,
        const util::ShPtr< math::MutateInterface< FragmentComplete> > &MUTATE,
        const util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > &SCORE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class FragmentMutateMCM

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_MCM_H_
