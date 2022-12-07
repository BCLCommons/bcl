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
#include "chemistry/bcl_chemistry_fragment_mutate_mcm.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_mutate_add_med_chem.h"
#include "chemistry/bcl_chemistry_fragment_mutate_alchemy.h"
#include "chemistry/bcl_chemistry_fragment_mutate_cyclize.h"
#include "chemistry/bcl_chemistry_fragment_mutate_extend_with_linker.h"
#include "chemistry/bcl_chemistry_fragment_mutate_fluorinate.h"
#include "chemistry/bcl_chemistry_fragment_mutate_halogenate.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_atom.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"
#include "chemistry/bcl_chemistry_fragment_mutate_ring_swap.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_score_function_generic.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_base.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "math/bcl_math_const_function.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_repeat.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_convergence_result.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_skipped_steps.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "random/bcl_random_uniform_distribution.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateMCM::FragmentMutateMCM() :
        m_OptiGoal( opti::e_SmallerIsBetter),
        m_Splitter( util::Implementation< FragmentSplitInterface>()),
        m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
        m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
        m_MaxSequentialMutates( 1),
        m_RingSwapProb( 0.1),
        m_CyclizeProb( 0.1),
        m_AlchemyProb( 0.1),
        m_RemoveAtomProb( 0.1),
        m_RemoveBondProb( 0.1),
        m_AddMedChemProb( 0.1),
        m_FluorinateProb( 0.1),
        m_HalogenateProb( 0.1),
        m_ExtendWithLinkerProb( 0.1)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateMCM::FragmentMutateMCM
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
    ) :
      m_OptiGoal( OPTI_GOAL),
      m_Splitter( SPLITTER),
      m_RotamerLibrarySearcher( TREE_SEARCH),
      m_FragmentPool( FRAGMENT_POOL),
      m_MaxSequentialMutates( 1),
      m_RingSwapProb( 0.1),
      m_CyclizeProb( 0.1),
      m_AlchemyProb( 0.1),
      m_RemoveAtomProb( 0.1),
      m_RemoveBondProb( 0.1),
      m_AddMedChemProb( 0.1),
      m_FluorinateProb( 0.1),
      m_HalogenateProb( 0.1),
      m_ExtendWithLinkerProb( 0.1)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_PropertyScorer = PROPERTY_SCORER;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor with MCM options
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateMCM::FragmentMutateMCM
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
    ) :
      m_OptiGoal( OPTI_GOAL),
      m_Splitter( SPLITTER),
      m_RotamerLibrarySearcher( TREE_SEARCH),
      m_FragmentPool( FRAGMENT_POOL),
      m_MaxSequentialMutates( MAX_SEQUENTIAL_MUTATES),
      m_RingSwapProb( RING_SWAP_PROB),
      m_CyclizeProb( CYCLIZE_PROB),
      m_AlchemyProb( ALCHEMY_PROB),
      m_RemoveAtomProb( REMOVE_ATOM_PROB),
      m_RemoveBondProb( REMOVE_BOND_PROB),
      m_AddMedChemProb( ADD_MEDCHEM_PROB),
      m_FluorinateProb( FLUORINATE_PROB),
      m_HalogenateProb( HALOGENATE_PROB),
      m_ExtendWithLinkerProb( EXTEND_WITH_LINKER_PROB)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_PropertyScorer = PROPERTY_SCORER;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate pose-sensitive constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateMCM::FragmentMutateMCM
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
    ) :
      m_OptiGoal( OPTI_GOAL),
      m_Splitter( SPLITTER),
      m_RotamerLibrarySearcher( TREE_SEARCH),
      m_FragmentPool( FRAGMENT_POOL),
      m_MaxSequentialMutates( 1),
      m_RingSwapProb( 0.1),
      m_CyclizeProb( 0.1),
      m_AlchemyProb( 0.1),
      m_RemoveAtomProb( 0.1),
      m_RemoveBondProb( 0.1),
      m_AddMedChemProb( 0.1),
      m_FluorinateProb( 0.1),
      m_HalogenateProb( 0.1),
      m_ExtendWithLinkerProb( 0.1)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_PropertyScorer = PROPERTY_SCORER;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate pose-sensitive constructor with MCM options
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateMCM::FragmentMutateMCM
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
    ) :
      m_OptiGoal( OPTI_GOAL),
      m_Splitter( SPLITTER),
      m_RotamerLibrarySearcher( TREE_SEARCH),
      m_FragmentPool( FRAGMENT_POOL),
      m_MaxSequentialMutates( MAX_SEQUENTIAL_MUTATES),
      m_RingSwapProb( RING_SWAP_PROB),
      m_CyclizeProb( CYCLIZE_PROB),
      m_AlchemyProb( ALCHEMY_PROB),
      m_RemoveAtomProb( REMOVE_ATOM_PROB),
      m_RemoveBondProb( REMOVE_BOND_PROB),
      m_AddMedChemProb( ADD_MEDCHEM_PROB),
      m_FluorinateProb( FLUORINATE_PROB),
      m_HalogenateProb( HALOGENATE_PROB),
      m_ExtendWithLinkerProb( EXTEND_WITH_LINKER_PROB)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_PropertyScorer = PROPERTY_SCORER;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief clone constructor
    FragmentMutateMCM *FragmentMutateMCM::Clone() const
    {
      return new FragmentMutateMCM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateMCM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateMCM::GetAlias() const
    {
      static const std::string s_name( "FragmentMutateMCM");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateMCM::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "FragmentMutateMCM!");

      // get the starting molecule minus the mutable region
      util::ShPtr< FragmentComplete> scaffold_fragment( new FragmentComplete());
      static FragmentTrackMutableAtoms atom_tracker;
      if( m_MutableAtomIndices.GetSize())
      {
        scaffold_fragment =
            util::ShPtr< FragmentComplete>(
              new FragmentComplete
              (
                atom_tracker.GetBaseFragment
                (
                  FRAGMENT,
                  FragmentComplete(),
                  m_MutableAtomIndices
                )
              )
            );
        BCL_Assert( scaffold_fragment->GetSize(), "Invalid atom indices specified");
      }
      else
      {
        scaffold_fragment = util::ShPtr< FragmentComplete>(
          new FragmentComplete(
            atom_tracker.GetRandomBaseFragment
            (
              FRAGMENT,
              m_Splitter
            )
          )
        );
      }

      // setup MCM mutate
      util::ShPtr< math::MutateInterface< FragmentComplete> > mutate
      (
        SetupMutate
        (
          FRAGMENT,
          FragmentEnsemble( storage::List< FragmentComplete>( 1, *scaffold_fragment)),
          m_MutableAtomIndices,
          m_MDL,
          m_PropertyScorer,
          m_ResolveClashes,
          m_Corina,
          m_FragmentPool,
          m_MaxSequentialMutates,
          m_RingSwapProb,
          m_CyclizeProb,
          m_AlchemyProb,
          m_RemoveAtomProb,
          m_RemoveBondProb,
          m_AddMedChemProb,
          m_FluorinateProb,
          m_HalogenateProb,
          m_ExtendWithLinkerProb
        )
      );

      // setup MCM score
      util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > score
      (
        SetupScore
        (
          m_PropertyScorer
        )
      );

      // setup MXM approximator
      mc::Approximator< FragmentComplete, double> appx
      (
        SetupMCM
        (
          m_OptiGoal,
          FRAGMENT,
          100, // iterations
          100, // consecutive unimproved
          10, // consecutive skipped
          0.90, // starting accept fraction
          0.10, // ending accept fraction
          mutate,
          score
        )
      );

      // run approximator
      appx.Approximate();

      // done
      if( appx.GetTracker().GetBest().IsDefined())
      {
        BCL_MessageStd( "Best this internal round: " + util::Format()( appx.GetTracker().GetBest()->Second()));
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>( new FragmentComplete( appx.GetTracker().GetBest()->First())), *this);
      }
      else
      {
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set up mutate object
    util::ShPtr< math::MutateInterface< FragmentComplete> > FragmentMutateMCM::SetupMutate
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
    ) const
    {

      // tree search for RingSwap
      util::ShPtr< SearchFragmentLibraryFromTree> tree_search
      (
        new SearchFragmentLibraryFromTree
        (
          *util::Implementation< RotamerLibraryInterface>( RotamerLibraryInterface::GetDefault())
        )
      );

      util::ShPtr< math::MutateDecisionNode< FragmentComplete> > mutater
      (
        new math::MutateDecisionNode< FragmentComplete>()
      );

      // POSE-DEPENDENT CONSTRUCTION OF MUTATES //
      if( !MDL.empty())
      {
        mutater->AddMutate( FragmentMutateRingSwap( tree_search, m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS, true, false, 0.1, true, true), RING_SWAP_PROB);
        mutater->AddMutate( FragmentMutateCyclize( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), CYCLIZE_PROB);
        mutater->AddMutate( FragmentMutateAlchemy( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), ALCHEMY_PROB);
        mutater->AddMutate( FragmentMutateRemoveAtom( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), REMOVE_ATOM_PROB);
        mutater->AddMutate( FragmentMutateRemoveBond( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), REMOVE_BOND_PROB);
        mutater->AddMutate( FragmentMutateExtendWithLinker( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), EXTEND_WITH_LINKER_PROB);
        mutater->AddMutate( FragmentMutateAddMedChem( FRAGMENT_POOL, m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), ADD_MEDCHEM_PROB);
        mutater->AddMutate( FragmentMutateFluorinate( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), FLUORINATE_PROB);
        mutater->AddMutate( FragmentMutateHalogenate( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, MDL, PROPERTY_SCORER, RESOLVE_CLASHES, storage::Vector< float>(), CORINA_CONFS), HALOGENATE_PROB);
      }
      // POSE-INDEPENDENT CONSTRUCTION OF MUTATES //
      else
      {
        mutater->AddMutate( FragmentMutateRingSwap( tree_search, m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS, true, false, 0.1, true, true), RING_SWAP_PROB);
        mutater->AddMutate( FragmentMutateCyclize( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), CYCLIZE_PROB);
        mutater->AddMutate( FragmentMutateAlchemy(  m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), ALCHEMY_PROB);
        mutater->AddMutate( FragmentMutateRemoveAtom( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), REMOVE_ATOM_PROB);
        mutater->AddMutate( FragmentMutateRemoveBond( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), REMOVE_BOND_PROB);
        mutater->AddMutate( FragmentMutateExtendWithLinker( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), EXTEND_WITH_LINKER_PROB);
        mutater->AddMutate( FragmentMutateAddMedChem( FRAGMENT_POOL, m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), ADD_MEDCHEM_PROB);
        mutater->AddMutate( FragmentMutateFluorinate( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), FLUORINATE_PROB);
        mutater->AddMutate( FragmentMutateHalogenate( m_DrugLikenessType, START_FRAGMENT, MUTABLE_FRAGMENTS, MUTABLE_ATOM_INDICES, CORINA_CONFS), HALOGENATE_PROB);
      }

      // set up sequential mutate to perform 1 to N mutates in a row prior to scoring (does not bypass druglikeness filtering)
      util::ShPtr< math::MutateInterface< FragmentComplete> > mutate_repeater
      (
        new math::MutateRepeat< FragmentComplete>
        (
          mutater,
          1,
          MAX_SEQUENTIAL_MUTATES
        )
      );
      return mutate_repeater;
    }

    //! @brief set up score object
    util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > FragmentMutateMCM::SetupScore
    (
      const descriptor::CheminfoProperty &PROPERTY
    ) const
    {
      util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> > scorer =
          util::ShPtr< math::FunctionInterfaceSerializable< FragmentComplete, double> >
      (
        new ScoreFunctionGeneric
        (
          PROPERTY
        )
      );
      return scorer;
    }

    //! @brief set up an MCM approximator object
    //! @return returns an approximator
    mc::Approximator< FragmentComplete, double> FragmentMutateMCM::SetupMCM
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
    ) const
    {
       // create the temperature control
       util::ShPtr< mc::TemperatureInterface> sp_temperature
       (
         new mc::TemperatureAccepted
         (
           TEMP_ACCEPT_START,
           TEMP_ACCEPT_END,
           ITERATIONS,
           1.0,
           std::max< size_t>( 1, 100 / 10)
         )
       );

       // create the metropolis
       mc::Metropolis< double> metropolis( sp_temperature, true);

       // create the termination criterion
       opti::CriterionCombine< FragmentComplete, double> criterion_combine;

       // insert termination criteria that depends on the total number of MC iterations
       opti::CriterionNumberIterations< FragmentComplete, double> maximum_number_iterations( ITERATIONS);
       criterion_combine.InsertCriteria( maximum_number_iterations);

       // insert termination criteria that depends on the total number of unimproved MC iterations
       opti::CriterionUnimproved< FragmentComplete, double> maximum_number_unimproved_iterations( MAX_UNIMPROVED);
       criterion_combine.InsertCriteria( maximum_number_unimproved_iterations);

       // insert termination criteria that depends on the total number of skipped MC iterations
       opti::CriterionSkippedSteps< FragmentComplete, double> maximum_number_skipped_iterations( MAX_SKIPPED);
       criterion_combine.InsertCriteria( maximum_number_skipped_iterations);

       // insert termination criteria that depends on convergence tolerance
//       opti::CriterionConvergenceResult< chemistry::FragmentComplete, double> convergence_result( ITERATIONS * 0.10, 0.10);
//       criterion_combine.InsertCriteria( convergence_result);

       mc::Approximator< FragmentComplete, double> approximator
       (
         SCORE,
         MUTATE,
         metropolis,
         criterion_combine,
         START_MOL,
         OPTI_GOAL
       );

       return approximator;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateMCM::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Monte Carlo - Metropolis wrapper for fragment mutate objects"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateMCM::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // call RISH function of the base class
      if( !FragmentMutateInterface::ReadInitializerSuccessHook( LABEL, ERROR_STREAM))
      {
        return false;
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
