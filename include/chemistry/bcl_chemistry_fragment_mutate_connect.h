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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_CONNECT_H_
#define BCL_CHEMISTRY_FRAGMENT_MUTATE_CONNECT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_map_conformer.h"
#include "bcl_chemistry_fragment_mutate_extend_with_linker.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
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
    //! @class FragmentMutateConnect
    //! @brief Used to connect disconnected atoms in fragment(s) with complex linkers
    //!
    //! @see @link example_chemistry_fragment_mutate_connect.cpp @endlink
    //! @author brownbp1
    //! @date May 19, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateConnect :
      public FragmentMutateInterface
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! blah
      util::ShPtr< SearchFragmentLibraryFromTree> m_RotamerLibrarySearcher;

      //! ordered list of strings corresponding to types of mutates with FragmentExtendWithLinker
      storage::Vector< std::string> m_Linkers;
      std::string m_LinkersString;

      //! fragment to which the input molecule will be connected
      FragmentComplete m_TerminalFragment;
      std::string m_TerminalFragmentFilename;

      //////// Begin terminal fragment atom selection options ////////

      //! atom indices that can be mutated
      storage::Vector< size_t> m_TerminalFragmentMutableAtomIndices;
      std::string m_TerminalFragmentMutableAtoms;

      //! explicit mutable subgraph of the current molecule
      FragmentEnsemble m_TerminalFragmentMutableFragments;
      std::string m_TerminalFragmentMutableFragmentsFilename;

      //! element types that can be mutated
      storage::Vector< ElementType> m_TerminalFragmentMutableElements;
      std::string m_TerminalFragmentMutableElementsString;

      //! atom indices that can be mutated
      storage::Vector< size_t> m_TerminalFragmentFixedAtomindices;
      std::string m_TerminalFragmentFixedAtoms;

      //! explicit fixed subgraph of the current molecule
      FragmentEnsemble m_TerminalFragmentFixedFragments;
      std::string m_TerminalFragmentFixedFragmentsFilename;

      //! element types that can be mutated
      storage::Vector< ElementType> m_TerminalFragmentFixedElements;
      std::string m_TerminalFragmentFixedElementsString;

      //! options to invert (find complement) the mutable and fixed fragment subgraphs with target molecule
      bool m_TerminalFragmentComplementMutableFragments;
      bool m_TerminalFragmentComplementFixedFragments;

      //! Scheme used for comparing whether two bonds are equivalent during mutable and fixed atom selection
      ConfigurationalBondTypeData::DataEnum m_TerminalFragmentMutableBondComparisonType;
      ConfigurationalBondTypeData::DataEnum m_TerminalFragmentFixedBondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent during mutable and fixed atom selection
      ConformationGraphConverter::AtomComparisonTypeEnum m_TerminalFragmentMutableAtomComparisonType;
      ConformationGraphConverter::AtomComparisonTypeEnum m_TerminalFragmentFixedAtomComparisonType;

      //////// End terminal fragment atom selection options ////////

      //! minimum distance cutoff between half-linker terminal atoms determining if join will occur
      float m_JoinDistanceCutoffMin;

      //! maximum distance cutoff between half-linker terminal atoms determining if join will occur
      float m_JoinDistanceCutoffMax;

      //! save an ensemble to this filename
      std::string m_OutputEnsembleFilename;

      //! starting fragmentRMSD cutoff below which conformers can be saved to output ensemble if desired
      float m_StartingFragmentRMSDCutoff;

      //! mutate to extend linker
      mutable FragmentMutateExtendWithLinker m_ExtendWithLinker;

      // Note on SampleConformations options -
      // some SampleConfs options are made accessible through member data
      // getters/setters and serializer because they may be useful to change
      // when experimenting with different linkers.
      // The global options refer to the ensembles that we generate of the
      // half-linker subgraph of the fragment_a+halflinker_a and
      // fragment_b+halflinker_b molecules. The local options refer to the
      // ensembles that we generate of the fully joined linker. The global
      // and local designations refer to the SetSamplingPreferences flag
      // in SampleConformations.

      //! conformation comparer for global ensembles
      std::string m_ConfComparerGlobal;

      //! unique conformer tolerance for global ensembles
      double m_ConfToleranceGlobal;

      //! maximum number of conformers to make for global ensembles
      size_t m_NMaxConfsGlobal;

      //! maximum number of iterations for global ensembles
      size_t m_NMaxItersGlobal;

      //! if performing clustering on global ensembles
      bool m_ClusterGlobal;

      //! clash resolution for global ensembles
      double m_ClashResolutionGlobal;

      //! clash score for global ensembles
      double m_ClashScoreGlobal;

      //! conformation comparer for local ensembles
      std::string m_ConfComparerLocal;

      //! unique conformer tolerance for local ensembles
      double m_ConfToleranceLocal;

      //! maximum number of conformers to make for local ensembles
      size_t m_NMaxConfsLocal;

      //! maximum number of iterations for local ensembles
      size_t m_NMaxItersLocal;

      //! if performing clustering on local ensembles
      bool m_ClusterLocal;

      //! clash resolution for local ensembles
      double m_ClashResolutionLocal;

      //! clash score for local ensembles
      double m_ClashScoreLocal;

      //! rings from fragment database
      util::ShPtr< FragmentEnsemble> m_Rings;
      std::string m_RingsFilename;

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
      FragmentMutateConnect();

      //! @brief druglikeness constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      FragmentMutateConnect
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const bool &CORINA_CONFS
      );

      //! @brief full constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      FragmentMutateConnect
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const bool &CORINA_CONFS
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
      FragmentMutateConnect
      (
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

      //! @brief local clash resolver constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMutateConnect
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA_CONFS
      );

      //! @brief clone constructor
      FragmentMutateConnect *Clone() const;

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

      //! @brief set options for the global conformational ensembles made of each fragment + half linker
      void InitializeSampleConfsGlobal( SampleConformations &SAMPLER) const;

      //! @brief set options for the local conformational ensembles made of the final joined molecule
      void InitializeSampleConfsLocal( SampleConformations &SAMPLER) const;

      //! @brief return an ensemble of conformers for the starting fragment extended with half the linker
      //! @param FRAGMENT_A molecule from which extension will start
      //! @param LINKER_COMPOSITION vector of strings indicating type of linker to add
      //! @param LINK_INDEX_A atom in FRAGMENT_A at which to begin extension
      //! @return an ensemble of conformers for the first half of the extended linker (for odd numbers, half+1),
      //! the index of the link point to the next half of the molecule (in hydrogenated numbering), and
      //! the samplebyparts indices
      storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> > ExtendHalfLinkerForward
      (
        const FragmentComplete FRAGMENT_A,
        const storage::Vector< std::string> &LINKER_COMPOSITION,
        const storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > &RINGS,
        const size_t LINK_AINDEX_A,
        const FragmentMapConformer &CLEANER,
        SampleConformations &CONFORMATOR
      ) const;

      //! @brief return an ensemble of conformers for the terminal fragment extended with half the linker
      //! @param FRAGMENT_B terminal fragment from which extension will start
      //! @param LINKER_COMPOSITION vector of strings indicating type of linker to add
      //! @param LINK_INDEX_B atom in FRAGMENT_B at which to begin extension
      //! @return an ensemble of conformers for the terminal half of the extended linker (for odd numbers, half-1)
      //! the index of the link point to the other half of the molecule (in hydrogenated numbering),
      //! and the samplebyparts indices
      storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> > ExtendHalfLinkerReverse
      (
        const FragmentComplete &FRAGMENT_B,
        const storage::Vector< std::string> &LINKER_COMPOSITION,
        const storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > &RINGS,
        const size_t LINK_AINDEX_B,
        const FragmentMapConformer &CLEANER,
        SampleConformations &CONFORMATOR
      ) const;

      //! @brief propose solutions for the fully linked molecule given two extended halves
      //! @param FWD_ENS fragment that was extended with half of the linker in the forward direction
      //! @param FWD_LINK_ATOM the atom at which this fragment will be joined
      //! @param FWD_SBP_INDICES the atoms that compose conformationally flexible dihedrals
      //! @param REV_ENS fragment that was extended with half of the linker in the reverse direction
      //! @param REV_LINK_ATOM the atom at which this fragment will be joined
      //! @param REV_SBP_INDICES the atoms that compose conformationally flexible dihedrals
      //! @return a conformational ensemble of the fully linked molecule and the final linker indices
      storage::Pair< FragmentEnsemble, storage::Vector< size_t> > JoinHalfExtendedFragments
      (
        const FragmentEnsemble &FWD_ENS,
        const size_t FWD_LINK_ATOM,
        const storage::Vector< size_t> &FWD_SBP_INDICES,
        const FragmentEnsemble &REV_ENS,
        const size_t REV_LINK_ATOM,
        const storage::Vector< size_t> &REV_SBP_INDICES,
        const RotamerLibraryFile &ROTLIB
      ) const;

      //! @brief select rings from the library and identify attachment atoms
      //! @params LINKER_COMPONENTS the components to check for rings
      //! @return a collection of rings and corresponding attachment indices
      storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > ChooseRings
      (
        const storage::Vector< std::string> &LINKER_COMPONENTS
      ) const;

      //! @brief modify input molecule by adding a link fragment
      //! @param FRAGMENT the fragment to be modified
      //! @param INDEX the index of FRAGMENT to which the linker fragment will be appended
      //! @param LINK_TYPE the string indicating the link type
      FragmentComplete AddLinkFragment
      (
        const FragmentComplete FRAGMENT,
        const size_t INDEX,
        const std::string &LINK_TYPE,
        const FragmentComplete &RING = FragmentComplete(),
        const size_t ATTACH_INDEX = util::GetUndefinedSize_t()
      ) const;

      //! @brief create supermolecule of the two fragments without a linker
      FragmentComplete CreateSuperMolecule
      (
        const FragmentComplete &FRAGMENT_A,
        const FragmentComplete &FRAGMENT_B
      ) const;

      //! @brief compute RMSD of linked molecule to the unlinked supermolecule of the two fragments
      double ComputeRMSDToSuperMolecule
      (
        const FragmentComplete &MOLECULE,
        const storage::Vector< size_t> &COMMON_INDICES,
        const FragmentComplete &SUPERMOLECULE
      ) const;

      //! @brief filter out linker pairs whose link atoms are too far apart
      //! @param FRAGMENT_A first fragment to be linked
      //! @param DIST_ATOM_A atom at which FRAGMENT_A will be linked to FRAGMENT_B
      //! @param FRAGMENT_B second fragment to be linked
      //! @param DIST_ATOM_B atom at which FRAGMENT_B will be linked to FRAGMENT_A
      //! @param DIST_THRESHOLD maximum allowed distance between DIST_ATOM_A and DIST_ATOM_B
      //! @return true if equal or below DIST_THRESHOLD, false otherwise
      bool DistanceFilter
      (
        const FragmentComplete &FRAGMENT_A,
        const size_t DIST_ATOM_A,
        const FragmentComplete &FRAGMENT_B,
        const size_t DIST_ATOM_B,
        const float MIN_DIST_THRESHOLD,
        const float MAX_DIST_THRESHOLD
      ) const;

      //! @brief try to remove clashes from the linker after joining the two halves
      // TODO deprecate with minimizer
      void ResolveClashes
      (
        FragmentComplete &MOLECULE,
        const storage::Vector< size_t> &LINKER_ATOM_INDICES,
        const RotamerLibraryFile &ROTLIB
      ) const;

      //! @brief try to correct bond lengths from the linker after joining the two halves
      // TODO deprecate with minimizer
      void CorrectBadBondLengths
      (
        FragmentComplete &MOLECULE,
        const storage::Vector< size_t> &LINKER_ATOM_INDICES
      ) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class FragmentMutateConnect

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_CONNECT_H_
