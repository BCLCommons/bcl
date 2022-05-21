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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_INTERFACE_H_
#define BCL_CHEMISTRY_FRAGMENT_MUTATE_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "sched/bcl_sched_mutex.h"
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
    //! @class FragmentMutateInterface
    //! @brief Interface for the drug design mutates. I guess a mix of an interface and an abstract base rather than a
    //! true interface class. Most importantly, no class object of this base should ever need construction and all
    //! derived classes adhere to the format
    //!
    //! @see @link example_chemistry_fragment_mutate_interface.cpp @endlink
    //! @author brownbp1
    //! @date Jan 23, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateInterface :
      public math::MutateInterface< FragmentComplete>
    {

    protected:

    //////////
    // data //
    //////////

      //////// General miscellaneous ////////

      // TODO: replace string with a CheminfoProperty object?
      //! type of drug likeness filter to apply during molecule cleaning
      std::string m_DrugLikenessType = "None";

      //! MDL property label specifying path to protein binding pocket
      std::string m_MDL = std::string();

      //! scoring property to be used during the clean phase
      descriptor::CheminfoProperty m_PropertyScorer = descriptor::CheminfoProperty();

      //! if true, sample conformations to resolve clashes with protein following mutate
      bool m_ResolveClashes = false;

      //! per-residue flexibility (lower numbers less rigid, higher numbers more rigid)
      storage::Vector< float> m_BFactors = storage::Vector< float>();

      //! enables Corina to generate the starting conformer during cleaning
      bool m_Corina = false;

      //! max number of tries for a single mutate
      size_t m_NumberMaxAttempts = 10;

      //! shuffle hydrogen atoms bonded to a target heavy atom prior to removal for opening valences
      bool m_OVShuffleH = true;

      //! reverse the direction of valence opening such that higher index hydrogen atoms are removed first
      bool m_OVReverse = false;

      //! reference molecule for substructure-based alignment during 3D conformer construction
      FragmentComplete m_ScaffoldFragment = FragmentComplete();
      std::string m_ScaffoldFragmentFilename;

      //////// Input fragment atom selection options ////////

      //! atom indices that can be mutated
      storage::Vector< size_t> m_MutableAtomIndices = storage::Vector< size_t>();
      std::string m_MutableAtoms;

      //! explicit mutable subgraph of the current molecule
      FragmentEnsemble m_MutableFragments = FragmentEnsemble();
      std::string m_MutableFragmentsFilename;
      storage::Vector< size_t> m_MutableFragmentsAtomindices = storage::Vector< size_t>();

      //! element types that can be mutated
      storage::Vector< ElementType> m_MutableElements = storage::Vector< ElementType>();
      std::string m_MutableElementsString;
      storage::Vector< size_t> m_MutableElementsAtomindices = storage::Vector< size_t>();

      //! atom indices that can be mutated
      storage::Vector< size_t> m_FixedAtomindices = storage::Vector< size_t>();
      std::string m_FixedAtoms;

      //! explicit fixed subgraph of the current molecule
      FragmentEnsemble m_FixedFragments = FragmentEnsemble();
      std::string m_FixedFragmentsFilename;
      storage::Vector< size_t> m_FixedFragmentsAtomIndices = storage::Vector< size_t>();

      //! element types that can be mutated
      storage::Vector< ElementType> m_FixedElements = storage::Vector< ElementType>();
      std::string m_FixedElementsString;
      storage::Vector< size_t> m_FixedElementsAtomindices = storage::Vector< size_t>();

      //! options to invert (find complement) the mutable and fixed fragment subgraphs with target molecule
      bool m_ComplementMutableFragments;
      bool m_ComplementFixedFragments;

      //! Scheme used for comparing whether two bonds are equivalent during mutable and fixed atom selection
      ConfigurationalBondTypeData::DataEnum m_MutableBondComparisonType;
      ConfigurationalBondTypeData::DataEnum m_FixedBondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent during mutable and fixed atom selection
      ConformationGraphConverter::AtomComparisonTypeEnum m_MutableAtomComparisonType;
      ConformationGraphConverter::AtomComparisonTypeEnum m_FixedAtomComparisonType;

      //////// Ring library atom selection options ////////

      //! atom indices that can be mutated
      storage::Vector< size_t> m_RingLibraryMutableAtomindices = storage::Vector< size_t>();
      std::string m_RingLibraryMutableAtoms;

      //! explicit mutable subgraph of the current molecule
      FragmentEnsemble m_RingLibraryMutableFragments = FragmentEnsemble();
      std::string m_RingLibraryMutableFragmentsFilename;
      storage::Vector< size_t> m_RingLibraryMutableFragmentsAtomindices = storage::Vector< size_t>();

      //! element types that can be mutated
      storage::Vector< ElementType> m_RingLibraryMutableElements = storage::Vector< ElementType>();
      std::string m_RingLibraryMutableElementsString;
      storage::Vector< size_t> m_RingLibraryMutableElementsAtomindices = storage::Vector< size_t>();

      //! atom indices that can be mutated
      storage::Vector< size_t> m_RingLibraryFixedAtomindices = storage::Vector< size_t>();
      std::string m_RingLibraryFixedAtoms;

      //! explicit fixed subgraph of the current molecule
      FragmentEnsemble m_RingLibraryFixedFragments = FragmentEnsemble();
      std::string m_RingLibraryFixedFragmentsFilename;
      storage::Vector< size_t> m_RingLibraryFixedFragmentsAtomIndices = storage::Vector< size_t>();

      //! element types that can be mutated
      storage::Vector< ElementType> m_RingLibraryFixedElements = storage::Vector< ElementType>();
      std::string m_RingLibraryFixedElementsString;
      storage::Vector< size_t> m_RingLibraryFixedElementsAtomindices = storage::Vector< size_t>();

      //! options to invert (find complement) the mutable and fixed fragment subgraphs with target molecule
      bool m_RingLibraryComplementMutableFragments;
      bool m_RingLibraryComplementFixedFragments;

      //! Scheme used for comparing whether two bonds are equivalent during mutable and fixed atom selection
      ConfigurationalBondTypeData::DataEnum m_RingLibraryMutableBondComparisonType;
      ConfigurationalBondTypeData::DataEnum m_RingLibraryFixedBondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent during mutable and fixed atom selection
      ConformationGraphConverter::AtomComparisonTypeEnum m_RingLibraryMutableAtomComparisonType;
      ConformationGraphConverter::AtomComparisonTypeEnum m_RingLibraryFixedAtomComparisonType;

      //////// Fragment libraries for chemical design ////////

      //! static libraries for certain derived classes
      static storage::Map< std::string, util::ShPtr< FragmentEnsemble> > s_RingLibraries;
      static storage::Map< std::string, util::ShPtr< FragmentEnsemble> > s_MedChemFragmentLibraries;

      //////// Thread safety ////////

      //! mutex for read protection
      static sched::Mutex s_Mutex;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual clone constructor
      virtual FragmentMutateInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const = 0;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      virtual const std::string &GetAlias() const = 0;

//      //! @brief returns all of the mutable atoms from all atom selection methods
//      //! @return the mutable atoms
//      const storage::Vector< size_t> &GetAllMutableAtomIndices() const;

      //! @brief returns passed mutable atoms
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetMutableAtomIndices() const;

      //! @brief returns mutable elements
      //! @return the mutable elements
      const storage::Vector< ElementType> &GetMutableElements() const;

      //! @brief returns atom indices of mutable elements
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetMutableElementsAtomIndices() const;

      //! @brief returns mutable fragments
      //! @return the mutable fragments
      const FragmentEnsemble &GetMutableFragments() const;

      //! @brief returns atom indices of mutable fragments
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetMutableFragmentsAtomIndices() const;

      //! @brief returns passed mutable atoms
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetFixedAtomIndices() const;

      //! @brief returns mutable elements
      //! @return the mutable elements
      const storage::Vector< ElementType> &GetFixedElements() const;

      //! @brief returns atom indices of mutable elements
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetFixedElementsAtomIndices() const;

      //! @brief returns mutable fragments
      //! @return the mutable fragments
      const FragmentEnsemble &GetFixedFragments() const;

      //! @brief returns atom indices of mutable fragments
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetFixedFragmentsAtomIndices() const;

      //! @brief returns passed mutable atoms
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryMutableAtomIndices() const;

      //! @brief returns mutable elements
      //! @return the mutable elements
      const storage::Vector< ElementType> &GetRingLibraryMutableElements() const;

      //! @brief returns atom indices of mutable elements
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryMutableElementsAtomIndices() const;

      //! @brief returns mutable fragments
      //! @return the mutable fragments
      const FragmentEnsemble &GetRingLibraryMutableFragments() const;

      //! @brief returns atom indices of mutable fragments
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryMutableFragmentsAtomIndices() const;

      //! @brief returns passed mutable atoms
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryFixedAtomIndices() const;

      //! @brief returns mutable elements
      //! @return the mutable elements
      const storage::Vector< ElementType> &GetRingLibraryFixedElements() const;

      //! @brief returns atom indices of mutable elements
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryFixedElementsAtomIndices() const;

      //! @brief returns mutable fragments
      //! @return the mutable fragments
      const FragmentEnsemble &GetRingLibraryFixedFragments() const;

      //! @brief returns atom indices of mutable fragments
      //! @return the mutable atoms
      const storage::Vector< size_t> &GetRingLibraryFixedFragmentsAtomIndices() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a fragment and returning a mutated fragment as a mutation result
      //! @param FRAGMENT small molecule of interest
      //! @return 3D conformer of the mutated fragment
      virtual math::MutateResult< FragmentComplete> operator()( const FragmentComplete &FRAGMENT) const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief set druglikeness filter type
      void SetDruglikenessType( const std::string &DRUGLIKENESS_TYPE);

      //! @brief set MDL SDF property label for the receptor path for BCL structure-based scoring
      void SetMDLReceptorLabel( const std::string &MDL);

      //! @brief set the scoring filter type
      void SetPropertyScorer( const descriptor::CheminfoProperty &PROPERTY_SCORER);

      //! @brief set maximum number of mutate attempts
      void SetNumberMaxMutates( const size_t N_MAX_MUTATES);

      //! @brief set bool to shuffle sequence of hydrogen atom removal when opening valences
      void SetOVShuffleH( const bool SHUFFLE_H);

      //! @brief set bool to reverse the order of hydrogen atom removal when opening valences
      void SetOVReverseH( const bool REVERSE_H);

      //! @brief set generic mutable atoms
      void SetMutableAtomIndices( const storage::Vector< size_t> &MUTABLE_ATOM_INDICES);

      //! @brief set the mutable atoms based on mutable elements
      void SetMutableElementsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< ElementType> &MUTABLE_ELEMENTS
      );

      //! @brief set the mutable atoms based on mutable fragments
      void SetMutableFragmentsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const bool COMPLEMENT
      );

      //! @brief set generic mutable atoms
      void SetFixedAtomIndices( const storage::Vector< size_t> &FIXED_ATOM_INDICES);

      //! @brief set the mutable atoms based on mutable elements
      void SetFixedElementsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< ElementType> &FIXED_ELEMENTS
      );

      //! @brief set the mutable atoms based on mutable fragments
      void SetFixedFragmentsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &FIXED_FRAGMENTS,
        const bool COMPLEMENT
      );

      //! @brief set generic mutable atoms
      void SetRingLibraryMutableAtomIndices( const storage::Vector< size_t> &MUTABLE_ATOM_INDICES);

      //! @brief set the mutable atoms based on mutable elements
      void SetRingLibraryMutableElementsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< ElementType> &MUTABLE_ELEMENTS
      );

      //! @brief set the mutable atoms based on mutable fragments
      void SetRingLibraryMutableFragmentsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const bool COMPLEMENT
      );

      //! @brief set generic mutable atoms
      void SetRingLibraryFixedAtomIndices( const storage::Vector< size_t> &FIXED_ATOM_INDICES);

      //! @brief set the mutable atoms based on mutable elements
      void SetRingLibraryFixedElementsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< ElementType> &FIXED_ELEMENTS
      );

      //! @brief set the mutable atoms based on mutable fragments
      void SetRingLibraryFixedFragmentsAtomIndices
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &FIXED_FRAGMENTS,
        const bool COMPLEMENT
      );

      //! @brief set the properties from one molecule onto the clean molecule;
      //! useful when properties computed during mutate are desirable to keep
      void SetPropertiesFromOther
      (
        FragmentComplete &TARGET_FRAGMENT,
        const FragmentComplete &REFERENCE_FRAGMENT
      ) const;

      //! @brief set the scaffold fragment
      void SetScaffoldFragment( const FragmentComplete &SCAFFOLD_FRAGMENT);

    //////////////////////
    // helper functions //
    //////////////////////

      // TODO make static
      //! @brief remove a hydrogen atom from a target atom
      //! @param FRAGMENT the molecule of interest
      //! @param ATOM_INDEX the index of the atom in the molecule of interest
      //! @param SHUFFLE_H if true, randomly select a hydrogen atom to remove
      //! @param REVERSE_H if true and not SHUFFLE_H, begin removal with the highest index hydrogen atom
      //! @return the new molecule, the index of the desired atom, and the original index of the removed hydrogen atom
      storage::Triplet< FragmentComplete, size_t, size_t> OpenValence
      (
        const FragmentComplete &FRAGMENT,
        const size_t &ATOM_INDEX,
        const bool SHUFFLE_H = true,
        const bool REVERSE_H = false
      ) const;

      // TODO make static
      //! @brief checks whether substitution at this atom is ortho, meta, or para directed
      //! @param MOLECULE the small molecule of interest
      //! @param ATOM simple pointer to the atom of interest in the molecule
      //! @return return true if the substitution is directed correctly
      bool IsRingSubstitutionDirected
      (
        const FragmentComplete &MOLECULE,
        util::SiPtr< const AtomConformationalInterface> &ATOM
      ) const;

      //! @brief select an atom from the target fragment
      //! @brief MOLECULE molecule from which to choose atom
      //! @param GLOBAL randomly select an atom from all atom indices
      //! @param RING_LIBRARY randomly select an atom from a ring (previously chosen
      //! from the user-passed ring library)
      //! @return return the chosen atom; base class default chooses global random atom
      //! using the input molecule atom selection specifications (not ring library data)
      virtual util::SiPtr< const AtomConformationalInterface> PickAtom
      (
        const FragmentComplete &MOLECULE,
        const bool GLOBAL = true,
        const bool RING_LIBRARY = false
      ) const;

      //! @brief selects a connection atom from the ring chosen from the ring library
      //! @param FRAGMENT molecule from which atom will be chosen
      //! @param RING_LIBRARY randomly select an atom from a ring (previously chosen
      //! from the user-passed ring library)
      //! @param H_BONDED_HEAVY if a hydrogen atom is selected, use its bonded heavy atom
      //! as the chosen atom instead
      //! @return the index of the chosen atom in the parent ring
      size_t RunPickAtom
      (
        const FragmentComplete &FRAGMENT,
        const bool RING_LIBRARY = false,
        const bool H_BONDED_HEAVY = false
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class FragmentMutateInterface

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_INTERFACE_H_
