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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_RING_SWAP_H_
#define BCL_CHEMISTRY_FRAGMENT_MUTATE_RING_SWAP_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "bcl_chemistry_search_fragment_library_from_tree.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_vector_nd.h"
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
    //! @class FragmentMutateRingSwap
    //! @brief swaps a ring system from the fragment library for one in the library. At most one ring may have different
    //!        size between swapped systems (which allows for removing/adding rings to a ring system if the rest of the
    //!        rings have the same size
    //!
    //! @see @link example_chemistry_fragment_mutate_ring_swap.cpp @endlink
    //! @author mendenjl, brownbp1
    //! @date Sep 10, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateRingSwap :
      public FragmentMutateInterface
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! rotamer library information
      util::ShPtr< SearchFragmentLibraryFromTree> m_RotamerLibrarySearcher;
      std::string m_RingLibraryFilename;

      //! pool of fragments to be picked from, indexed by # of double bond valences
      storage::Vector< storage::Vector< FragmentComplete> > m_FragmentPool;

      //! fragment pool scaffold sums
      storage::Vector< size_t> m_FragmentPoolScaffoldSums;

      //! probability of expanding a single atom into a ring (1-this probability is probability of swapping rings)
      //! This also defines the probability of removing a ring entirely
      double m_RingInitiationProbability;

      //! Try to fix the conformation after the ring swap (bond angles and lengths only. Dihedrals are not preserved)
      //! This is *very* slow and not recommended if you're only using 2d descriptors
      bool m_FixGeometry;

      //! Neutralize the molecule after ring swap
      //! Recommended unless explicitly using a model trained on formally charged molecules
      bool m_Neutralize;

      //! if true, at most one ring may change size. For an additional ring to be added to the system then, all other
      //! rings must have the same size
      bool m_RestrictToNoMoreThanOneRingSizeChange;

      //! if false, do not allow collapse of rings greater than 4 atoms (prevents implosion of molecules whose cores are rings)
      bool m_AllowLargeRingCollapse;

      //! if true, align new ring to current ring prior to substitution to preserve topological distances between substituents as best as possible
      bool m_AlignRings;

      //! extend atoms included in conformational sampling this many bonds out from any perturbed atom
      size_t m_ExtendAdjacentAtoms;

      //! perform a quick substructure-based ensemble align and choose best conformer based on ChargeRMSD
      bool m_ChooseBestAlignedConf;

      //! Scheme used for comparing whether two bonds are equivalent
      ConfigurationalBondTypeData::DataEnum m_BondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparisonType;

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
      FragmentMutateRingSwap();

      //! @brief pose-dependent constructor
      explicit FragmentMutateRingSwap
      (
        const util::ShPtr< SearchFragmentLibraryFromTree> &FRAGMENT_LIBRARY,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA,
        const bool &FIX_GEOMETRY,
        const bool &NEUTRALIZE,
        const double &RING_INITIATION_PROBABILITY,
        const bool &PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE,
        const bool &ALLOW_LARGE_RING_COLLAPSE
      );

      //! @brief pose-independent constructor
      explicit FragmentMutateRingSwap
      (
        const util::ShPtr< SearchFragmentLibraryFromTree> &FRAGMENT_LIBRARY,
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const bool &CORINA,
        const bool &FIX_GEOMETRY,
        const bool &NEUTRALIZE,
        const double &RING_INITIATION_PROBABILITY,
        const bool &PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE,
        const bool &ALLOW_LARGE_RING_COLLAPSE
      );

      //! @brief clone constructor
      FragmentMutateRingSwap *Clone() const;

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

      //! @brief compute whether two vectors differ by at most one element
      static bool SortedVectorsDifferByNoMoreThanOneElement
      (
        const linal::Vector< float> &VEC_A,
        const linal::Vector< float> &VEC_B
      );

      //! @brief add ring size information as a property to the molecule
      static typename descriptor::CacheMap::value_type GetRingSizeInformation( const FragmentComplete &FRAG);

      //! @brief add ring size information as a property to the molecule
      static typename descriptor::CacheMap::value_type GetRingSizeInformation
      (
        const FragmentComplete &FRAG,
        const storage::Vector< size_t> &RING
      );

      //! @brief add ring size information as a property to the molecule
      static size_t GetNumberDoubleBondValences( const FragmentComplete &FRAG);

      //! @brief add ring size information as a property to the molecule
      static size_t GetNumberDoubleBondValences( const FragmentComplete &FRAG, const storage::Vector< size_t> &RING);

      //! @brief Get the counts for the given fragment
      static size_t GetCounts( const FragmentComplete &FRAG);

      //! @brief get a random ring system of the given fragment provided an atom index
      static storage::Vector< size_t> GetRandomRingSystem( const FragmentComplete &FRAG, const size_t ATOM_INDEX, bool ONLY_SELECT_COLLAPSABLE_RINGS = false);

      //! @brief get a random atom that is not in a ring system for the given fragment provided an atom index (not really random anymore)
      static size_t GetRandomNonringAtom( const FragmentComplete &FRAG, const size_t ATOM_INDEX);

      //! @brief get a random ring system of the given fragment
      static storage::Vector< size_t> GetGlobalRandomRingSystem( const FragmentComplete &FRAG, bool ALLOW_LARGE_RING_COLLAPSE = false);

      //! @brief get a random atom that is not in a ring system for the given fragment
      static size_t GetGlobalRandomNonringAtom( const FragmentComplete &FRAG);

      //! @brief return the immediate substituent indices for all atoms of a given ring
      static storage::Vector< storage::Vector< storage::Pair< size_t, ConfigurationalBondType> > > GetSubstitutents
      (
        const FragmentComplete &FRAG,
        const storage::Vector< size_t> &RING,
        const bool &COLLAPSE
      );

      //! @brief return whether the ring can be collapsed to a single atom
      static bool IsCollapsible
      (
        const FragmentComplete &FRAG,
        const storage::Vector< size_t> &RING
      );

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

    }; // class FragmentMutateRingSwap

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_RING_SWAP_H_
