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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_EXTEND_WITH_LINKER_H_
#define BCL_CHEMISTRY_FRAGMENT_MUTATE_EXTEND_WITH_LINKER_H_

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
    //! @class FragmentMutateExtendWithLinker
    //! @brief Used to extend a molecule either by breaking a bond and inserting a motif before reconnecting,
    //! or by linking to a new ring system
    //!
    //! @see @link example_chemistry_fragment_mutate_extend_with_linker.cpp @endlink
    //! @author brownbp1
    //! @date April 20, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateExtendWithLinker :
      public FragmentMutateInterface
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! rings from fragment database
      util::ShPtr< FragmentEnsemble> m_Rings;
      std::string m_RingsFilename;

      //! probability to insert the extension within the fragment;
      //! remainder of the probability goes toward extending out from the fragment
      float m_ExtendWithinProb;

      //! relative probability for an amide link
      float m_AmideLinkProb;

      //! probability amide link will attach nitrogen to fragment A
      float m_AmideNToAProb;

      //! relative probability for an ester link
      float m_EsterLinkProb;

      //! probability ester link will attach oxygen to fragment A
      float m_EsterOToAProb;

      //! relative probability for a methoxy link
      float m_MethoxyLinkProb;

      //! probability methoxy link will attach oxygen to fragment A
      float m_MethoxyOToAProb;

      //! relative probability for an ethoxy link
      float m_EthoxyLinkProb;

      //! probability ethoxy link will attach oxygen to fragment A
      float m_EthoxyOToAProb;

      //! relative probability for a single element link
      float m_SingleElementLinkProb;

      //! probabilities for selecting different elements for the single element link
      float m_B; // boron
      float m_C; // carbon
      float m_O; // oxygen
      float m_N; // nitrogen
      float m_P; // phosphorous
      float m_S; // sulfur
      float m_Se; // selenium

      //! relative probability of direct link
      float m_DirectLinkProb;

      //! relative probability of an alkyl link
      float m_AlkylLinkProb;

      //! probability that the alkyl linker will have 3 carbons; only other option is 2
      float m_ThreeC;

      //! probability of converting alkane to alkyne
      float m_Alkyne;

      //! relative probability for extending a fragment internally with a ring
      float m_RingLinkProb;

      //! fragment min size
      size_t m_FragmentMinSize;

      //! allow duplication of fragment during intramolecular extension
      bool m_AllowFragmentDuplication;

      //! atom indices that can be mutated
      storage::Vector< size_t> m_PairedAtomIndices;
      std::string m_PairedAtoms;

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
      FragmentMutateExtendWithLinker();

      //! @brief constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      FragmentMutateExtendWithLinker
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
      FragmentMutateExtendWithLinker
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
      FragmentMutateExtendWithLinker
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
      FragmentMutateExtendWithLinker *Clone() const;

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

      //! @brief set the absolute probability of extending the fragment internally
      void SetExtendWithinProb( const float EXTEND_WITHIN_PROB);

      //! @brief set the size of the smallest fragment to be isolated
      void SetMinFragmentSize( const size_t MIN_FRAGMENT_SIZE);

      //! @brief set whether to allow duplicates if only one fragment meets min size requirement
      void SetAllowDuplicates( const bool ALLOW_DUPLICATES);

      //! @brief set the relative probability of extending with an amide linker
      void SetAmideLinkProb( const float AMIDE_LINK_PROB);

      //! @brief set the probability of orienting the amide N toward fragment A with the amide linker
      void SetAmideNToAProb( const float AMIDE_N_TO_A_PROB);

      //! @brief set the relative probability of extending with an ester linker
      void SetEsterLinkProb( const float ESTER_LINK_PROB);

      //! @brief set the probability of orienting the ester O toward fragment A with the ester linker
      void SetEsterOToAProb( const float ESTER_O_TO_A_PROB);

      //! @brief set the relative probability of extending with a methoxy linker
      void SetMethoxyLinkProb( const float METHOXY_LINK_PROB);

      //! @brief set the probability of orienting the methoxy O toward fragment A with the methoxy linker
      void SetMethoxyOToAProb( const float METHOXY_O_TO_A_PROB);

      //! @brief set the relative probability of extending with an ethoxy linker
      void SetEthoxyLinkProb( const float ETHOXY_LINK_PROB);

      //! @brief set the probability of orienting the ethoxy O toward fragment A with the ethoxy linker
      void SetEthoxyOToAProb( const float ETHOXY_O_TO_A_PROB);

      //! @brief set the relative probability of extending with a single element linker
      void SetSingleElementLinkProb( const float SINGLE_ELEMENT_LINK_PROB);

      //! @brief set the relative probability of linking with a single B
      void SetBProb( const float B_PROB);

      //! @brief set the relative probability of linking with a single C
      void SetCProb( const float C_PROB);

      //! @brief set the relative probability of linking with a single O
      void SetOProb( const float O_PROB);

      //! @brief set the relative probability of linking with a single N
      void SetNProb( const float N_PROB);

      //! @brief set the relative probability of linking with a single P
      void SetPProb( const float P_PROB);

      //! @brief set the relative probability of linking with a single S
      void SetSProb( const float S_PROB);

      //! @brief set the relative probability of linking with a single Se
      void SetSeProb( const float SE_PROB);

      //! @brief set the relative probability of extending with a direct link
      void SetDirectLinkProb( const float DIRECT_LINK_PROB);

      //! @brief set the relative probability of extending with an alkyl linker
      void SetAlkylLinkProb( const float ALKYL_LINK_PROB);

      //! @brief set the relative probability that the alkyl linker will contain a triple bond
      void SetAlkyneProb( const float ALKYNE_PROB);

      //! @brief set the relative probability of extending with a ring
      void SetRingLinkProb( const float RING_LINK_PROB);

      //! @brief directly link two fragments via one of the linker strategies below
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDEX_A link indices in first molecule
      //! @param LINK_INDEX_B link indices in second molecule
      //! @param RING_LIBRARY optional ring library if not provided through map
      //! @param RING_LINK_INDEX_A optional specific index by which to connect fragment A to the ring
      //! @param RING_LINK_INDEX_B optional specific index by which to connect fragment B to the ring
      //! @return the newly generated molecules
      FragmentComplete RingLink
      (
        const FragmentComplete &FRAGMENT_A,
        const FragmentComplete &FRAGMENT_B,
        const size_t &LINK_INDEX_A,
        const size_t &LINK_INDEX_B,
        const FragmentEnsemble &RING_LIBRARY = FragmentEnsemble(),
        const size_t RING_LINK_INDEX_A = util::GetUndefinedSize_t(),
        const size_t RING_LINK_INDEX_B = util::GetUndefinedSize_t()
      ) const;

      //! @brief link two fragments via a ring
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDEX_A link indices in first molecule
      //! @param LINK_INDEX_B link indices in second molecule
      //! @param RING the ring with which to bridge the two fragments
      //! @param RING_LINK_INDEX_A optional specific index by which to connect fragment A to the ring
      //! @param RING_LINK_INDEX_B optional specific index by which to connect fragment B to the ring
      //! @return the newly generated molecules
      FragmentComplete RingLink
      (
        const FragmentComplete &FRAGMENT_A,
        const FragmentComplete &FRAGMENT_B,
        const size_t &LINK_INDEX_A,
        const size_t &LINK_INDEX_B,
        const FragmentComplete &RING,
        const size_t RING_LINK_INDEX_A = util::GetUndefinedSize_t(),
        const size_t RING_LINK_INDEX_B = util::GetUndefinedSize_t()
      ) const;

      //! @brief link two fragments via a ring
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDEX_A link indices in first molecule
      //! @param LINK_INDEX_B link indices in second molecule
      //! @return the newly generated molecules
      FragmentComplete DirectLink
      (
        const FragmentComplete &FRAGMENT_A,
        const FragmentComplete &FRAGMENT_B,
        const size_t&LINK_INDEX_A,
        const size_t&LINK_INDEX_B
      ) const;

      //! @brief link two fragments via a single element
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDEX_A link indices in first molecule
      //! @param LINK_INDEX_B link indices in second molecule
      //! @param ELEMENT_TYPE element type serving as a link
      //! @return the newly generated molecules
      FragmentComplete SingleElementLink
      (
        const FragmentComplete &FRAGMENT_A,
        const FragmentComplete &FRAGMENT_B,
        const size_t&LINK_INDEX_A,
        const size_t&LINK_INDEX_B,
        const std::string &ELEMENT_TYPE
      ) const;

       //! @brief link two fragments with an alkyl chain
       //! @param FRAGMENT_A first molecule
       //! @param FRAGMENT_B second molecule
       //! @param LINK_INDEX_A link indices in first molecule
       //! @param LINK_INDEX_B link indices in second molecule
       //! @param REPEATS the number of linker repeats
       //! @return the newly generated molecules
       FragmentComplete AlkylLink
       (
         const FragmentComplete &FRAGMENT_A,
         const FragmentComplete &FRAGMENT_B,
         const size_t&LINK_INDEX_A,
         const size_t&LINK_INDEX_B,
         const size_t &REPEATS
       ) const;

       //! @brief link two fragments with a methoxy repeat
       //! @param FRAGMENT_A first molecule
       //! @param FRAGMENT_B second molecule
       //! @param LINK_INDEX_A link indices in first molecule
       //! @param LINK_INDEX_B link indices in second molecule
       //! @param REPEATS the number of linker repeats
       //! @return the newly generated molecules
       FragmentComplete MethoxyLink
       (
         const FragmentComplete &FRAGMENT_A,
         const FragmentComplete &FRAGMENT_B,
         const size_t&LINK_INDEX_A,
         const size_t&LINK_INDEX_B,
         const size_t &REPEATS
       ) const;

       //! @brief link two fragments with an ethoxy repeat
       //! @param FRAGMENT_A first molecule
       //! @param FRAGMENT_B second molecule
       //! @param LINK_INDEX_A link indices in first molecule
       //! @param LINK_INDEX_B link indices in second molecule
       //! @param REPEATS the number of linker repeats
       //! @return the newly generated molecules
       FragmentComplete EthoxyLink
       (
         const FragmentComplete &FRAGMENT_A,
         const FragmentComplete &FRAGMENT_B,
         const size_t&LINK_INDEX_A,
         const size_t&LINK_INDEX_B,
         const size_t &REPEATS
       ) const;

       //! @brief link two fragments with an amide repeat
       //! @param FRAGMENT_A first molecule
       //! @param FRAGMENT_B second molecule
       //! @param LINK_INDEX_A link indices in first molecule
       //! @param LINK_INDEX_B link indices in second molecule
       //! @param REPEATS the number of linker repeats
       //! @return the newly generated molecules
       FragmentComplete AmideLink
       (
         const FragmentComplete &FRAGMENT_A,
         const FragmentComplete &FRAGMENT_B,
         const size_t&LINK_INDEX_A,
         const size_t&LINK_INDEX_B,
         const size_t &REPEATS
       ) const;

       //! @brief link two fragments with an ester repeat
       //! @param FRAGMENT_A first molecule
       //! @param FRAGMENT_B second molecule
       //! @param LINK_INDEX_A link indices in first molecule
       //! @param LINK_INDEX_B link indices in second molecule
       //! @param REPEATS the number of linker repeats
       //! @return the newly generated molecules
       FragmentComplete EsterLink
       (
         const FragmentComplete &FRAGMENT_A,
         const FragmentComplete &FRAGMENT_B,
         const size_t&LINK_INDEX_A,
         const size_t&LINK_INDEX_B,
         const size_t &REPEATS
       ) const;

      //! @brief link two fragments with an amide repeat
      //! @param CONNECTION whether to connect the amide to molecule A via C or N or both
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      FragmentEnsemble GenerateAmideLinker
      (
        const std::string &CONNECTION,
        const size_t &REPEATS
      ) const;

      //! @brief extends a new ring from a current ring with one of the above linkers
      //! @param FRAGMENT_A first molecule
      //! @param LINK_INDEX_A link indices in first molecule
      //! @param REPEATS the number of linker repeats
      //! @param ELEMENT_TYPE element type serving as a link of SingleElementLinker
      //! @return the newly generated molecules
      FragmentComplete RingExtension
      (
        const FragmentComplete &FRAGMENT_A,
        const size_t &LINK_INDEX_A,
        const size_t &REPEATS,
        const std::string &ELEMENT_TYPE
      ) const;

    private:

      //! @brief selects the method with which to link the fragments
      //! @return a string indicating the method to use
      std::string ChooseLinkMethod( const bool &EXTEND_WITHIN) const;

      //! @brief selects an element to use with the SingleElementLink function
      //! @return string corresponding to chosen element
      std::string ChooseLinkElement() const;

      //! @brief checks if two atoms in rings are in the same ring or different rings
      //! @param ATOM_INDEX_A first atom
      //! @param ATOM_INDEX_B second atom
      //! @param PARENT_MOL the fragment containing both atoms
      //! @return true if the atoms are in the same ring, false otherwise
      bool IfAtomsInSameRing
      (
        const size_t &ATOM_INDEX_A,
        const size_t &ATOM_INDEX_B,
        const FragmentComplete &PARENT_MOL
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

    }; // class FragmentMutateExtendWithLinker

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_EXTEND_WITH_LINKER_H_
