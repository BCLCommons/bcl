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

// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)
#ifndef BCL_DESCRIPTOR_MOLECULE_TOTAL_TOXIC_FRAGMENTS_H_
#define BCL_DESCRIPTOR_MOLECULE_TOTAL_TOXIC_FRAGMENTS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeTotalToxicFragments
    //! Identifies common mutagenic fragments within molecules. These are fragments that would be expected
    //! to result in a failed Ames test. Based on https://pubs.acs.org/doi/10.1021/jm040835a
    //!
    //! @see @link example_descriptor_molecule_total_toxic_fragments.cpp @endlink
    //! @author chenz50, brownbp1
    //! @date Sep 15, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeTotalToxicFragments :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    private:

      //! pool of fragments to be picked from
      chemistry::FragmentEnsemble m_StabilizingFragments;
      std::string m_StabilizingFragmentsFilename;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeTotalToxicFragments();

      //! @brief Clone function
      //! @return pointer to new MoleculeTotalToxicFragments
      MoleculeTotalToxicFragments *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

      //! @brief calculate the total number of aromatic amines
      //! @param MOLECULE the molecule of interest
      size_t CountAromaticAmines( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of aromatic nitro groups
      //! @param MOLECULE the molecule of interest
      size_t CountAromaticNitro( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of aliphatic halides
      //! @param MOLECULE the molecule of interest
      size_t CountAliphaticHalide( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of three-membered heterocycles
      //! @param MOLECULE the molecule of interest
      size_t CountThreeMemberHeterocycle( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of nitroso groups
      //! @param MOLECULE the molecule of interest
      size_t CountNitroso( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of unsubstituted heteroatoms
      //! @param MOLECULE the molecule of interest
      size_t CountUnsubstitutedHeteroAtoms( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of azo groups
      //! @param MOLECULE the molecule of interest
      size_t CountAzoType( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate and return the total number of ring structures that contains polycyclic aromatic substructure
      //! @param MOLECULE the molecule of interest
      size_t CountPolycyclicAromatic( const chemistry::FragmentComplete &MOLECULE) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Count the number of heavy atoms connected to a given atom
      //! @param ATOM the atom of interest
      //! @return the number of heavy atoms bonded to the atom of interest
      size_t CountNumConnectedHeavyAtoms( const chemistry::AtomConformationalInterface &ATOM) const;

      //! @brief Helper function to recursively count the number of interface atoms between 2 aromatic rings in a given connected ring substructure
      //! @param sub_ring the connected ring substructure of a molecule
      //! @param visited: a list of boolean keeping track of if the atoms in the structure has been visited in the recursive calls or not
      //! @return the number of interface atom in the Ring substructure
      size_t CountAromaticInterfaceAtoms
      (
        const chemistry:: FragmentComplete &SUB_RING,
        const chemistry::AtomConformationalInterface &CURR_ATOM,
        std::vector< bool> &VISITED
      ) const;

      //! @brief Helper function to recursively visited neighbors of atoms in a subring structure to search for 3 member ring
      //! @param sub_ring the connected ring substructure of a molecule
      //! @param curr_idx, the index of current atom within the sub_ring fragment
      //! @param target_idx, the index of the original hetero atom from which we started the search
      //! @param num_bond_visited, number of bonded having stepped through since start, should be <= 3 for a 3 member ring
      //! @param visited: a list of boolean keeping track of if the atoms in the structure has been visited in the recursive calls or not
      //! @return 1 if the 3rd atom visited is the target atom (the original heteroatom), 0 otherwise

      size_t SearchThreeMemberCycle
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &SUB_RING,
        const chemistry::AtomConformationalInterface &CURR_ATOM,
        const chemistry::AtomConformationalInterface &ORIG_ATOM,
        const size_t NUM_BOND_VISITED
      ) const;

      //! @brief Searches for stabilizing groups nearby mutagenic substructures
      //! @return a vector (corresponding to each stabilizing group) mapping
      //! the bond distance between the stabilizing group and mutagenic group
      //! connecting atoms
      storage::Vector< storage::Map< storage::Pair< size_t, size_t>, size_t> >
      FindStabilizingGroups
      (
        const chemistry::FragmentComplete &MOLECULE
      ) const;

      //! @brief print out the result for each mutagenic fragment filter
      void PrintOutput( const storage::Vector< size_t> RESULTS) const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MoleculeTotalToxicFragments

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_TOTAL_TOXIC_FRAGMENTS_H_
