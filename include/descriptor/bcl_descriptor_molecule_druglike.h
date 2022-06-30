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

#ifndef BCL_DESCRIPTOR_MOLECULE_DRUGLIKE_H_
#define BCL_DESCRIPTOR_MOLECULE_DRUGLIKE_H_

// include the namespace header
#include "bcl_descriptor.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "sched/bcl_sched_mutex.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeDruglike
    //! @brief checks multiple parameters to determine if a molecule is druglike
    //! @details Rotatable bond definition is from:
    //!
    //! @see @link example_descriptor_molecule_druglike.cpp @endlink
    //! @author brownbp1
    //! @date Dec 13, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeDruglike :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

      float m_MinWeight;
      float m_MaxWeight;
      size_t m_MaxRingSize;
      float m_MinLogP;
      float m_MaxLogP;
      float m_MaxTPSA;
      size_t m_MaxHBDA;
      size_t m_MaxNRotBonds;
      float m_MaxTotalBondEnergy;
      size_t m_MaxNF;
      size_t m_MaxNCl;
      size_t m_MaxNBr;
      size_t m_MaxNI;
      size_t m_MaxNHalogens;
      size_t m_MaxRingHalogens;
      size_t m_MaxNonAromaticClBrI;
      float m_MaxComplexity;
      bool m_EnforceHitlike;
      bool m_EnforceDruglikeRings;
      bool m_EnforceDruglikeECFP;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeDruglike();

      //! @brief hitlike constructor
      explicit MoleculeDruglike
      (
        const bool &ENFORCE_HITLIKE = false
      );

      //! @brief constructor
      explicit MoleculeDruglike
      (
        const float &MIN_WEIGHT,
        const float &MAX_WEIGHT,
        const size_t &MAX_RING_SIZE,
        const float &MIN_LOGP,
        const float &MAX_LOGP,
        const float &MAX_TPSA,
        const size_t &MAX_HBDA,
        const size_t &MAX_N_ROT_BONDS,
        const float &MAX_BOND_PROPENSITY,
        const size_t &MAX_N_F,
        const size_t &MAX_N_CL,
        const size_t &MAX_N_BR,
        const size_t &MAX_N_I,
        const size_t &MAX_TOTAL_HALOGENS,
        const size_t &MAX_RING_HALOGENS,
        const size_t &MAX_NON_AROMATIC_CLBRI,
        const float &MAX_COMPLEXITY,
        const bool &ENFORCE_HITLIKE,
        const bool &ENFORCE_DRUGLIKE_RINGS,
        const bool &ENFORCE_DRUGLIKE_ECFP
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeDruglike
      MoleculeDruglike *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate the max number of halogen atoms on a ring in a molecule
      //! @param MOLECULE the molecule of interest
      size_t CountRingHalogens( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief determines if the ratio of ring heavy atoms to substituted halogens is druglike
      //! @param MOLECULE the molecule of interest
      bool HasGoodRingAtomHalogenRatio( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the max number of atoms in a ring in a molecule
      //! @param MOLECULE the molecule of interest
      size_t MaxRingFragmentSize( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of non-aromatic ring Cl, Br, I halogen atoms
      //! @param MOLECULE the molecule of interest
      size_t CountNonAroClBrI( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief determine whether molecule rings are drug-like
      //! @param MOLECULE the molecule of interest
      bool ContainsOnlyDruglikeRings( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief determine if undesirable nitrogen covalent bonds exist in a molecule
      //! @param MOLECULE the molecule of interest
      bool ContainsBadNitrogenLinkages( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief determine if molecule contains aromatic amine
      //! @param MOLECULE the molecule of interest
      bool ContainsAromaticAmine( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief determine if molecules contains a reactive alkene
      //! note that this will need to be disabled if we choose to include certain covalent warheads, e.g. acrylamide
      //! @param MOLECULE the molecule of interest
      bool ContainsReactiveAlkene( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class MoleculeDruglike

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_DRUGLIKE_H_
