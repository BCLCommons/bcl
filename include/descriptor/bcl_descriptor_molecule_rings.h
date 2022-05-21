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

#ifndef BCL_DESCRIPTOR_MOLECULE_RINGS_H_
#define BCL_DESCRIPTOR_MOLECULE_RINGS_H_

// include the namespace header
#include "bcl_descriptor.h"

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
    //! @class MoleculeRings
    //! @brief calculates the # of rings with a given conjugation from the molecule
    //!
    //! @see @link example_descriptor_molecule_rings.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeRings :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    private:

    //////////
    // data //
    //////////

      chemistry::ConstitutionalBondTypeData::ConjugationEnum m_Conjugation; //!< bond property to retrieve from the small molecule
      // Aromatic rings have all aromatic bonds
      // Conjugated rings have no non-conjugated single bonds but are not aromatic
      // Non-conjugated rings have at least one bond that is a non-conjugated single bond
      // Any rings include all of the above

      // Whether to count all macrocycles (>8 atoms)
      bool m_MacrocyclesOnly;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a property
      //! @param PROPERTY the property to retrieve
      explicit MoleculeRings(
        const chemistry::ConstitutionalBondTypeData::Conjugation &CONJUGATION = chemistry::ConstitutionalBondTypeData::e_Any,
        const bool &MACROCYCLES_ONLY = false
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeRings
      MoleculeRings *Clone() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class MoleculeRings

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_RINGS_H_
