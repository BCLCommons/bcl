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

#ifndef BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPES_H_
#define BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configurational_bond_type_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConfigurationalBondTypes
    //! @brief enumeration class for bond types
    //!
    //! @see @link example_chemistry_configurational_bond_types.cpp @endlink
    //! @author mendenjl
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationalBondTypes :
      public util::Enumerate< ConfigurationalBondTypeData, ConfigurationalBondTypes>
    {
      friend class util::Enumerate< ConfigurationalBondTypeData, ConfigurationalBondTypes>;

    public:

    //////////
    // data //
    //////////

      // Bond types can be have one of four different orders: single/double/triple/unknown
      // one of three different conjugations non-conjugated/conjugated/aromatic
      // be either in a ring or not in a ring (aromatic bonds are always in a ring)
      const ConfigurationalBondType e_NonConjugatedSingleBond;
      const ConfigurationalBondType e_ConjugatedSingleBond;            //!< Single bond between two atoms with conjugated atom type
      const ConfigurationalBondType e_ConjugatedDoubleBond;            //!< Conjugated double bond; non-isometric
      const ConfigurationalBondType e_ConjugatedTripleBond;            //!< Conjugated triple bond
      const ConfigurationalBondType e_AromaticSingleBond;
      const ConfigurationalBondType e_AromaticDoubleBond;
      const ConfigurationalBondType e_AromaticTripleBond;              //!< Occurs in some large macrocycles
      const ConfigurationalBondType e_NonConjugatedSingleBondInRing;
      const ConfigurationalBondType e_ConjugatedSingleBondInRing;      //!< Conjugated single bond in ring
      const ConfigurationalBondType e_ConjugatedDoubleBondInRing;      //!< Conjugated double bond
      const ConfigurationalBondType e_ConjugatedTripleBondInRing;      //!< Conjugated triple bond

      // bond types with unknown order
      const ConfigurationalBondType e_ConjugatedBond;       //!< Conjugated bond of unknown order
      const ConfigurationalBondType e_AromaticBond;         //!< Aromatic bond of unknown order
      const ConfigurationalBondType e_ConjugatedBondInRing; //!< Conjugated bond of unknown order in ring

      // amide bond types
      const ConfigurationalBondType e_AmideSingleBond;       //!< Single amide bond

      // Stereo-specific bonds

      // a double bond with unknown isometry
      const ConfigurationalBondType e_ConjugatedDoubleBond_X;

      // E = entgegen = highest priority substituent on opposite side of bond
      const ConfigurationalBondType e_ConjugatedDoubleBond_E;

      // Z = zusammen = highest priority substituent on same side of bond
      const ConfigurationalBondType e_ConjugatedDoubleBond_Z;

      // E = entgegen = highest priority substituent on opposite side of bond
      const ConfigurationalBondType e_ConjugatedDoubleBondInRing_X;

      // E = entgegen = highest priority substituent on opposite side of bond
      const ConfigurationalBondType e_ConjugatedDoubleBondInRing_E;

      // Z = zusammen = highest priority substituent on same side of bond
      const ConfigurationalBondType e_ConjugatedDoubleBondInRing_Z;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct ConfigurationalBondTypes with all instances of the enums
      ConfigurationalBondTypes();

    ////////////////
    // operations //
    ////////////////

      //! Add a bond type to the enum
      //! @param BASE_BOND_TYPE the constitutional bond type
      //! @param ISOMETRY isometry of the bond type
      ConfigurationalBondType &AddBond
      (
        ConstitutionalBondType BASE_BOND_TYPE,
        const BondIsometry &ISOMETRY
      );

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class ConfigurationalBondTypes

    BCL_API
    const ConfigurationalBondTypes &GetConfigurationalBondTypes();

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< chemistry::ConfigurationalBondTypeData, chemistry::ConfigurationalBondTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPES_H_

