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

#ifndef BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPES_H_
#define BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configurational_bond_types.h"
#include "bcl_chemistry_constitutional_bond_type_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstitutionalBondTypes
    //! @brief enumeration class for bond types
    //!
    //! @see @link example_chemistry_constitutional_bond_types.cpp @endlink
    //! @author mendenjl
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConstitutionalBondTypes :
      public util::Enumerate< ConstitutionalBondTypeData, ConstitutionalBondTypes>
    {
      friend class util::Enumerate< ConstitutionalBondTypeData, ConstitutionalBondTypes>;

    public:

    //////////
    // data //
    //////////

      // Bond types can be have one of four different orders: single/double/triple/unknown
      // one of three different conjugations non-conjugated/conjugated/aromatic
      // be either in a ring or not in a ring (aromatic bonds are always in a ring)
      const ConstitutionalBondType e_NonConjugatedSingleBond;
      const ConstitutionalBondType e_ConjugatedSingleBond;            //!< Conjugated single bond (single bond between two atoms with conjugated atom type)
      const ConstitutionalBondType e_ConjugatedDoubleBond;            //!< Conjugated double bond
      const ConstitutionalBondType e_ConjugatedTripleBond;            //!< Conjugated triple bond
      const ConstitutionalBondType e_ConjugatedBond;                  //!< Conjugated bond of unknown order
      const ConstitutionalBondType e_AromaticSingleBond;
      const ConstitutionalBondType e_AromaticDoubleBond;
      const ConstitutionalBondType e_AromaticTripleBond;              //!< Occurs in some large macrocycles
      const ConstitutionalBondType e_AromaticBond;                    //!< Aromatic bond of unknown order
      const ConstitutionalBondType e_NonConjugatedSingleBondInRing;
      const ConstitutionalBondType e_ConjugatedSingleBondInRing;      //!< Conjugated single bond in ring
      const ConstitutionalBondType e_ConjugatedDoubleBondInRing;      //!< Conjugated double bond
      const ConstitutionalBondType e_ConjugatedTripleBondInRing;      //!< Conjugated triple bond
      const ConstitutionalBondType e_ConjugatedBondInRing;            //!< Conjugated bond of unknown order
      const ConstitutionalBondType e_AmideSingleBond;                 //!< Single amide bond (not in ring)

    private:

    ///////////////
    // constants //
    ///////////////
      enum
      {
        s_MaxSdfID = 5,            //!< maximum legitimate value for the bond order in an sdf file
        s_MaxValenceBonds = 4,     //!< maximum number of valence bonds
        s_MaxValenceElectrons = 6  //!< maximum electrons in valence bonds
      };

      // array of constitutional bond types based on sdf id
      util::SiPtr< const ConstitutionalBondType> m_BondTypesBySdfId[ s_MaxSdfID + 1];

      // vectors of valence bonds returned depending on
      storage::Vector< storage::Vector< ConstitutionalBondType> > m_ValenceBonds[ s_MaxValenceBonds + 1];

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct ConstitutionalBondTypes with all instances of the enums
      ConstitutionalBondTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief retrieve valence bonds given the # of valence bonds and e- in valence bonds
      //! @param VALENCE_BONDS # of valence bonds
      //! @param ELECTRONS_IN_VALENCE_BONDS # of electrons in valence bonds
      //! @return the valence bond types, if known
      const storage::Vector< ConstitutionalBondType> &GetValenceBonds
      (
        const size_t &VALENCE_BONDS,
        const size_t &ELECTRONS_IN_VALENCE_BONDS
      ) const;

      //! @brief retrieve Bond type from sdf-id, sdf-alternative id, and bond annotation id
      //! @param SDF_ID the id of the bond type in the sdf file (never 0)
      //! @return the corresponding ConstitutionalBondType
      ConstitutionalBondType FindBondTypeFromSDFInfo( const size_t &SDF_ID) const;

    }; // class ConstitutionalBondTypes

    BCL_API
    const ConstitutionalBondTypes &GetConstitutionalBondTypes();

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< chemistry::ConstitutionalBondTypeData, chemistry::ConstitutionalBondTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPES_H_

