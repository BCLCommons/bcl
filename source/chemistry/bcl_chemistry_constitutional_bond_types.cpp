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
#include "chemistry/bcl_chemistry_constitutional_bond_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct all bond types
    ConstitutionalBondTypes::ConstitutionalBondTypes() :
      util::Enumerate< ConstitutionalBondTypeData, ConstitutionalBondTypes>( false),
      e_NonConjugatedSingleBond(         AddEnum( "NonConjugatedSingleBond"        , ConstitutionalBondTypeData( 2, 1, 1, ConstitutionalBondTypeData::e_Nonconjugated, false))),
      e_ConjugatedSingleBond(            AddEnum( "ConjugatedSingleBond"           , ConstitutionalBondTypeData( 2, 1, 5, ConstitutionalBondTypeData::e_Conjugated   , false))),
      e_ConjugatedDoubleBond(            AddEnum( "ConjugatedDoubleBond"           , ConstitutionalBondTypeData( 4, 2, 5, ConstitutionalBondTypeData::e_Conjugated   , false))),
      e_ConjugatedTripleBond(            AddEnum( "ConjugatedTripleBond"           , ConstitutionalBondTypeData( 6, 3, 3, ConstitutionalBondTypeData::e_Conjugated   , false))),
      e_ConjugatedBond(                  AddEnum( "ConjugatedBond"                 , ConstitutionalBondTypeData( 3, 5, 5, ConstitutionalBondTypeData::e_Conjugated   , false))),
      e_AromaticSingleBond(              AddEnum( "AromaticSingleBond"             , ConstitutionalBondTypeData( 2, 1, 6, ConstitutionalBondTypeData::e_Aromatic     , true))),
      e_AromaticDoubleBond(              AddEnum( "AromaticDoubleBond"             , ConstitutionalBondTypeData( 4, 2, 7, ConstitutionalBondTypeData::e_Aromatic     , true))),
      e_AromaticTripleBond(              AddEnum( "AromaticTripleBond"             , ConstitutionalBondTypeData( 6, 3, 3, ConstitutionalBondTypeData::e_Aromatic     , true))),
      e_AromaticBond(                    AddEnum( "AromaticBond"                   , ConstitutionalBondTypeData( 3, 4, 4, ConstitutionalBondTypeData::e_Aromatic     , true))),
      e_NonConjugatedSingleBondInRing(   AddEnum( "NonConjugatedSingleBondInRing"  , ConstitutionalBondTypeData( 2, 1, 1, ConstitutionalBondTypeData::e_Nonconjugated, true))),
      e_ConjugatedSingleBondInRing(      AddEnum( "ConjugatedSingleBondInRing"     , ConstitutionalBondTypeData( 2, 1, 5, ConstitutionalBondTypeData::e_Conjugated   , true))),
      e_ConjugatedDoubleBondInRing(      AddEnum( "ConjugatedDoubleBondInRing"     , ConstitutionalBondTypeData( 4, 2, 5, ConstitutionalBondTypeData::e_Conjugated   , true))),
      e_ConjugatedTripleBondInRing(      AddEnum( "ConjugatedTripleBondInRing"     , ConstitutionalBondTypeData( 6, 3, 3, ConstitutionalBondTypeData::e_Conjugated   , true))),
      e_ConjugatedBondInRing(            AddEnum( "ConjugatedBondInRing"           , ConstitutionalBondTypeData( 3, 5, 5, ConstitutionalBondTypeData::e_Conjugated   , true))),
      e_AmideSingleBond(                 AddEnum( "AmideSingleBond"                , ConstitutionalBondTypeData( 2, 1, 1, ConstitutionalBondTypeData::e_Amide        , false)))
    {
      m_BondTypesBySdfId[ 0] = &e_Undefined;
      m_BondTypesBySdfId[ 1] = &e_NonConjugatedSingleBond;
      m_BondTypesBySdfId[ 2] = &e_ConjugatedDoubleBond;
      m_BondTypesBySdfId[ 3] = &e_ConjugatedTripleBond;
      m_BondTypesBySdfId[ 4] = &e_AromaticBond;
      m_BondTypesBySdfId[ 5] = &e_ConjugatedBond;

      // attach each bond type to the related types
      for( size_t index( 0), num_bond_types( GetEnumCount()); index < num_bond_types; ++index)
      {
        ConstitutionalBondType &type( GetEnumFromIndex( index));
        type->AttachToRelatedBondTypes( *this);
      }

      // setup relations for the undefined bond type; never used unless there are undefined atom types
      e_Undefined->AttachToRelatedBondTypes( *this);

      const size_t max_bond_order( 3);

      // handle cases with valences
      for( size_t valence_bond_count( 0); valence_bond_count <= s_MaxValenceBonds; ++valence_bond_count)
      {
        // determine the max excess electrons for this type
        const size_t max_excess_electrons
        (
          std::min( size_t( s_MaxValenceElectrons), max_bond_order * valence_bond_count) - valence_bond_count
        );
        m_ValenceBonds[ valence_bond_count].Resize( max_excess_electrons + 1);
        for( size_t excess_electrons( 0); excess_electrons <= max_excess_electrons; ++excess_electrons)
        {
          size_t remaining_electrons( valence_bond_count + excess_electrons);
          size_t remaining_bonds( valence_bond_count);

          // make all single bonds into double bonds first
          const size_t single_bonds( std::max( valence_bond_count, excess_electrons) - excess_electrons);
          remaining_electrons -= single_bonds;
          remaining_bonds -= single_bonds;

          // constraint: # of bonds:
          // # of double bonds (x) + # of triple bonds (y) = # of remaining bonds (B)
          // constraint: # of electrons in bonds
          // 2 * x + 3 * y = # of remaining electrons in bonds (E)
          // Solving for x: 3 * B - E
          const size_t double_bonds( 3 * remaining_bonds - remaining_electrons);
          const size_t triple_bonds( remaining_bonds - double_bonds);

          m_ValenceBonds[ valence_bond_count]( excess_electrons).Append
          (
            storage::Vector< ConstitutionalBondType>( single_bonds, e_NonConjugatedSingleBond)
          );
          m_ValenceBonds[ valence_bond_count]( excess_electrons).Append
          (
            storage::Vector< ConstitutionalBondType>( double_bonds, e_ConjugatedDoubleBond)
          );
          m_ValenceBonds[ valence_bond_count]( excess_electrons).Append
          (
            storage::Vector< ConstitutionalBondType>( triple_bonds, e_ConjugatedTripleBond)
          );
        }
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConstitutionalBondTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief retrieve Bond type from sdf-id, sdf-alternative id, and bond annotation id
    //! @param SDF_ID the id of the bond type in the sdf file (never undefined)
    //! @return the corresponding ConstitutionalBondType
    ConstitutionalBondType ConstitutionalBondTypes::FindBondTypeFromSDFInfo( const size_t &SDF_ID) const
    {
      if( !SDF_ID || SDF_ID > s_MaxSdfID)
      {
        BCL_MessageCrt
        (
          "Warning: Found bond line containing invalid bond type = " + util::Format()( SDF_ID)
          + ". Bond will be omitted"
        );
      }
      // if the sdf file id is legitimate, return the corresponding bond type by sdf id
      if( SDF_ID <= s_MaxSdfID)
      {
        // returns undefined for 0
        return *m_BondTypesBySdfId[ SDF_ID];
      }

      // sdf file id out of range; return undefined
      return e_Undefined;
    }

    //! @brief retrieve valence bonds given the # of valence bonds and e- in valence bonds
    //! @param VALENCE_BONDS # of valence bonds
    //! @param ELECTRONS_IN_VALENCE_BONDS # of electrons in valence bonds
    //! @return the valence bond types, if known
    const storage::Vector< ConstitutionalBondType> &ConstitutionalBondTypes::GetValenceBonds
    (
      const size_t &VALENCE_BONDS,
      const size_t &ELECTRONS_IN_VALENCE_BONDS
    ) const
    {
      if
      (
        VALENCE_BONDS == 0
        || VALENCE_BONDS > s_MaxValenceBonds
        || ELECTRONS_IN_VALENCE_BONDS > s_MaxValenceElectrons
        || VALENCE_BONDS > ELECTRONS_IN_VALENCE_BONDS
        || ELECTRONS_IN_VALENCE_BONDS > VALENCE_BONDS * 3
      )
      {
        // none, or undefine valence bonds
        return m_ValenceBonds[ 0]( 0);
      }

      // return the defined vector
      return m_ValenceBonds[ VALENCE_BONDS]( ELECTRONS_IN_VALENCE_BONDS - VALENCE_BONDS);
    }

    const ConstitutionalBondTypes &GetConstitutionalBondTypes()
    {
      return ConstitutionalBondTypes::GetEnums();
    }

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< chemistry::ConstitutionalBondTypeData, chemistry::ConstitutionalBondTypes>;

  } // namespace util
} // namespace bcl

