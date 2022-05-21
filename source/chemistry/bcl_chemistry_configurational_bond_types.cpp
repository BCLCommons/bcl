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
#include "chemistry/bcl_chemistry_configurational_bond_types.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitutional_bond_types.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct all bond types
    ConfigurationalBondTypes::ConfigurationalBondTypes() :
      util::Enumerate< ConfigurationalBondTypeData, ConfigurationalBondTypes>( false),
      e_NonConjugatedSingleBond(       AddBond( GetConstitutionalBondTypes().e_NonConjugatedSingleBond      , e_NonIsometric)),
      e_ConjugatedSingleBond(          AddBond( GetConstitutionalBondTypes().e_ConjugatedSingleBond         , e_NonIsometric)),
      e_ConjugatedDoubleBond(          AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBond         , e_NonIsometric)),
      e_ConjugatedTripleBond(          AddBond( GetConstitutionalBondTypes().e_ConjugatedTripleBond         , e_NonIsometric)),
      e_AromaticSingleBond(            AddBond( GetConstitutionalBondTypes().e_AromaticSingleBond           , e_NonIsometric)),
      e_AromaticDoubleBond(            AddBond( GetConstitutionalBondTypes().e_AromaticDoubleBond           , e_NonIsometric)),
      e_AromaticTripleBond(            AddBond( GetConstitutionalBondTypes().e_AromaticTripleBond           , e_NonIsometric)),
      e_NonConjugatedSingleBondInRing( AddBond( GetConstitutionalBondTypes().e_NonConjugatedSingleBondInRing, e_NonIsometric)),
      e_ConjugatedSingleBondInRing(    AddBond( GetConstitutionalBondTypes().e_ConjugatedSingleBondInRing   , e_NonIsometric)),
      e_ConjugatedDoubleBondInRing(    AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBondInRing   , e_NonIsometric)),
      e_ConjugatedTripleBondInRing(    AddBond( GetConstitutionalBondTypes().e_ConjugatedTripleBondInRing   , e_NonIsometric)),
      e_ConjugatedBond(                AddBond( GetConstitutionalBondTypes().e_ConjugatedBond               , e_NonIsometric)),
      e_AromaticBond(                  AddBond( GetConstitutionalBondTypes().e_AromaticBond                 , e_NonIsometric)),
      e_ConjugatedBondInRing(          AddBond( GetConstitutionalBondTypes().e_ConjugatedBondInRing         , e_NonIsometric)),
      e_AmideSingleBond(               AddBond( GetConstitutionalBondTypes().e_AmideSingleBond              , e_NonIsometric)),
      e_ConjugatedDoubleBond_X(        AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBond         , e_UnknownIsometry)),
      e_ConjugatedDoubleBond_E(        AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBond         , e_EIsometry)),
      e_ConjugatedDoubleBond_Z(        AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBond         , e_ZIsometry)),
      e_ConjugatedDoubleBondInRing_X(  AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBondInRing   , e_UnknownIsometry)),
      e_ConjugatedDoubleBondInRing_E(  AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBondInRing   , e_EIsometry)),
      e_ConjugatedDoubleBondInRing_Z(  AddBond( GetConstitutionalBondTypes().e_ConjugatedDoubleBondInRing   , e_ZIsometry))
    {
      // all enums have been created

      // initialize isometries of all consitutional bond types
      ConstitutionalBondType undefined_constitutional_bond_type( e_Undefined->GetConstitutionalBondType());
      undefined_constitutional_bond_type->InitializeIsometries( e_Undefined);
      for( size_t i( 0), number_bond_types( GetEnumCount()); i < number_bond_types; ++i)
      {
        ConfigurationalBondType &bond_type( GetEnumFromIndex( i));
        ConstitutionalBondType constitutional_bond_type( bond_type->GetConstitutionalBondType());
        constitutional_bond_type->InitializeIsometries( e_Undefined);
      }

      // attach all configurational bond types to the related constitutional types
      for( size_t i( 0), number_bond_types( GetEnumCount()); i < number_bond_types; ++i)
      {
        ConfigurationalBondType &bond_type( GetEnumFromIndex( i));
        ConstitutionalBondType constitutional_bond_type( bond_type->GetConstitutionalBondType());
        constitutional_bond_type->SetConfigurationalBondType( bond_type);
      }

      // attach each configurational bond type to its related types
      // This must be done after the linking of the constitutional types to the corresponding configurational types,
      // because ConfigurationalBondTypeData::AttachToRelatedTypes depends on the link between constitutional bond types
      // and configurational bond types
      e_Undefined->AttachToRelatedBondTypes( e_Undefined);
      for( size_t i( 0), number_bond_types( GetEnumCount()); i < number_bond_types; ++i)
      {
        ConfigurationalBondType &bond_type( GetEnumFromIndex( i));
        bond_type->AttachToRelatedBondTypes( bond_type);
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! Add a bond type to the enum
    //! @param BASE_BOND_TYPE the constitutional bond type
    //! @param ISOMETRY isometry of the bond type
    ConfigurationalBondType &ConfigurationalBondTypes::AddBond
    (
      ConstitutionalBondType BASE_BOND_TYPE,
      const BondIsometry &ISOMETRY
    )
    {
      // create the suffix of the name based on the suffix
      const std::string isometry_name( ISOMETRY == e_NonIsometric ? std::string() : "_" + GetIsometryName( ISOMETRY));

      return AddEnum( BASE_BOND_TYPE.GetName() + isometry_name, ConfigurationalBondTypeData( BASE_BOND_TYPE, ISOMETRY));
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConfigurationalBondTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    const ConfigurationalBondTypes &GetConfigurationalBondTypes()
    {
      return ConfigurationalBondTypes::GetEnums();
    }

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< chemistry::ConfigurationalBondTypeData, chemistry::ConfigurationalBondTypes>;

  } // namespace util
} // namespace bcl

