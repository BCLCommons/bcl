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
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // enums //
  ///////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

    //! @brief return the number of valence bonds (this is faster than asking for valence bonds)
    //! @return a reference the vector of valence bonds
    size_t AtomConformationalInterface::GetNumberValenceBonds() const
    {
      // bad or unknown atom type, so return undefined valence bonds-
      if( GetAtomType()->GetNumberBonds() < GetBonds().GetSize())
      {
        return util::GetUndefined< size_t>();
      }
      return GetAtomType()->GetNumberBonds() - GetBonds().GetSize();
    }

    //! @brief return the number of electrons in valence bonds (this is faster than asking for valence bonds)
    //! @return a reference the vector of valence bonds
    size_t AtomConformationalInterface::GetNumberElectronsInValenceBonds() const
    {
      // bad or unknown atom type, so return undefined valence bond e-
      if( GetAtomType()->GetNumberBonds() < GetBonds().GetSize())
      {
        return util::GetUndefined< size_t>();
      }
      else if( GetAtomType()->GetNumberBonds() == GetBonds().GetSize())
      {
        return size_t( 0);
      }

      size_t number_valence_electrons( 2 * GetAtomType()->GetNumberElectronsInBonds());

      // subtract electrons in all bonds
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( GetBonds().Begin()), itr_end( GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        number_valence_electrons -= itr->GetBondType()->GetNumberOfElectrons();
      }

      // return # of electrons over 2, since e- are shared
      return number_valence_electrons / 2;
    }

    //! @brief return the number of valence bonds of a particular bond order
    //! @param BOND_ORDER desired bond order for the valences
    //! @return a the number of valence bonds of a particular bond order
    size_t AtomConformationalInterface::GetNumberofValenceBondsWithOrder
    (
      const size_t BOND_ORDER
    ) const
    {
      size_t number_bonds( 0);
      for
      (
        storage::Vector< ConstitutionalBondType>::const_iterator
          itr( GetValenceBonds().Begin()), itr_end( GetValenceBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetBondData( ConstitutionalBondTypeData::e_BondOrder) == BOND_ORDER)
        {
          // single bond to be saturated
          ++number_bonds;
        }
      }
      // return # of electrons over 2, since e- are shared
      return number_bonds;
    }

    //! @brief return the number of covalently bonded hydrogens
    //! @return the number of covalently bonded hydrogens
    size_t AtomConformationalInterface::GetNumberCovalentlyBoundHydrogens() const
    {
      size_t h_count( 0);

      // subtract electrons in all bonds
      for
      (
        storage::Vector< BondConformational>::const_iterator itr( GetBonds().Begin()), itr_end( GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
        {
          ++h_count;
        }
      }
      return h_count;
    }

    //! @brief return the number of bonds with a particular property (e.g. is in ring, is aromatic, is isometric)
    //! @param DATA the data to retrieve for each bond
    //! @param VALUE the value to count, e.g. if ConfigurationalBondType->GetData( DATA) == VALUE, ++return
    //! @return the number of bonds with the property of interest
    size_t AtomConformationalInterface::CountNonValenceBondsWithProperty
    (
      const ConfigurationalBondTypeData::Data &DATA,
      const size_t &VALUE
    ) const
    {
      size_t number_bonds_satisfy_query( 0);

      // subtract electrons in all bonds
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr( GetBonds().Begin()), itr_end( GetBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        number_bonds_satisfy_query += size_t( itr->GetBondType()->GetBondData( DATA) == VALUE);
      }

      return number_bonds_satisfy_query;
    }

    //! @brief return a reference to the list of valence bonds
    //! @return a reference to a vector of valence bonds
    const storage::Vector< ConstitutionalBondType> &AtomConformationalInterface::GetValenceBonds() const
    {
      return GetConstitutionalBondTypes().GetValenceBonds( GetNumberValenceBonds(), GetNumberElectronsInValenceBonds());
    }

    //! @brief convert the constitutional interface into an sdf::AtomInfo
    //! @return the constitutional interface converted into an sdf::AtomInfo
    sdf::AtomInfo AtomConformationalInterface::GetAtomInfo() const
    {
      return sdf::AtomInfo( GetAtomType(), GetChirality(), GetPosition(), GetNumberValenceBonds() > 0);
    }

    //! @brief find the bond from this atom to to ATOM
    //! @param ATOM the atom to locate
    //! @return an iterator to the bond connected to ATOM
    storage::Vector< BondConformational>::const_iterator
      AtomConformationalInterface::FindBondTo( const AtomConformationalInterface &ATOM) const
    {
      const AtomConformationalInterface *target_atom_ptr( &ATOM);

      storage::Vector< BondConformational>::const_iterator
        itr_bond( GetBonds().Begin()), itr_bond_end( GetBonds().End());
      // iterate over bonds until the target atom is found
      for( ; itr_bond != itr_bond_end && &itr_bond->GetTargetAtom() != target_atom_ptr; ++itr_bond)
      {
      }

      // return the iterator
      return itr_bond;
    }

  } // namespace chemistry
} // namespace bcl
