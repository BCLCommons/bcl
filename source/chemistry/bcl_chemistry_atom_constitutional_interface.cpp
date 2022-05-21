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
#include "chemistry/bcl_chemistry_atom_constitutional_interface.h"

// includes from bcl - sorted alphabetically

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
    size_t AtomConstitutionalInterface::GetNumberValenceBonds() const
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
    size_t AtomConstitutionalInterface::GetNumberElectronsInValenceBonds() const
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
        storage::Vector< BondConstitutional>::const_iterator
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

    //! @brief return a reference to the list of valence bonds
    //! @return a reference to a vector of valence bonds
    const storage::Vector< ConstitutionalBondType> &AtomConstitutionalInterface::GetValenceBonds() const
    {
      return GetConstitutionalBondTypes().GetValenceBonds( GetNumberValenceBonds(), GetNumberElectronsInValenceBonds());
    }

    //! @brief convert the constitutional interface into an sdf::AtomInfo
    //! @return the constitutional interface converted into an sdf::AtomInfo
    sdf::AtomInfo AtomConstitutionalInterface::GetAtomInfo() const
    {
      return sdf::AtomInfo( GetAtomType(), e_UnknownChirality, linal::Vector3D( 0.0), GetNumberValenceBonds() > 0);
    }

  } // namespace chemistry
} // namespace bcl
