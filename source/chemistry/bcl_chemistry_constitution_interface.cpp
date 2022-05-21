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
#include "chemistry/bcl_chemistry_constitution_interface.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief return a number of hydrogens in the conformation
    //! @return the number of hydrogens in the conformation
    size_t ConstitutionInterface::GetNumberHydrogens() const
    {
      size_t number_h( 0);
      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        number_h += itr_atoms->GetElementType() == GetElementTypes().e_Hydrogen;
      }
      return number_h;
    }

    //! @brief get the number of valences for the entire molecule or fragment
    //! @return the number of valences on the entire molecule or fragment
    size_t ConstitutionInterface::GetNumberValences() const
    {
      size_t number_valences( 0);

      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        number_valences += itr_atoms->GetNumberValenceBonds();
      }
      return number_valences;
    }

    //! @brief get a string with each atom type in the molecule
    //! @return a string with each atom type in the molecule
    std::string ConstitutionInterface::GetAtomTypesString() const
    {
      std::string atom_types;

      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        atom_types += itr_atoms->GetAtomType().GetName() + " ";
      }
      return atom_types;
    }

    //! @brief get a string with each bond type in the molecule
    //! @return a string with each bond type in the molecule
    std::string ConstitutionInterface::GetBondTypesString() const
    {
      std::ostringstream bond_types;
      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConstitutional> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConstitutional>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( &*itr < &itr_connected->GetTargetAtom())
          {
            bond_types << itr_connected->GetBondType().GetName() << ' ';
          }
        }
      }

      return bond_types.str();
    }

    //! @brief returns a string with the sum formula
    //! @return a string with the sum formula
    std::string ConstitutionInterface::GetSumFormula() const
    {
      // vector of positive integers that count the number of atoms that appear in the molecule
      storage::Vector< size_t> formula( GetElementTypes().GetEnumCount(), size_t( 0));

      // number of unknown elements
      size_t number_unknown( 0);

      // iterate through all atoms of the molecule and count their appearance
      for( iterate::Generic< const AtomConstitutionalInterface> itr( GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( itr->GetElementType().IsDefined())
        {
          ++formula( itr->GetElementType().GetIndex());
        }
        else
        {
          ++number_unknown;
        }
      }

      // stringstream for the final sum formula
      std::stringstream sum_formula( "");

      // sum formula begins always with a carbon if there is one
      if( formula( GetElementTypes().e_Carbon) != 0)
      {
        sum_formula << GetElementTypes().e_Carbon->GetChemicalSymbol();

        // only print the number of atoms if it is greater than 1
        if( formula( GetElementTypes().e_Carbon) > 1)
        {
          sum_formula << formula( GetElementTypes().e_Carbon);
        }
      }

      // assemble sum formula with delimiter "_"
      for
      (
        ElementTypes::const_iterator itr( GetElementTypes().Begin()), itr_end( GetElementTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->IsDefined() && formula( *itr) != 0 && ( *itr) != GetElementTypes().e_Carbon)
        {
          sum_formula << ( *itr)->GetChemicalSymbol();

          // only print the number of atoms if it is greater than 1
          if( formula( *itr) > 1)
          {
            sum_formula << formula( *itr);
          }
        }
      }

      if( number_unknown != 0)
      {
        sum_formula << "X";

        // only print the number of atoms if it is greater than 1
        if( number_unknown > 1)
        {
          sum_formula << number_unknown;
        }
      }

      // return sum formula
      return sum_formula.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ConstitutionInterface::WriteMDL( std::ostream &OSTREAM) const
    {
      return
        sdf::MdlHandler::WriteToSDF
        (
          OSTREAM,
          GetName(),
          GetAtomInfo(),
          GetBondInfo(),
          GetStoredProperties().GetMDLProperties()
        );
    }

  } // namespace chemistry
} // namespace bcl
