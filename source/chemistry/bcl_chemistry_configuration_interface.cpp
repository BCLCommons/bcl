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
#include "chemistry/bcl_chemistry_configuration_interface.h"
#include "sdf/bcl_sdf_mdl_handler.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief get the number of valences for the entire molecule or fragment
    //! @return the number of valences on the entire molecule or fragment
    size_t ConfigurationInterface::GetNumberValences() const
    {
      size_t number_valences( 0);

      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr_atoms( GetAtomsIterator());
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
    std::string ConfigurationInterface::GetAtomTypesString() const
    {
      std::string atom_types;

      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr_atoms( GetAtomsIterator());
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
    std::string ConfigurationInterface::GetBondTypesString() const
    {
      std::ostringstream bond_types;
      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConfigurational> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConfigurational>::const_iterator
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

    //! @brief get a string containing stereochemistry of atoms
    //! @brief return a string containing the stereochemistry of the atoms of this molecule as a string
    std::string ConfigurationInterface::GetChiralityString() const
    {
      std::ostringstream stereochem;
      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        stereochem << GetChiralityName( itr->GetChirality()) << ' ';
      }

      return stereochem.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ConfigurationInterface::WriteMDL( std::ostream &OSTREAM) const
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
