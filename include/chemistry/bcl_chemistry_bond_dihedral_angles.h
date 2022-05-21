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

#ifndef BCL_CHEMISTRY_BOND_DIHEDRAL_ANGLES_H_
#define BCL_CHEMISTRY_BOND_DIHEDRAL_ANGLES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_type_data.h"
#include "bcl_chemistry_configurational_bond_types.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondDihedralAngles
    //! @brief calculates and stores bond length information
    //! Determines bond lengths of atoms based on element type & bond type (Single/Double/Triple/Aromatic)
    //! by using the covalent radius around an atom of that type
    //!
    //! @see @link example_chemistry_bond_dihedral_angles.cpp @endlink
    //! @author mendenjl
    //! @date Sep 17, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondDihedralAngles
    {

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief estimate the standard deviation for 0 degree-centered bond for a dihedral given only the two atoms
      //! @param ATOM_B, ATOM_C Atoms in the conformation that are bonded
      //! @return estimated standard deviation for the dihedral about this bond
      static double GetEstimatedStdForDihedralBondAngleBin
      (
        const AtomConformationalInterface &ATOM_B,
        const AtomConformationalInterface &ATOM_C
      );

      //! @brief Get the average and standard deviation of a given dihedral angle
      //! @param ATOM_TYPE)A,ATOM_TYPE_B,ATOM_TYPE_C, ATOM_TYPE_D the atom types in the bond
      //! @param BOND_TYPE_AB, BOND_TYPE_BC, BOND_TYPE_CD the bond types in the bond
      //! @return return the estimated ave/std dihedral deviation for a centered bin w/ 30 degree window
      //!         Values will be undefined (nan) if unavailable for the diehdral
      static storage::Pair< double, double> GetAveStdDihedralDeviation
      (
        const AtomType &ATOM_TYPE_A,
        const ConfigurationalBondType &BOND_TYPE_AB,
        const AtomType &ATOM_TYPE_B,
        const ConfigurationalBondType &BOND_TYPE_BC,
        const AtomType &ATOM_TYPE_C,
        const ConfigurationalBondType &BOND_TYPE_CD,
        const AtomType &ATOM_TYPE_D
      );

    }; // class BondDihedralAngles

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_DIHEDRAL_ANGLES_H_
