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

#ifndef BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_INTERFACE_H_
#define BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_types.h"
#include "bcl_chemistry_bond_constitutional.h"
#include "bcl_chemistry_constitutional_bond_types.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomConstitutionalInterface
    //! @brief Basic methods implemented in the AtomConstitutional class
    //! @details Contains methods that return the AtomType, BondsTypes, valence bonds, element type and charge of atom
    //!
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomConstitutionalInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns a reference of atom type data attribute
      //! @return a reference of atom type data attribute
      virtual const AtomType &GetAtomType() const = 0;

      //! @brief returns bonds that the atom is involved in
      //! @return bonds that the atom is involved in
      virtual const storage::Vector< BondConstitutional> &GetBonds() const = 0;

      //! @brief return the number of valence bonds (this is faster than asking for valence bonds)
      //! @return a reference the vector of valence bonds
      virtual size_t GetNumberValenceBonds() const;

      //! @brief return the number of electrons in valence bonds (this is faster than asking for valence bonds)
      //! @return a reference the vector of valence bonds
      virtual size_t GetNumberElectronsInValenceBonds() const;

      //! @brief return a reference to the list of valence bonds
      //! @return a reference to a vector of valence bonds
      virtual const storage::Vector< ConstitutionalBondType> &GetValenceBonds() const;

      //! @brief return a reference to the element type of the atom
      //! @return the element type of the atom
      const ElementType &GetElementType() const
      {
        return GetAtomType()->GetElementType();
      }

      //! @brief get the charge
      //! @return the charge of the atom
      const short &GetCharge() const
      {
        return GetAtomType()->GetFormalCharge();
      }

      //! @brief convert the constitutional interface into an sdf::AtomInfo
      //! @return the constitutional interface converted into an sdf::AtomInfo
      sdf::AtomInfo GetAtomInfo() const;

    ////////////////
    // operations //
    ////////////////

    }; // class AtomConstitutionalInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_INTERFACE_H_ 
