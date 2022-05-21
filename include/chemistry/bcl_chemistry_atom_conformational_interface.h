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

#ifndef BCL_CHEMISTRY_ATOM_CONFORMATIONAL_INTERFACE_H_
#define BCL_CHEMISTRY_ATOM_CONFORMATIONAL_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
//#include "bcl_chemistry_atom_configurational_interface.h"
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "coord/bcl_coord_orientation_interface.h"
#include "sdf/bcl_sdf_atom_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomConformationalInterface
    //! @brief Basic methods implemented in the AtomConformational class
    //! @details Contains methods that return the bond information, valence bond information, atom conformation
    //! attributes of atom
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl
    //! @date Dec 05, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomConformationalInterface :
      public coord::OrientationInterface
    {

    /////////////
    // friends //
    /////////////

      friend class ConformationInterface; //!< Calls set position, but only to perform translations, rotations, etc.

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns bonds that the atom is involved in
      //! @return bonds that the atom is involved in
      virtual const storage::Vector< BondConformational> &GetBonds() const = 0;

      //! @brief return a reference to the list of valence bonds
      //! @return a reference to a vector of valence bonds
      virtual const storage::Vector< ConstitutionalBondType> &GetValenceBonds() const;

      //! @brief returns chirality attribute of atom
      //! @return chirality attribute of atom
      virtual const Chirality &GetChirality() const = 0;

      //! @brief return the number of valence bonds (this is faster than asking for valence bonds)
      //! @return # of valence bonds
      size_t GetNumberValenceBonds() const;

      //! @brief return the number of valence bonds of a particular bond order
      //! @param BOND_ORDER desired bond order for the valences
      //! @return a the number of valence bonds of a particular bond order
      size_t GetNumberofValenceBondsWithOrder( const size_t BOND_ORDER) const;

      //! @brief return the number of electrons in valence bonds (this is faster than asking for valence bonds)
      //! @return # of electrons in valence bonds
      size_t GetNumberElectronsInValenceBonds() const;

      //! @brief return the number of covalently bonded hydrogens
      //! @return the number of covalently bonded hydrogens
      size_t GetNumberCovalentlyBoundHydrogens() const;

      //! @brief return the number of bonds with a particular property (e.g. is in ring, is aromatic, is isometric)
      //! @param DATA the data to retrieve for each bond
      //! @param VALUE the value to count, e.g. if ConfigurationalBondType->GetData( DATA) == VALUE, ++return
      //! @return the number of bonds with the property of interest
      size_t CountNonValenceBondsWithProperty
      (
        const ConfigurationalBondTypeData::Data &DATA,
        const size_t &VALUE
      ) const;

      //! @brief convert the constitutional interface into an sdf::AtomInfo
      //! @return the constitutional interface converted into an sdf::AtomInfo
      sdf::AtomInfo GetAtomInfo() const;

      //! @brief find the bond from this atom to to ATOM
      //! @param ATOM the atom to locate
      //! @return an iterator to the bond connected to ATOM
      storage::Vector< BondConformational>::const_iterator FindBondTo( const AtomConformationalInterface &ATOM) const;

      //! @brief returns a reference of atom type data attribute
      //! @return a reference of atom type data attribute
      virtual const AtomType &GetAtomType() const = 0;

      //! @brief return a reference to the position vector
      //! @return a reference to the position vector
      virtual const linal::Vector3D &GetPosition() const = 0;

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

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const
      {
        return GetPosition();
      }

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief set the position
      //! @param POSITION the new coordinates
      virtual void SetPosition( const linal::Vector3D &POSITION) = 0;

    }; // class AtomConformationalInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CONFORMATIONAL_INTERFACE_H_
