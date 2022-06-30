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

#ifndef BCL_CHEMISTRY_ATOM_COMPLETE_H_
#define BCL_CHEMISTRY_ATOM_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_bond_conformational.h"
#include "bcl_chemistry_configurational_bond_type_data.h"

// external includes - sorted alphabetically
#include <cstddef>

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomComplete
    //! @brief This class contains atom conformation information read directly from ISTREAM
    //! @details stores atom configuration, 3D coordinates of atoms and conformation of bonds attached to atom.
    //!
    //! @see @link example_chemistry_atom_complete.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Jan 11, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomComplete :
      public AtomConformationalInterface
    {
    private:

    /////////////
    // friends //
    /////////////

      // atom vector is the only class that can create or contain AtomComplete
      friend class AtomVector< AtomComplete>;   //!< changes the atom/bond types during standardization
      friend class AtomsCompleteStandardizer;   //!< changes the atom/bond types during standardization
      friend class AtomPropertiesStereocenters; //!< changes chirality
      friend class StereocentersHandler;        //!< changes chirality
      friend class BondIsometryHandler;         //!< changes bond isometry
      friend class MutateChirality;             //!< changes chirality and bond isometry
      friend class MutateBondLengths;           //!< changes bond lengths
      friend class FragmentAssemble;            //!< fixes bond orders (usually aromatic bond orders) during assembly
      friend class FragmentComplete;            //!< changes atom positions
      friend class MergeFragmentComplete;       //!< changes atom positions
      friend class FragmentMutateRingSwap;      //!< changes connectivity
      friend class FragmentMutateAlchemy;       //!< changes atom types and bond types to alter the chemical identity
      friend class FragmentMapConformer;        //!< maps positions of atoms from one fragment to another

    //////////////
    // typedefs //
    //////////////

      // these typedefs are needed so that AtomVector can be templated on a single type (AtomComplete), rather than
      // having to pass all these types in as well
      typedef BondConformational                t_Bond;      //!< Type of bonds this class contains
      typedef ConfigurationalBondTypeData::Data t_BondData;  //!< enum of data that can be obtained from a bond
      typedef AtomConformationalInterface       t_Interface; //!< atom interface class used by the bond

    //////////
    // data //
    //////////

      AtomType                             m_AtomType;    //!< Atom type of atom
      Chirality                            m_Chirality;   //!< chirality of atom
      linal::Vector3D                      m_Coordinates; //!< 3D coordinates of atom
      storage::Vector< BondConformational> m_Bonds;       //!< bonds of the object

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AtomComplete();

      //! @brief constructor from initializer, needed by atom vector
      //! @param ATOM_INFO constitution, configuration, and conformational info for this atom
      AtomComplete( const sdf::AtomInfo &ATOM_INFO);

      //! virtual copy constructor
      AtomComplete *Clone() const;

      //! @brief copy constructor
      //! @param ATOM atom conformation object of atom of interest
      AtomComplete( const AtomComplete &ATOM);

      //! @brief assignment constructor
      //! @param ATOM atom conformation object of atom of interest
      AtomComplete &operator =( const AtomComplete &ATOM);

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns bonds
      //! @return the bond conformations of the atom
      const storage::Vector< BondConformational> &GetBonds() const;

      //! @brief returns chirality
      //! @return reference of chirality attribute
      const Chirality &GetChirality() const;

      //! @brief returns atom type
      //! @return a reference of stom type data attribute
      const AtomType &GetAtomType() const;

      //! @brief returns coordinates
      //! @return a reference to the position vector
      const linal::Vector3D &GetPosition() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief get the bond type of a bond between two atoms
      //! @param ATOM the atom to which bond has to be reset
      ConfigurationalBondType GetBondTypeTo( const AtomConformationalInterface &ATOM) const;

      //! @brief set position
      //! @params POSITION set the position of atom
      void SetPosition( const linal::Vector3D &POSITION);

      //! @brief set atom type of the atom
      //! @param ATOM_TYPE the desired atom type
      void SetAtomType( const AtomType &ATOM_TYPE);

      //! @brief set charge on atom
      //! @param CHARGE the charge on atom
      void SetCharge( const short &CHARGE);

      //! @brief set chirality
      //! @params CHIRALITY chirality attribute to be set
      void SetChirality( const Chirality &CHIRALITY);

    private:

      //! @brief set bonds
      //! @param BONDS the bonds that the atom connects
      void SetBonds( const storage::Vector< BondConformational> &BONDS);

      //! @brief reset bond between atoms
      //! @param ATOM the atom to which bond has to be reset
      //! @param BOND_TYPE set the new bond as the given type
      void SetBondTypeMonoDirectional
      (
        const AtomConformationalInterface &ATOM,
        const ConfigurationalBondType &BOND_TYPE
      );

      //! @brief reset bond between atoms
      //! @param ATOM the atom to which bond has to be reset
      //! @param BOND_TYPE set the new bond as the given type
      void SetBondTypeTo( AtomComplete &ATOM, const ConfigurationalBondType &BOND_TYPE);

      //! @brief set bond between unbonded atoms
      //! @param ATOM the target atom with which bond has to be made
      //! @param BOND set the bond as the given type
      //! @param REPLACE_OPPOSITE_ATOM_VALENCE true if the the target atom connectivity has to be updated else false
      void ReplaceValenceWithBond
      (
        AtomComplete &ATOM,
        const ConfigurationalBondType &BOND,
        const bool &REPLACE_OPPOSITE_ATOM_VALENCE = true
      );

      //! @brief remove bond between bonded atoms
      //! @param ATOM the target atom with which bond has to be removed
      //! @param REMOVE_RETURN_BOND true if the the target atom connectivity has to be updated else false
      void AddValenceBondByRemovingBondTo( AtomComplete &ATOM, const bool &REMOVE_RETURN_BOND = true);

      //! @brief copy the atom and bonds using the specified difference in pointer address
      //! @param ATOM the atom to copy
      //! @param DIFF difference in begin of vector that atoms are located in (in bytes)
      void Copy( const AtomComplete &ATOM, const ptrdiff_t &DIFF);

    }; // class AtomComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_COMPLETE_H_
