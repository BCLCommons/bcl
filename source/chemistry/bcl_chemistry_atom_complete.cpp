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
#include "chemistry/bcl_chemistry_atom_complete.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomComplete::AtomComplete()
    {
    }

    //! @brief constructor from initializer, needed by atom vector
    //! @param ATOM_INFO constitution, configuration, and conformational info for this atom
    AtomComplete::AtomComplete( const sdf::AtomInfo &ATOM_INFO) :
      m_AtomType( ATOM_INFO.GetAtomType()),
      m_Chirality( ATOM_INFO.GetChirality()),
      m_Coordinates( ATOM_INFO.GetCoordinates())
    {
    }

    //! virtual copy constructor
    AtomComplete *AtomComplete::Clone() const
    {
      return new AtomComplete( *this);
    }

    //! @brief copy constructor
    //! @param ATOM atom conformation object of atom of interest
    AtomComplete::AtomComplete( const AtomComplete &ATOM)
    {
      *this = ATOM;
    }

    //! @brief assignment constructor
    //! @param ATOM atom conformation object of atom of interest
    AtomComplete &AtomComplete::operator =( const AtomComplete &ATOM)
    {
      if( this != &ATOM)
      {
        // this class can only be constructed by AtomVector
        // likewise, we know the atom is in a vector, so the addresses of the other atoms in the bonds must be offset
        // by the difference between the address of this atom and the address of ATOM
        // to get this difference in bytes, it is necessary to first cast the this pointer and atom pointer to char *
        Copy( ATOM, ( char *)( this) - ( const char *)( &ATOM));
      }

      return *this;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns bonds
    //! @return the bond conformations of the atom
    const storage::Vector< BondConformational> &AtomComplete::GetBonds() const
    {
      return m_Bonds;
    }

    //! @brief returns chirality
    //! @return reference of chirality attribute
    const Chirality &AtomComplete::GetChirality() const
    {
      return m_Chirality;
    }

    //! @brief returns atom type
    //! @return a reference of stom type data attribute
    const AtomType &AtomComplete::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns coordinates
    //! @return a reference to the position vector
    const linal::Vector3D &AtomComplete::GetPosition() const
    {
      return m_Coordinates;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomComplete::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Bonds, ISTREAM);
      io::Serialize::Read( m_Coordinates, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Bonds, OSTREAM, INDENT);
      io::Serialize::Write( m_Coordinates, OSTREAM, INDENT);
      // end
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the bond type of a bond between two atoms
    //! @param ATOM the atom to which bond has to be reset
    ConfigurationalBondType AtomComplete::GetBondTypeTo( const AtomConformationalInterface &ATOM) const
    {
      storage::Vector< BondConformational>::const_iterator itr( FindBondTo( ATOM));
      if( itr != m_Bonds.End()) // bond found
      {
        // return the type
        return itr->GetBondType();
      }

      // bond not found
      return GetConfigurationalBondTypes().e_Undefined;
    }

    //! @brief set bonds
    //! @param BONDS the bonds that the atom connects
    void AtomComplete::SetBonds( const storage::Vector< BondConformational> &BONDS)
    {
      m_Bonds = BONDS;
    }

    //! @brief set atom type of the atom
    //! @param ATOM_TYPE the desired atom type
    void AtomComplete::SetAtomType( const AtomType &ATOM_TYPE)
    {
      m_AtomType = ATOM_TYPE;
    }

    //! @brief set charge on atom
    //! @param CHARGE the charge on atom
    void AtomComplete::SetCharge( const short &CHARGE)
    {
      m_AtomType = GetAtomTypes().GetAtomType( GetElementType(), CHARGE);
    }

    //! @brief set chirality
    //! @params CHIRALITY chirality attribute to be set
    void AtomComplete::SetChirality( const Chirality &CHIRALITY)
    {
      m_Chirality = CHIRALITY;
    }

    //! @brief set position
    //! @params POSITION set the position of atom
    void AtomComplete::SetPosition( const linal::Vector3D &POSITION)
    {
      m_Coordinates = POSITION;
    }

    //! @brief reset bond between atoms
    //! @param ATOM the atom to which bond has to be reset
    //! @param BOND_TYPE set the new bond as the given type
    void AtomComplete::SetBondTypeMonoDirectional
    (
      const AtomConformationalInterface &ATOM,
      const ConfigurationalBondType &BOND_TYPE
    )
    {
      const AtomConformationalInterface *target_atom_ptr( &ATOM);

      // iterate over bonds
      for
      (
        storage::Vector< BondConformational>::iterator
          itr_bond( m_Bonds.Begin()), itr_bond_end( m_Bonds.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // find the target atom
        if( &itr_bond->GetTargetAtom() == target_atom_ptr)
        {
          // set the bond type for this atom
          itr_bond->SetBondType( BOND_TYPE);
          // return
          return;
        }
      }

      BCL_Exit( "Tried to change bond type; but no bond existed between the atoms!", 1);
    }

    //! @brief reset bond between atoms
    //! @param ATOM the atom to which bond has to be reset
    //! @param BOND_TYPE set the new bond as the given type
    void AtomComplete::SetBondTypeTo
    (
      AtomComplete &ATOM,
      const ConfigurationalBondType &BOND_TYPE
    )
    {
      const AtomConformationalInterface *target_atom_ptr( &ATOM);

      // iterate over bonds
      for
      (
        storage::Vector< BondConformational>::iterator
          itr_bond( m_Bonds.Begin()), itr_bond_end( m_Bonds.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // find the target atom
        if( &itr_bond->GetTargetAtom() == target_atom_ptr)
        {
          // if the bond type isn't already set properly, set it for this atom and the reverse bond on the other atom
          if( itr_bond->GetBondType() != BOND_TYPE)
          {
            itr_bond->SetBondType( BOND_TYPE);
            ATOM.SetBondTypeMonoDirectional( *this, BOND_TYPE);
          }
          return;
        }
      }

      BCL_Exit( "Tried to change bond type; but no bond existed between the atoms!", 1);
    }

    //! @brief set bond between unbonded atoms
    //! @param ATOM the target atom with which bond has to be made
    //! @param BOND set the bond as the given type
    //! @param REPLACE_OPPOSITE_ATOM_VALENCE true if the the target atom connectivity has to be updated else false
    void AtomComplete::ReplaceValenceWithBond
    (
      AtomComplete &ATOM,
      const ConfigurationalBondType &BOND,
      const bool &REPLACE_OPPOSITE_ATOM_VALENCE
    )
    {
      bool have_bond_to_atom( false);
      const AtomConformationalInterface *target_atom_ptr( &ATOM);

      // first ensure that this atom does not already have a bond to atom
      for
      (
        storage::Vector< BondConformational>::iterator
          itr_bond( m_Bonds.Begin()), itr_bond_end( m_Bonds.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // find the target atom
        if( &itr_bond->GetTargetAtom() == target_atom_ptr)
        {
          have_bond_to_atom = true;
          break;
        }
      }
      BCL_Assert( !have_bond_to_atom, "Cannot add a bond to an atom that was already bonded to");

      const size_t valence_pos
      (
        GetValenceBonds().Find( BOND->GetConstitutionalBondType())
      );

      BCL_Assert
      (
        valence_pos < GetValenceBonds().GetSize(),
        "ReplaceValenceWithBond called with bond type that was not a valid valence"
      );

      m_Bonds.PushBack( BondConformational( ATOM, BOND));

      if( REPLACE_OPPOSITE_ATOM_VALENCE)
      {
        ATOM.ReplaceValenceWithBond( *this, BOND, false);
      }
    }

    //! @brief remove bond between bonded atoms
    //! @param ATOM the target atom with which bond has to be removed
    //! @param REMOVE_RETURN_BOND true if the the target atom connectivity has to be updated else false
    void AtomComplete::AddValenceBondByRemovingBondTo( AtomComplete &ATOM, const bool &REMOVE_RETURN_BOND)
    {
      const AtomConformationalInterface *target_atom_ptr( &ATOM);

      // keep track of whether the bond was found and could be set
      bool bond_was_set( false);

      // iterate over bonds
      for
      (
        storage::Vector< BondConformational>::iterator
          itr_bond( m_Bonds.Begin()), itr_bond_end( m_Bonds.End());
        itr_bond != itr_bond_end;
      )
      {
        // find the target atom
        if( &itr_bond->GetTargetAtom() == target_atom_ptr)
        {
          if( REMOVE_RETURN_BOND)
          {
            ATOM.AddValenceBondByRemovingBondTo( *this, false);
          }
          itr_bond = storage::Vector< BondConformational>::iterator( m_Bonds.Remove( itr_bond));
          bond_was_set = true;
          break;
        }
        ++itr_bond;
      }
      BCL_Assert( bond_was_set, "Tried to set a bond to an atom that this atom was not bonded to!");
    }

    //! @brief copy the atom and bonds using the specified difference in pointer address
    //! @param ATOM the atom to copy
    //! @param DIFF difference in begin of vector that atoms are located in
    void AtomComplete::Copy( const AtomComplete &ATOM, const ptrdiff_t &DIFF)
    {
      m_AtomType = ATOM.m_AtomType;
      m_Chirality = ATOM.m_Chirality;
      m_Coordinates = ATOM.m_Coordinates;

      const size_t number_bonds( ATOM.m_Bonds.GetSize());
      m_Bonds.Reset();
      m_Bonds.AllocateMemory( number_bonds);

      // copy the bonds, changing the targeted atoms by DIFF
      for( size_t bond_number( 0); bond_number < number_bonds; ++bond_number)
      {
        const BondConformational &bond_to_copy( ATOM.m_Bonds( bond_number));
        m_Bonds.PushBack
        (
          BondConformational
          (
            *( const AtomConformationalInterface *)( ( const char *)&bond_to_copy.GetTargetAtom() + DIFF),
            bond_to_copy.GetBondType()
          )
        );
      }
    }

  } // namespace chemistry
} // namespace bcl

