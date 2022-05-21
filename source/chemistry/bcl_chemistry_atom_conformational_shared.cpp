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
#include "chemistry/bcl_chemistry_atom_conformational_shared.h"

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
    AtomConformationalShared::AtomConformationalShared()
    {
    }

    //! @brief constructor from initializer, needed by atom vector
    //! @param ATOM_INFO all info about the atom, only coordinates will be used
    AtomConformationalShared::AtomConformationalShared( const sdf::AtomInfo &ATOM_INFO) :
      m_Configuration(),
      m_Coordinates( ATOM_INFO.GetCoordinates())
    {
    }

    //! virtual copy constructor
    AtomConformationalShared *AtomConformationalShared::Clone() const
    {
      return new AtomConformationalShared( *this);
    }

    //! @brief copy constructor
    //! @param ATOM atom conformation object of atom of interest
    AtomConformationalShared::AtomConformationalShared( const AtomConformationalShared &ATOM)
    {
      *this = ATOM;
    }

    //! @brief assignment constructor
    //! @param ATOM atom conformation object of atom of interest
    AtomConformationalShared &AtomConformationalShared::operator =( const AtomConformationalShared &ATOM)
    {
      if( this != &ATOM)
      {
        // this class can only be created by AtomVector
        // likewise, we know the atom is in a vector, so the addresses of the other atoms in the bonds must be offset
        // by the difference between the address of this atom and the address of ATOM
        Copy( ATOM, ( char *)( this) - ( const char *)( &ATOM));
      }

      return *this;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomConformationalShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns bonds
    //! @return the bond conformations of bonds connected to the atom
    const storage::Vector< BondConformational> &AtomConformationalShared::GetBonds() const
    {
      return m_Bonds;
    }

    //! @brief returns chirality
    //! @return reference of chirality attribute
    const Chirality &AtomConformationalShared::GetChirality() const
    {
      return m_Configuration->GetChirality();
    }

    //! @brief returns atom type
    //! @return a reference of stom type data attribute
    const AtomType &AtomConformationalShared::GetAtomType() const
    {
      return m_Configuration->GetAtomType();
    }

    //! @brief returns coordinates
    //! @return a reference to the position vector
    const linal::Vector3D &AtomConformationalShared::GetPosition() const
    {
      return m_Coordinates;
    }

    //! @brief return a simple pointer to atom configurational object
    //! @return pointer to atom configurational object
    const util::SiPtr< const AtomConfigurationalShared> AtomConformationalShared::GetConfiguration() const
    {
      return m_Configuration;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomConformationalShared::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Configuration, ISTREAM);
      io::Serialize::Read( m_Bonds, ISTREAM);
      io::Serialize::Read( m_Coordinates, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomConformationalShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Configuration, OSTREAM, INDENT);
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

    //! @brief set bonds
    //! @param BONDS the bonds that the atom connects
    void AtomConformationalShared::SetBonds( const storage::Vector< BondConformational> &BONDS)
    {
      m_Bonds = BONDS;
    }

    //! @brief copy the atom and bonds using the specified difference in pointer address
    //! @param ATOM the atom to copy
    //! @param DIFF difference in begin of vector that atoms are located in
    void AtomConformationalShared::Copy( const AtomConformationalShared &ATOM, const ptrdiff_t &DIFF)
    {
      m_Configuration = ATOM.m_Configuration;
      m_Coordinates = ATOM.m_Coordinates;

      const size_t number_bonds( ATOM.m_Bonds.GetSize());
      m_Bonds.Reset();
      m_Bonds.AllocateMemory( number_bonds);

      // copy the bonds, changing the target atoms by DIFF bytes
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

