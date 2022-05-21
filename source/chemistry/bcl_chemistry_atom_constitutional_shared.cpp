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
#include "chemistry/bcl_chemistry_atom_constitutional_shared.h"

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
    AtomConstitutionalShared::AtomConstitutionalShared()
    {
    }

    //! @brief constructor from initializer, needed by atom vector
    //! @param ATOM_INFO atom information
    AtomConstitutionalShared::AtomConstitutionalShared( const sdf::AtomInfo &ATOM_INFO) :
      m_AtomType( ATOM_INFO.GetAtomType()),
      m_Bonds()
    {
    }

    //! virtual copy constructor
    AtomConstitutionalShared *AtomConstitutionalShared::Clone() const
    {
      return new AtomConstitutionalShared( *this);
    }

    //! @brief copy constructor
    //! @param ATOM atom constitution object of atom of interest
    AtomConstitutionalShared::AtomConstitutionalShared( const AtomConstitutionalShared &ATOM)
    {
      *this = ATOM;
    }

    //! @brief copy constructor
    //! @param ATOM atom constitution object of atom of interest
    AtomConstitutionalShared &AtomConstitutionalShared::operator =( const AtomConstitutionalShared &ATOM)
    {
      if( this != &ATOM)
      {
        // this class can only be created by AtomVector
        // likewise, we know the atom is in a vector, so the addresses of the other atoms in the bonds must be offset
        // by the difference between the address of this atom and the address of ATOM
        Copy( ATOM, ( const char *)( this) - ( const char *)( &ATOM));
      }

      return *this;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomConstitutionalShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns atom attribute
    //! @return a reference to atom attribute
    const AtomType &AtomConstitutionalShared::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns bonds
    //! @return a reference to bond configuration objects
    const storage::Vector< BondConstitutional> &AtomConstitutionalShared::GetBonds() const
    {
      return m_Bonds;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomConstitutionalShared::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AtomType, ISTREAM);
      io::Serialize::Read( m_Bonds, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomConstitutionalShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_AtomType, OSTREAM, INDENT);
      io::Serialize::Write( m_Bonds, OSTREAM, INDENT);
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

    //! @brief copy the atom and bonds using the specified difference in pointer address
    //! @param ATOM the atom to copy
    //! @param DIFF difference in begin of vector that atoms are located in
    void AtomConstitutionalShared::Copy( const AtomConstitutionalShared &ATOM, const ptrdiff_t &DIFF)
    {
      m_AtomType = ATOM.m_AtomType;

      const size_t number_bonds( ATOM.m_Bonds.GetSize());
      m_Bonds.Reset();
      m_Bonds.AllocateMemory( number_bonds);

      // copy the bonds, changing the target atoms by DIFF bytes
      for( size_t bond_number( 0); bond_number < number_bonds; ++bond_number)
      {
        const BondConstitutional &bond_to_copy( ATOM.m_Bonds( bond_number));
        m_Bonds.PushBack
        (
          BondConstitutional
          (
            *( const AtomConstitutionalInterface *)( ( const char *)&bond_to_copy.GetTargetAtom() + DIFF),
            bond_to_copy.GetBondType()
          )
        );
      }
    }

    //! @brief set bonds constitution of atom
    //! @param BONDS the bonds that the atom connects
    void AtomConstitutionalShared::SetBonds( const storage::Vector< BondConstitutional> &BONDS)
    {
      m_Bonds = BONDS;
    }

  } // namespace chemistry
} // namespace bcl

