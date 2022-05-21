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
#include "chemistry/bcl_chemistry_atom_configurational_shared.h"

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
    AtomConfigurationalShared::AtomConfigurationalShared()
    {
    }

    //! @brief constructor from initializer, needed by atom vector
    //! @param ATOM_INFO all info about the atom; only chirality will be used
    AtomConfigurationalShared::AtomConfigurationalShared( const sdf::AtomInfo &ATOM_INFO) :
      m_Chirality( ATOM_INFO.GetChirality()),
      m_Constitution(),
      m_Bonds()
    {
    }

    //! virtual copy constructor
    AtomConfigurationalShared *AtomConfigurationalShared::Clone() const
    {
      return new AtomConfigurationalShared( *this);
    }

    //! @brief copy constructor
    //! @param ATOM atom configuration object of atom of interest
    AtomConfigurationalShared::AtomConfigurationalShared( const AtomConfigurationalShared &ATOM)
    {
      // call operator = to avoid code redundancy
      *this = ATOM;
    }

    //! @brief assignment constructor
    //! @param ATOM atom configuration object of atom of interest
    AtomConfigurationalShared &AtomConfigurationalShared::operator =( const AtomConfigurationalShared &ATOM)
    {
      if( this != &ATOM)
      {
        // AtomConfigurationalShared must always be in an AtomVector
        // likewise, we know the atom is in a vector, so the addresses of the other atoms in the bonds must be offset
        // by the difference between the address of this atom and the address of ATOM
        //Copy( ATOM, ( const AtomConfigurationalShared * const)( this) - ( const AtomConfigurationalShared * const)( &ATOM));
        Copy( ATOM, ( char *)( this) - ( const char *)( &ATOM));
      }

      return *this;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomConfigurationalShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns chirality
    //! @return a reference to chirality attribute
    const Chirality &AtomConfigurationalShared::GetChirality() const
    {
      return m_Chirality;
    }

    //! @brief returns atom type
    //! @return a reference of stom type data attribute
    const AtomType &AtomConfigurationalShared::GetAtomType() const
    {
      return m_Constitution->GetAtomType();
    }

    //! @brief returns bonds
    //! @return the bond configuration of bond connected the atom
    const storage::Vector< BondConfigurational> &AtomConfigurationalShared::GetBonds() const
    {
      return m_Bonds;
    }

    //! @brief return a simple pointer to atom constitutional object
    //! @return pointer to atom constitutional object
    const util::SiPtr< const AtomConstitutionalShared> AtomConfigurationalShared::GetConstitution() const
    {
      return m_Constitution;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomConfigurationalShared::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Chirality, ISTREAM);
      io::Serialize::Read( m_Constitution, ISTREAM);
      io::Serialize::Read( m_Bonds, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomConfigurationalShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Chirality, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Constitution, OSTREAM, INDENT) << '\n';
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

    //! @brief set bonds
    //! @param BONDS the bonds that the atom connects
    void AtomConfigurationalShared::SetBonds( const storage::Vector< BondConfigurational> &BONDS)
    {
      m_Bonds = BONDS;
    }

    //! @brief copy the atom and bonds using the specified difference in pointer address
    //! @param ATOM the atom to copy
    //! @param DIFF difference in begin of vector that atoms are located in
    void AtomConfigurationalShared::Copy( const AtomConfigurationalShared &ATOM, const ptrdiff_t &DIFF)
    {
      m_Constitution = ATOM.m_Constitution;
      m_Chirality = ATOM.m_Chirality;

      const size_t number_bonds( ATOM.m_Bonds.GetSize());
      m_Bonds.Reset();
      m_Bonds.AllocateMemory( number_bonds);

      // update the bonds accordingly
      for( size_t bond_number( 0); bond_number < number_bonds; ++bond_number)
      {
        const BondConfigurational &bond_to_copy( ATOM.m_Bonds( bond_number));
        const AtomConfigurationalInterface &target_atom( bond_to_copy.GetTargetAtom());

        m_Bonds.PushBack
        (
          BondConfigurational
          (
            *( const AtomConfigurationalInterface *)( ( const char *)( &target_atom) + DIFF),
            bond_to_copy.GetBondType()
          )
        );
      }
    }

  } // namespace chemistry
} // namespace bcl

