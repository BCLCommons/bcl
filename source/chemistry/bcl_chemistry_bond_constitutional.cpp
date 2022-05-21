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
#include "chemistry/bcl_chemistry_bond_constitutional.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BondConstitutional::s_Instance
    (
      GetObjectInstances().AddInstance( new BondConstitutional())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BondConstitutional::BondConstitutional()
    {
    }

    //! @brief constructor from AtomConstitutional and ConfigurationalBondType charge
    //! @param ATOM_CONSTITUTIONAL atom involved in the bond
    //! @param BASE_BOND bondtype of the bond
    BondConstitutional::BondConstitutional
    (
      const AtomConstitutionalInterface &ATOM_CONSTITUTIONAL,
      const ConstitutionalBondType &BASE_BOND
    ) :
      m_AtomConstitutional( &ATOM_CONSTITUTIONAL),
      m_BondType( BASE_BOND)
    {
    }

    //! @brief constructor from AtomConstitutional and configurational BondType
    //! @param ATOM_CONSTITUTIONAL atom involved in the bond
    //! @param BASE_BOND bondtype of the bond
    BondConstitutional::BondConstitutional
    (
      const AtomConstitutionalInterface &ATOM_CONSTITUTIONAL,
      const ConfigurationalBondType &BASE_BOND
    ) :
      m_AtomConstitutional( &ATOM_CONSTITUTIONAL),
      m_BondType( BASE_BOND->GetConstitutionalBondType())
    {
    }

    //! virtual copy constructor
    BondConstitutional *BondConstitutional::Clone() const
    {
      return new BondConstitutional( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BondConstitutional::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a reference to bond constitutional attribute
    //! @return a reference to bond constitutional attribute
    const ConstitutionalBondType &BondConstitutional::GetBondType() const
    {
      return m_BondType;
    }

    //! @brief get a reference to atom constitution to which bond belongs
    //! @return a reference to atom constitution to which bond belongs
    const AtomConstitutionalInterface &BondConstitutional::GetTargetAtom() const
    {
      return *m_AtomConstitutional;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BondConstitutional::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_BondType, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &BondConstitutional::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_BondType, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

  } // namespace chemistry
} // namespace bcl

