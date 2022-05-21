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
#include "chemistry/bcl_chemistry_bond_conformational.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BondConformational::s_Instance
    (
      GetObjectInstances().AddInstance( new BondConformational())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BondConformational::BondConformational() :
      m_AtomConformational(),
      m_BondType( GetConfigurationalBondTypes().e_Undefined)
    {
    }

    //! @brief constructor from AtomConformational and ConfigurationalBondType charge
    //! @param ATOM_CONFORMATIONAL atom involved in the bond
    //! @param BASE_BOND bondtype of the bond
    BondConformational::BondConformational
    (
      const AtomConformationalInterface &ATOM_CONFORMATIONAL,
      const ConfigurationalBondType &BASE_BOND
    ) :
      m_AtomConformational( &ATOM_CONFORMATIONAL),
      m_BondType( BASE_BOND)
    {
    }

    //! virtual copy constructor
    BondConformational *BondConformational::Clone() const
    {
      return new BondConformational( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BondConformational::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a reference to bond configurational attribute
    //! @return a reference to bond configurational attribute
    const ConfigurationalBondType &BondConformational::GetBondType() const
    {
      return m_BondType;
    }

    //! @brief get a reference to atom conformational to which bond belongs
    //! @return a reference to atom conformational to which bond belongs
    const AtomConformationalInterface &BondConformational::GetTargetAtom() const
    {
      return *m_AtomConformational;
    }

    //! @brief test whether the pointer is to a defined atom
    //! @return true if this bond is defined
    bool BondConformational::IsDefined() const
    {
      return m_AtomConformational.IsDefined() && m_BondType.IsDefined();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BondConformational::Read( std::istream &ISTREAM)
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
    std::ostream &BondConformational::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set bond type of the bond
    //! @param BOND_TYPE the desired bond type
    void BondConformational::SetBondType( const ConfigurationalBondType &BOND_TYPE)
    {
      m_BondType = BOND_TYPE;
    }

  } // namespace chemistry
} // namespace bcl

