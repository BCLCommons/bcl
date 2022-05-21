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
#include "sdf/bcl_sdf_bond_info.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitutional_bond_type_data.h"
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BondInfo::BondInfo() :
      m_AtomIndexLow( util::GetUndefined< size_t>()),
      m_AtomIndexHigh( util::GetUndefined< size_t>()),
      m_BondType()
    {
    }

    //! @brief constructor from members
    BondInfo::BondInfo
    (
      const size_t &ATOM_INDEX_A,
      const size_t &ATOM_INDEX_B,
      const chemistry::ConfigurationalBondType &BOND_TYPE
    ) :
      m_AtomIndexLow( std::min( ATOM_INDEX_A, ATOM_INDEX_B)),
      m_AtomIndexHigh( std::max( ATOM_INDEX_A, ATOM_INDEX_B)),
      m_BondType( BOND_TYPE)
    {
    }

    //! @brief constructor from members
    BondInfo::BondInfo
    (
      const size_t &ATOM_INDEX_A,
      const size_t &ATOM_INDEX_B,
      const chemistry::ConstitutionalBondType &BOND_TYPE
    ) :
      m_AtomIndexLow( std::min( ATOM_INDEX_A, ATOM_INDEX_B)),
      m_AtomIndexHigh( std::max( ATOM_INDEX_A, ATOM_INDEX_B)),
      m_BondType( BOND_TYPE->WithIsometry( chemistry::e_UnknownIsometry))
    {
    }

    //! @brief virtual copy constructor
    BondInfo *BondInfo::Clone() const
    {
      return new BondInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set the chirality from a given property string
    //! @param ISOMETRY the desired isometry
    void BondInfo::SetIsometry( const chemistry::BondIsometry &ISOMETRY)
    {
      m_BondType = m_BondType->WithIsometry( ISOMETRY);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine if the specified line fits the formatting for a bond line
    //! @param LINE the line to check
    //! @return true if the line fits the formatting
    bool BondInfo::FormattedAsMdlBondLine( const std::string &LINE)
    {
      return GetMdlEntryTypes().Bond_FirstAtomIndex->IsUnsignedInt( LINE) &&
        GetMdlEntryTypes().Bond_SecondAtomIndex->IsUnsignedInt( LINE) &&
        GetMdlEntryTypes().Bond_Type->IsUnsignedInt( LINE);
    }

    //! @brief extract information from a line of the sdf believed to be an mdl bond line
    //! @param LINE line from an sdf file that is believed to be a bond line
    //! @return reference to this
    //! This only extracts the bond order or whether the bond is aromatic
    //! the bond type can be refined by later calls to set functions
    BondInfo &BondInfo::ExtractMdlBondLineInfo( const std::string &LINE)
    {
      // load in the atom indices from the line
      m_AtomIndexLow = GetMdlEntryTypes().Bond_FirstAtomIndex->GetUnsignedInt( LINE) - 1;
      m_AtomIndexHigh = GetMdlEntryTypes().Bond_SecondAtomIndex->GetUnsignedInt( LINE) - 1;

      // reorder atom indices
      if( m_AtomIndexHigh < m_AtomIndexLow)
      {
        std::swap( m_AtomIndexHigh, m_AtomIndexLow);
      }

      // get the basic bond type
      m_BondType =
        chemistry::GetConstitutionalBondTypes().FindBondTypeFromSDFInfo
        (
          GetMdlEntryTypes().Bond_Type->GetUnsignedInt( LINE)
        )->WithIsometry( chemistry::e_UnknownIsometry);

      return *this;
    }

    //! @brief create an mdl bond line string out of the bond info
    //! @return an mdl bond line string created from the bond info
    std::string BondInfo::ToMdlBondLine() const
    {
      std::string mdl_bond_line( GetDefaultLine( e_BondLine));
      GetMdlEntryTypes().Bond_FirstAtomIndex->Set( mdl_bond_line, m_AtomIndexLow + 1);
      GetMdlEntryTypes().Bond_SecondAtomIndex->Set( mdl_bond_line, m_AtomIndexHigh + 1);
      m_BondType->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic
          && GetExplicitAromaticBondsPref() ?
              GetMdlEntryTypes().Bond_Type->Set( mdl_bond_line, size_t( 4)) : // aromatic
              GetMdlEntryTypes().Bond_Type->Set( mdl_bond_line, m_BondType->GetSDFileID()); // kekulized
      return mdl_bond_line;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare this bond info to other bond info
    //! @param BOND_INFO the object to compare this bond info to
    //! @return true if the bond info is the same
    bool BondInfo::operator ==( const BondInfo &BOND_INFO) const
    {
      // data type to use to compare bonds.  This type is chosen so that aromatic bonds compare equal regardless of
      // arbitrary bond order
      static const chemistry::ConfigurationalBondTypeData::Data complete
      (
        chemistry::ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness
      );
      return
          m_AtomIndexLow == BOND_INFO.m_AtomIndexLow
          && m_AtomIndexHigh == BOND_INFO.m_AtomIndexHigh
          && m_BondType->GetBondData( complete) == BOND_INFO.m_BondType->GetBondData( complete);
    }

    //! @brief compare this bond info to other bond info
    //! @param BOND_INFO the object to compare this bond info to
    //! @return true if the bond info is less than BOND_INFO
    bool BondInfo::operator <( const BondInfo &BOND_INFO) const
    {
      if( m_AtomIndexLow != BOND_INFO.m_AtomIndexLow)
      {
        return m_AtomIndexLow < BOND_INFO.m_AtomIndexLow;
      }
      if( m_AtomIndexHigh != BOND_INFO.m_AtomIndexHigh)
      {
        return m_AtomIndexHigh < BOND_INFO.m_AtomIndexHigh;
      }
      return m_BondType->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness)
             < BOND_INFO.m_BondType->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read BondInfo object from std::istream
    //! @param ISTREAM istream that contains BondInfo object
    //! @return istream after BondInfo object was extracted
    std::istream &BondInfo::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AtomIndexLow,  ISTREAM);
      io::Serialize::Read( m_AtomIndexHigh, ISTREAM);
      std::string bond_type_name;
      io::Serialize::Read( bond_type_name, ISTREAM);
      m_BondType = chemistry::ConfigurationalBondType( bond_type_name);

      // return
      return ISTREAM;
    }

    //! @brief write BondInfo into std::ostream
    //! @param OSTREAM ostream that gets BondInfo object
    //! @param INDENT indentation
    //! @return ostream after BondInfo object was inserted
    std::ostream &BondInfo::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AtomIndexLow,  OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_AtomIndexHigh, OSTREAM, 0)      << '\t';
      io::Serialize::Write( m_BondType.GetName(), OSTREAM, 0);

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BondInfo::s_Instance
    (
      GetObjectInstances().AddInstance( new BondInfo())
    );

  } // namespace sdf
} // namespace bcl

