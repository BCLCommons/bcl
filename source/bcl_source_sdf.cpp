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
#include "sdf/bcl_sdf_atom_info.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief default constructor
    AtomInfo::AtomInfo() :
      m_AtomType(),
      m_Chirality( chemistry::e_UnknownChirality),
      m_Coordinates(),
      m_CanAddH( true)
    {
    }

    //! @brief constructor from all known atom data
    AtomInfo::AtomInfo
    (
      const chemistry::AtomType  &ATOM_TYPE,
      const chemistry::Chirality &CHIRALITY,
      const linal::Vector3D      &COORDINATES,
      const bool                 &CAN_ADD_H
    ) :
      m_AtomType( ATOM_TYPE),
      m_Chirality( CHIRALITY),
      m_Coordinates( COORDINATES),
      m_CanAddH( CAN_ADD_H)
    {
    }

    //! @brief virtual copy constructor
    AtomInfo *AtomInfo::Clone() const
    {
      return new AtomInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AtomInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief change the atom type
    //! @param ATOM_TYPE the new atom type
    //! This function is necessary if the mdl file contained a property called AtomType, with one entry for each non-H
    void AtomInfo::SetAtomType( const chemistry::AtomType &ATOM_TYPE)
    {
      m_AtomType = ATOM_TYPE;
    }

    //! @brief change the chirality of the atom type
    //! @param CHIRALITY the new chirality to give the atom
    //! This function is necessary if the mdl file contained a property called Chirality, one entry for each atom with
    //! 4 bonds
    void AtomInfo::SetChirality( const chemistry::Chirality &CHIRALITY)
    {
      m_Chirality = CHIRALITY;
    }

    //! @brief change the charge on the atom type
    //! @param NEW_CHARGE the new charge to give the atom type
    //! This function is necessary if the mdl file contained M  CHG lines that reference this atom
    void AtomInfo::SetCharge( const short &CHARGE)
    {
      if( m_AtomType->GetFormalCharge() != CHARGE)
      {
        m_AtomType = chemistry::AtomTypes::GetAtomType( m_AtomType->GetElementType(), CHARGE);
      }
    }

    //! @brief change the coordinates
    //! @param COORDINATES the new coordinates
    void AtomInfo::SetCoordinates( const linal::Vector3D &COORDINATES)
    {
      m_Coordinates = COORDINATES;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare this atom info to other atom info
    //! @param ATOM_INFO the object to compare this atom info to
    //! @return true if the atom info is the same
    bool AtomInfo::operator ==( const AtomInfo &ATOM_INFO) const
    {
      // compare all data members
      return
           m_AtomType    == ATOM_INFO.m_AtomType
        && m_Chirality   == ATOM_INFO.m_Chirality
        && m_Coordinates == ATOM_INFO.m_Coordinates
        && m_CanAddH     == ATOM_INFO.m_CanAddH;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief create an mdl atom line string out of the atom info
    //! @return an mdl atom line string created from the atom info
    std::string AtomInfo::ToMdlAtomLine() const
    {
      std::string line( GetDefaultLine( e_AtomLine));
      GetMdlEntryTypes().Atom_CoordinateX->Set( line, m_Coordinates.X());
      GetMdlEntryTypes().Atom_CoordinateY->Set( line, m_Coordinates.Y());
      GetMdlEntryTypes().Atom_CoordinateZ->Set( line, m_Coordinates.Z());
      GetMdlEntryTypes().Atom_Symbol->Set( line, m_AtomType->GetElementType()->GetChemicalSymbol());
      if( m_AtomType->GetFormalCharge() && std::abs( m_AtomType->GetFormalCharge()) <= 3)
      {
        GetMdlEntryTypes().Atom_Charge->Set( line, 4 - m_AtomType->GetFormalCharge());
      }
      if( !m_CanAddH) // if no additional hydrogens can be added, set the Atom_HydrogenCount to 1 to indicate it, as per the mdl spec
      {
        GetMdlEntryTypes().Atom_HydrogenCount->Set( line, 1);
      }

      // Note: Valence/Number Bonds column is interpreted by programs differently, so it is best to just not write it
      // out at all @see @link http://molmatinf.com/whynotmolsdf.html @endlink
      //if( m_AtomType->IsGasteigerAtomType())
      //{
      //  GetMdlEntryTypes().Atom_Valence->Set( line, m_AtomType->GetNumberBonds());
      //}
      return line;
    }

    //! @brief extract information from the mdl atom line string
    //! @param LINE line from an sdf file that is believed to be an atom line
    //! @return reference to this
    //! This does not extract chirality or the gasteiger atom type (only element type and possibly charge)
    bool AtomInfo::FormattedAsMdlAtomLine( const std::string &LINE)
    {
      return GetMdlEntryTypes().Atom_CoordinateX->IsDouble( LINE) &&
        GetMdlEntryTypes().Atom_CoordinateY->IsDouble( LINE) &&
        GetMdlEntryTypes().Atom_CoordinateZ->IsDouble( LINE) &&
        !GetMdlEntryTypes().Atom_Symbol->GetTrimmedString( LINE).empty() &&
        !GetMdlEntryTypes().Atom_Charge->GetTrimmedString( LINE).empty();
    }

    //! @brief extract information from the mdl atom line string
    //! @param LINE line from an sdf file that is believed to be an atom line
    //! @return reference to this
    //! This does not extract chirality or the gasteiger atom type (only element type and possibly charge)
    AtomInfo &AtomInfo::ExtractMdlAtomLineInfo( const std::string &LINE)
    {
      // get the coordinates
      m_Coordinates.X() = GetMdlEntryTypes().Atom_CoordinateX->GetDouble( LINE);
      m_Coordinates.Y() = GetMdlEntryTypes().Atom_CoordinateY->GetDouble( LINE);
      m_Coordinates.Z() = GetMdlEntryTypes().Atom_CoordinateZ->GetDouble( LINE);

      // tell whether H can be added
      m_CanAddH = GetMdlEntryTypes().Atom_HydrogenCount->GetTrimmedString( LINE) != "1";

      // get the element type
      const chemistry::ElementType element_type
      (
        chemistry::GetElementTypes().ElementTypeLookup( GetMdlEntryTypes().Atom_Symbol->GetTrimmedString( LINE))
      );

      // get the mdl charge, which is 4 - the real charge, except for neutral atoms
      const std::string mdl_charge_str( GetMdlEntryTypes().Atom_Charge->GetTrimmedString( LINE));

      // valid MDL charges are always 1 character between 0 - 7
      if( mdl_charge_str.size() != size_t( 1) || !isdigit( mdl_charge_str[ 0]) || mdl_charge_str[ 0] > '7')
      {
        BCL_MessageCrt
        (
          "Atom with non-mdl charge id (0-7 is the accepted range): " + mdl_charge_str
        );
      }
      else
      {
        // one character charge; normal case
        const short mdl_charge( mdl_charge_str[ 0] - '0');

        // from the mdl documentation on charges:
        // 0 = uncharged, 1 = +3, 2 = +2, 3 = +1,
        // 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3, which in general means
        // charge = 4 - charge_id

        // map the charge id to an actual value
        short charge( 0);
        if( mdl_charge != 0)
        {
          charge = 4 - mdl_charge;
        }
        m_AtomType = chemistry::AtomTypes::GetAtomType( element_type, charge);
      }

      // chirality info in the sdf block is ignored, as per the mdl spec
      return *this;
    }

    //! @brief read AtomInfo object from std::istream
    //! @param ISTREAM istream that contains AtomInfo object
    //! @return istream after AtomInfo object was extracted
    std::istream &AtomInfo::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AtomType,    ISTREAM);
      io::Serialize::Read( m_Chirality,   ISTREAM);
      io::Serialize::Read( m_CanAddH,     ISTREAM);
      //io::Serialize::Read( m_AtomMapping, ISTREAM);
      io::Serialize::Read( m_Coordinates, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write AtomInfo into std::ostream
    //! @param OSTREAM ostream that gets AtomInfo object
    //! @param INDENT indentation
    //! @return ostream after AtomInfo object was inserted
    std::ostream &AtomInfo::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AtomType,    OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Chirality,   OSTREAM,      0) << '\t';
      io::Serialize::Write( m_CanAddH,     OSTREAM,      0) << '\t';
      //io::Serialize::Write( m_AtomMapping, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Coordinates, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AtomInfo::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomInfo())
    );

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_ctab_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief standard constructor
    CTabHandler::CTabHandler() :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
    }

    //! @brief constructor from an input stream
    CTabHandler::CTabHandler( std::istream &ISTREAM) :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      ReadCTab( ISTREAM);
    }

    //! @brief constructor from a pre-read set of lines
    CTabHandler::CTabHandler( const storage::List< std::string> &LINES) :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      ReadCTab( LINES.Begin(), LINES.End());
    }

    //! @brief constructor from a list of AtomInfo and BondInfos
    CTabHandler::CTabHandler
    (
      const storage::Vector< AtomInfo> &ATOM_INFOS,
      const storage::Vector< BondInfo> &BOND_INFOS
    ) :
      m_IsValid( true),
      m_Header( ATOM_INFOS.GetSize(), BOND_INFOS.GetSize()),
      m_AtomInfos( ATOM_INFOS),
      m_BondInfos( BOND_INFOS),
      m_AtomMapping( ATOM_INFOS.GetSize(), size_t( 0)),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      size_t n_atoms( ATOM_INFOS.GetSize());

      // check sanity of the bond data
      for
      (
        storage::Vector< BondInfo>::const_iterator itr_bond( m_BondInfos.Begin()),
          itr_bond_end( m_BondInfos.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        if( itr_bond->GetAtomIndexLow() >= n_atoms || itr_bond->GetAtomIndexHigh() >= n_atoms)
        {
          BCL_MessageCrt
          (
            "Could not make a connection table, invalid data given.  Bond tried to join atoms "
            + util::Format()( itr_bond->GetAtomIndexLow()) + " and " + util::Format()( itr_bond->GetAtomIndexHigh()) +
            " but maximum index is " + util::Format()( n_atoms - 1)
          );
          m_AtomInfos.Reset();
          m_BondInfos.Reset();
          m_AtomMapping.Reset(),
          m_Header = MdlHeader();
          m_IsValid = false;
          break;
        }
      }
    }

    //! @brief virtual copy constructor
    CTabHandler *CTabHandler::Clone() const
    {
      return new CTabHandler( *this);
    }

    //! @brief Read the CTab using iterators to strings
    //! @param LINE_BEGIN a line that represents a header/counts line
    //! @param LINE_END one-past-end of possible lines
    //! @details if the CTab ends before LINE_END, not all lines will be read
    //! @return an iterator to the first line that was not read
    storage::List< std::string>::const_iterator CTabHandler::ReadCTab
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END
    )
    {

      Reset();

      storage::List< std::string>::const_iterator itr( LINE_BEGIN);
      if( LINE_BEGIN == LINE_END)
      {
        BCL_MessageStd( "CTabHandler: nothing to read");
        return LINE_BEGIN;
      }

      // Read the header, should be the first line
      m_Header.SetFromMdlLine( *itr, 0);

      if( !m_Header.IsValid())
      {
        BCL_MessageStd( "CTab header \"" + *itr + "\" is invalid, not reading any further");
        return ++itr;
      }
      ++itr;

      const size_t n_atoms( m_Header.GetNumberAtoms());
      const size_t n_bonds( m_Header.GetNumberBonds());

      m_AtomInfos.AllocateMemory( n_atoms);
      m_BondInfos.AllocateMemory( n_bonds);
      m_AtomMapping.AllocateMemory( n_atoms);

      // read atom information
      size_t found_atoms = 0;
      for( ; itr != LINE_END && found_atoms < n_atoms; ++itr)
      {
        if( !ContainsNonspaceCharacters( *itr))
        {
          BCL_MessageCrt( "Should not find an empty line in the atom line section! - skipping line");
          continue;
        }

        if( !AtomInfo::FormattedAsMdlAtomLine( *itr))
        {
          BCL_MessageCrt( "Line \"" + *itr + "\" was not formatted as an MDL atom line; not reading any further in this molecule");
          break;
        }

        // get the basic atom information and mapping
        m_AtomInfos.PushBack( AtomInfo().ExtractMdlAtomLineInfo( *itr));
        if( GetMdlEntryTypes().Atom_AtomMappingNumber->IsUnsignedInt( *itr))
        {
          m_AtomMapping.PushBack( GetMdlEntryTypes().Atom_AtomMappingNumber->GetUnsignedInt( *itr));
        }
        else
        {
          m_AtomMapping.PushBack( 0);
        }
        ++found_atoms;
      }

      // check we read enough atoms
      if( found_atoms != n_atoms)
      {
        BCL_MessageCrt
        (
          "Did not find the appropriate number of atoms!  The MDL header line specified " +
          util::Format()( n_atoms) + " atom lines but only found " + util::Format()( found_atoms)
        );
        return itr; // this will be one after the last atom line that was read
      }

      // read bonds
      size_t found_bonds = 0;
      for( ; itr != LINE_END && found_bonds < n_bonds; ++itr)
      {
        if( !ContainsNonspaceCharacters( *itr))
        {
          BCL_MessageCrt( "Should not find an empty line in the bond line section! - skipping line");
          continue;
        }

        // check that formatting is correct
        if( !BondInfo::FormattedAsMdlBondLine( *itr))
        {
          BCL_MessageCrt( "Line \"" + *itr + "\" was not formatted as an MDL bond line; not reading any further in this molecule");
          break;
        }

        // create the bond line
        BondInfo bond;
        bond.ExtractMdlBondLineInfo( *itr);

        // check the high atom.  If < n_atoms we're fine, if > n_atoms then at least this one is incorrect
        if( bond.GetAtomIndexHigh() >= n_atoms)
        {
          BCL_MessageCrt
          (
            "Warning: line \"" + *itr + "\" references an invalid atom (" +
            util::Format()( bond.GetAtomIndexHigh() + 1) + ")."
            "  Maximum atom number for this molecule is " + util::Format()( n_atoms) + ". Skipping this bond"
          );
          continue;
        }

        // check that the bond type is recognized, otherwise skip it
        if( !bond.GetConfigurationalBondType().IsDefined())
        {
          BCL_MessageCrt
          (
            "Warning: bond line \"" + *itr + "\" contains an unrecognized bond order; skipping this bond"
          );

          // increment this because we still read everything properly, we just don't know what to do with it
          ++found_bonds;
          continue;
        }
        m_BondInfos.PushBack( bond);
        ++found_bonds;
      }

      // check that enough bonds were read
      if( found_bonds != n_bonds)
      {
        BCL_MessageCrt
        (
          "Did not find the appropriate number of bonds!  The MDL header specified " +
          util::Format()( n_bonds) + " bond lines but only found " + util::Format()( found_bonds)
        );
        return itr;
      }

      // read in MDL properties
      for( ; itr != LINE_END; ++itr)
      {
        MdlProperty prop( *itr);

        //if( !MdlProperty::IsBCLPropertyLine( *itr))
        if( !prop.IsValid())
        {
          BCL_MessageStd( "Unrecognized MDL property line \"" + *itr + "\"; skipping line.");
          continue;
        }

        // found the end of the block, so everything was read correctly
        if( prop.GetProperty() == MdlProperty::e_BlockTerminator)
        {
          ++itr; // go to the next line, i.e. last unread line
          m_IsValid = true; // we are only valid if we actually read a terminal line
          break;
        }

        // some, but not all, properties need to be stored with the associated ctab.  check for this
        if( prop.ShouldCache())
        {
          m_MDLProperties[ prop.GetLabel()] = prop.GetPropertyStrings();
        }

        ApplyMDLProperty( prop);
      }
      return itr;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CTabHandler::ReadCTab( std::istream &ISTREAM)
    {
      // buffer for the stream
      storage::List< std::string> lines;

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        BCL_MessageCrt( "Error reading CTab: passed bad istream");
        return ISTREAM;
      }

      while( !ISTREAM.eof())
      {
        lines.PushBack( std::string());
        std::string &last_line( lines.LastElement());
        std::getline( ISTREAM, lines.LastElement());
        if( IsTerminalLine( last_line))
        {
          ReadCTab( lines.Begin(), lines.End());
          break;
        }
      }
      return ISTREAM;
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param FORCE_WRITE_ATOM_TYPES whether to force writing BCL atom types (and other info)
    std::ostream &CTabHandler::WriteCTab
    (
      std::ostream &OSTREAM,
      const bool &FORCE_WRITE_ATOM_TYPES
    ) const
    {

      // don't write anything if our CTab isn't valid
      if( !IsValid())
      {
        return OSTREAM;
      }

      // write to a stringstream first and output at the end of writing
      // if anything bad happens (can't convert int etc) then we don't want to write a fragmented file
      std::stringstream tmp_ostream;

      // v2000 format allots 3 characters for # of atoms / bonds, so it cannot handle more than 999 atoms or bonds
      if( m_AtomInfos.GetSize() > size_t( 999) || m_BondInfos.GetSize() > size_t( 999))
      {
        BCL_MessageCrt
        (
          "Ignoring request to write a connection table (CTab);  "
          "CTab contains " + util::Format()( m_AtomInfos.GetSize()) + " atoms and "
          + util::Format()( m_BondInfos.GetSize()) + " bonds, but the V2000 format can only"
          "support up to 999 of each"
        );
        return OSTREAM;
      }
      //std::ostringstream tmp_ostream;

      // write out the header
      tmp_ostream << m_Header.ToMdlLine() << '\n';

      // write out all atom info
      size_t atom_no( 0);
      for
      (
        storage::Vector< AtomInfo>::const_iterator
          itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++atom_no
      )
      {
        // basic atom information
        std::string line( itr_atom->ToMdlAtomLine());

        // add mapping
        GetMdlEntryTypes().Atom_AtomMappingNumber->Set( line, m_AtomMapping( atom_no));

        tmp_ostream << line << '\n';
      }

      // write out all bond lines
      for
      (
        storage::Vector< BondInfo>::const_iterator
          itr_bond( m_BondInfos.Begin()), itr_bond_end( m_BondInfos.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        tmp_ostream << itr_bond->ToMdlBondLine() << '\n';
      }

      // determine which MDL properties need to be added
      storage::Vector< MdlProperty::PropertyEnum> properties_to_add( GetNecessaryMDLProperties( FORCE_WRITE_ATOM_TYPES));

      for
      (
        storage::Vector< MdlProperty::PropertyEnum>::const_iterator
          itr_prop( properties_to_add.Begin()), itr_prop_end( properties_to_add.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // create the string from the property
        const std::string property_str( MdlProperty( *itr_prop, m_AtomInfos, m_BondInfos).GetString());

        if( !property_str.empty())
        {
          // write the given property
          tmp_ostream << property_str << '\n';
        }
        else
        {
          BCL_MessageVrb( "Warning: Property \"" + MdlProperty::GetPropertyName( *itr_prop) + "\" was needed but no info was written");
        }
      }

      // write properties block terminator
      tmp_ostream << MdlProperty( MdlProperty::e_BlockTerminator).GetString() << "\n";

      OSTREAM << tmp_ostream.str();
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains 'M  END'
    bool CTabHandler::IsTerminalLine( const std::string &LINE)
    {
      return LINE == MdlProperty( MdlProperty::e_BlockTerminator).GetString();
    }

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool CTabHandler::ContainsNonspaceCharacters( const std::string &STRING)
    {
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( !isspace( *itr))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief apply instructions from an MDL property to this class
    //! @param PROP the property to apply
    void CTabHandler::ApplyMDLProperty( const MdlProperty &PROP)
    {
      if( !PROP.IsValid())
      {
        BCL_MessageStd( "Warning: Property \"" + PROP.GetString() + "\" is not a valid property, not reading it");
      }
      if( PROP.GetProperty() < MdlProperty::s_NumberProperties && PROP.IsValid())
      {
        PROP.ApplyProperty( m_AtomInfos, m_BondInfos);
        if( PROP.GetProperty() == MdlProperty::e_BclAtomType)
        {
          m_AtomTypesRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclBondType)
        {
          m_BondTypesRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclChirality)
        {
          m_ChiralityWasRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclDoubleBondIsometry)
        {
          m_DoubleBondIsometryWasRead = true;
        }
      }
    }

    //! @brief determine the MDL properties that must be used when writing CTab information
    //! @return a vector of MDL PropertyEnums that denote which properties must be computed
    storage::Vector< MdlProperty::PropertyEnum> CTabHandler::GetNecessaryMDLProperties( const bool &FORCE_ATOM_TYPES) const
    {
      storage::Vector< MdlProperty::PropertyEnum> properties_to_add;

      // determine whether more than 1 coordinate is undefined / 0, in which case all
      // isometry information must be written as well
      // If the positions are all defined, and no more than 1 atom is at the origin, then
      // the isometry can be determined later on
      bool positions_invalid( false);
      {
        size_t zero_position_count( 0);
        for
        (
          storage::Vector< AtomInfo>::const_iterator
            itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // keep track of whether the coordinates are all defined and that no more than one of them is at the origin
          if( !itr_atom->GetCoordinates().IsDefined())
          {
            positions_invalid = true;
            break;
          }
          else if( math::EqualWithinAbsoluteTolerance( itr_atom->GetCoordinates().Norm(), double( 0.0), 1.0e-4))
          {
            if( ++zero_position_count == 2)
            {
              positions_invalid = true;
              break;
            }
          }
        }
      }

      // check whether charges must be written out as a separate property
      bool charges_outside_of_mdl_range( false);
      {
        for
        (
          storage::Vector< AtomInfo>::const_iterator
            itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // if the absolute value of a charge is greater than 3 it must be put in the "M  CHG' field
          if( std::abs( itr_atom->GetAtomType()->GetFormalCharge()) > short( 3))
          {
            charges_outside_of_mdl_range = true;
            break;
          }
        }
      }

      // if valences are present, the bcl atom type must be written out, otherwise
      // the atom type may be determined incorrectly when the file is reloaded
      if( FORCE_ATOM_TYPES || TestMustWriteAtomTypes( m_AtomInfos, m_BondInfos))
      {
        properties_to_add.PushBack( MdlProperty::e_BclAtomType);
        properties_to_add.PushBack( MdlProperty::e_BclBondType);
      }

      // if positions were invalid, or atom types are being written, isometry info must be written as well
      if( positions_invalid || properties_to_add.GetSize())
      {
        properties_to_add.PushBack( MdlProperty::e_BclChirality);
        properties_to_add.PushBack( MdlProperty::e_BclDoubleBondIsometry);
      }

      // if any charges were outside the normal mdl range, add them as a property
      if( charges_outside_of_mdl_range)
      {
        properties_to_add.PushBack( MdlProperty::e_Charge);
      }

      return properties_to_add;
    }

    //! @brief from a set of atoms and bonds, determine whether atom types must be written
    //! @param ATOMS vector of atom info
    //! @param BONDS vector of bond info
    //! @return true if atom types must be written to ensure that they will be the same upon rereading the molecule
    bool CTabHandler::TestMustWriteAtomTypes
    (
      const storage::Vector< AtomInfo> &ATOMS,
      const storage::Vector< BondInfo> &BONDS
    )
    {
      storage::Vector< size_t> bond_counts( ATOMS.GetSize(), size_t( 0));
      // add up the bond counts for each atom
      for( storage::Vector< BondInfo>::const_iterator itr( BONDS.Begin()), itr_end( BONDS.End()); itr != itr_end; ++itr)
      {
        ++bond_counts( itr->GetAtomIndexLow());
        ++bond_counts( itr->GetAtomIndexHigh());
      }

      // test whether any atom has a valence
      bool had_valence( false);
      size_t atom_index( 0);
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr, ++atom_index
      )
      {
        // only consider well-defined atom types
        if( itr->GetAtomType().IsDefined() && itr->GetAtomType()->IsGasteigerAtomType())
        {
          if( bond_counts( atom_index) < itr->GetAtomType()->GetNumberBonds())
          {
            had_valence = true;
            break;
          }
        }
      }

      // return true if there were any valences
      return had_valence;
    }

    //! @brief read CTabHandler object from std::istream
    //! @param ISTREAM istream that contains CTabHandler object
    //! @return istream after CTabHandler object was extracted
    std::istream &CTabHandler::Read( std::istream &ISTREAM)
    {
      return ReadCTab( ISTREAM);
    }

    //! @brief write CTabHandler into std::ostream
    //! @param OSTREAM ostream that gets CTabHandler object
    //! @param INDENT indentation
    //! @return ostream after CTabHandler object was inserted
    std::ostream &CTabHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CTabHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new CTabHandler())
    );

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecular_configuration_shared.h"
#include "chemistry/bcl_chemistry_molecular_conformation_shared.h"
#include "chemistry/bcl_chemistry_molecule_complete.h"
#include "sdf/bcl_sdf_fragment_factory.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief read complete molecule from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a molecule complete
    chemistry::MoleculeComplete Factory::MakeMolecule( const MdlHandler &HANDLER)
    {
      return FragmentFactory::MakeFragment( HANDLER, e_Saturate);
    }

    //! @brief read molecule conformation from MDL file
    //! @param MDL_HANDLER handler that has the conformation information
    //! @return a molecule conformation shared
    chemistry::MolecularConformationShared Factory::MakeConformation( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConformationShared( MakeMolecule( HANDLER));
    }

    //! @brief read molecule configuration from MDL file
    //! @param MDL_HANDLER handler that has the configuration information
    //! @return a molecule configuration shared
    chemistry::MolecularConfigurationShared Factory::MakeConfiguration( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConfigurationShared( MakeMolecule( HANDLER));
    }

    //! @brief read molecule constitution from MDL file
    //! @param MDL_HANDLER handler that has the constitution information
    //! @return a molecule constitution shared
    chemistry::MolecularConstitutionShared Factory::MakeConstitution( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConstitutionShared( MakeMolecule( HANDLER));
    }

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_fragment_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const CTabHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      // create atoms
      chemistry::AtomVector< chemistry::AtomComplete> atoms( HANDLER.GetAtomInfo(), HANDLER.GetBondInfo());

      if( !HANDLER.WereBondTypesRead() || !HANDLER.WereAtomTypesRead())
      {
        // determine atom and bond types
        chemistry::AtomsCompleteStandardizer standardizer( atoms, MOLECULE_ID, false);
      }
      else
      {
        chemistry::AtomsCompleteStandardizer::SetConjugationOfBondTypes( atoms);
      }
      chemistry::AtomsCompleteStandardizer::TryNeutralize( atoms, NEUTRALIZATION_PREF);

      if( !HANDLER.WasDoubleBondIsometryRead())
      {
        // add bond isometry information
        chemistry::BondIsometryHandler::AddIsometryInformation( atoms, true);
      }
      if( !HANDLER.WasChiralityRead())
      {
        // add stereocenter information
        chemistry::StereocentersHandler::AddChiralityFromConformation( atoms);
      }

      // for each atom, put the preference as to whether to add H into a vector
      storage::Vector< size_t> add_hydrogen;
      add_hydrogen.AllocateMemory( HANDLER.GetAtomInfo().GetSize());
      // iterate over all atoms
      for
      (
        storage::Vector< AtomInfo>::const_iterator
          itr( HANDLER.GetAtomInfo().Begin()), itr_end( HANDLER.GetAtomInfo().End());
        itr != itr_end;
        ++itr
      )
      {
        add_hydrogen.PushBack( itr->CanAddH());
      }
      chemistry::HydrogensHandler::HandleHydrogenPref( atoms, H_PREF, add_hydrogen);

      // Give this molecule a blank name for the time being
      return chemistry::FragmentComplete( atoms, MOLECULE_ID);
    }

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const MolfileHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      chemistry::FragmentComplete new_fragment( MakeFragment( static_cast< CTabHandler>( HANDLER), H_PREF, NEUTRALIZATION_PREF, MOLECULE_ID));
      new_fragment.SetName( HANDLER.GetDescription());
      return new_fragment;
    }

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const MdlHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      chemistry::FragmentComplete new_fragment( MakeFragment( HANDLER.GetMolfile(), H_PREF, NEUTRALIZATION_PREF, MOLECULE_ID));
      const storage::Map< std::string, std::string> &prop_map( HANDLER.GetMiscProperties());
      for
      (
        storage::Map< std::string, std::string>::const_iterator itr_map( prop_map.Begin()), itr_map_end( prop_map.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        new_fragment.StoreProperty( itr_map->first, itr_map->second);
      }
      new_fragment.CacheNumeric( prop_map);
      return new_fragment;
    }

    //! @brief read fragment conformation from MDL file
    //! @param MDL_HANDLER handler that has the conformation information
    //! @return a fragment conformation shared
    chemistry::FragmentConformationShared
      FragmentFactory::MakeConformation( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConformationShared( MakeFragment( HANDLER, H_PREF));
    }

    //! @brief read fragment configuration from MDL file
    //! @param MDL_HANDLER handler that has the configuration information
    //! @return a fragment configuration shared
    chemistry::FragmentConfigurationShared
      FragmentFactory::MakeConfiguration( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConfigurationShared( MakeFragment( HANDLER, H_PREF));
    }

    //! @brief read fragment constitution from MDL file
    //! @param MDL_HANDLER handler that has the constitution information
    //! @return a fragment constitution shared
    chemistry::FragmentConstitutionShared
      FragmentFactory::MakeConstitution( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConstitutionShared( MakeFragment( HANDLER, H_PREF));
    }

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_mdl_entry_type_data.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MdlEntryTypeData::MdlEntryTypeData() :
      m_MdlLineType(),
      m_Start( util::GetUndefined< size_t>()),
      m_Length( util::GetUndefined< size_t>()),
      m_Default(),
      m_Format(),
      m_DataType( util::CPPDataTypes::e_Unknown)
    {
    }

    //! @brief construct from data
    //! @param LINE_TYPE type of line the entry is found in
    //! @param START starting position in line
    //! @param LENGTH length of entry type in line
    //! @param DEFAULT default string
    //! @param FORMAT format for the data field
    //! @param DATA_TYPE datatype of this entry
    MdlEntryTypeData::MdlEntryTypeData
    (
      const MdlLineType LINE_TYPE,
      const size_t START,
      const size_t LENGTH,
      const std::string &DEFAULT,
      const util::Format &FORMAT,
      const util::CPPDataTypes::Types DATA_TYPE
    ) :
      m_MdlLineType( LINE_TYPE),
      m_Start( START),
      m_Length( LENGTH),
      m_Default( DEFAULT),
      m_Format( util::Format( FORMAT).W( LENGTH)),
      m_DataType( DATA_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MdlEntryTypeData
    MdlEntryTypeData *MdlEntryTypeData::Clone() const
    {
      return new MdlEntryTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MdlEntryTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the string from a particular line, without any beginning or end spaces
    //! @param LINE the line to retrieve this entry from
    //! @return the entry from that line
    std::string MdlEntryTypeData::GetTrimmedString( const std::string &LINE) const
    {
      if( m_Start >= LINE.size())
      {
        return std::string();
      }
      size_t start_pos( m_Start);
      size_t end_pos( std::min( m_Start + m_Length, LINE.size()));
      BCL_Assert( LINE.size() >= end_pos, "Line was too small to get string from");
      while( start_pos < end_pos && isspace( LINE[ start_pos]))
      {
        ++start_pos;
      }
      if( start_pos < end_pos)
      {
        while( isspace( LINE[ --end_pos]))
        {
        }
        ++end_pos;
      }
      return std::string( &LINE[ start_pos], end_pos - start_pos);
    }

    //! @brief get this entry as a size_t from a particular line
    //! @param LINE the line to retrieve this entry from
    //! @return the size_t given by this entry from a particular line
    bool MdlEntryTypeData::IsUnsignedInt( const std::string &LINE) const
    {
      unsigned int value( 0);
      // some of the integral fields for MDL types are sometimes ommitted
      // so just return 0 if an empty string is retrieved
      const std::string trimmed_string( GetTrimmedString( LINE));
      if( trimmed_string.empty())
      {
        return 0;
      }
      
      return util::TryConvertFromString( value, trimmed_string, util::GetLogger());
    }

    //! @brief get this entry as a size_t from a particular line
    //! @param LINE the line to retrieve this entry from
    //! @return the size_t given by this entry from a particular line
    unsigned int MdlEntryTypeData::GetUnsignedInt( const std::string &LINE) const
    {
      unsigned int value( 0);
      // some of the integral fields for MDL types are sometimes ommitted
      // so just return 0 if an empty string is retrieved
      const std::string trimmed_string( GetTrimmedString( LINE));
      if( trimmed_string.empty())
      {
        return 0;
      }
      BCL_Assert
      (
        util::TryConvertFromString( value, trimmed_string, util::GetLogger()),
        "Could not convert " + trimmed_string + " to a size_t"
      );
      return value;
    }

    //! @brief get this entry as a double from a particular line
    //! @param LINE the line to retrieve this entry from
    //! @return the double given by this entry from a particular line
    bool MdlEntryTypeData::IsDouble( const std::string &LINE) const
    {
      double value;
      return util::TryConvertFromString( value, GetTrimmedString( LINE), util::GetLogger());
    }

    //! @brief get this entry as a double from a particular line
    //! @param LINE the line to retrieve this entry from
    //! @return the double given by this entry from a particular line
    double MdlEntryTypeData::GetDouble( const std::string &LINE) const
    {
      double value;
      BCL_Assert
      (
        util::TryConvertFromString( value, GetTrimmedString( LINE), util::GetLogger()),
        "Could not convert " + GetTrimmedString( LINE) + " to a float"
      );
      return value;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MdlEntryTypeData::Read( std::istream &ISTREAM)
    {
      // read all member
      io::Serialize::Read( m_MdlLineType, ISTREAM);
      io::Serialize::Read( m_Start,       ISTREAM);
      io::Serialize::Read( m_Length,      ISTREAM);
      io::Serialize::Read( m_Default,     ISTREAM);
      m_Format.Read(                      ISTREAM);
      io::Serialize::Read( m_DataType,    ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MdlEntryTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write all member
      io::Serialize::Write( m_MdlLineType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Start,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Length,      OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Default,     OSTREAM, INDENT) << '\n';
      m_Format.Write(                      OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DataType,    OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_mdl_entry_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief default constructor
    MdlEntryTypes::MdlEntryTypes() :
      Header_NumberAtomLines(        AddEnum( "Header_NumberAtomLines",        MdlEntryTypeData( e_HeaderLine,       0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberBondLines(        AddEnum( "Header_NumberBondLines",        MdlEntryTypeData( e_HeaderLine,       3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberAtomLists(        AddEnum( "Header_NumberAtomLists",        MdlEntryTypeData( e_HeaderLine,       6,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete1(              AddEnum( "Header_Obsolete1",              MdlEntryTypeData( e_HeaderLine,       9,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_ChiralFlag(             AddEnum( "Header_ChiralFlag",             MdlEntryTypeData( e_HeaderLine,      12,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberSTextEntries(     AddEnum( "Header_NumberSTextEntries",     MdlEntryTypeData( e_HeaderLine,      15,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete2(              AddEnum( "Header_Obsolete2",              MdlEntryTypeData( e_HeaderLine,      18,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete3(              AddEnum( "Header_Obsolete3",              MdlEntryTypeData( e_HeaderLine,      21,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete4(              AddEnum( "Header_Obsolete4",              MdlEntryTypeData( e_HeaderLine,      24,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete5(              AddEnum( "Header_Obsolete5",              MdlEntryTypeData( e_HeaderLine,      27,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete6(              AddEnum( "Header_Obsolete6",              MdlEntryTypeData( e_HeaderLine,      30,   3,    "999", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_CtabVersion(            AddEnum( "Header_CtabVersion",            MdlEntryTypeData( e_HeaderLine,      33,   6,  "V2000", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXNHeader_NumberReactantLines( AddEnum( "RXNHeader_NumberReactantLines", MdlEntryTypeData( e_RXNHeaderLine,    0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      RXNHeader_NumberProductLines(  AddEnum( "RXNHeader_NumberProductLines",  MdlEntryTypeData( e_RXNHeaderLine,    3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Atom_CoordinateX(              AddEnum( "Atom_CoordinateX",              MdlEntryTypeData( e_AtomLine,         0,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_CoordinateY(              AddEnum( "Atom_CoordinateY",              MdlEntryTypeData( e_AtomLine,        10,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_CoordinateZ(              AddEnum( "Atom_CoordinateZ",              MdlEntryTypeData( e_AtomLine,        20,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_Symbol(                   AddEnum( "Atom_Symbol",                   MdlEntryTypeData( e_AtomLine,        31,   3,      "X", util::Format().L(),         util::CPPDataTypes::e_String))),
      Atom_MassDifference(           AddEnum( "Atom_MassDifference",           MdlEntryTypeData( e_AtomLine,        34,   2,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_Charge(                   AddEnum( "Atom_Charge",                   MdlEntryTypeData( e_AtomLine,        36,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_StereoParity(             AddEnum( "Atom_StereoParity",             MdlEntryTypeData( e_AtomLine,        39,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_HydrogenCount(            AddEnum( "Atom_HydrogenCount",            MdlEntryTypeData( e_AtomLine,        42,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_StereoCareBox(            AddEnum( "Atom_StereoCareBox",            MdlEntryTypeData( e_AtomLine,        45,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_Valence(                  AddEnum( "Atom_Valence",                  MdlEntryTypeData( e_AtomLine,        48,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_HODesignator(             AddEnum( "Atom_HODesignator",             MdlEntryTypeData( e_AtomLine,        51,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_NotUsed1(                 AddEnum( "Atom_NotUsed1",                 MdlEntryTypeData( e_AtomLine,        54,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_NotUsed2(                 AddEnum( "Atom_NotUsed2",                 MdlEntryTypeData( e_AtomLine,        57,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_AtomMappingNumber(        AddEnum( "Atom_AtomMappingNumber",        MdlEntryTypeData( e_AtomLine,        60,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_InversionFlag(            AddEnum( "Atom_InversionFlag",            MdlEntryTypeData( e_AtomLine,        63,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_ExactChangeFlag(          AddEnum( "Atom_ExactChangeFlag",          MdlEntryTypeData( e_AtomLine,        66,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_FirstAtomIndex(           AddEnum( "Bond_FirstAtomIndex",           MdlEntryTypeData( e_BondLine,         0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Bond_SecondAtomIndex(          AddEnum( "Bond_SecondAtomIndex",          MdlEntryTypeData( e_BondLine,         3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Bond_Type(                     AddEnum( "Bond_Type",                     MdlEntryTypeData( e_BondLine,         6,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_Stereo(                   AddEnum( "Bond_Stereo",                   MdlEntryTypeData( e_BondLine,         9,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_NotUsed(                  AddEnum( "Bond_NotUsed",                  MdlEntryTypeData( e_BondLine,        12,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_Topology(                 AddEnum( "Bond_Topology",                 MdlEntryTypeData( e_BondLine,        15,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_ReactingCenterStatus(     AddEnum( "Bond_ReactingCenterStatus",     MdlEntryTypeData( e_BondLine,        18,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      RXN_RXNStartLine(              AddEnum( "RXN_RXNStartLine",              MdlEntryTypeData( e_RXNStartLine,     0,   4,   "$RXN", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXN_MolStartLine(              AddEnum( "RXN_MolStartLine",              MdlEntryTypeData( e_RXNMolStartLine,  0,   4,   "$MOL", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXN_RXNTerminationLine(        AddEnum( "RXN_TerminationLine",           MdlEntryTypeData( e_RXNTerminationLine, 0, 0,       "", util::Format().R(),         util::CPPDataTypes::e_String))),
      Terminator(                    AddEnum( "Terminator",                    MdlEntryTypeData( e_TerminationLine,  0,   4,   "$$$$", util::Format().R(),         util::CPPDataTypes::e_String)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MdlEntryTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief global access to all entry types
    const MdlEntryTypes &GetMdlEntryTypes()
    {
      return MdlEntryTypes::GetEnums();
    }

  } // namespace sdf

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< sdf::MdlEntryTypeData, sdf::MdlEntryTypes>;

  } // namespace util
} // namespace bcl
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
#include "sdf/bcl_sdf_mdl_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief standard constructor
    MdlHandler::MdlHandler() :
      m_IsValid( false),
      m_Molfile(),
      m_MiscProperties(),
      m_WasParsed( false),
      m_Lines()
    {
    }

    //! @brief standard constructor
    MdlHandler::MdlHandler( std::istream &ISTREAM) :
      m_IsValid( false),
      m_Molfile(),
      m_MiscProperties(),
      m_WasParsed( false),
      m_Lines()
    {
      ReadFromSDF( ISTREAM);
    }

    MdlHandler::MdlHandler( const storage::List< std::string> &LINES) :
      m_IsValid( false),
      m_Molfile(),
      m_MiscProperties(),
      m_WasParsed( false),
      m_Lines( LINES)
    {
    }

    //! @brief virtual copy constructor
    MdlHandler *MdlHandler::Clone() const
    {
      return new MdlHandler( *this);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MdlHandler::ReadFromSDF( std::istream &ISTREAM)
    {
      Reset();

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        BCL_MessageCrt( "Error reading SDFile: passed bad istream");
        return ISTREAM;
      }

      // buffer a single molecule, up to a $$$$ delimiter
      bool contains_data( false);
      bool found_term( false);
      while( !ISTREAM.eof())
      {
        m_Lines.PushBack( std::string());
        std::string &last_line( m_Lines.LastElement());
        std::getline( ISTREAM, m_Lines.LastElement());
        if( !contains_data && ContainsNonspaceCharacters( last_line))
        {
          contains_data = true;
        }
        if( IsTerminalLine( last_line))
        {
          found_term = true;
          break;
        }
      }

      if( !contains_data)
      {
        // everything is blank, just reset everything and exit silently
        m_Lines.Reset();
        m_WasParsed = true;
        m_IsValid = false;
      }
      else
      {
        if( !found_term)
        {
          // no terminator was found, print out an error and reset everything because the data is incomplete
          BCL_MessageCrt( "Unexpected end of SDF, expected $$$$; not parsing this molecule");
          m_IsValid = false;
          m_WasParsed = true;
          m_Lines.Reset();
        }
      }

      return ISTREAM;
    }

    void MdlHandler::FinalizeParsing() const
    {
      if( m_WasParsed || m_Lines.IsEmpty())
      {
        return;
      }

      // set this immediately so that we don't do it again through for some strange reason
      m_WasParsed = true;

      storage::List< std::string>::const_iterator itr_lines( m_Lines.Begin()), itr_lines_end( m_Lines.End());
      
      // Read the molfile portion of the file
      itr_lines = m_Molfile.ReadMolfile( itr_lines, itr_lines_end);
      if( !m_Molfile.IsValid())
      {
        BCL_MessageCrt( "Molfile portion of SDF was faulty, not reading this molecule");
        return;
      }

      // read misc properties
      for( ; itr_lines != itr_lines_end; ++itr_lines)
      {
        // skip empty lines
        if( !ContainsNonspaceCharacters( *itr_lines))
        {
          continue;
        }

        // if end of mdl block is reached
        if( IsTerminalLine( *itr_lines))
        {
          m_IsValid = true;
          break;
        }

        const std::string data_label( GetMDLDataLabel( *itr_lines));

        // check that label is defined
        if( data_label.empty())
        {
          BCL_MessageCrt
          (
            "Warning: Blank line did not contain a data label (i.e. >  <LABEL>), but one was expected; "
          );

          // next line
          continue;
        }

        // get a reference on the value string from the misc properties
        std::string &value( m_MiscProperties[ data_label]);

        ++itr_lines;

        // save all lines associated with this mdl property
        while( itr_lines != itr_lines_end && ContainsNonspaceCharacters( *itr_lines))
        {
          // check for terminal line in case someone forgot a blank line after the property
          if( IsTerminalLine( *itr_lines))
          {
            m_IsValid = true;
            break;
          }

          // check if another data element was given, i.e. '>  <something>'
          // if so it should be separated from the previous one with a blank line, so print 
          // a message to the user.  
          if( !GetMDLDataLabel( *itr_lines).empty())
          {
            BCL_MessageCrt
            (
              "Property with label \"" + GetMDLDataLabel( *itr_lines) +
              "\" should be separated from property \"" + data_label + "\" with a new line" 
            );
            break; // break out of this loop so the next property can be read
          }

          // add the data to the MDL property
          value += *( itr_lines++);
          value += '\n';
        }
      } // misc properties reading

      // *itr_lines should be '$$$$' here
      if( itr_lines == itr_lines_end || !IsTerminalLine( *itr_lines))
      {
        BCL_MessageStd( "Unexpected end of SDF file.  No '$$$$' found");
        return;
      }
    }

    //! @brief write to std::ostream in mdl format
    std::ostream &MdlHandler::WriteToSDF( std::ostream &OSTREAM) const
    {
      return this->WriteToSDF( OSTREAM, m_Molfile, m_MiscProperties);
    }

    //! @brief write to std::ostream in mdl format
    std::ostream &MdlHandler::WriteToSDF
    (
      std::ostream &OSTREAM,
      const MolfileHandler &MOLFILE,
      const storage::Map< std::string, std::string> &MISC_PROPERTIES
    )
    {
      if( !MOLFILE.IsValid())
      {
        BCL_MessageStd( "Cannot write SDF formatted molecule, invalid molfile provided");
        return OSTREAM;
      }

      // write into a temporary buffer.  if something goes wrong we don't want to output
      // a fragmented file
      std::ostringstream ostream;
      MOLFILE.WriteMolfile( ostream, GetAddAtomMdlLineFlag()->GetFlag());
      WriteMiscProperties( ostream, MISC_PROPERTIES);
      ostream << GetDefaultLine( e_TerminationLine) << '\n'; // ending $$$$ delimiter
      OSTREAM << ostream.str();
      return OSTREAM;
    }

    //! @brief write to std::ostream in mdl format
    std::ostream &MdlHandler::WriteToSDF
    (
      std::ostream &OSTREAM,
      const std::string &DESCRIPTION,
      const storage::Vector< AtomInfo> &ATOM_INFO,
      const storage::Vector< BondInfo> &BOND_INFO,
      const storage::Map< std::string, std::string> &MISC_PROPERTIES
    )
    {
      MolfileHandler molfile( DESCRIPTION, ATOM_INFO, BOND_INFO);
      if( !molfile.IsValid())
      {
        BCL_MessageCrt( "Could not make a valid molfile to write to SDF file");
        return OSTREAM;
      }

      return WriteToSDF( OSTREAM, molfile, MISC_PROPERTIES);
    }

    //! @brief write a simple string that should be unique for a constitution, independent of H
    //! @note hash is dependent on ordering of atoms; use accordingly
    std::string MdlHandler::CreateConstitutionalHashString
    (
      const storage::Vector< AtomInfo> &ATOM_INFO,
      const storage::Vector< BondInfo> &BOND_INFO,
      chemistry::ConfigurationalBondTypeData::Data BOND_TYPE
    )
    {
      std::ostringstream stream;

      storage::Vector< size_t> heavy_atom_index( ATOM_INFO.GetSize(), util::GetUndefined< size_t>());
      size_t heavy_atom_counter( 0);
      for( size_t atom_number( 0), number_atoms( ATOM_INFO.GetSize()); atom_number < number_atoms; ++atom_number)
      {
        const chemistry::AtomType &atom_type( ATOM_INFO( atom_number).GetAtomType());
        // skip hydrogens
        if( atom_type->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        stream << atom_type->GetNumberBonds();
        stream << atom_type->GetElementType()->GetChemicalSymbol();
        stream << atom_type->GetFormalCharge() << ',';

        // map the heavy atom index so that bonds can be remapped accordingly
        heavy_atom_index( atom_number) = heavy_atom_counter++;
      }

      // write out the bonds
      util::Format index_formatter;
      index_formatter.W( 3);
      for
      (
        storage::Vector< BondInfo>::const_iterator
          itr_bond( BOND_INFO.Begin()), itr_bond_end( BOND_INFO.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // skip bonds to H
        if
        (
          !util::IsDefined( heavy_atom_index( itr_bond->GetAtomIndexLow()))
          || !util::IsDefined( heavy_atom_index( itr_bond->GetAtomIndexHigh()))
        )
        {
          continue;
        }
        stream << index_formatter( heavy_atom_index( itr_bond->GetAtomIndexLow()))
               << index_formatter( heavy_atom_index( itr_bond->GetAtomIndexHigh()))
               << itr_bond->GetConfigurationalBondType()->GetBondData( BOND_TYPE);
      }
      return stream.str();
    }

    //! @brief write a simple string that should be unique for this molecular configuration
    //! @note hash is dependent on ordering of atoms, etc.; use accordingly
    std::string MdlHandler::CreateConfigurationalHashString
    (
      const storage::Vector< AtomInfo> &ATOM_INFO,
      const storage::Vector< BondInfo> &BOND_INFO,
      chemistry::ConfigurationalBondTypeData::Data BOND_TYPE
    )
    {
      // create a stream, initially containing the constitutional information
      std::ostringstream stream;

      util::Format index_formatter;
      index_formatter.W( 3);

      size_t heavy_atom_index( 0);

      // add on the chirality info for the relevant atoms
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_INFO.Begin()), itr_end( ATOM_INFO.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }
        if( itr->GetChirality() != chemistry::e_NonChiral) // only write out chiral centers
        {
          stream << index_formatter( heavy_atom_index);
          stream << chemistry::ChiralityEnum( itr->GetChirality());
        }
        ++heavy_atom_index;
      }

      // add the double bond isometry for the relevant bonds
      stream << ' ';
      for
      (
        storage::Vector< BondInfo>::const_iterator itr_bond( BOND_INFO.Begin()), itr_bond_end( BOND_INFO.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // skip bonds that could not have isometry
        if( itr_bond->GetConfigurationalBondType()->GetIsometry() != chemistry::e_NonIsometric)
        {
          // push back the next bond isometry
          if( !itr_bond->GetConfigurationalBondType()->IsBondInRing())
          {
            stream << chemistry::BondIsometryEnum( itr_bond->GetConfigurationalBondType()->GetIsometry());
          }
        }
      }

      stream << " " << CreateConstitutionalHashString( ATOM_INFO, BOND_INFO, BOND_TYPE);
      return stream.str();
    }

    //! @brief write a simple string that should be unique for this molecular conformation
    //! @note hash is dependent on orientation of molecule, ordering of atoms, etc.; use accordingly
    std::string MdlHandler::CreateConformationalHashString
    (
      const storage::Vector< AtomInfo> &ATOM_INFO,
      const storage::Vector< BondInfo> &BOND_INFO
    )
    {
      // determine whether the conformation is actually present
      // determine whether more than 1 coordinate is undefined / 0, in which case all
      // isometry information must be written as well
      // If the positions are all defined, and no more than 1 atom is at the origin, then
      // the isometry is given by the positions, so there is no need to include configurational information
      // in the hash
      {
        size_t zero_position_count( 0);
        for
        (
          storage::Vector< AtomInfo>::const_iterator
            itr_atom( ATOM_INFO.Begin()), itr_atom_end( ATOM_INFO.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // keep track of whether the coordinates are all defined and that no more than one of them is at the origin
          if( !itr_atom->GetCoordinates().IsDefined())
          {
            break;
          }
          else if( math::EqualWithinAbsoluteTolerance( itr_atom->GetCoordinates().Norm(), double( 0.0), 1.0e-4))
          {
            if( ++zero_position_count == 2)
            {
              break;
            }
          }
        }
        if( zero_position_count > 1)
        {
          return CreateConfigurationalHashString( ATOM_INFO, BOND_INFO);
        }
      }

      // create a stream, initially containing the constitutional information
      // it is not necessary to add the configurational information, since that is given by the positions
      std::ostringstream stream;
      stream << CreateConstitutionalHashString( ATOM_INFO, BOND_INFO) << ' ';

      static const util::Format coordinates_formatter( util::Format().FFP( 3).R()); // only print out to 3 precision
      for
      (
        storage::Vector< AtomInfo>::const_iterator
          itr_atom( ATOM_INFO.Begin()), itr_atom_end( ATOM_INFO.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // skip hydrogens
        if( itr_atom->GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        // only write out coordinates to three positions to avoid numerical roundoff issues
        stream << coordinates_formatter( itr_atom->GetCoordinates().X()) << ' ';
        stream << coordinates_formatter( itr_atom->GetCoordinates().Y()) << ' ';
        stream << coordinates_formatter( itr_atom->GetCoordinates().Z()) << ' ';
      }
      return stream.str();
    }

    //! @brief read MdlHandler object from std::istream
    //! @param ISTREAM istream that contains MdlHandler object
    //! @return istream after MdlHandler object was extracted
    std::istream &MdlHandler::Read( std::istream &ISTREAM)
    {
      return ReadFromSDF( ISTREAM);
    }

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    util::ShPtr< command::FlagInterface> &MdlHandler::GetAddAtomMdlLineFlag()
    {
      static util::ShPtr< command::FlagInterface> s_atom_type_add
        (
          new command::FlagStatic( "add_atom_type", "add atom types in mdl line while writing out")
        );

      return s_atom_type_add;
    }

    //! @brief add hydrogen handling preferences to the command line flag
    //! @param CMD command to add the hydrogen handling preference flags to
    void MdlHandler::AddAtomMdlLineFlag( command::Command &CMD)
    {
      CMD.AddFlag( GetAddAtomMdlLineFlag());
    }

    //! @brief write MdlHandler into std::ostream
    //! @param OSTREAM ostream that gets MdlHandler object
    //! @param INDENT indentation
    //! @return ostream after MdlHandler object was inserted
    std::ostream &MdlHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool MdlHandler::ContainsNonspaceCharacters( const std::string &STRING)
    {
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( !isspace( *itr))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief write mdl lines into std::ostream
    //! @param OSTREAM ostream that gets MdlHandler object
    void MdlHandler::WriteMiscProperties
    (
      std::ostream &OSTREAM,
      const storage::Map< std::string, std::string> &MISC_PROPERTIES
    )
    {
      static const std::string pre_property_name_str( "> <");  // goes before each property's name
      static const std::string post_property_name_str( ">\n"); // goes after each property's name

      // iterate over all MDL misc property lines
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          itr_map( MISC_PROPERTIES.Begin()),
          itr_map_end( MISC_PROPERTIES.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->first.size() != 0 && itr_map->second.size() != 0) // the property has a name and value
        {
          OSTREAM << pre_property_name_str        // property name start delimiter
                  << itr_map->first               // property name
                  << post_property_name_str;      // property name deliminater

          const std::string &value( itr_map->second);

          // write the value of the given property line
          OSTREAM << value;

          // if the last character of the string was not a new-line, add one here
          if( value[ value.size() - 1] != '\n')
          {
            OSTREAM << '\n'; // followed by a newline
          }
          OSTREAM << '\n';     // put a blank line follows the last line of the property value
        }
      }
    }

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains $$$$
    bool MdlHandler::IsTerminalLine( const std::string &LINE)
    {
      static const std::string s_default_terminal_line( GetDefaultLine( e_TerminationLine));
      return util::StartsWith( LINE, s_default_terminal_line);
    }

    //! @brief return the datalable, if line is a datalabel line
    //! @param LINE line from mdl section
    //! @return string that constains datalable, string will be empty for non-data lable lines
    std::string MdlHandler::GetMDLDataLabel( const std::string &LINE)
    {
      // data label delimiter left
      static const char s_misc_property_delimiter_left( '<');

      // data label delimiter right
      static const char s_misc_property_delimiter_right( '>');

      // handle empty lines and lines that do not start with <
      if( LINE.empty() || LINE[ 0] != s_misc_property_delimiter_right)
      {
        return std::string();
      }

      // find label start and end
      const std::string::size_type label_start( LINE.find( s_misc_property_delimiter_left, 1));
      const std::string::size_type label_end( LINE.rfind( s_misc_property_delimiter_right));

      // return an empty string for non data label line
      if( label_start == std::string::npos || label_end == std::string::npos)
      {
        return std::string();
      }

      // misc property name
      return LINE.substr( label_start + 1, label_end - label_start - 1);
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MdlHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new MdlHandler())
    );

  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_mdl_header.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief default constructor
    MdlHeader::MdlHeader() :
      m_NumberAtoms( 0),
      m_NumberBonds( 0),
      m_IsValid( false)
    {
    }

    //! @brief constructor from # atoms, # bonds
    MdlHeader::MdlHeader
    (
      const size_t &NUMBER_ATOMS,
      const size_t &NUMBER_BONDS
    ) :
      m_NumberAtoms( NUMBER_ATOMS),
      m_NumberBonds( NUMBER_BONDS),
      m_IsValid( true)
    {
    }

    //! @brief virtual copy constructor
    MdlHeader *MdlHeader::Clone() const
    {
      return new MdlHeader( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MdlHeader::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the string of the header
    //! @return the string of the header
    std::string MdlHeader::ToMdlLine() const
    {
      std::string mdl_header( GetDefaultLine( e_HeaderLine));
      GetMdlEntryTypes().Header_NumberAtomLines->Set( mdl_header, m_NumberAtoms);
      GetMdlEntryTypes().Header_NumberBondLines->Set( mdl_header, m_NumberBonds);
      return mdl_header;
    }

    //! @brief return the string of the header
    //! @return the string of the header
    void MdlHeader::SetFromMdlLine( const std::string &MDL_HEADER, const size_t NUMBER_DESC_LINES)
    {
      // reinitialize
      m_IsValid = false;
      m_NumberAtoms = m_NumberBonds = 0;

      // first check that the size is correct
      static const size_t s_min_size( GetDefaultLine( e_HeaderLine).size());

      if( MDL_HEADER.size() >= s_min_size)
      {
        // mdl V3000 format is not supported, so check that the header is == the default (V2000)
        m_IsValid =
          GetMdlEntryTypes().Header_CtabVersion->GetTrimmedString( MDL_HEADER)
          == GetMdlEntryTypes().Header_CtabVersion->GetDefault();
        if( m_IsValid)
        {
          m_NumberAtoms = GetMdlEntryTypes().Header_NumberAtomLines->GetUnsignedInt( MDL_HEADER);
          m_NumberBonds = GetMdlEntryTypes().Header_NumberBondLines->GetUnsignedInt( MDL_HEADER);
        }
      }
      else if
      (
        NUMBER_DESC_LINES == size_t( 3) && MDL_HEADER.size() >= 6
        && isdigit( MDL_HEADER[ 2]) && isdigit( MDL_HEADER[ 5])
        && ( isdigit( MDL_HEADER[ 1]) || ( MDL_HEADER[ 1] == ' ' && MDL_HEADER[ 0] == ' '))
        && ( isdigit( MDL_HEADER[ 0]) || MDL_HEADER[ 0] == ' ')
        && ( isdigit( MDL_HEADER[ 4]) || ( MDL_HEADER[ 4] == ' ' && MDL_HEADER[ 3] == ' '))
        && ( isdigit( MDL_HEADER[ 3]) || MDL_HEADER[ 3] == ' ')
      )
      {
        // we've already read 3 description lines, test whether it looks enough like a header line to proceed
        m_IsValid = true;
        m_NumberAtoms = GetMdlEntryTypes().Header_NumberAtomLines->GetUnsignedInt( MDL_HEADER);
        m_NumberBonds = GetMdlEntryTypes().Header_NumberBondLines->GetUnsignedInt( MDL_HEADER);
      }
    }

    //! @brief read MdlHeader object from std::istream
    //! @param ISTREAM istream that contains MdlHeader object
    //! @return istream after MdlHeader object was extracted
    std::istream &MdlHeader::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_NumberAtoms, ISTREAM);
      io::Serialize::Read( m_NumberBonds, ISTREAM);
      io::Serialize::Read( m_IsValid,     ISTREAM);
      // return
      return ISTREAM;
    }

    //! @brief write MdlHeader into std::ostream
    //! @param OSTREAM ostream that gets MdlHeader object
    //! @param INDENT indentation
    //! @return ostream after MdlHeader object was inserted
    std::ostream &MdlHeader::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_NumberAtoms, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NumberBonds, OSTREAM, 0) << '\t';
      io::Serialize::Write( m_IsValid,     OSTREAM, 0);
      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MdlHeader::s_Instance
    (
      GetObjectInstances().AddInstance( new MdlHeader())
    );

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_mdl_line_types.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief LineType as string
    //! @param LINE_TYPE the line type for which a name is desired
    //! @return the line type as string
    const std::string &GetLineTypeName( const MdlLineType &LINE_TYPE)
    {
      static const std::string s_names[] =
      {
        "Header",
        "Atom",
        "Bond",
        "Terminal",
        "RXNHeader",
        "RXNStart",
        "RXNMolStart",
        "RXNTermination",
        GetStaticClassName< MdlLineType>()
      };
      return s_names[ LINE_TYPE];
    }

    //! @brief create a vector containing the default lines for each line type
    //! @return a vector containing the default lines for each line type
    storage::Vector< std::string> GetDefaultLines()
    {
      // length of each line
      storage::Vector< size_t> line_length( s_NumberLineTypes, size_t( 0));

      // iterate over all entry types
      for
      (
        MdlEntryTypes::const_iterator itr( GetMdlEntryTypes().Begin()),
        itr_end( GetMdlEntryTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        const MdlLineType type( ( *itr)->GetMdlLineType());

        // get the position this entry ends at
        const size_t end_pos( ( *itr)->GetStart() + ( *itr)->GetLength());

        // set line length to the max of this end_pos and the existing end_pos
        line_length( size_t( type)) = std::max( line_length( size_t( type)), end_pos);
      }

      storage::Vector< std::string> lines( s_NumberLineTypes);

      // initialize all line types with their length
      for( size_t line_type_number( 0); line_type_number < s_NumberLineTypes; ++line_type_number)
      {
        // create a default string full of spaces
        lines( line_type_number) = std::string( line_length( line_type_number), ' ');
      }

      // iterate over all entry types; set the default value in the corresponding string
      for
      (
        MdlEntryTypes::const_iterator itr( GetMdlEntryTypes().Begin()), itr_end( GetMdlEntryTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        const MdlEntryType &type( *itr);
        type->Set( lines( type->GetMdlLineType()), type->GetDefault());
      }
      return lines;
    }

    //! @brief gather all entry types for given line type
    //! @param LINE_TYPE line type the entry types should be collected for
    //! @return set of entry types for that line type
    const std::string &GetDefaultLine( const MdlLineType &LINE_TYPE)
    {
      static const storage::Vector< std::string> default_lines( GetDefaultLines());
      if( LINE_TYPE == s_NumberLineTypes)
      {
        static const std::string undefined_line;
        return undefined_line;
      }
      return default_lines( size_t( LINE_TYPE));
    }

  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_mdl_property.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief property as string
    //! @param PROPERTY the property desired
    //! @return the property as string
    const std::string &MdlProperty::GetPropertyName( const MdlProperty::Property &PROPERTY)
    {
      static const std::string s_names[] =
      {
        "CHG",
        "RBC",
        "ALS",
        "BCL ATM",
        "BCL BND",
        "BCL CHI",
        "BCL DBI",
        "BCL ARO", // atom aromaticity
        "BCL RXN",
        "END",
        GetStaticClassName< Property>()
      };
      return s_names[ PROPERTY];
    }

    //! @brief get the possible property prefixes
    const storage::Vector< std::string> &MdlProperty::GetAllPropertyPrefixes()
    {
      static const storage::Vector< std::string> s_prefixes = storage::Vector< std::string>::Create( "M  ", "G  ", "A  ", "V  ");
      return s_prefixes;
    }

    //! @brief property as string, with M  prefix
    //! @param PROPERTY the property desired
    //! @return the property as string with the M  prefix
    const std::string &MdlProperty::GetPropertyNameWithPrefix( const MdlProperty::Property &PROPERTY)
    {
      static const std::string s_names[] =
      {
        "M  " + GetPropertyName( e_Charge),
        "M  " + GetPropertyName( e_RingBondCount),
        "M  " + GetPropertyName( e_AtomList),
        "M  " + GetPropertyName( e_BclAtomType),
        "M  " + GetPropertyName( e_BclBondType),
        "M  " + GetPropertyName( e_BclChirality),
        "M  " + GetPropertyName( e_BclDoubleBondIsometry),
        "M  " + GetPropertyName( e_BclAtomAromaticity),
        "M  " + GetPropertyName( e_BclRXNProperty),
        "M  " + GetPropertyName( e_BlockTerminator),
        GetStaticClassName< Property>()
      };
      return s_names[ PROPERTY];
    }

    //! @brief fixed width of property
    //! @param PROPERTY the property
    //! @return the properties fixed width (undefined if property is not fixed width)
    const size_t &MdlProperty::GetFixedWidthSize( const MdlProperty::PropertyEnum &PROPERTY)
    {
      static const size_t s_fixed_widths[] =
      {
        3, // CHG
        3, // RBC
        util::GetUndefined< size_t>(), // ALS
        util::GetUndefined< size_t>(), // BCL ATM
        util::GetUndefined< size_t>(), // BCL BND
        3, // BCL CHI
        1, // BCL DBI
        util::GetUndefined< size_t>(), // atom aromaticity
        util::GetUndefined< size_t>(), // BCL RXN
        util::GetUndefined< size_t>() // END
      };
      return s_fixed_widths[ PROPERTY];
    }

    //! @brief get the multiplicity of the property (e.g. the number of values per property)
    //! @param PROPERTY the property desired
    //! @return the multiplicity
    const size_t &MdlProperty::GetMultiplicity( const MdlProperty::PropertyEnum &PROPERTY)
    {
      static const size_t s_multiplities[] =
      {
        2, // atom index (1 offset) and charge
        2, // atom index (1 offset) and number of ring bonds
        1, // atom index (1 offset) followed by list of element types
        1, // atom type
        1, // bond type
        2, // atom index (1 offset) and chirality
        1, // just isomorphism string
        2, // list of aromatic atoms
        1, // rxn properties
        0, // no optional strings
        0
      };
      return s_multiplities[ PROPERTY];
    }

    //! @brief determine whether the # of values should be the first number on the line
    //! @param PROPERTY the property desired
    //! @return whether to print the size in the first field
    bool MdlProperty::FirstFieldIsSize( const PropertyEnum &PROPERTY)
    {
      static const bool s_first_field_is_size[] =
      {
        true,  // e_Charge
        true,  // e_RingBondCount
        false, // e_AtomList
        false, // e_BclAtomType
        false, // e_BclBondType
        false, // e_BclChirality
        false, // e_BclDoubleBondIsometry
        false, // e_BclAtomAromaticity
        false, // e_BclRXNProperty
        false, // e_BlockTerminator
        false  // s_NumberProperties
      };
      return s_first_field_is_size[ PROPERTY];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MdlProperty::MdlProperty( const PropertyEnum &PROPERTY) :
      m_Property( PROPERTY)
    {
      SetLabel( PROPERTY);
    }

    //! @brief constructor given atom lines, bond lines, and the property desired
    MdlProperty::MdlProperty
    (
      const PropertyEnum &PROPERTY,
      const storage::Vector< AtomInfo> &ATOM_LINES,
      const storage::Vector< BondInfo> &BOND_LINES
    ) :
      m_Property( PROPERTY)
    {
      SetLabel( PROPERTY);

      // choose the function based on the property
      switch( m_Property)
      {
        case e_Charge:
          RetrieveCharges( ATOM_LINES);
          break;
        case e_RingBondCount:
          break;
        case e_AtomList:
          break;
        case e_BclAtomType:
          RetrieveAtomTypes( ATOM_LINES);
          break;
        case e_BclBondType:
          RetrieveBondTypes( BOND_LINES);
          break;
        case e_BclChirality:
          RetrieveChiralities( ATOM_LINES);
          break;
        case e_BclDoubleBondIsometry:
          RetrieveDoubleBondIsometries( BOND_LINES);
          break;
        case e_BclAtomAromaticity:
          break; // read-only from a Molfile now, can't do anything with atom/bond lines
        case e_BclRXNProperty:
          //RetrieveRXNProperties( ATOM_LINES, BOND_LINES);
          //break;
        case e_BlockTerminator:
        case s_NumberProperties:
        default:
          break;
      }
    }

    //! @brief construct from string containing property
    MdlProperty::MdlProperty( const std::string &STRING) :
      m_Property( s_NumberProperties)
    {
      SetFromString( STRING);
    }

    //! @brief Clone function
    //! @return pointer to new MdlProperty
    MdlProperty *MdlProperty::Clone() const
    {
      return new MdlProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MdlProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the string of the property
    //! @return the string of the property
    std::string MdlProperty::GetString() const
    {

      // unknown/generic property
      if( m_Property == s_NumberProperties)
      {
        //return std::string();

        // write out the properties separated by a space
        std::stringstream strm;
        strm << m_PropertyLabel;
        if( !m_PropertyLabel.empty() && !m_PropertyStrings.IsEmpty())
        {
          strm << " ";
          for( size_t p( 0), end_p( m_PropertyStrings.GetSize() - 1); p < end_p; ++p)
          {
            strm << m_PropertyStrings( p) << " ";
          }
          strm << m_PropertyStrings( m_PropertyStrings.GetSize() - 1);
        }
        return strm.str();
      }

      // if the property string requires values (e.g. multiplicity != 0), but none are available, then
      // just return
      const size_t multiplicity( GetMultiplicity( m_Property));
      if( multiplicity != size_t( 0) && m_PropertyStrings.IsEmpty())
      {
        return std::string();
      }

      std::string property_string( GetPropertyNameWithPrefix( m_Property));

      const size_t fixed_width( GetFixedWidthSize( m_Property));

      util::Format string_formatter;
      if( util::IsDefined( fixed_width))
      {
        string_formatter.W( fixed_width);
      }

      if( FirstFieldIsSize( m_Property))
      {
        property_string += string_formatter( m_PropertyStrings.GetSize() / multiplicity);
      }
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr
      )
      {
        property_string += ' ';
        property_string += string_formatter( *itr);
      }
      return property_string;
    }

    //! @brief report if the property is one recognized by the BCL
    //! @return true if it is a property explicitly handled by the BCL
    bool MdlProperty::IsBCLProperty() const
    {
      return m_Property < s_NumberProperties;
    }

    //! @brief determine if an MDL property is one recognized by the BCL
    //! @param LINE the line to examine
    //! @return true if the property line contains information the BCL can explicitly recognize
    bool MdlProperty::IsBCLPropertyLine( const std::string &LINE)
    {
      bool is_bcl_prop( false);
      if( util::StartsWith( LINE, "M  ")) // all BCL properties start with 'M  ' currently
      {

        // determine whether the property is one of the known properties
        Property test_property;
        for( size_t property_id( 0); property_id < s_NumberProperties; ++property_id)
        {
          test_property = Property( property_id);
          if( util::StartsWith( LINE, GetPropertyNameWithPrefix( test_property)))
          {
            // set variable if the property matches one known to the BCL
            is_bcl_prop = true;
            break;
          }
        }
      }

      return is_bcl_prop;
    }

    //! @brief helper function to determine whether a line is an MDL property line without copying or fully parsing it
    //! @param LINE the line of interest
    //! @return true if the line appears to be an MDL property line
    bool MdlProperty::IsMDLPropertyLine( const std::string &LINE)
    {
      // ensure that the property line looks like a normal starting property line
      const storage::Vector< std::string> &prefixes( GetAllPropertyPrefixes());
      bool has_good_prefix( false);
      for( size_t p( 0), end_p( prefixes.GetSize()); p < end_p; ++p)
      {
        //if( !util::StartsWith( LINE, "M  "))
        if( util::StartsWith( LINE, prefixes( p).c_str()))
        {
          has_good_prefix = true;
          break;
        }
      }
      return has_good_prefix;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief apply the property to a vector of atom and bond lines
    //! @param ATOM_LINES Mdl Atoms on which to apply the property
    //! @param BOND_LINES Mdl Bonds on which to apply the property
    void MdlProperty::ApplyProperty( storage::Vector< AtomInfo> &ATOM_LINES, storage::Vector< BondInfo> &BOND_LINES) const
    {
      // choose the function based on the property
      switch( m_Property)
      {
        case e_Charge:
          ApplyCharges( ATOM_LINES);
          break;
        case e_BclAtomType:
          ApplyAtomTypes( ATOM_LINES);
          break;
        case e_BclBondType:
          ApplyBondTypes( BOND_LINES);
          break;
        case e_BclChirality:
          ApplyChiralities( ATOM_LINES);
          break;
        case e_BclDoubleBondIsometry:
          ApplyDoubleBondIsometries( BOND_LINES);
          break;
        case e_BclRXNProperty:
        case e_BlockTerminator:
        case s_NumberProperties:
        default:
          break;
      }
    }

    //! @brief read MdlProperty object from std::istream
    //! @param ISTREAM istream that contains MdlProperty object t
    //! @return istream after MdlProperty object was extracted
    std::istream &MdlProperty::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Property, ISTREAM);
      io::Serialize::Read( m_PropertyStrings, ISTREAM);
      // return
      return ISTREAM;
    }

    //! @brief write MdlProperty into std::ostream
    //! @param OSTREAM ostream that gets MdlProperty object
    //! @param INDENT indentation
    //! @return ostream after MdlProperty object was inserted
    std::ostream &MdlProperty::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Property, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PropertyStrings, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the label according to the property number
    //! @details looks up the full label if a known property, e.g. "M  END", or leaves it blank
    void MdlProperty::SetLabel( const PropertyEnum &PROPERTY)
    {
      if( PROPERTY < s_NumberProperties)
      {
        m_PropertyLabel = GetPropertyNameWithPrefix( PROPERTY);
      }
    }

    //! @brief return the string of the header
    //! @return the string of the header
    void MdlProperty::SetFromString( const std::string &MDL_STRING)
    {
      m_Property = s_NumberProperties;
      m_PropertyLabel = std::string();
      m_PropertyStrings.Reset();

      // ensure that the property line looks like a normal starting property line
      // test all possible prefixes to see if we match any of them
      // all prefixes are 3 characters long, e.g. 'M  '
      const storage::Vector< std::string> &possible_prefixes( GetAllPropertyPrefixes());
      std::string str_prefix( MDL_STRING.substr( 0, 3));
      bool found_prefix( false);
      for( size_t pre_no( 0), end_no( possible_prefixes.GetSize()); pre_no < end_no; ++pre_no)
      {
        if( str_prefix == possible_prefixes( pre_no))
        {
          found_prefix = true;
          break;
        }
      }
      if( !found_prefix)
      {
        BCL_MessageStd( "String \"" + MDL_STRING + "\" does not contain a valid MDL property prefix; ignoring");
        return;
      }

      const size_t property_name_start( 3);
      size_t property_values_start( 0);

      // determine whether the property is one of the known properties
      bool is_known_property( false);
      Property test_property;
      for( size_t property_id( 0); property_id < s_NumberProperties; ++property_id)
      {
        test_property = Property( property_id);
        const std::string &property_string( GetPropertyName( test_property));
        if( MDL_STRING.substr( property_name_start, property_string.size()) == property_string)
        {
          property_values_start = property_name_start + property_string.size();
          m_Property = test_property;
          is_known_property = true;
          break;
        }
      }

      if( is_known_property)
      {

        // store the label, including the 'M  ' portion
        m_PropertyLabel = MDL_STRING.substr( 0, property_values_start);

        // if the multiplicity is 0, then there should be nothing else to read
        if( GetMultiplicity( m_Property) == size_t( 0))
        {
          return;
        }

        // determine width of this type
        const size_t fixed_width( GetFixedWidthSize( m_Property));

        if( !FirstFieldIsSize( m_Property)) // non-size fields have a space between property name and first value
        {
          ++property_values_start;
        }

        if( util::IsDefined( fixed_width))
        {
          for
          (
            size_t token_start( property_values_start), string_size( MDL_STRING.size());
            token_start < string_size;
            token_start += fixed_width + 1
          )
          {
            m_PropertyStrings.PushBack( util::TrimString( MDL_STRING.substr( token_start, fixed_width)));
          }
        }
        else
        {
          m_PropertyStrings = util::SplitString( MDL_STRING.substr( property_values_start));
        }

        if( FirstFieldIsSize( m_Property))
        {
          m_PropertyStrings.Remove( m_PropertyStrings.Begin());
        }

        // check multiplicity
        BCL_Assert
        (
          !( m_PropertyStrings.GetSize() % GetMultiplicity( m_Property)),
          util::Format()( m_PropertyStrings.GetSize()) + "values given for property: " + GetPropertyName( m_Property)
          + "; should have been a multiple of " + util::Format()( GetMultiplicity( m_Property))
        );
      }
      else
      {
        // unknown property, but we can keep the information anyway.
        // Assume the property label is space-separated from information, i.e. of the format 'X  YYY ..."

        // find the first space character in the string
        size_t label_pos( 3);
        size_t str_len( MDL_STRING.length());

        // skip immediate blank characters
        for( ; label_pos < str_len && isspace( MDL_STRING[ label_pos]); ++label_pos);
        size_t prop_name_start( label_pos);

        // find the first blank space character
        for( ; label_pos < str_len && !isspace( MDL_STRING[ label_pos]); ++label_pos);

        // if label is empty, stop now
        if( label_pos == prop_name_start)
        {
          return;
        }

        // set the property label and then split the remaining part of the line on spaces for the property strings
        m_PropertyLabel = MDL_STRING.substr( 0, label_pos);
        m_PropertyStrings = util::SplitString( MDL_STRING.substr( label_pos));
      }
    }

    //! @brief apply charges to atom lines
    //! @param ATOM_LINES mdl atoms to apply charges to
    void MdlProperty::ApplyCharges( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      // set the charges on each atom
      for
      (
        size_t charge_number( 0), number_charges( m_PropertyStrings.GetSize());
        charge_number < number_charges;
        charge_number += 2
      )
      {
        const size_t atom_index( util::ConvertStringToNumericalValue< size_t>( m_PropertyStrings( charge_number)) - 1);
        const short charge( util::ConvertStringToNumericalValue< short>( m_PropertyStrings( charge_number + 1)));
        BCL_Assert
        (
          atom_index < ATOM_LINES.GetSize(),
          "Tried to set charge on atom index " + util::Format()( atom_index) + " beyond end"
        );

        ATOM_LINES( atom_index).SetCharge( charge);
      }
    }

    //! @brief apply atom types to atom lines
    //! @param ATOM_LINES mdl atoms to apply atom types to
    void MdlProperty::ApplyAtomTypes( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      if( ATOM_LINES.GetSize() != m_PropertyStrings.GetSize())
      {
        BCL_MessageCrt
        (
          "Incorrect # of atom types specified!  Needed " + util::Format()( ATOM_LINES.GetSize()) + " but only " +
          util::Format()( m_PropertyStrings.GetSize()) + " were provided. Will leave atom types undefined for this molecule"
        );
        return;
      }
      size_t counter( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        chemistry::AtomType type( chemistry::GetAtomTypes().GetAtomType( *itr));
        BCL_Assert( type.IsDefined(), "Undefined atom type in atom type property: " + *itr);
        BCL_Assert
        (
          type->GetElementType() == ATOM_LINES( counter).GetAtomType()->GetElementType(),
          "Tried to change element type when setting atom type!"
        );
        ATOM_LINES( counter).SetAtomType( type);
      }
    }

    //! @brief apply bond types to bond lines
    //! @param BOND_LINES mdl bond to apply bond types to
    void MdlProperty::ApplyBondTypes( storage::Vector< BondInfo> &BOND_LINES) const
    {
      if( BOND_LINES.GetSize() != m_PropertyStrings.GetSize())
      {
        BCL_MessageCrt
        (
          "Incorrect # of bond types! Needed " + util::Format()( BOND_LINES.GetSize()) + " but " + util::Format()( m_PropertyStrings.GetSize()) +
          " were provided. Will leave bond types on this molecule undefined"
        );
        return;
      }
      size_t counter( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        chemistry::ConstitutionalBondType type( *itr);
        BCL_Assert( type.IsDefined(), "Undefined bond type in bond type property: " + *itr);
        BondInfo old( BOND_LINES( counter));
        BOND_LINES( counter) = BondInfo( old.GetAtomIndexLow(), old.GetAtomIndexHigh(), type);
        BOND_LINES( counter).SetIsometry( old.GetConfigurationalBondType()->GetIsometry());
      }
    }

    //! @brief apply chiralities to mdl atoms
    //! @param ATOM_LINES mdl atoms to apply chirality property to
    void MdlProperty::ApplyChiralities( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      // reset all chiralities to non-chiral
      for
      (
        storage::Vector< AtomInfo>::iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->SetChirality( chemistry::e_NonChiral);
      }
      // apply chiralities to each atom
      for
      (
        size_t chirality_number( 0), number_chiralities( m_PropertyStrings.GetSize());
        chirality_number < number_chiralities;
        chirality_number += 2
      )
      {
        // written index is 1-offset (for consistency with MDL format), so convert to 0-offset
        const size_t atom_index
        (
          util::ConvertStringToNumericalValue< size_t>( m_PropertyStrings( chirality_number)) - 1
        );
        const chemistry::Chirality chirality
        (
          chemistry::ChiralityEnum( m_PropertyStrings( chirality_number + 1))
        );
        BCL_Assert
        (
          atom_index < ATOM_LINES.GetSize(),
          "Tried to set chirality on atom index " + util::Format()( atom_index) + " beyond end"
        );
        BCL_Assert
        (
          chirality != chemistry::s_NumberChiralities,
          "Bad chirality name " + m_PropertyStrings( chirality_number + 1)
        );
        if( chirality == chemistry::e_UnknownChirality)
        {
          if( ATOM_LINES( atom_index).GetAtomType()->GetNumberBonds() == size_t( 4)
              || !ATOM_LINES( atom_index).GetAtomType()->IsGasteigerAtomType())
          {
            ATOM_LINES( atom_index).SetChirality( chirality);
          }
            continue;
        }
        BCL_Assert
        (
          ATOM_LINES( atom_index).GetAtomType()->GetNumberBonds() == size_t( 4)
          || !ATOM_LINES( atom_index).GetAtomType()->IsGasteigerAtomType(),
          "Tried to set chirality on atom that does not have 4 bonds! "
        );

        ATOM_LINES( atom_index).SetChirality( chirality);
      }
    }

    //! @brief apply double bond isometries to MdlBond with type ConfigurationalBondTypes().e_DoubleBond
    //! @param BOND_LINES mdl bonds to apply isometry property to
    void MdlProperty::ApplyDoubleBondIsometries( storage::Vector< BondInfo> &BOND_LINES) const
    {
      storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
      storage::Vector< BondInfo>::iterator itr_bond( BOND_LINES.Begin()), itr_bond_end( BOND_LINES.End());

      for( ; itr != itr_end && itr_bond != itr_bond_end; ++itr, ++itr_bond)
      {
        // get the next bond isometry
        const chemistry::BondIsometry isometry = chemistry::BondIsometryEnum( *itr);

        // ensure it is defined
        BCL_Assert( isometry != chemistry::s_NumberOfIsometries, "Undefined double bond isometry: " + *itr);

        // skip over non-double bonds
        while( itr_bond != itr_bond_end && itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() != size_t( 4))
        {
          ++itr_bond;
        }

        // stop if there are no further double bonds
        if( itr_bond == itr_bond_end)
        {
          break;
        }

        itr_bond->SetIsometry( isometry);
      }

      if( itr != itr_end)
      {
        BCL_MessageCrt( "More double bond isometries given than double bonds in molecule! File corruption likely.");
        return;
      }

      // skip any remaining non-double bonds
      while( itr_bond != itr_bond_end && itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() != size_t( 4))
      {
        ++itr_bond;
      }

      if( itr != itr_end)
      {
        BCL_MessageCrt( "Double bond isometries not given for all double bonds! File corruption likely.");
        return;
      }
    }

    //! @brief retrieve charges from atom lines
    //! @param ATOM_LINES mdl atoms to retrieve charges from
    void MdlProperty::RetrieveCharges( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();

      size_t index( 1); // charge indices are 1-offset
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        if( itr->GetAtomType()->GetFormalCharge()) // only write non-zero charges
        {
          m_PropertyStrings.PushBack( util::Format()( index));
          m_PropertyStrings.PushBack( util::Format()( itr->GetAtomType()->GetFormalCharge()));
        }
      }
    }

    //! @brief retrieve atom types to atom lines
    //! @param ATOM_LINES mdl atoms to retrieve atom types from
    void MdlProperty::RetrieveAtomTypes( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();
      m_PropertyStrings.AllocateMemory( ATOM_LINES.GetSize());
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        m_PropertyStrings.PushBack( itr->GetAtomType().GetName());
      }
    }

    //! @brief retrieve bond types to bond lines
    //! @param BOND_LINES mdl bonds to retrieve bond types from
    void MdlProperty::RetrieveBondTypes( const storage::Vector< BondInfo> &BOND_LINES)
    {
      m_PropertyStrings.Reset();
      m_PropertyStrings.AllocateMemory( BOND_LINES.GetSize());
      for
      (
        storage::Vector< BondInfo>::const_iterator itr( BOND_LINES.Begin()), itr_end( BOND_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        m_PropertyStrings.PushBack( itr->GetConstitutionalBondType().GetName());
      }
    }

    //! @brief retrieve chiralities to mdl atoms
    //! @param ATOM_LINES mdl atoms to retrieve chirality property from
    void MdlProperty::RetrieveChiralities( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();

      size_t index( 1); // chirality indices are 1-offset
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        if( itr->GetChirality() != chemistry::e_NonChiral && itr->GetChirality() != chemistry::e_UnknownChirality) // only write out chiral centers
        {
          m_PropertyStrings.PushBack( util::Format()( index));
          m_PropertyStrings.PushBack( chemistry::ChiralityEnum( itr->GetChirality()));
        }
      }
    }

    //! @brief retrieve double bond isometries from MdlBonds
    //! @param BOND_LINES mdl bonds to apply isometry property from
    void MdlProperty::RetrieveDoubleBondIsometries( const storage::Vector< BondInfo> &BOND_LINES)
    {
      m_PropertyStrings.Reset();

      for
      (
        storage::Vector< BondInfo>::const_iterator itr_bond( BOND_LINES.Begin()), itr_bond_end( BOND_LINES.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // skip bonds that could not have isometry
        if( itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() == size_t( 4))
        {
          // push back the next bond isometry
          m_PropertyStrings.PushBack
          (
            chemistry::BondIsometryEnum( itr_bond->GetConfigurationalBondType()->GetIsometry())
          );
        }
      }
    }
  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief MoleculeReadingPref as string
    //! @param PREF the MoleculeReadingPref desired
    //! @return the MoleculeReadingPref as string
    const std::string &GetPrefName( const HydrogenHandlingPref &PREF)
    {
      static const std::string s_names[] =
      {
        "SelectiveSaturate",
        "Saturate",
        "Maintain",
        "Remove",
        GetStaticClassName< HydrogenHandlingPref>()
      };

      return s_names[ PREF];
    }

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    const util::ShPtr< command::FlagInterface> &GetAddHydrogensFlag()
    {
      static util::ShPtr< command::FlagInterface> s_h_add
        (
          new command::FlagStatic( "add_h", "add hydrogens to molecules when loaded")
        );

      return s_h_add;
    }

    //! @brief get access to a global flag defining whether hydrogens should be removed
    //! @return access to a global flag defining whether hydrogens should be removed
    const util::ShPtr< command::FlagInterface> &GetRemoveHydrogensFlag()
    {
      static util::ShPtr< command::FlagInterface> s_h_del
        (
          new command::FlagStatic( "remove_h", "remove hydrogens from molecules when loaded")
        );

      return s_h_del;
    }

    //! @brief get access to a global flag defining whether charges should be neutralized wherever possible
    //! @return access to a global flag defining whether charges should be neutralized wherever possible
    const util::ShPtr< command::FlagInterface> &GetNeutralizeChargesFlag()
    {
      static util::ShPtr< command::FlagInterface> s_neutralize
      (
        new command::FlagStatic
        (
          "neutralize",
          "neutralize charges; by default, if the flag is specified but no neutralization type is given, "
          "BondOrderAndpH will be used, otherwise, no neutralization is used. All neutralization algorithms preserve "
          "aromaticity except BondOrderAndpHAromaticityLossOk",
          command::Parameter
          (
            "method",
            "method used to neutralize charged atoms in the molecule",
            command::ParameterCheckSerializable( NeutralizationPrefEnum()),
            "BondOrderAndpH"
          )
        )
      );

      return s_neutralize;
    }

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    const util::ShPtr< command::FlagInterface> &GetExplicitAromaticityFlag()
    {
      static util::ShPtr< command::FlagInterface> s_explicit_aro
        (
          new command::FlagStatic
          (
            "explicit_aromaticity",
            "write MDL bonds with aromatic bonds specified explicitly (as 4); "
            "alternatively, the default behavior during MDL writing is to "
            "kekulize aromatic rings (write alternating single (1) / double (2) "
            "bonds)"
            )
        );

      return s_explicit_aro;
    }

    //! @brief add molecule reading and writing preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeIOPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetAddHydrogensFlag());
      CMD.AddFlag( GetRemoveHydrogensFlag());
      CMD.AddFlag( GetNeutralizeChargesFlag());
      CMD.AddFlag( GetExplicitAromaticityFlag());
    }

    //! @brief add molecule reading  preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeReadingPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetAddHydrogensFlag());
      CMD.AddFlag( GetRemoveHydrogensFlag());
      CMD.AddFlag( GetNeutralizeChargesFlag());
    }

    //! @brief add molecule writing  preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeWritingPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetExplicitAromaticityFlag());
    }

    //! @brief get access to the desired hydrogen handling based on the add/remove flags
    //! @return access to the desired hydrogen handling based on the add/remove flags
    HydrogenHandlingPref GetCommandLineHydrogensPref()
    {
      const bool add_h( GetAddHydrogensFlag()->GetFlag());
      const bool remove_h( GetRemoveHydrogensFlag()->GetFlag());
      BCL_Assert( !add_h || !remove_h, "Cannot add and remove hydrogens; choose one or the other");
      return add_h ? e_Saturate : remove_h ? e_Remove : e_Maintain;
    }

    //! @brief NeutralizationPref as string
    //! @param PREF the NeutralizationPref desired
    //! @return the NeutralizationPref as string
    const std::string &GetNeutralizationPrefName( const NeutralizationPref &PREF)
    {
      static const std::string s_names[] =
      {
        "None",
        "BondOrder",
        "pH",
        "BondOrderAndpH",
        "BondOrderAndpHAromaticityLossOk",
        GetStaticClassName< NeutralizationPref>()
      };
      return s_names[ PREF];
    }

    //! @brief determine whether the preference allows alteration of the bond order
    //! @return true if the neutralization pref allows alteration of the bond order
    bool GetNeutralizationPrefAllowsBondOrderChange( const NeutralizationPref &PREF)
    {
      static const bool s_bond_order_chng_allowed[] = { false, true, false, true, true, false};
      return s_bond_order_chng_allowed[ size_t( PREF)];
    }

    //! @brief determine whether the preference allows alteration of the pH
    //! @return true if the neutralization pref allows alteration of the pH
    bool GetNeutralizationPrefAllowspHChange( const NeutralizationPref &PREF)
    {
      static const bool s_chng_allowed[] = { false, false, true, true, true, false};
      return s_chng_allowed[ size_t( PREF)];
    }

    //! @brief determine whether the preference allows alteration of the pH
    //! @return true if the neutralization pref allows alteration of the pH
    bool GetNeutralizationPrefAllowsAromaticityChange( const NeutralizationPref &PREF)
    {
      static const bool s_chng_allowed[] = { false, false, false, false, true, false};
      return s_chng_allowed[ size_t( PREF)];
    }

    //! @brief whether to explicitly write aromatic bonds instead of default kekulized bonds
    //! @return true if writing aromatic bonds explicitly
    BCL_API bool GetExplicitAromaticBondsPref()
    {
      return GetExplicitAromaticityFlag()->GetFlag();
    }

    //! @brief get whether to try to neutralize charges on read in molecules
    //! @return whether to try to neutralize charges on read in molecules
    NeutralizationPrefEnum GetCommandLineNeutralizationPref()
    {
      if( !GetNeutralizeChargesFlag()->GetFlag())
      {
        return e_None;
      }
      NeutralizationPrefEnum type( GetNeutralizeChargesFlag()->GetFirstParameter()->GetValue());
      return type;
    }

  } // namespace sdf
} // namespace bcl
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
#include "sdf/bcl_sdf_molfile_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! number molecule descriptor lines
    const size_t MolfileHandler::s_NumberDescriptionLines = 3;

    //! @brief standard constructor
    MolfileHandler::MolfileHandler() :
      m_Description()
    {
    }

    //! @brief constructor from a description and atom/bond info
    MolfileHandler::MolfileHandler
    (
      const std::string &DESCRIPTION,
      const storage::Vector< AtomInfo> &ATOM_INFOS,
      const storage::Vector< BondInfo> &BOND_INFOS
    ) :
      CTabHandler( ATOM_INFOS, BOND_INFOS),
      m_Description( DESCRIPTION)
    {
      m_Description = StandardizeDescription( m_Description);
    }

    //! @brief constructor from input stream
    MolfileHandler::MolfileHandler( std::istream &ISTREAM) :
      m_Description()
    {
      ReadMolfile( ISTREAM);
    }

    //! @brief constructor from a pre-read set of lines
    MolfileHandler::MolfileHandler( const storage::List< std::string> &LINES)
    {
      ReadMolfile( LINES.Begin(), LINES.End());
    }

    //! @brief virtual copy constructor
    MolfileHandler *MolfileHandler::Clone() const
    {
      return new MolfileHandler( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MolfileHandler::ReadMolfile( std::istream &ISTREAM)
    {
      // reset all members
      m_Description.clear();

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        BCL_MessageCrt( "Error reading Molfile: passed bad istream");
        return ISTREAM;
      }

      // read descriptor line till the header
      int description_lines = 0;
      std::string line;
      for( std::getline( ISTREAM, line); !ISTREAM.eof(); std::getline( ISTREAM, line))
      {
        // skip empty lines
        if( !ContainsNonspaceCharacters( line))
        {
          m_Description += '\n';
          ++description_lines;
          continue;
        }

        // if we found the header line, stop reading in the description
        MdlHeader header;
        header.SetFromMdlLine( line, description_lines);
        if( header.IsValid())
        {
          break;
        }

        // insert the descriptor line
        m_Description += line;
        m_Description += '\n';
        ++description_lines;
      } // search for header line and inserting descriptor lines

      if( ISTREAM.eof())
      {
        BCL_MessageStd( "Unexpected end of input, cannot read molfile");
        m_Description.clear();
        return ISTREAM;
      }

      if( description_lines != s_NumberDescriptionLines)
      {
        BCL_MessageStd
        (
          "Warning: description contains a non-standard number of lines (" +
          util::Format()( description_lines) + " provided, standard is " + util::Format()( s_NumberDescriptionLines) +
          "); proceeding anyway."
        );
      }

      // standardize the description
      m_Description = StandardizeDescription( m_Description);

      // read the CTab from the stream
      return ReadCTab( ISTREAM);
    }

    //! @brief read a molfile from a set of iterators
    //! @param LINE_BEGIN where to begin reading the molfile from
    //! @param LINE_END one-past-last iterator of where reading should happen
    //! @return an iterator to the first unread line
    storage::List< std::string>::const_iterator MolfileHandler::ReadMolfile
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END
    )
    {
      Reset();

      if( LINE_BEGIN == LINE_END)
      {
        BCL_MessageStd( "MolfileHandler: Nothing to read.");
        return LINE_BEGIN;
      }

      storage::List< std::string>::const_iterator itr_lines( LINE_BEGIN);

      // read descriptor line till the header
      int description_lines = 0;
      for( ; itr_lines != LINE_END; ++itr_lines)
      {
        // skip empty lines
        if( !ContainsNonspaceCharacters( *itr_lines))
        {
          m_Description += '\n';
          ++description_lines;
          continue;
        }

        // if we found the header line, stop reading in the description
        MdlHeader header;
        header.SetFromMdlLine( *itr_lines, description_lines);
        if( header.IsValid())
        {
          break;
        }

        // insert the descriptor line
        m_Description += *itr_lines;
        m_Description += '\n';
        ++description_lines;
      } // search for header line and inserting descriptor lines

      if( itr_lines == LINE_END)
      {
        BCL_MessageStd( "Unexpected end of molfile, cannot read molfile");
        m_Description.clear();
        return itr_lines;
      }

      if( description_lines != s_NumberDescriptionLines)
      {
        BCL_MessageStd
        (
          "Warning: molfile description contains a non-standard number of lines (" +
          util::Format()( description_lines) + ", standard is " + util::Format()( s_NumberDescriptionLines) +
          "); proceeding anyway."
        );
      }

      // standardize the description
      m_Description = StandardizeDescription( m_Description);

      itr_lines = ReadCTab( itr_lines, LINE_END);
      if( !IsValid())
      {
        BCL_MessageCrt( "Connection table was faulty, not reading the rest of this molecule: " + m_Description);
      }
      return itr_lines;
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param FORCE_WRITE_ATOM_TYPES whether to force writing atom types and other BCL info
    //! @return the output stream that was written to
    std::ostream &MolfileHandler::WriteMolfile( std::ostream &OSTREAM, const bool &FORCE_WRITE_ATOM_TYPES) const
    {
      if( !IsValid())
      {
        return OSTREAM;
      }
      return WriteMolfile( OSTREAM, m_Description, *this, FORCE_WRITE_ATOM_TYPES);
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param DESCRIPTION the three-line description to use
    //! @param CTAB the connection table data for this molfile
    //! @param FORCE_WRITE_ATOM_TYPES whether to force writing atom types and other BCL info
    //! @return the output stream that was written to
    std::ostream &MolfileHandler::WriteMolfile
    (
      std::ostream &OSTREAM,
      const std::string &DESCRIPTION,
      const CTabHandler &CTAB,
      const bool &FORCE_WRITE_ATOM_TYPES
    )
    {
      if( CTAB.IsValid()) // from CTabHandler
      {
        // if something goes wrong, don't write a fragmented file
        std::ostringstream ostream;
        ostream << DESCRIPTION << '\n';
        CTAB.WriteCTab( ostream, FORCE_WRITE_ATOM_TYPES);
        OSTREAM << ostream.str();
      }
      else
      {
        BCL_MessageCrt( "Ignoring request to write molecule with name \"" + DESCRIPTION + "\"; an invalid connection table was given");
      }
      return OSTREAM;
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param DESCRIPTION the three-line description to use
    //! @param ATOM_INFO the atom infos to write
    //! @param BOND_INFO the bond infos to write
    //! @return the ostream that was written to
    std::ostream &MolfileHandler::WriteMolfile
    (
      std::ostream &OSTREAM,
      const std::string &DESCRIPTION,
      const storage::Vector< AtomInfo> &ATOM_INFO,
      const storage::Vector< BondInfo> &BOND_INFO
    )
    {
      CTabHandler ctab( ATOM_INFO, BOND_INFO);
      if( !ctab.IsValid())
      {
        BCL_MessageCrt( "Could not make a valid CTab to write to SDF file");
        return OSTREAM;
      }

      return WriteMolfile( OSTREAM, DESCRIPTION, ctab);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the description; checks # of newlines, ensures that it is exactly s_NumberDescriptionLines
    //! If the # of newlines is < s_NumberDescriptionLines - 1, adds new lines as necessary
    //! Any newlines >= s_NumberDescriptionLines are replaced with spaces
    std::string MolfileHandler::StandardizeDescription( const std::string &DESCRIPTION)
    {
      const size_t last_non_space( DESCRIPTION.find_last_not_of( " \n\t\r"));
      const size_t description_size( last_non_space + 1);
      std::string description;
      description.reserve( description_size);

      // keep track of the number of new lines seen so far
      size_t number_new_lines( 0);

      for( size_t i( 0); i < description_size; ++i)
      {
        if( DESCRIPTION[ i] == '\n' && ++number_new_lines >= s_NumberDescriptionLines)
        {
          description += ' ';
          continue;
        }
        else if( DESCRIPTION[ i] == '\r') // skip carriage returns (occur when reading sdfs from windows on non-windows machine)
        {
          continue;
        }
        description += DESCRIPTION[ i];
      }
      // add new lines until there are s_NumberDescriptionLines - 1 of them
      while( ++number_new_lines < s_NumberDescriptionLines)
      {
        description += '\n';
      }

      return description;
    }

    //! @brief read MolfileHandler object from std::istream
    //! @param ISTREAM istream that contains MolfileHandler object
    //! @return istream after MolfileHandler object was extracted
    std::istream &MolfileHandler::Read( std::istream &ISTREAM)
    {
      return ReadMolfile( ISTREAM);
    }

    //! @brief write MolfileHandler into std::ostream
    //! @param OSTREAM ostream that gets MolfileHandler object
    //! @param INDENT indentation
    //! @return ostream after MolfileHandler object was inserted
    std::ostream &MolfileHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MolfileHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new MolfileHandler())
    );

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_rxn_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "sdf/bcl_sdf_fragment_factory.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief construct a ReactionComplete from a rxn file
    //! @param HANDLER the handler to use
    //! @param H_PREF hydrogen handling preference
    chemistry::ReactionComplete RXNFactory::MakeReactionComplete( const RXNHandler &HANDLER)
    {
      // The MDL handlers for products and reactants
      //const storage::Vector< MdlHandler> &reactant_handlers( HANDLER.GetReactantMdlHandlers());
      //const storage::Vector< MdlHandler> &product_handlers( HANDLER.GetProductMdlHandlers());
      const storage::Vector< MolfileHandler> &reactant_handlers( HANDLER.GetReactantHandlers());
      const storage::Vector< MolfileHandler> &product_handlers( HANDLER.GetProductHandlers());
      storage::Vector< storage::Vector< AtomInfo> > reactant_atom_info;
      storage::Vector< storage::Vector< BondInfo> > reactant_bond_info;

      storage::Vector< storage::Vector< AtomInfo> > product_atom_info;
      storage::Vector< storage::Vector< BondInfo> > product_bond_info;

      //storage::Vector< chemistry::FragmentComplete> reactants;
      //storage::Vector< chemistry::FragmentComplete> products;

      // make the reactants
      storage::Vector< chemistry::ReactionStructure> reactant_structs;
      for( size_t r( 0); r < reactant_handlers.GetSize(); ++r)
      {
        reactant_structs.PushBack( MakeReactionStructure( reactant_handlers( r)));
      }

      // Make the products
      storage::Vector< chemistry::ReactionStructure> product_structs;
      for( size_t p( 0); p < product_handlers.GetSize(); ++p)
      {
        product_structs.PushBack( MakeReactionStructure( product_handlers( p)));
      }

      // Make the reaction complete
      return chemistry::ReactionComplete
      (
        reactant_structs,
        product_structs,
        HANDLER.GetReactiveAtomsInReactants(),
        HANDLER.GetReactiveAtomsInProducts(),
        HANDLER.GetDescription()
      );

    }

    //! @brief construct a ReactionStructure from a molfile
    chemistry::ReactionStructure RXNFactory::MakeReactionStructure( const CTabHandler &CTAB)
    {
      chemistry::ReactionStructure rxn_struct;
      if( !CTAB.IsValid())
      {
        return rxn_struct;
      }

      const storage::Map< std::string, storage::Vector< std::string> > &mdl_props( CTAB.GetCachedProperties());

      chemistry::AtomVector< chemistry::AtomComplete> atom_v( CTAB.GetAtomInfo(), CTAB.GetBondInfo());
      chemistry::FragmentComplete new_frag( atom_v, std::string());
      
      std::vector< chemistry::ReactionStructure::Aromaticity> allowed_aromaticities( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);

      const std::string arom_prop( MdlProperty( MdlProperty::e_BclAtomAromaticity).GetLabel());
      if( mdl_props.Has( arom_prop))
      {
        const storage::Vector< std::string> &arom_str( mdl_props.GetValue( arom_prop));
        std::vector< size_t> arom_atoms( atom_v.GetSize(), size_t( 0));

        // store the aromaticity of the molecule
        for( size_t i( 0); i < arom_str.GetSize(); i += 2)
        {
          size_t atom_no( 0);
          if( !util::TryConvertFromString( atom_no, arom_str( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse aromaticity information from reaction");
            allowed_aromaticities = std::vector< chemistry::ReactionStructure::Aromaticity>( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);
            break;
          }
          if( atom_no == 0 || atom_no > atom_v.GetSize())
          {
            BCL_MessageStd( "Could not parse aromaticity information from reaction");
            if( !atom_no)
            {
              BCL_MessageStd( "Atom number was 0, minimum value is 1 as per the SDF numbering standard");
            }
            allowed_aromaticities = std::vector< chemistry::ReactionStructure::Aromaticity>( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);
            break;
          }
          if( arom_str( i + 1) == "ARO")
          {
            allowed_aromaticities[ atom_no - 1] = chemistry::ReactionStructure::e_Aromatic;
          }
          else if( arom_str( i + 1) == "ALI")
          {
            allowed_aromaticities[ atom_no - 1] = chemistry::ReactionStructure::e_Aliphatic;
          }
        }
      }
      
      rxn_struct = chemistry::ReactionStructure( new_frag, allowed_aromaticities);

      return rxn_struct;
    }

  } // namespace sdf
} // namespace bcl

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
#include "sdf/bcl_sdf_rxn_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"
#include "storage/bcl_storage_vector.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    const size_t RXNHandler::s_NumberDescriptionLines = 3;

    //! @brief standard constructor
    RXNHandler::RXNHandler() :
      m_Description(),
      m_NumberReactants( 0),
      m_NumberProducts( 0),
      m_ReactantMolfiles(),
      m_ProductMolfiles(),
      m_IsValid( false)
    {
    }

    //! @brief constructor with initial input
    RXNHandler::RXNHandler( std::istream &ISTREAM) :
      m_Description(),
      m_NumberReactants( 0),
      m_NumberProducts( 0),
      m_ReactantMolfiles(),
      m_ProductMolfiles(),
      m_IsValid( false)
    {
      ReadFromRXN( ISTREAM);
    }

    //! @brief virtual copy constructor
    RXNHandler *RXNHandler::Clone() const
    {
      return new RXNHandler( *this);
    }

    //! @brief read from std::istream
    //! @note THIS READS TO THE END OF THE STREAM.  RXN files are designed to contain only one
    //!       reaction, which is an assumption made here.  All of the lines from the stream may not
    //!       be processed but the whole stream WILL be processed
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RXNHandler::ReadFromRXN( std::istream &ISTREAM)
    {
      // RXN files are only supposed to contain one reaction
      // assume ISTREAM contains only this and buffer everything appropriately

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        if( !ISTREAM.good())
        {
          BCL_MessageCrt( "Error reading RXN: passed bad stream");
        }
        return ISTREAM;
      }

      storage::List< std::string> rxn_buff;
      std::string buf;
      while( std::getline( ISTREAM, buf))
      {
        rxn_buff.PushBack( buf);
      }

      // read information from the buffer
      ReadFromRXN( rxn_buff.Begin(), rxn_buff.End());
      
      return ISTREAM;
    }

    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM,
      const std::string &DESCRIPTION,
      const storage::Vector< MolfileHandler> &REACTANT_MOLFILES,
      const storage::Vector< MolfileHandler> &PRODUCT_MOLFILES
    )
    {
      size_t n_reactants( REACTANT_MOLFILES.GetSize());
      size_t n_products( REACTANT_MOLFILES.GetSize());
      if( n_reactants < 1000 && n_products < 1000)
      {
        std::stringstream oss;
        
        // Add the $RXN label and description
        oss << GetDefaultLine( e_RXNStartLine) << std::endl;
        oss << MolfileHandler::StandardizeDescription( DESCRIPTION) << std::endl;

        // Add the reaction header, i.e. (n reactants)(n products)
        std::string rxn_header( GetDefaultLine( e_RXNHeaderLine));
        GetMdlEntryTypes().RXNHeader_NumberReactantLines->Set( rxn_header, n_reactants);
        GetMdlEntryTypes().RXNHeader_NumberProductLines->Set( rxn_header, n_products);
        oss << rxn_header << std::endl;

        size_t mol_no( 0);
        for
        (
          storage::Vector< MolfileHandler>::const_iterator itr_mf( REACTANT_MOLFILES.Begin()),
            itr_mf_end( REACTANT_MOLFILES.End());
          itr_mf != itr_mf_end;
          ++itr_mf, ++mol_no
        )
        {
          if( !itr_mf->IsValid())
          {
            BCL_MessageCrt
            ( 
              "Cannot write RXN file: invalid molfile was provided for reactant " + util::Format()( mol_no)
              + " with description \"" + itr_mf->GetDescription() + "\""
            );
            return OSTREAM;
          }

          // write $MOL
          oss << GetDefaultLine( e_RXNMolStartLine) << std::endl;

          // write out the molfile
          itr_mf->WriteMolfile( oss);
        }

        mol_no = 0;
        for
        (
          storage::Vector< MolfileHandler>::const_iterator itr_mf( PRODUCT_MOLFILES.Begin()),
            itr_mf_end( PRODUCT_MOLFILES.End());
          itr_mf != itr_mf_end;
          ++itr_mf, ++mol_no
        )
        {
          if( !itr_mf->IsValid())
          {
            BCL_MessageCrt
            ( 
              "Cannot write RXN file: invalid molfile was provided for product " + util::Format()( mol_no)
              + " with description \"" + itr_mf->GetDescription() + "\""
            );
            return OSTREAM;
          }

          // write $MOL
          oss << GetDefaultLine( e_RXNMolStartLine) << std::endl;

          // write the molfile out
          itr_mf->WriteMolfile( oss);
        }

        OSTREAM << oss.str() << "\n";
      }
      else
      {
        BCL_MessageCrt
        ( 
          "Cannot write RXN file with >999 reactants or products (" + util::Format()( n_reactants)
          + " and " + util::Format()( n_products) + " were given, respectively)"
        );
        BCL_MessageCrt( "Unable to write RXN with description \"" + DESCRIPTION + "\"");
      }
      return OSTREAM;
    }

    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM, 
      const std::string &DESCRIPTION,
      const storage::Vector< chemistry::FragmentComplete> &REACTANTS,
      const storage::Vector< chemistry::FragmentComplete> &PRODUCTS,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP
    )
    {
      storage::Vector< std::string> reactant_descriptions;
      reactant_descriptions.AllocateMemory( REACTANTS.GetSize());

      storage::Vector< std::string> product_descriptions;
      product_descriptions.AllocateMemory( PRODUCTS.GetSize());

      storage::Vector< storage::Vector< AtomInfo> > reactant_atom_info;
      reactant_atom_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< BondInfo> > reactant_bond_info;
      reactant_bond_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< AtomInfo> > product_atom_info;
      product_atom_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< BondInfo> > product_bond_info;
      product_bond_info.AllocateMemory( REACTANTS.GetSize());

      for( size_t rno( 0), end_rno( REACTANTS.GetSize()); rno < end_rno; ++rno)
      {
        reactant_descriptions.PushBack( REACTANTS( rno).GetName());
        reactant_atom_info.PushBack( REACTANTS( rno).GetAtomInfo());
        reactant_bond_info.PushBack( REACTANTS( rno).GetBondInfo());
      }

      for( size_t pno( 0), end_pno( PRODUCTS.GetSize()); pno < end_pno; ++pno)
      {
        product_descriptions.PushBack( PRODUCTS( pno).GetName());
        product_atom_info.PushBack( PRODUCTS( pno).GetAtomInfo());
        product_bond_info.PushBack( PRODUCTS( pno).GetBondInfo());
      }

      return WriteToRXN
             ( 
               OSTREAM, DESCRIPTION, reactant_atom_info, 
               reactant_bond_info, product_atom_info, product_bond_info, 
               REACTANT_ATOM_MAP, PRODUCT_ATOM_MAP,
               reactant_descriptions, product_descriptions
             );
    } 

    //! @brief write a RXN-formatted reaction to an output stream 
    //! @param OSTREAM the output stream
    //! @param DESCRIPTION the description of the reaction
    //! @param REACTANTS the reactants to write
    //! @param PRODUCTS the products to write
    //! @param REACTANT_ATOM_MAP reactive atom map for reactants
    //! @param PRODUCT_ATOM_MAP reactive atom map for products
    //! @param TERMINATION_LINE the line to terminate the RXN with
    //! @return the output stream that was written
    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM, 
      const std::string &DESCRIPTION,
      const storage::Vector< storage::Vector< AtomInfo> > &REACTANT_ATOM_INFO,
      const storage::Vector< storage::Vector< BondInfo> > &REACTANT_BOND_INFO,
      const storage::Vector< storage::Vector< AtomInfo> > &PRODUCT_ATOM_INFO,
      const storage::Vector< storage::Vector< BondInfo> > &PRODUCT_BOND_INFO,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP,
      const storage::Vector< std::string> &REACTANT_DESCRIPTIONS,
      const storage::Vector< std::string> &PRODUCT_DESCRIPTIONS
    ) 
    {
      size_t n_reactants( REACTANT_ATOM_INFO.GetSize());
      size_t n_products( PRODUCT_ATOM_INFO.GetSize());
      
      if( REACTANT_BOND_INFO.GetSize() != n_reactants)
      {
        BCL_MessageStd
        ( 
          "Cannot write reaction, atom infos were specified for " + util::Format()( n_reactants) + " reactant(s), "
          "but bond infos were specified for " + util::Format()( REACTANT_BOND_INFO.GetSize()) + " reactant(s).  "
          "These should be equivalent"
        );
        return OSTREAM;
      }

      if( PRODUCT_BOND_INFO.GetSize() != n_reactants)
      {
        BCL_MessageStd
        ( 
          "Cannot write reaction, atom infos were specified for " + util::Format()( n_reactants) + " product(s), "
          "but bond infos were specified for " + util::Format()( PRODUCT_BOND_INFO.GetSize()) + " product(s).  "
          "These should be equivalent"
        );
        return OSTREAM;
      }
      
      if( n_reactants > 999 || n_products > 999)
      {
        BCL_MessageCrt
        ( 
          "Ignoring request to write a reaction with " + util::Format()( n_reactants) 
          + " reactants and " + util::Format()( n_products) + " products"
        );
        BCL_MessageCrt( "The MDL specification does not allow for reactions with more than 999 molecules");
        return OSTREAM;
      }

      // translate the reactive atom maps to something that maps atom index (key) to mapped value (value)
      // for each reactant/product (vector index)
      storage::Vector< storage::Map< size_t, size_t> > reactive_atoms_reactants( n_reactants);
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( REACTANT_ATOM_MAP.Begin()),
          itr_map_end( REACTANT_ATOM_MAP.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() >= n_reactants)
        {
          BCL_MessageStd
          ( 
            "Cannot output reaction: malformed atom mapping.  Tried to add a map value of " + 
            util::Format()( itr_map->first) + " to reactant " + util::Format()( itr_map->second.First()) +
            " but maximum reactant index is " + util::Format()( n_reactants - 1)
          ); 
          return OSTREAM;
        }
        reactive_atoms_reactants( itr_map->second.First())[ itr_map->second.Second()] = itr_map->first;
      }

      storage::Vector< storage::Map< size_t, size_t> > reactive_atoms_products( n_products);
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( PRODUCT_ATOM_MAP.Begin()),
          itr_map_end( PRODUCT_ATOM_MAP.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() >= n_products)
        {
          BCL_MessageStd
          ( 
            "Cannot output reaction: malformed atom mapping.  Tried to add a map value of " + 
            util::Format()( itr_map->first) + " to product " + util::Format()( itr_map->second.First()) +
            " but maximum product index is " + util::Format()( n_products - 1)
          ); 
          return OSTREAM;
        }
        reactive_atoms_products( itr_map->second.First())[ itr_map->second.Second()] = itr_map->first;
      }

      storage::Vector< MolfileHandler> reactant_molfiles;
      reactant_molfiles.AllocateMemory( n_reactants);
      storage::Vector< MolfileHandler> product_molfiles;
      product_molfiles.AllocateMemory( n_products);

      bool use_reactant_descriptions( REACTANT_DESCRIPTIONS.GetSize() == n_reactants);
      bool use_product_descriptions( PRODUCT_DESCRIPTIONS.GetSize() == n_products);

      // generate reactant molfile handlers
      for( size_t rno( 0); rno < n_reactants; ++rno)
      {
        std::string desc;
        if( use_reactant_descriptions)
        {
          desc = REACTANT_DESCRIPTIONS( rno);
        }
        else
        {
          desc = std::string( "Reactant " + util::Format()( rno) + "\n\n\n");
        }
        desc = MolfileHandler::StandardizeDescription( desc);
        reactant_molfiles.PushBack( MolfileHandler( desc, REACTANT_ATOM_INFO( rno), REACTANT_BOND_INFO( rno)));

        MolfileHandler &r_molfile( reactant_molfiles.LastElement());

        // set the atom mapping
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_map( reactive_atoms_reactants( rno).Begin()), 
            itr_map_end( reactive_atoms_reactants( rno).End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          r_molfile.SetAtomMapping( itr_map->first, itr_map->second);
        }

        if( !r_molfile.IsValid())
        {
          BCL_MessageStd( "Cannot write reaction, reactant " + util::Format()( rno) + " could not be formatted properly");
          return OSTREAM;
        }
      }

      // generate product molfile handlers
      for( size_t pno( 0); pno < n_products; ++pno)
      {
        std::string desc;
        if( use_product_descriptions)
        {
          desc = PRODUCT_DESCRIPTIONS( pno);
        }
        else
        {
          desc = std::string( "Product " + util::Format()( pno) + "\n\n\n");
        }
        desc = StandardizeDescription( desc);
        product_molfiles.PushBack( MolfileHandler( desc, PRODUCT_ATOM_INFO( pno), PRODUCT_BOND_INFO( pno)));

        MolfileHandler &p_molfile( product_molfiles.LastElement());

        // set the atom mapping
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_map( reactive_atoms_products( pno).Begin()), 
            itr_map_end( reactive_atoms_products( pno).End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          p_molfile.SetAtomMapping( itr_map->first, itr_map->second);
        }

        if( !p_molfile.IsValid())
        {
          BCL_MessageStd( "Cannot write reaction, product " + util::Format()( pno) + " could not be formatted properly");
          return OSTREAM;
        }
      }

      return WriteToRXN( OSTREAM, DESCRIPTION, reactant_molfiles, product_molfiles);
    }

    //! @brief gets the class name
    //! @return the class name
    const std::string &RXNHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get all description lines stored in handler
    //! @return list of mdl description lines stored in handler
    const std::string &RXNHandler::GetDescription() const
    {
      return m_Description;
    }

    //! @brief check if the handler is valid
    //! @return true if the mdl handler is in a valid state
    bool RXNHandler::IsValid() const
    {
      return m_IsValid;
    }

    //! @brief Get the number of Reactants
    const size_t &RXNHandler::GetNumberReactants() const
    {
      return m_NumberReactants;
    }

    //! @brief Get the number of Reactants
    const size_t &RXNHandler::GetNumberProducts() const
    {
      return m_NumberProducts;
    }

    //! @brief Get reacting atom mapping on the reactants
    const storage::Map< size_t, storage::Pair< size_t, size_t> > &RXNHandler::GetReactiveAtomsInReactants() const
    {
      return m_ReactiveAtomsReactants;
    }

    //! @brief Get reacting atom mapping on the products
    const storage::Map< size_t, storage::Pair< size_t, size_t> > &RXNHandler::GetReactiveAtomsInProducts() const
    {
      return m_ReactiveAtomsProducts;
    }

    //! @brief read RXN from a file buffer as iterators to strings
    //! @param LINE_BEGIN iterator to first line to begin reading (should be "$RXN" but may be blank/whitespace only)
    //! @param LINE_END iterator to one past the last line that should be read
    //! @return iterator to the first line that was not read.  If all were read, should be LINE_END
    storage::List< std::string>::const_iterator RXNHandler::ReadFromRXN
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END
    )
    {
      // reset all members
      m_Description.erase();
      m_NumberReactants = m_NumberProducts = 0;
      m_ReactiveAtomsReactants.Reset();
      m_ReactiveAtomsProducts.Reset();
      m_IsValid = false;

      if( LINE_BEGIN == LINE_END)
      {
        return LINE_BEGIN;
      }

      storage::List< std::string>::const_iterator itr_line( LINE_BEGIN);

      // Skip lines until the header line is reached, then skip the header line
      size_t n_skipped_lines( 0);
      for( ; itr_line != LINE_END && !IsRXNDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
      ++itr_line; // move past $RXN line

      if( n_skipped_lines)
      {
        BCL_MessageCrt
        ( 
          "Warning: there were " + util::Format()( n_skipped_lines) + 
          " blank lines before $RXN.  This is non-standard (should be 0); continuing anyway"
        );
      }

      // Read 3 lines in as the description
      for( size_t d( 0); d < s_NumberDescriptionLines && itr_line != LINE_END; ++d, ++itr_line)
      {
        m_Description += *itr_line + "\n";
      }
      m_Description = MolfileHandler::StandardizeDescription( m_Description);

      // Next two lines should be the header line and the first $MOL line
      if
      ( 
        !GetMdlEntryTypes().RXNHeader_NumberReactantLines->IsUnsignedInt( *itr_line) ||
        !GetMdlEntryTypes().RXNHeader_NumberProductLines->IsUnsignedInt( *itr_line)
      )
      {
        BCL_MessageCrt
        ( 
          "Line \"" + *itr_line + "\" should have been an RXN header line, but it was not parseable as one.  "
          "Not reading remaining information for RXN with description \"" + m_Description + "\""
        );
        return ++itr_line;
      }

      // read in reactant and product information
      if( itr_line != LINE_END)
      {
        m_NumberReactants = GetMdlEntryTypes().RXNHeader_NumberReactantLines->GetUnsignedInt( *itr_line);
        m_NumberProducts = GetMdlEntryTypes().RXNHeader_NumberProductLines->GetUnsignedInt( *itr_line);
      }
      else
      {
        BCL_MessageCrt( "Unexpected end of RXN file: no header was found for RXN with description \"" + m_Description + "\"");
        return itr_line;
      }
      ++itr_line;

      // if there are no reactants or products then we have technically read in a legitimate RXN file
      if( !m_NumberReactants && !m_NumberProducts)
      {
        BCL_MessageStd( "Note: RXN with description " + m_Description + " contains no molecules");
        m_IsValid = true;
        return itr_line;
      }

      while( itr_line != LINE_END && m_ReactantMolfiles.GetSize() < m_NumberReactants)
      {
        // find the $MOL
        // TODO more error checking here
        n_skipped_lines = 0;
        for( ; itr_line != LINE_END && !IsRXNMolDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
        ++itr_line; // skip past the $MOL
        
        // if we're at the end then bail out
        if( itr_line == LINE_END)
        {
          BCL_MessageCrt( "Unexpected end of RXN: reactants were not all read");
          return itr_line;
        }
        
        // post a warning about skipped lines, files really shouldn't have any blanks but we'll tolerate it
        if( n_skipped_lines)
        {
          BCL_MessageStd( "Warning: reactant block of RXN: there were " + util::Format()( n_skipped_lines) + " blank lines before $MOL, this is non-standard");
        }

        // read the reactant molfile
        m_ReactantMolfiles.PushBack( MolfileHandler());
        itr_line = m_ReactantMolfiles.LastElement().ReadMolfile( itr_line, LINE_END);

        // if there was a problem, stop reading
        if( !m_ReactantMolfiles.LastElement().IsValid())
        {
          BCL_MessageCrt( "Error reading RXN file: could not read reactant #" + util::Format()( m_ReactantMolfiles.GetSize()) + "; exitting");
          return itr_line;
        }
      }

      // make sure we read all reactants
      if( m_ReactantMolfiles.GetSize() != m_NumberReactants)
      {
        BCL_MessageCrt
        ( 
          "Unexpected end of RXN file: could not read the required number of reactants.  "
          "Read " + util::Format()( m_ReactantMolfiles.GetSize()) + " but file specified " + util::Format()( m_NumberReactants)
        );
        return itr_line;
      }

      // now do the same thing for products
      while( itr_line != LINE_END && m_ProductMolfiles.GetSize() < m_NumberProducts)
      {
        // find the $MOL
        // TODO more error checking here
        n_skipped_lines = 0;
        for( ; itr_line != LINE_END && !IsRXNMolDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
        ++itr_line; // skip past the $MOL
        
        // if we're at the end then bail out
        if( itr_line == LINE_END)
        {
          BCL_MessageCrt( "Unexpected end of RXN: products were not all read");
          return itr_line;
        }
        
        // post a warning about skipped lines, files really shouldn't have any blanks but we'll tolerate it
        if( n_skipped_lines)
        {
          BCL_MessageStd( "Warning: product block in RXN: there were " + util::Format()( n_skipped_lines) + " blank lines before $MOL, this is non-standard");
        }

        // read the reactant molfile
        m_ProductMolfiles.PushBack( MolfileHandler());
        itr_line = m_ProductMolfiles.LastElement().ReadMolfile( itr_line, LINE_END);

        // if there was a problem, stop reading
        if( !m_ProductMolfiles.LastElement().IsValid())
        {
          BCL_MessageCrt( "Error reading RXN file: could not read product #" + util::Format()( m_ProductMolfiles.GetSize()) + "; exitting");
          return itr_line;
        }
      }

      // make sure we read all reactants
      if( m_ProductMolfiles.GetSize() != m_NumberProducts)
      {
        BCL_MessageCrt
        ( 
          "Unexpected end of RXN file: could not read the required number of products.  "
          "Read " + util::Format()( m_ReactantMolfiles.GetSize()) + " but file specified " + util::Format()( m_NumberProducts)
        );
        return itr_line;
      }

      // everything is good, itr_line should now point to either LINE_END or whatever one line after the last 'M  END' line is
      // finish off the parsing and set up necessary internal structures
      m_IsValid = ValidateInput(); 

      return itr_line;
    }

    //! @brief finalize parsing; translates m_Stream into internal members, sets m_Parsed to true, cleans m_Stream
    bool RXNHandler::ValidateInput() const
    {
      bool did_succeed( false);

      // set up mappings for reactants
      for( size_t r_no( 0); r_no < m_NumberReactants; ++r_no)
      {
        storage::Vector< size_t> r_mapping( m_ReactantMolfiles( r_no).GetAtomMapping());
        if( r_mapping.GetSize() != m_ReactantMolfiles( r_no).GetAtomInfo().GetSize())
        {
          BCL_MessageCrt( "Error parsing RXN file: reactant " + util::Format()( r_no) + " atoms did not have a proper mapping");
          return did_succeed;
        }
        
        // if atom mappings are non-zero then add them to the reactant atom map
        for( size_t a_no( 0), end_a( r_mapping.GetSize()); a_no < end_a; ++a_no)
        {
          const size_t &map_val( r_mapping( a_no));

          // if the mapping value is 0 then skip it
          if( !map_val)
          {
            continue;
          }

          if( m_ReactiveAtomsReactants.Has( map_val))
          {
            const storage::Pair< size_t, size_t> &old_val( m_ReactiveAtomsReactants.GetValue( map_val));
            BCL_MessageCrt( "Error parsing RXN file: multiple atoms were mapped to value " + util::Format()( map_val));
            BCL_MessageCrt
            ( 
              "  Reactant " + util::Format()( r_no) + " atom index " + util::Format()( a_no) + " and reactant "
              + util::Format()( old_val.First()) + " atom index " + util::Format()( old_val.Second()) + " both "
              "contained this value"
            );
            return did_succeed;
          }

          m_ReactiveAtomsReactants[ map_val].First() = r_no;
          m_ReactiveAtomsReactants[ map_val].Second() = a_no;
        }
      }

      // set up mappings for products
      for( size_t p_no( 0); p_no < m_NumberProducts; ++p_no)
      {
        storage::Vector< size_t> p_mapping( m_ProductMolfiles( p_no).GetAtomMapping());
        if( p_mapping.GetSize() != m_ProductMolfiles( p_no).GetAtomInfo().GetSize())
        {
          BCL_MessageCrt( "Error parsing RXN file: product " + util::Format()( p_no) + " atoms did not have a proper mapping");
          return did_succeed;
        }
        
        // if atom mappings are non-zero then add them to the product atom map
        for( size_t a_no( 0), end_a( p_mapping.GetSize()); a_no < end_a; ++a_no)
        {
          const size_t &map_val( p_mapping( a_no));

          // if the mapping value is 0 then skip it
          if( !map_val)
          {
            continue;
          }

          if( m_ReactiveAtomsProducts.Has( map_val))
          {
            const storage::Pair< size_t, size_t> &old_val( m_ReactiveAtomsProducts.GetValue( map_val));
            BCL_MessageCrt( "Error parsing RXN file: multiple atoms were mapped to value " + util::Format()( map_val));
            BCL_MessageCrt
            ( 
              "  Product " + util::Format()( p_no) + " atom index " + util::Format()( a_no) + " and product "
              + util::Format()( old_val.First()) + " atom index " + util::Format()( old_val.Second()) + " both "
              "contained this value"
            );
            return did_succeed;
          }

          m_ReactiveAtomsProducts[ map_val].First() = p_no;
          m_ReactiveAtomsProducts[ map_val].Second() = a_no;
        }
      }

      did_succeed = true;
      return did_succeed;
    }

    //! @brief read RXNHandler object from std::istream
    //! @param ISTREAM istream that contains RXNHandler object
    //! @return istream after RXNHandler object was extracted
    std::istream &RXNHandler::Read( std::istream &ISTREAM)
    {
      return ReadFromRXN( ISTREAM);
    }

    //! @brief write RXNHandler into std::ostream
    //! @param OSTREAM ostream that gets RXNHandler object
    //! @param INDENT indentation
    //! @return ostream after RXNHandler object was inserted
    std::ostream &RXNHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the description; checks # of newlines, ensures that it is exactly s_NumberDescriptionLines
    //! If the # of newlines is < s_NumberDescriptionLines - 1, adds new lines as necessary
    //! Any newlines >= s_NumberDescriptionLines are replaced with spaces
    std::string RXNHandler::StandardizeDescription( const std::string &DESCRIPTION)
    {
      const size_t last_non_space( DESCRIPTION.find_last_not_of( " \n\t\r"));
      const size_t description_size( last_non_space + 1);
      std::string description;
      description.reserve( description_size);

      // keep track of the number of new lines seen so far
      size_t number_new_lines( 0);

      for( size_t i( 0); i < description_size; ++i)
      {
        if( DESCRIPTION[ i] == '\n' && ++number_new_lines >= s_NumberDescriptionLines)
        {
          description += ' ';
          continue;
        }
        // skip carriage returns (occur when reading sdfs from windows on non-windows machine)
        else if( DESCRIPTION[ i] == '\r')
        {
          continue;
        }
        description += DESCRIPTION[ i];
      }
      // add new lines until there are s_NumberDescriptionLines - 1 of them
      while( ++number_new_lines < s_NumberDescriptionLines)
      {
        description += '\n';
      }
      return description;
    }

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool RXNHandler::ContainsNonspaceCharacters( const std::string &STRING)
    {
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( !isspace( *itr))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief write mdl lines into std::ostream
    //! @param OSTREAM ostream that gets RXNHandler object
    void RXNHandler::WriteMiscProperties
    (
      std::ostream &OSTREAM,
      const storage::Map< std::string, std::string> &MISC_PROPERTIES
    )
    {
      static const std::string pre_property_name_str( "> <");  // goes before each property's name
      static const std::string post_property_name_str( ">\n"); // goes after each property's name

      // iterate over all MDL misc property lines
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          itr_map( MISC_PROPERTIES.Begin()),
          itr_map_end( MISC_PROPERTIES.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->first.size() != 0 && itr_map->second.size() != 0) // the property has a name and value
        {
          OSTREAM << pre_property_name_str        // property name start delimiter
                  << itr_map->first               // property name
                  << post_property_name_str;      // property name deliminater

          const std::string &value( itr_map->second);

          // write the value of the given property line
          OSTREAM << value;

          // if the last character of the string was not a new-line, add one here
          if( value[ value.size() - 1] != '\n')
          {
            OSTREAM << '\n'; // followed by a newline
          }
          OSTREAM << '\n';     // put a blank line follows the last line of the property value
        }
      }
    }

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains $RXN
    bool RXNHandler::IsRXNDelimiter( const std::string &LINE)
    {
      static const std::string s_default_terminal_line( GetDefaultLine( e_RXNStartLine));
      return util::StartsWith( LINE, s_default_terminal_line);
    }

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains $MOL
    bool RXNHandler::IsRXNMolDelimiter( const std::string &LINE)
    {
      static const std::string s_default_terminal_line( GetDefaultLine( e_RXNMolStartLine));
      return util::StartsWith( LINE, s_default_terminal_line);
    }

    //! @brief return the datalable, if line is a datalabel line
    //! @param LINE line from mdl section
    //! @return string that constains datalable, string will be empty for non-data lable lines
    std::string RXNHandler::GetMDLDataLabel( const std::string &LINE)
    {
      // data label delimiter left
      static const char s_misc_property_delimiter_left( '<');

      // data label delimiter right
      static const char s_misc_property_delimiter_right( '>');

      // handle empty lines and lines that do not start with <
      if( LINE.empty() || LINE[ 0] != s_misc_property_delimiter_right)
      {
        return std::string();
      }

      // find label start and end
      const std::string::size_type label_start( LINE.find( s_misc_property_delimiter_left, 1));
      const std::string::size_type label_end( LINE.rfind( s_misc_property_delimiter_right));

      // return an empty string for non data label line
      if( label_start == std::string::npos || label_end == std::string::npos)
      {
        return std::string();
      }

      // misc property name
      return LINE.substr( label_start + 1, label_end - label_start - 1);
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RXNHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new RXNHandler())
    );

  } // namespace sdf
} // namespace bcl

