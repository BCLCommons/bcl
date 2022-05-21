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

