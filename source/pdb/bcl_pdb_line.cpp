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
#include "pdb/bcl_pdb_line.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Line::s_Instance
    (
      GetObjectInstances().AddInstance( new Line())
    );

    //! @brief record format
    const util::Format Line::s_RecordFormat
    (
      util::Format().W( 6).Fill( ' ').L()
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Line::Line() :
      m_LineType( GetLineTypes().e_Undefined),
      m_String( 128, ' ')
    {
    }

    //! construct Line from LineType and it will write the description
    Line::Line( const LineType &LINE_TYPE) :
      m_LineType( LINE_TYPE),
      m_String( s_LineLength, ' ')
    {
      m_String.replace( s_RecordStart, s_RecordLength, s_RecordFormat( LINE_TYPE.GetName()), s_RecordStart, s_RecordLength);
    }

    //! construct Line from std::string (complete line from pdb) and determines linetype
    Line::Line( const std::string &STRING) :
      m_LineType( GetLineTypes().LineTypeFromPDBLine( STRING)),
      m_String( STRING)
    {
      //ensure that line is long enough
      m_String.resize( s_LineLength, ' ');
    }

    //! copy constructor
    Line *Line::Clone() const
    {
      return new Line( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Line::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns coordinates for ATOM or HETATM line
    //! @return Vector3D with x, Y and Z coordinate
    linal::Vector3D Line::RetrieveCoordinates() const
    {
      //function only works for coordinates in ATOM or HETATM line
      if( m_LineType == GetLineTypes().ATOM)
      {
        //return
        return linal::Vector3D
               (
                 GetNumericalValue< double>( GetEntryTypes().ATOMX),
                 GetNumericalValue< double>( GetEntryTypes().ATOMY),
                 GetNumericalValue< double>( GetEntryTypes().ATOMZ)
               );
      }
      else if( m_LineType == GetLineTypes().HETATM)
      {
        //return
        return linal::Vector3D
               (
                 GetNumericalValue< double>( GetEntryTypes().HETATMX),
                 GetNumericalValue< double>( GetEntryTypes().HETATMY),
                 GetNumericalValue< double>( GetEntryTypes().HETATMZ)
               );
      }
      BCL_MessageCrt( "Function Position called for a non-ATOM line!");
      return linal::Vector3D( util::GetUndefined< double>());
    }

    //! @brief sets x, y and z coordinates of an ATOM or HETATM line
    //! @param COORDINATES the coordinates of that atom
    void Line::PutCoordinates( const linal::Vector3D &COORDINATES)
    {
      if( m_LineType == GetLineTypes().ATOM)
      {
        Put( GetEntryTypes().ATOMX, COORDINATES.X());
        Put( GetEntryTypes().ATOMY, COORDINATES.Y());
        Put( GetEntryTypes().ATOMZ, COORDINATES.Z());
      }
      else if( m_LineType == GetLineTypes().HETATM)
      {
        Put( GetEntryTypes().HETATMX, COORDINATES.X());
        Put( GetEntryTypes().HETATMY, COORDINATES.Y());
        Put( GetEntryTypes().HETATMZ, COORDINATES.Z());
      }
      else
      {
        BCL_MessageCrt( "cannot write Coordinates in non ATOM line");
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief equal operator
    Line &Line::operator =( const Line &LINE)
    {
      // copy member
      m_LineType = LINE.m_LineType;
      m_String   = LINE.m_String;

      // end
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns bool whether Line matches certain criterium
    //! @param CRITERIUM pair of entry type and the string that should be present
    //! @return true if the line type matches the entry type and if the entry is equal to the string
    bool Line::MatchesCriteria( const storage::Pair< EntryType, std::string> &CRITERIUM) const
    {
      //return if LineType matches Entries and string comparison between criteriastring and trimmed string within ENTRIES
      return m_LineType                    == CRITERIUM.First()->GetLineType() &&
             GetString( CRITERIUM.First()) == CRITERIUM.Second();
    }

    //! get certain entry as char
    char Line::GetChar( const EntryType &ENTRY) const
    {
      BCL_Assert
      (
        m_LineType == ENTRY->GetLineType(),
        "EntryType: " + ENTRY.GetName() + " does not match LineType " + m_LineType.GetName() +
        " for pdb line:\n" + m_String
      );
      BCL_Assert( ENTRY->GetDataType() == util::CPPDataTypes::e_Char, "cannot get non-char entry");

      return m_String[ ENTRY->GetStart()];
    }

    //! get certain entry as string
    std::string Line::GetString( const EntryType &ENTRY) const
    {
      BCL_Assert
      (
        m_LineType == ENTRY->GetLineType(),
        "EntryType: " + ENTRY.GetName() + " does not match LineType " + m_LineType.GetName() +
        " for pdb line:\n" + m_String
      );

      return m_String.substr( ENTRY->GetStart(), ENTRY->GetLength());
    }

    //! @brief clear the line except for the line type and line type record
    void Line::Clear()
    {
      // replace all positions behind record locator with ' '
      std::fill( m_String.begin() + LineTypes::s_TypeRecordStart + LineTypes::s_TypeRecordLength, m_String.end(), ' ');
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Line into std::ostream
    std::ostream &Line::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LineType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_String, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Line from std::istream
    std::istream &Line::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_LineType, ISTREAM);
      io::Serialize::Read( m_String, ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief copy a HETATM to ATOM line
    //! @param HETATM_LINE HETATM line
    //! @return ATOM line with entries from given HETATM line
    Line Line::CopyHetatmToAtomLine( const Line &HETATM_LINE)
    {
      // check that line is HETATM line
      if( HETATM_LINE.GetType() != GetLineTypes().HETATM)
      {
        BCL_MessageCrt( "passed non HETATM line but: " + HETATM_LINE.GetString());

        return Line();
      }

      // copy line, set type and change record string
      Line atom_line( HETATM_LINE);
      atom_line.m_LineType = GetLineTypes().ATOM;
      atom_line.m_String.replace( s_RecordStart, s_RecordLength, s_RecordFormat( GetLineTypes().ATOM.GetName()), s_RecordStart, s_RecordLength);

      // end
      return atom_line;
    }

  } // namespace pdb
} // namespace bcl
