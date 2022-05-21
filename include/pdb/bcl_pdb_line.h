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

#ifndef BCL_PDB_LINE_H_
#define BCL_PDB_LINE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_entry_types.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Line
    //! @brief A helper class, that wraps the information contained in each pdb line
    //! @details For each line in pdb, it stores the contents of the line as a string and an associated line type which
    //! is dependent on the first keyword on the pdb line
    //!
    //! @see @link example_pdb_line.cpp @endlink
    //! @author staritrd, meilerj, woetzen
    //! @date 11.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Line :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      LineType    m_LineType; //!< the LineType associated with each pdb line
      std::string m_String;   //!< the content

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief length of pdb line
      static const size_t s_LineLength = 80;

      //! @brief record start
      static const size_t s_RecordStart = 0;

      //! @brief record length
      static const size_t s_RecordLength = 6;

      //! @brief record format
      static const util::Format s_RecordFormat;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Line();

      //! construct Line from LineType and it will write the description
      Line( const LineType &LINE_TYPE);

      //! construct Line from std::string (complete line from pdb) and determines linetype
      explicit Line( const std::string &STRING);

      //! copy constructor
      Line *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! data access to LineType
      const LineType &GetType() const
      {
        return m_LineType;
      }

      //! data access to the line
      const std::string &GetString() const
      {
        return m_String;
      }

      //! @brief returns coordinates for ATOM or HETATM line
      //! @return Vector3D with x, Y and Z coordinate
      linal::Vector3D RetrieveCoordinates() const;

      //! @brief sets x, y and z coordinates of an ATOM or HETATM line
      //! @param COORDINATES the coordinates of that atom
      void PutCoordinates( const linal::Vector3D &COORDINATES);

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      Line &operator =( const Line &LINE);

    ////////////////
    // operations //
    ////////////////

      //! returns bool whether Line matches certain criterium
      //! @param CRITERIUM pair of entry type and the string that should be present
      //! @return true if the line type matches the entry type and if the entry is equal to the string
      bool MatchesCriteria( const storage::Pair< EntryType, std::string> &CRITERIUM) const;

      //! write certain entry to pdb line
      template< typename t_DataType>
      void Put( const EntryType &ENTRY, const t_DataType &DATA)
      {
        BCL_Assert( m_LineType == ENTRY->GetLineType(), "Entry does not match LineType");
        m_String.replace( ENTRY->GetStart(), ENTRY->GetLength(), ENTRY->GetFormat()( DATA), 0, ENTRY->GetLength());
      }

      //! get certain entry as char
      char GetChar( const EntryType &ENTRY) const;

      //! @brief get certain entry as string
      //! @param ENTRY the entry of interest
      //! @return string wich is the substring of the line for that entry
      std::string GetString( const EntryType &ENTRY) const;

      //! get certain entry as numerical value of template type
      //! @param ENTRY the ntry requested
      //! @return a numerical value converted from the string entry
      template< typename t_DataType>
      t_DataType GetNumericalValue( const EntryType &ENTRY) const
      {
        BCL_Assert( ENTRY->IsNumeric(), "cannot get numerical data from non-numerical entry");

        //convert the entry to a numerical value
        return util::ConvertStringToNumericalValue< t_DataType>( GetString( ENTRY));
      }

      //! @brief clear the line except for the line type and line type record
      void Clear();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write Line into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read Line from std::istream
      std::istream &Read( std::istream &ISTREAM);

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief copy a HETATM to ATOM line
      //! @param HETATM_LINE HETATM line
      //! @return ATOM line with entries from given HETATM line
      static Line CopyHetatmToAtomLine( const Line &HETATM_LINE);

    }; // class PdbLine

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_LINE_H_
