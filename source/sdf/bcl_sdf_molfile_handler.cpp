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

