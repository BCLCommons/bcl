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
