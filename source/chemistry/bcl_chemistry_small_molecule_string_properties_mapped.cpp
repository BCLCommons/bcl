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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_mapped.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // single instance of this class
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesMapped::s_Instance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesMapped()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! virtual copy constructor
    SmallMoleculeStringPropertiesMapped *SmallMoleculeStringPropertiesMapped::Clone() const
    {
      return new SmallMoleculeStringPropertiesMapped( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SmallMoleculeStringPropertiesMapped::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SmallMoleculeStringPropertiesMapped::GetAlias() const
    {
      static const std::string s_name( "MappedString");
      return s_name;
    }

    //! @brief operator the implements the assignment operation on the two arguments returning a result
    //! @param MOLECULE the molecule to calculate the property for
    //! @return the property as a string
    std::string SmallMoleculeStringPropertiesMapped::operator()( const ConformationInterface &MOLECULE) const
    {
      // compute the key property
      const std::string key( m_KeyProperty->operator()( MOLECULE));

      // check for the key in the map
      storage::Map< std::string, std::string>::const_iterator itr( m_Map.Find( key));

      // handle the null case; property is not in the map
      if( itr == m_Map.End())
      {
        return std::string();
      }

      // return the string
      return itr->second;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SmallMoleculeStringPropertiesMapped::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
         "uses one property as a key to return the corresponding value from a 2-column table in a file"
         "\n    Input files have rows with a key, followed by spaces or tabs, followed by a value string"
         "\n    2   5.0 6.0 12.0"
         "\n    1   3.0 14.0 18.0"
      );

      parameters.AddInitializer
      (
        "file",
        "file with every line containing a key followed by a value",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );

      parameters.AddInitializer
      (
        "property",
        "string property to calculate the key for a given molecule",
        io::Serialization::GetAgent( &m_KeyProperty)
      );

      parameters.AddInitializer
      (
        "delimiters",
        "delimiters between key and value. Any of these may be used to delimit keys and values, so set this if your "
        "file contains any of the default delimiters in a non-delimiter context",
        io::Serialization::GetAgent( &m_Delimiters),
        "\" \t,;=\""
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool SmallMoleculeStringPropertiesMapped::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      // reset the map
      m_Map.Reset();

      // no errors, so continue

      // open the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // create a key and value string
      std::string key, value;

      // check that the file was not empty
      if( input.eof())
      {
        ERROR_STREAM << "Mapped file was empty!\n";
        return false;
      }

      // read in all the line of "ISTREAM"
      while( !input.eof())
      {
        // create std::string "current_line" to hold the current line of the accessibility file
        std::string current_line;

        // get the next line in "ISTREAM" and put it into "current_line"
        std::getline( input, current_line);

        // skip empty lines; particularly important if the file was generated on Windows and read on Linux/Mac,
        // in which case every other call to getline returns an empty line
        if( util::TrimString( current_line).empty())
        {
          continue;
        }

        // read the key and value
        if( !ReadKeyAndValue( current_line, key, value))
        {
          ERROR_STREAM << "Mapped file was must not contain any rows that do not have a tab-delimited key-value pair, but found line:\n";
          ERROR_STREAM << current_line;
          return false;
        }
        m_Map[ key] = value;
      }

      return true;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief given a line, read the next key and value into passed-in strings
    //! @param LINE the next line
    //! @param KEY storage for the key
    //! @param VALUE storage for the value
    //! @return true iff a value and key were each read
    bool SmallMoleculeStringPropertiesMapped::ReadKeyAndValue
    (
      const std::string &LINE,
      std::string &KEY,
      std::string &VALUE
    ) const
    {
      // find the first tab; this separates the key from the value
      const size_t delimiter_pos( LINE.find_first_of( m_Delimiters));

      KEY.erase();
      VALUE.erase();

      // return false if the delimiter was absent or no value was given
      if( delimiter_pos == std::string::npos)
      {
        return false;
      }

      // get the key
      KEY = util::TrimString( LINE.substr( 0, delimiter_pos));

      // find the first non-delimiter
      const size_t first_non_delim( LINE.find_first_not_of( m_Delimiters, delimiter_pos));

      // return false if no value was given
      if( first_non_delim == std::string::npos)
      {
        return false;
      }

      // get the value, which is everything after the delimiter on the line
      VALUE = util::TrimString( LINE.substr( first_non_delim));
      return true;
    }

  } // namespace chemistry
} // namespace bcl
