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
#include "opencl/bcl_opencl_extension_data.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////
  // data //
  //////////

    //! @brief prefix to all extensions "cl"
    //! @return prefix string for opencl extensions
    const std::string &ExtensionData::GetPrefix()
    {
      static const std::string s_prefix( "cl");
      return s_prefix;
    }

    //! @brief separator for prefix, vendor and name
    //! @return separator '_'
    char ExtensionData::GetSeparator()
    {
      return '_';
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ExtensionData::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExtensionData::ExtensionData()
    {
    }

    //! @brief construct from extension string
    //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
    ExtensionData::ExtensionData( const std::string &EXTENSION_STRING)
    {
      BCL_Assert( SetMembersFromExtensionString( EXTENSION_STRING), "incorrect extension string: " + EXTENSION_STRING);
    }

    //! @brief Clone function
    //! @return pointer to new Extension
    ExtensionData *ExtensionData::Clone() const
    {
      return new ExtensionData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ExtensionData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the vendor, eg. "amd"
    //! @return name of vendor
    const std::string &ExtensionData::GetVendor() const
    {
      return m_Vendor;
    }

    //! @brief return the name of extension, eg. "nv"
    //! @return name of extension
    const std::string &ExtensionData::GetName() const
    {
      return m_Name;
    }

    //! @brief get the extension string made up of prefix, vendor and name
    //! @return full extension string, e.g. "cl_amd_fp64"
    std::string ExtensionData::GetString() const
    {
      return GetPrefix() + GetSeparator() + m_Vendor + GetSeparator() + m_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExtensionData::Read( std::istream &ISTREAM)
    {
      // read extension string
      std::string extension_string;
      io::Serialize::Read( extension_string, ISTREAM);
      BCL_Assert( SetMembersFromExtensionString( extension_string), "incorrect extension string: " + extension_string);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ExtensionData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write extension string
      io::Serialize::Write( GetString(), OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set vendor name and extension name form extension string
    //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
    //! @return true, if successful
    bool ExtensionData::SetMembersFromExtensionString( const std::string &EXTENSION_STRING)
    {
      // check prefix
      static const std::string prefix( GetPrefix() + GetSeparator());
      if( EXTENSION_STRING.compare( 0, prefix.length(), prefix) != 0)
      {
        return false;
      }

      // position of second separator
      const std::string::size_type pos_2nd_separator( EXTENSION_STRING.find( GetSeparator(), prefix.length()));
      if( pos_2nd_separator == std::string::npos)
      {
        return false;
      }

      // extract vendor name between prefix and second separator
      const std::string::const_iterator itr_beg( EXTENSION_STRING.begin());
      m_Vendor = std::string( itr_beg + prefix.length(), itr_beg + pos_2nd_separator);

      // extract name
      m_Name = std::string( itr_beg + pos_2nd_separator + 1, EXTENSION_STRING.end());

      // end
      return true;
    }

  } // namespace opencl
} // namespace bcl
