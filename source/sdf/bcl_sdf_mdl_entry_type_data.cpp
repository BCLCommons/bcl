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
