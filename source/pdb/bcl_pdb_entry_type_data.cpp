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
#include "pdb/bcl_pdb_entry_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EntryTypeData::EntryTypeData() :
      m_LineType( GetLineTypes().e_Undefined),
      m_Start( util::GetUndefined< size_t>()),
      m_Length( util::GetUndefined< size_t>()),
      m_DataType( util::CPPDataTypes::e_Unknown),
      m_IsNumerical( false)
    {
    }

    //! @brief construct from information about entry
    //! @param LINE_TYPE
    //! @param START
    //! @param LENGTH
    //! @param PRECISION
    //! @param RIGHT_ALIGNED
    //! @param DATA_TYPE
    EntryTypeData::EntryTypeData
    (
      const LineType &LINE_TYPE,
      const size_t START,
      const size_t LENGTH,
      const size_t PRECISION,
      const bool RIGHT_ALIGNED,
      const util::CPPDataTypes::Types DATA_TYPE
    ) :
      m_LineType( LINE_TYPE),
      m_Start( START),
      m_Length( LENGTH),
      m_Precision( PRECISION),
      m_RightAligned( RIGHT_ALIGNED),
      m_DataType( DATA_TYPE),
      m_IsNumerical
      (
        m_DataType == util::CPPDataTypes::e_Double ||
        m_DataType == util::CPPDataTypes::e_Float ||
        m_DataType == util::CPPDataTypes::e_Int ||
        m_DataType == util::CPPDataTypes::e_SizeT
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new EntryTypeData
    EntryTypeData *EntryTypeData::Clone() const
    {
      return new EntryTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &EntryTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief construct a format object
    util::Format EntryTypeData::GetFormat() const
    {
      // initialize format object
      util::Format format;

      // set width
      format.W( m_Length);

      // set precision for double or float types
      if( util::IsDefined( m_Precision))
      {
        format.FFP( m_Precision);
      }

      // set alignment
      m_RightAligned ? format.R() : format.L();

      // return
      return format;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EntryTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_LineType    , ISTREAM);
      io::Serialize::Read( m_Start       , ISTREAM);
      io::Serialize::Read( m_Length      , ISTREAM);
      io::Serialize::Read( m_Precision   , ISTREAM);
      io::Serialize::Read( m_RightAligned, ISTREAM);
      io::Serialize::Read( m_DataType    , ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &EntryTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_LineType    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Start       , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Length      , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_Precision   , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_RightAligned, OSTREAM,      0) << '\t';
      io::Serialize::Write( m_DataType    , OSTREAM,      0);

      // return
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
