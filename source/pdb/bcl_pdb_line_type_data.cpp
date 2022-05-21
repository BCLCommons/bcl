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
#include "pdb/bcl_pdb_line_type_data.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LineTypeData::LineTypeData() :
      m_MultipleTimes( false),
      m_MultipleLines( false),
      m_First( util::GetUndefined< size_t>()),
      m_Last( util::GetUndefined< size_t>())
    {
    }

    //! @brief construct from record type
    //! @param MULTIPLE_TIMES occurs multiple times
    //! @param MULTIPLE_LINES over multiple lines
    LineTypeData::LineTypeData( const bool MULTIPLE_TIMES, const bool MULTIPLE_LINES) :
      m_MultipleTimes( MULTIPLE_TIMES),
      m_MultipleLines( MULTIPLE_LINES),
      m_First( util::GetUndefined< size_t>()),
      m_Last( util::GetUndefined< size_t>())
    {
    }

    //! @brief Clone function
    //! @return pointer to new LineTypeData
    LineTypeData *LineTypeData::Clone() const
    {
      return new LineTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LineTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief the first entry type for that linetype
    //! @return entry type
    EntryType LineTypeData::GetFirstEntryType() const
    {
      return EntryType( m_First);
    }

    //! @brief the last entry type for that linetype
    //! @return entry type
    EntryType LineTypeData::GetLastEntryType() const
    {
      return EntryType( m_Last);
    }

    //! @brief number fo entry types for tat linetype
    //! @return number of entries
    size_t LineTypeData::GetNumberOfEntryTypes() const
    {
      return m_Last - m_First + 1;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief consider an additional EntryType
    //! @param ENTRY_TYPE
    void LineTypeData::ConsiderNewEntryType( const EntryType &ENTRY_TYPE)
    {
      // update last
      m_Last = ENTRY_TYPE.GetIndex();

      // check if first is still undefined
      if( !util::IsDefined( m_First))
      {
        m_First = m_Last;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LineTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_First, ISTREAM);
      io::Serialize::Read( m_Last, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &LineTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_First, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Last, OSTREAM, 0);

      // return
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
