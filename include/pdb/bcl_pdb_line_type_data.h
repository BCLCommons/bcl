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

#ifndef BCL_PDB_LINE_TYPE_DATA_H_
#define BCL_PDB_LINE_TYPE_DATA_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineTypeData
    //! @brief Data class that stores the number of entries in a pdb line
    //! @details this class is data class to be used in the enumeration of LineTypes
    //!
    //! @see @link example_pdb_line_type_data.cpp @endlink
    //! @author woetzen
    //! @date Feb 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      // relevant for type
      //! @see http://www.wwpdb.org/documentation/format33/sect1.html#Types @endsee

      //! occurs multiple times
      bool m_MultipleTimes;

      //! over multiple lines
      bool m_MultipleLines;

      //! index of first EntryType
      size_t m_First;

      //! index of last EntryType
      size_t m_Last;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LineTypeData();

      //! @brief construct from record type
      //! @param MULTIPLE_TIMES occurs multiple times
      //! @param MULTIPLE_LINES over multiple lines
      LineTypeData( const bool MULTIPLE_TIMES, const bool MULTIPLE_LINES);

      //! @brief Clone function
      //! @return pointer to new LineTypeData
      LineTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief the first entry type for that linetype
      //! @return entry type
      EntryType GetFirstEntryType() const;

      //! @brief the last entry type for that linetype
      //! @return entry type
      EntryType GetLastEntryType() const;

      //! @brief number fo entry types for tat linetype
      //! @return number of entries
      size_t GetNumberOfEntryTypes() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief consider an additional EntryType
      //! @param ENTRY_TYPE
      void ConsiderNewEntryType( const EntryType &ENTRY_TYPE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LineTypeData

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_LINE_TYPE_DATA_H_
