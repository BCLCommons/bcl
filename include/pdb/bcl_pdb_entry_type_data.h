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

#ifndef BCL_PDB_ENTRY_TYPE_DATA_H_
#define BCL_PDB_ENTRY_TYPE_DATA_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line_types.h"
#include "util/bcl_util_cpp_data_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EntryTypeData
    //! @brief for organizing pdb entry type specific constants
    //! @details protein data bank files are organized in different line types. each line has different entries, which are
    //! written in specific places within that line. This class stores, for each possible entry in every line, the
    //! LineType this entry belongs to, the position of the entry within that line, the length of the entry and the data
    //! type (double, string..).
    //! It assists in retrieving data from pdb lines easily and also in writing data to a pdb line in the proper place.
    //! With the data type it is possible to perform checks on the data retrieved and written to a pdb line.
    //!
    //! @see @link example_pdb_entry_type_data.cpp @endlink
    //! @author woetzen, staritrd
    //! @date Feb 8, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EntryTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      LineType            m_LineType;     //!< type of line where this entry can be found in
      size_t              m_Start;        //!< the location within the string of the line, where the entry can be found
      size_t              m_Length;       //!< length of string
      size_t              m_Precision;    //!< Precision for a double data type
      bool                m_RightAligned; //!< is entry right bound
      util::CPPDataTypes::TypeEnum m_DataType; //!< the data type
      bool                m_IsNumerical;  //!< Is the data type numerical

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EntryTypeData();

      //! @brief construct from information abput entry
      //! @param LINE_TYPE
      //! @param START
      //! @param LENGTH
      //! @param PRECISION
      //! @param RIGHT_ALIGNED
      //! @param DATA_TYPE
      EntryTypeData
      (
        const LineType &LINE_TYPE,
        const size_t START,
        const size_t LENGTH,
        const size_t PRECISION,
        const bool RIGHT_ALIGNED,
        const util::CPPDataTypes::Types DATA_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new EntryTypeData
      EntryTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! returns the line type
      const LineType &GetLineType() const
      {
        return m_LineType;
      }

      //! return the start of the entry
      size_t GetStart() const
      {
        return m_Start;
      }

      //! return length of entry
      size_t GetLength() const
      {
        return m_Length;
      }

      //! returns the data type of an entry
      const util::CPPDataTypes::Types &GetDataType() const
      {
        return m_DataType;
      }

      //! @brief return true if the given entry type is numeric
      bool IsNumeric() const
      {
        return m_IsNumerical;
      }

      //! @brief construct a format object
      util::Format GetFormat() const;

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

    }; // EntryTypeData

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_ENTRY_TYPE_DATA_H_
