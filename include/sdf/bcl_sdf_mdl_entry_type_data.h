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

#ifndef BCL_SDF_MDL_ENTRY_TYPE_DATA_H_
#define BCL_SDF_MDL_ENTRY_TYPE_DATA_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_line_types.h"
#include "util/bcl_util_cpp_data_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MdlEntryTypeData
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_sdf_mdl_entry_type_data.cpp @endlink
    //! @author butkiem1
    //! @date Feb 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MdlEntryTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      MdlLineTypeEnum     m_MdlLineType; //!< line type this entry is associated with
      size_t              m_Start;       //!< position of entry
      size_t              m_Length;      //!< length of entry
      std::string         m_Default;     //!< default value
      util::Format        m_Format;      //!< Format of the entry type
      util::CPPDataTypes::TypeEnum m_DataType; //!< datatype of that entry

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MdlEntryTypeData();

      //! @brief construct from data
      //! @param LINE_TYPE type of line the entry is found in
      //! @param START starting position in line
      //! @param LENGTH length of entry type in line
      //! @param DEFAULT default string
      //! @param FORMAT format for the data field
      //! @param DATA_TYPE datatype of this entry
      MdlEntryTypeData
      (
        const MdlLineType LINE_TYPE,
        const size_t START,
        const size_t LENGTH,
        const std::string &DEFAULT,
        const util::Format &FORMAT,
        const util::CPPDataTypes::Types DATA_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new MdlEntryTypeData
      MdlEntryTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! returns the mdl line type
      const MdlLineTypeEnum &GetMdlLineType() const
      {
        return m_MdlLineType;
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

      //! @brief default value for that entry
      //! @return strign represeting the default
      const std::string &GetDefault() const
      {
        return m_Default;
      }

      //! returns the data type of an entry
      const util::CPPDataTypes::Types &GetDataType() const
      {
        return m_DataType;
      }

      //! @brief return the format object
      //! @return Format object to properly format the data to be inserted in a mdl line
      const util::Format &GetFormat() const
      {
        return m_Format;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get the string from a particular line, without any beginning or end spaces
      //! @param LINE the line to retrieve this entry from
      //! @return the entry from that line
      std::string GetTrimmedString( const std::string &LINE) const;

      //! @brief get this entry as an unsigned int from a particular line
      //! MDL integer types are limited to 999, so using a size_t is unnecessary
      //! @param LINE the line to retrieve this entry from
      //! @return the unsigned int given by this entry from a particular line
      bool IsUnsignedInt( const std::string &LINE) const;

      //! @brief get this entry as an unsigned int from a particular line
      //! MDL integer types are limited to 999, so using a size_t is unnecessary
      //! @param LINE the line to retrieve this entry from
      //! @return the unsigned int given by this entry from a particular line
      unsigned int GetUnsignedInt( const std::string &LINE) const;

      //! @brief get this entry as a double from a particular line
      //! @param LINE the line to retrieve this entry from
      //! @return the double given by this entry from a particular line
      bool IsDouble( const std::string &LINE) const;

      //! @brief get this entry as a double from a particular line
      //! @param LINE the line to retrieve this entry from
      //! @return the double given by this entry from a particular line
      double GetDouble( const std::string &LINE) const;

      //! @brief write this entry to mdl line as a t_DataType
      //! @param STRING string to accept the given entry
      //! @param DATA size_t to be inserted
      template< typename t_DataType>
      void Set( std::string &STRING, const t_DataType &DATA) const
      {
        BCL_Assert( STRING.size() >= m_Start + m_Length, "Not enough space in string to set entry");

        // width has to be enforced, otherwise line will shrink
        STRING.replace( m_Start, m_Length, m_Format( DATA), 0, m_Length);
      }

    ///////////////
    // operators //
    ///////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class MdlEntryTypeData

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_MDL_ENTRY_TYPE_DATA_H_
