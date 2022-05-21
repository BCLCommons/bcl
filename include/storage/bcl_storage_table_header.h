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

#ifndef BCL_STORAGE_TABLE_HEADER_H_
#define BCL_STORAGE_TABLE_HEADER_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TableHeader
    //! @brief Vector of std::strings to hold column names for a row and a table
    //!
    //! @see @link example_storage_table_header.cpp @endlink
    //! @author karakam
    //! @date 5/25/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TableHeader :
      public Vector< std::string>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      // character to be used to separate coloumns in tables
      static const char s_ColSeparator = ' ';

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TableHeader();

      //! @brief construct TableHeader from a vector of strings
      //! @param TABLE_HEADER vector of strings
      TableHeader( const Vector< std::string> &TABLE_HEADER);

      //! @brief copy constructor
      //! @param TABLE_HEADER for which a copy is desired
      TableHeader( const TableHeader &TABLE_HEADER);

      //! @brief virtual copy constructor, needed to make hard copies from pointer to it
      TableHeader *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief check if there is a col with given name
      //! @param COLUMN_NAME name of col
      //! @return true if col with that name exists
      bool HasColumn( const std::string &COLUMN_NAME) const;

      //! @brief check if this header includes all columns with the given names
      //! @param COLUMN_NAME_VECTOR vector of names of cols
      //! @return true if all column with the given names exists
      bool HasColumns( const Vector< std::string> &COLUMN_NAME_VECTOR) const;

      //! @brief create a vector that contains format for each column using the given template
      //!        if a column name is longer than the given format, then the format will be updated
      //! @param TEMPLATE_FORMAT format to be used as template
      //! @return a vector that contains format for each column
      std::vector< util::Format> GetColumnFormats( const util::Format &TEMPLATE_FORMAT) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return index for the given COLUMN_NAME
      //! @param COLUMN_NAME name of the column of interest
      //! @return index for the given COLUMN_NAME
      size_t operator[]( const std::string &COLUMN_NAME) const;

      //! @brief operator = for equating two TableHeaders
      //! @param TABLE_HEADER the right hand operand of the operator to equate the left hand operand to
      //! @return gives the TableHeader which has been equated to TABLE_HEADER
      TableHeader &operator =( const TableHeader &TABLE_HEADER);

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief Read inputs a TableHeader from a stream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief Write outputs TableHeader to a stream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write table header in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param FORMAT_COLS column width
      //! @return output stream which was written to
      std::ostream &WriteFormatted
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_COLS = util::Format().W( 10)
      ) const;

      //! @brief write table header in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param COLUMN_FORMATS vector of column formats
      //! @return output stream which was written to
      std::ostream &WriteFormatted
      (
        std::ostream &OSTREAM,
        const std::vector< util::Format> &COLUMN_FORMATS
      ) const;

    }; // end TableHeader class

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_TABLE_HEADER_H_
