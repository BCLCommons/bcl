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
#include "storage/bcl_storage_table_header.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> TableHeader::s_Instance
    (
      GetObjectInstances().AddInstance( new TableHeader())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TableHeader::TableHeader() :
      Vector< std::string>()
    {
    }

    //! @brief construct TableHeader from a vector of strings
    //! @param TABLE_HEADER vector of strings
    TableHeader::TableHeader( const Vector< std::string> &TABLE_HEADER) :
      Vector< std::string>( TABLE_HEADER)
    {
      // ensure that all column names are unique
      BCL_Assert
      (
        Set< std::string>( TABLE_HEADER.Begin(), TABLE_HEADER.End()).GetSize() == TABLE_HEADER.GetSize(),
        "The vector of column names provided are not unique : " + util::Format()( TABLE_HEADER)
      );
    }

    //! @brief copy constructor
    //! @param TABLE_HEADER for which a copy is desired
    TableHeader::TableHeader( const TableHeader &TABLE_HEADER) :
      Vector< std::string>( TABLE_HEADER)
    {
    }

    //! @brief virtual copy constructor, needed to make hard copies from pointer to it
    TableHeader *TableHeader::Clone() const
    {
      return new TableHeader( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &TableHeader::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief check if there is a col with given name
    //! @param COLUMN_NAME name of col
    //! @return true if col with that name exists
    bool TableHeader::HasColumn( const std::string &COLUMN_NAME) const
    {
      // check if col name can be found
      return std::find( Begin(), End(), COLUMN_NAME) != End();
    }

    //! @brief check if this header includes all columns with the given names
    //! @param COLUMN_NAME_VECTOR vector of names of cols
    //! @return true if all column with the given names exists
    bool TableHeader::HasColumns( const Vector< std::string> &COLUMN_NAME_VECTOR) const
    {
      // iterate over vector
      for
      (
        Vector< std::string>::const_iterator
          name_itr( COLUMN_NAME_VECTOR.Begin()), name_itr_end( COLUMN_NAME_VECTOR.End());
        name_itr != name_itr_end; ++name_itr
      )
      {
        // if any of the columns do not exist, then return false
        if( !HasColumn( *name_itr))
        {
          BCL_MessageStd( "No column found w/ name: " + *name_itr);
          return false;
        }
      }

      // return true, this is reached only if all columns existed
      return true;
    }

    //! @brief create a vector that contains format for each column using the given template
    //!        if a column name is longer than the given format, then the format will be updated
    //! @param TEMPLATE_FORMAT format to be used as template
    //! @return a vector that contains format for each column
    std::vector< util::Format> TableHeader::GetColumnFormats( const util::Format &TEMPLATE_FORMAT) const
    {
      // initialize vector of formats
      std::vector< util::Format> formats;

      // iterate over the columns
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        // insert the template
        formats.push_back( TEMPLATE_FORMAT);

        // if the template format is shorter than the length
        if( TEMPLATE_FORMAT.GetWidth() < itr->size())
        {
          // update the width
          formats.back().W( itr->size());
        }
      }

      // end
      return formats;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return index for the given COLUMN_NAME
    //! @param COLUMN_NAME name of the column of interest
    //! @return index for the given COLUMN_NAME
    size_t TableHeader::operator[]( const std::string &COLUMN_NAME) const
    {
      // calculate the index by using std::find to search for COLUMN_NAME
      const size_t index( std::find( Begin(), End(), COLUMN_NAME) - Begin());

      // ensure that the COLUMN NAME is found within this header
      BCL_Assert
      (
        index < GetSize(),
        "Could not find the provided column name \"" + COLUMN_NAME + "\" in this header " + util::Format()( *this)
      );

      // end
      return index;
    }

    //! @brief operator = for equating two TableHeaders
    //! @param TABLE_HEADER the right hand operand of the operator to equate the left hand operand to
    //! @return gives the TableHeader which has been equated to TABLE_HEADER
    TableHeader &TableHeader::operator =( const TableHeader &TABLE_HEADER)
    {
      // update members
      Vector< std::string>::operator =( TABLE_HEADER);

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Read inputs a TableHeader from a stream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &TableHeader::Read( std::istream &ISTREAM)
    {
      // read base class
      Vector< std::string>::Read( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief Write outputs TableHeader to a stream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &TableHeader::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      Vector< std::string>::Write( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write table header in a formatted way to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param FORMAT_COLS column formating flags
    //! @return output stream which was written to
    std::ostream &
    TableHeader::WriteFormatted
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_COLS
    ) const
    {
      // iterate over column names
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        OSTREAM << FORMAT_COLS( *itr) << s_ColSeparator;
      }

      // end
      return OSTREAM;
    }

    //! @brief write table header in a formatted way to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param COLUMN_FORMATS vector of column formats
    //! @return output stream which was written to
    std::ostream &TableHeader::WriteFormatted
    (
      std::ostream &OSTREAM,
      const std::vector< util::Format> &COLUMN_FORMATS
    ) const
    {
      // make sure the number of columns is same as the number of formats given
      BCL_Assert( GetSize() <= COLUMN_FORMATS.size(), "The provided column formats does not have enough entries");

      // initialize format iterators
      std::vector< util::Format>::const_iterator format_itr( COLUMN_FORMATS.begin());

      // iterate over the column names and the formats at the same time
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr, ++format_itr)
      {
        OSTREAM << format_itr->operator()( *itr) << s_ColSeparator;
      }

      // end
      return OSTREAM;

    }

  } // namespace storage
} // namespace bcl
