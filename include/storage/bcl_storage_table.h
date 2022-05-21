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

#ifndef BCL_STORAGE_TABLE_H_
#define BCL_STORAGE_TABLE_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_list.h"
#include "bcl_storage_pair.h"
#include "bcl_storage_row.h"
#include "bcl_storage_table_header.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Table
    //! @brief Table class allows storing information in rows and columns with descriptions for both
    //!
    //! @tparam t_DataType indicates the type of data that the Table will hold
    //!
    //! @see @link example_storage_table.cpp @endlink
    //! @author karakam, woetzen
    //! @date 05/15/07
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Table :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Table header
      util::ShPtr< TableHeader> m_Header;

      //! list of pairs of rows and their descriptors
      List< Pair< std::string, Row< t_DataType> > > m_Data;

    public:

      //! typedef for iterator
      typedef typename List< Pair< std::string, Row< t_DataType> > >::iterator
        iterator;

      //! typedef for const_iterator
      typedef typename List< Pair< std::string, Row< t_DataType> > >::const_iterator
        const_iterator;

      //! typedef for reverse_iterator
      typedef typename List< Pair< std::string, Row< t_DataType> > >::reverse_iterator
        reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef typename List< Pair< std::string, Row< t_DataType> > >::const_reverse_iterator
        const_reverse_iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! indentation to use for subheadings in the table
      static const std::string s_Indentation;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Table();

      //! @brief construct a Table from a TABLE_HEADER
      //! @param TABLE_HEADER header of the table
      Table( const TableHeader &TABLE_HEADER);

      //! @brief construct a Table from a SP_TABLE_HEADER
      //! @param SP_TABLE_HEADER ShPtr to header of the table
      Table( const util::ShPtr< TableHeader> &SP_TABLE_HEADER);

      //! @brief construct a Table from a range of rows
      //! @param ITR_BEGIN itr that points to begin of range of rows
      //! @param ITR_END itr that points to the end of a range of rows
      Table( const const_iterator &ITR_BEGIN, const const_iterator &ITR_END);

      //! @brief construct a Table from another Table
      //! @param TABLE Table< t_DataType> to copy
      Table( const Table< t_DataType> &TABLE);

      //! @brief virtual copy constructor
      //! @return pointer to a copy of the actual object
      Table< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns a const reference to TableHeader
      //! @return const reference to TableHeader
      const TableHeader &GetHeader() const;

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const;

      //! @brief returns the number of columns
      //! @return the number of columns
      size_t GetNumberColumns() const;

      //! @brief returns the number of rows
      //! @return the number of rows
      size_t GetNumberRows() const;

      //! @brief return the maximal size of the container
      //! @return maximum size, i.e. maximal number of elements to store
      size_t MaxSize() const;

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      const List< Pair< std::string, Row< t_DataType> > > &InternalData() const;

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin();

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const;

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End();

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const;

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin();

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const;

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd();

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const;

      //! @brief GetRowNames returns a list of the row names
      //! @return list of the row names
      const List< std::string> GetRowNames();

      //! @brief returns the longest row name
      //! @return name
      std::string GetLongestRowName() const;

      //! @brief checks if a row with given name exists
      //! @param ROW_NAME name of the row of interest
      //! @return true if a row with the given ROW_NAME exists, false if not
      bool HasRow( const std::string &ROW_NAME) const;

      //! @brief index of row with given name
      //! @param ROW_NAME name of the row of interest
      //! @return the index of the row with the given ROW_NAME - if not found return undefined
      size_t RowIndex( const std::string &ROW_NAME) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return row with given name
      //! @param ROW_NAME name of the row of interest
      //! @return row with given name
      Row< t_DataType> &operator[]( const std::string &ROW_NAME);

      //! @brief return const row with given name
      //! @param ROW_NAME name of the row of interest
      //! @return const row with given name
      const Row< t_DataType> &operator[]( const std::string &ROW_NAME) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether container is empty
      //! @return if the container is empty or not
      bool IsEmpty() const;

      //! @brief delete all elements, keep the header
      void Reset();

      //! @brief insert given DATA as a row with given ROW_NAME
      //! @param ROW_NAME name of the row to be created
      //! @param DATA data to be stored in the row to be created
      //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
      void InsertRow
      (
        const std::string &ROW_NAME,
        const Vector< t_DataType> &DATA,
        const bool INCLUDE_SIMILAR_ROW_NAMES = false
      );

      //! @brief insert an empty row with given ROW_NAME
      //! @param ROW_NAME name of the row to be created
      //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
      //! @return reference to inserted row that was created
      Row< t_DataType> &InsertRow
      (
        const std::string &ROW_NAME,
        const bool INCLUDE_SIMILAR_ROW_NAMES = false
      );

      //! @brief insert empty rows with given ROW_NAMES
      //! @param ROW_NAMES vector of names of the rows to be created
      void InsertRows( const Vector< std::string> &ROW_NAMES);

      //! @brief remove row with given ROW_NAME
      //! @param ROW_NAME name of the row to be removed
      void RemoveRow( const std::string &ROW_NAME);

      //! @brief delete single element at ITR
      //! @param ITR iterator pointing to the element that will be destroyed
      //! @return iterator pointing to the element immediately following the deleted one
      iterator Remove( iterator ITR);

      //! @brief delete a range of elements [FIRST, LAST)
      //! @param ITR_FIRST first element of range to be deleted
      //! @param ITR_LAST last element, in the range, will be kept
      void Remove( iterator ITR_FIRST, iterator ITR_LAST);

      //! @brief append a TABLE to the end of this Table
      //! @param TABLE Table to be appended of same format
      //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
      void Append( const Table< t_DataType> &TABLE, const bool INCLUDE_SIMILAR_ROW_NAMES = false);

      //! @brief insert pair of row with identical header
      void PushBack( const Pair< std::string, Row< t_DataType> > &NAME_ROW);

      //! @brief sort the table by given column
      //! @param COLUMN_NAME the name of the column according to which the data will be sorted
      //! @param BINARY_PREDICATE binary predicate to be used in sorting the table using a column's values
      void SortByColumn
      (
        const std::string &COLUMN_NAME,
        const util::BinaryFunctionInterface< t_DataType, t_DataType, bool> &BINARY_PREDICATE =
        ( **math::Comparisons< t_DataType>::GetEnums().e_Less)
      );

      //! @brief extract a new table composed of the columns with the names in the specified table header
      //! @param SUB_TABLE_HEADER header of the subtable to be extracted
      //! @return a new table that only contains the given column names
      Table< t_DataType> ExtractSubTableFromColumns( const TableHeader &SUB_TABLE_HEADER) const;

      //! @brief sort the table by given column
      //! @param BINARY_PREDICATE ???
      template< typename t_BinaryPredicate>
      void SortByRowName( const t_BinaryPredicate &BINARY_PREDICATE)
      {
        m_Data.Sort
        (
          PairBinaryPredicateFirst< std::string, Row< t_DataType> >
          (
            util::BinaryFunctionSTLWrapper< t_BinaryPredicate>
            (
              BINARY_PREDICATE
            )
          )
        );
      }

      //! @brief extract two columns from the table as a list of pairs
      //! @param COL_NAME_A first column of interest
      //! @param COL_NAME_B second column of interest
      List< Pair< t_DataType, t_DataType> > ExtractDataPairs
      (
        const std::string &COL_NAME_A,
        const std::string &COL_NAME_B
      ) const;

      //! @brief return a List of tables of identical sizes that are different subsets of that table
      //! @param NUMBER_TABLES number fo tables that are created
      //! @param NUMBER_ROWS number of rows per table
      List< Table< t_DataType> > CreateDifferentSubTables( const size_t NUMBER_TABLES, const size_t NUMBER_ROWS) const;

      //! @brief returns a transposed table where rows and columns are swapped
      //! @return transposed table where rows and columns are swapped
      Table< t_DataType> GetTransposedTable() const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read container from io::IFStream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief read container from io::IFStream
      //! @param ISTREAM input stream
      //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
      //! @return istream which was read from
      std::istream &ReadFormatted
      (
        std::istream &ISTREAM,
        const bool INCLUDE_SIMILAR_ROW_NAMES = false
      );

      //! @brief write table in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param FORMAT_ROW_NAMES the format of the row names
      //! @param FORMAT_COLS the format of the columns
      //! @param TITLE An optional title to be used instead of storage::Table< t_DataType>
      //! @return output stream which was written to
      std::ostream &
      WriteFormatted
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_COLS = util::Format().W( 10).R(),
        const std::string &TITLE = GetStaticClassName< Table< t_DataType> >()
      ) const;

      //! @brief write table in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param FORMAT_COLS column formating flags
      //! @return output stream which was written to
      std::ostream &
      WriteFormattedWithoutNames
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_COLS = util::Format().W( 10)
      ) const;

      //! @brief write sub table
      //! @param OSTREAM the stream to write to
      //! @param COL_NAMES the column names to be considered
      //! @param ROW_NAMES the row name to be considered
      //! @param FORMAT formatting to be used when writing the values in the subtable specified
      //! @param WRITE_HEADER whether or not to write header for that sub table
      std::ostream &WriteSubTable
      (
        std::ostream &OSTREAM,
        const Vector< std::string> &COL_NAMES,
        const Vector< std::string> &ROW_NAMES,
        const util::Format &FORMAT = util::Format().W( 10),
        const bool WRITE_HEADER = false
      ) const;

      //! @brief write sub table
      //! @param OSTREAM the stream to write to
      //! @param COL_NAMES the column names to be considered
      //! @param ROW_NAMES the row name to be considered
      //! @param FORMATS formatting for each col in the sub table
      //! @param WRITE_HEADER whether or not to write header for that sub table
      std::ostream &WriteSubTable
      (
        std::ostream &OSTREAM,
        const Vector< std::string> &COL_NAMES,
        const Vector< std::string> &ROW_NAMES,
        const std::vector< util::Format> &FORMATS,
        const bool WRITE_HEADER
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief report the widest row name ( also compares the static class name)
      //! @return the widest row name ( also compares the static class name)
      const std::string &GetWidestRowName() const;

    }; // template class Table

    //! @brief operator == checks if two Tables are the same
    //! @param TABLE_A first container
    //! @param TABLE_B second container
    //! @return true, if tables are identical
    template< typename t_DataType>
    inline bool operator ==
    (
      const Table< t_DataType> &TABLE_A,
      const Table< t_DataType> &TABLE_B
    )
    {
      return TABLE_A.InternalData() == TABLE_B.InternalData();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Table< char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Table< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Table< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Table< size_t>;

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_TABLE_H_
