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

#ifndef BCL_STORAGE_TABLE_HPP_
#define BCL_STORAGE_TABLE_HPP_

// include the header of this class
#include "bcl_storage_table.h"

// includes from bcl - sorted alphabetically
#include "bcl_storage.h"
#include "bcl_storage_row.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <functional>

namespace bcl
{
  namespace storage
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Table< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Table< t_DataType>())
    );

    template< typename t_DataType>
    const std::string Table< t_DataType>::s_Indentation( "  ");

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Table< t_DataType>::Table() :
      m_Header(),
      m_Data()
    {
    }

    //! @brief construct a Table from a TABLE_HEADER
    //! @param TABLE_HEADER header of the table
    template< typename t_DataType>
    Table< t_DataType>::Table( const TableHeader &TABLE_HEADER) :
      m_Header( TABLE_HEADER.Clone()),
      m_Data()
    {
    }

    //! @brief construct a Table from a SP_TABLE_HEADER
    //! @param SP_TABLE_HEADER ShPtr to header of the table
    template< typename t_DataType>
    Table< t_DataType>::Table( const util::ShPtr< TableHeader> &SP_TABLE_HEADER) :
      m_Header( SP_TABLE_HEADER),
      m_Data()
    {
    }

    //! @brief construct a Table from a range of rows
    //! @param ITR_BEGIN itr that points to begin of range of rows
    //! @param ITR_END itr that points to the end of a range of rows
    template< typename t_DataType>
    Table< t_DataType>::Table( const const_iterator &ITR_BEGIN, const const_iterator &ITR_END) :
      m_Header(),
      m_Data( ITR_BEGIN, ITR_END)
    {
      // initialize header if a range was actually given
      if( ITR_BEGIN != ITR_END)
      {
        m_Header = util::ShPtr< TableHeader>( ITR_BEGIN->Second().GetHeader().Clone());
      }
    }

    //! @brief construct a Table from another Table
    //! @param TABLE Table< t_DataType> to copy
    template< typename t_DataType>
    Table< t_DataType>::Table( const Table< t_DataType> &TABLE) :
      m_Header( TABLE.m_Header),
      m_Data( TABLE.m_Data)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a copy of the actual object
    template< typename t_DataType>
    Table< t_DataType> *Table< t_DataType>::Clone() const
    {
      return new Table< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief GetClassIdentifier returns class name of the object
    //! @return returns string with the class name
    template< typename t_DataType>
    const std::string &Table< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a const reference to TableHeader
    //! @return const reference to TableHeader
    template< typename t_DataType>
    const TableHeader &Table< t_DataType>::GetHeader() const
    {
      return *m_Header;
    }

    //! @brief returns size of the container
    //! @return size, i.e. number of elements stored
    template< typename t_DataType>
    size_t Table< t_DataType>::GetSize() const
    {
      return m_Data.GetSize();
    }

    //! @brief returns the number of columns
    //! @return the number of columns
    template< typename t_DataType>
    size_t Table< t_DataType>::GetNumberColumns() const
    {
      return m_Header->GetSize();
    }

    //! @brief returns the number of rows
    //! @return the number of rows
    template< typename t_DataType>
    size_t Table< t_DataType>::GetNumberRows() const
    {
      return m_Data.GetSize();
    }

    //! @brief return the maximal size of the container
    //! @return maximum size, i.e. maximal number of elements to store
    template< typename t_DataType>
    size_t Table< t_DataType>::MaxSize() const
    {
      return m_Data.MaxSize();
    }

    //! @brief returns a const reference to the used stl container
    //! @return const reference to the internal stl container
    template< typename t_DataType>
    const List< Pair< std::string, Row< t_DataType> > > &Table< t_DataType>::InternalData() const
    {
      return m_Data;
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    template< typename t_DataType>
    typename Table< t_DataType>::iterator Table< t_DataType>::Begin()
    {
      return m_Data.Begin();
    }

    //! @brief return const_iterator on begin
    //! @return const_iterator pointing to the beginning of the container, i.e. the first element
    template< typename t_DataType>
    typename Table< t_DataType>::const_iterator Table< t_DataType>::Begin() const
    {
      return m_Data.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    template< typename t_DataType>
    typename Table< t_DataType>::iterator Table< t_DataType>::End()
    {
      return m_Data.End();
    }

    //! @brief return const_iterator on end
    //! @return const_iterator pointing to the end of the container, i.e. behind the last element
    template< typename t_DataType>
    typename Table< t_DataType>::const_iterator Table< t_DataType>::End() const
    {
      return m_Data.End();
    }

    //! @brief return iterator to reverse begin
    //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
    template< typename t_DataType>
    typename Table< t_DataType>::reverse_iterator Table< t_DataType>::ReverseBegin()
    {
      return m_Data.ReverseBegin();
    }

    //! @brief return const_iterator to reverse begin
    //! @return const_reverse_iterator pointing to the beginning of the reversed container
    template< typename t_DataType>
    typename Table< t_DataType>::const_reverse_iterator Table< t_DataType>::ReverseBegin() const
    {
      return m_Data.ReverseBegin();
    }

    //! @brief return iterator to reverse end
    //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
    template< typename t_DataType>
    typename Table< t_DataType>::reverse_iterator Table< t_DataType>::ReverseEnd()
    {
      return m_Data.ReverseEnd();
    }

    //! @brief return const_iterator to reverse end
    //! @return const_reverse_iterator pointing to the end of the reversed container
    template< typename t_DataType>
    typename Table< t_DataType>::const_reverse_iterator Table< t_DataType>::ReverseEnd() const
    {
      return m_Data.ReverseEnd();
    }

    //! @brief GetRowNames returns a list of the row names
    //! @return list of the row names
    template< typename t_DataType>
    const List< std::string> Table< t_DataType>::GetRowNames()
    {
      // initialize list of row names
      List< std::string> row_list;

      // iterate over rows
      for
      (
        iterator itr( m_Data.Begin()), itr_end( m_Data.End());
        itr != itr_end;
        ++itr
      )
      {
        // store this row name
        row_list.InsertElement( itr->First());
      }
      // end
      return row_list;
    }

    //! @brief returns the longest row name
    //! @return longest row name
    template< typename t_DataType>
    std::string Table< t_DataType>::GetLongestRowName() const
    {
      // initialize longest row name
      std::string longest_row_name( "");

      // iterate over rows
      for
      (
        const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
        itr != itr_end;
        ++itr
      )
      {
        // if this row has a longer name
        if( itr->First().length() > longest_row_name.length())
        {
          longest_row_name = itr->First();
        }
      }
      // end
      return longest_row_name;
    }

    //! @brief checks if a row with given name exists
    //! @param ROW_NAME name of the row of interest
    //! @return true if a row with the given ROW_NAME exists, false if not
    template< typename t_DataType>
    bool Table< t_DataType>::HasRow( const std::string &ROW_NAME) const
    {
      // returns true if a row with the given ROW_NAME exists, false if not
      return std::find_if
      (
        m_Data.Begin(), m_Data.End(), PairEqualFirst< std::string>( ROW_NAME)
      ) != m_Data.End();
    }

    //! @brief index of row with given name
    //! @param ROW_NAME name of the row of interest
    //! @return the index of the row with the given ROW_NAME - if not found return undefined
    template< typename t_DataType>
    size_t Table< t_DataType>::RowIndex( const std::string &ROW_NAME) const
    {
      // start with index 0
      size_t index( 0);

      // iterate over all rows
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr, ++index)
      {
        // if row name matches, return index
        if( itr->First() == ROW_NAME)
        {
          return index;
        }
      }

      // return undefined is nothing was found
      return util::GetUndefined< size_t>();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return row with given name
    //! @param ROW_NAME name of the row of interest
    //! @return row with given name
    template< typename t_DataType>
    Row< t_DataType> &Table< t_DataType>::operator[]( const std::string &ROW_NAME)
    {
      // find the row with the given identifier and return it
      iterator itr( std::find_if( m_Data.Begin(), m_Data.End(), PairEqualFirst< std::string>( ROW_NAME)));

      // assert that it is valid
      BCL_Assert( itr != m_Data.End(), "The provided Row does not exist with given name " + ROW_NAME);

      // return
      return itr->Second();
    }

    //! @brief return const row with given name
    //! @param ROW_NAME name of the row of interest
    //! @return const row with given name
    template< typename t_DataType>
    const Row< t_DataType> &Table< t_DataType>::operator[]( const std::string &ROW_NAME) const
    {
      // find the row with the given identifier
      const_iterator itr
      (
        std::find_if( m_Data.Begin(), m_Data.End(), PairEqualFirst< std::string>( ROW_NAME))
      );

      // assert that it is valid
      BCL_Assert( itr != m_Data.End(), "The provided Row does not exist with given name " + ROW_NAME);

      // return
      return itr->Second();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief checks whether container is empty
    //! @return if the container is empty or not
    template< typename t_DataType>
    bool Table< t_DataType>::IsEmpty() const
    {
      return m_Data.IsEmpty();
    }

    //! @brief delete all elements, keep the header
    template< typename t_DataType>
    void Table< t_DataType>::Reset()
    {
      m_Data.Reset();
    }

    //! @brief insert given DATA as a row with given ROW_NAME
    //! @param ROW_NAME name of the row to be created
    //! @param DATA data to be stored in the row to be created
    //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
    template< typename t_DataType>
    void Table< t_DataType>::InsertRow
    (
      const std::string &ROW_NAME,
      const Vector< t_DataType> &DATA,
      const bool INCLUDE_SIMILAR_ROW_NAMES
    )
    {
      if( !INCLUDE_SIMILAR_ROW_NAMES)
      {
        // assert that there is no row in table with given name
        BCL_Assert
        (
          !HasRow( ROW_NAME),
          "This table already contains a row with given name: " + ROW_NAME +
          " with data: " + util::Format()( DATA) + "\n" + util::Format()( *this)
        );
      }

      // check that the data has the correct size
      if( DATA.GetSize() != m_Header->GetSize())
      {
        BCL_MessageCrt
        (
          "This table has " + util::Format()( m_Header->GetSize()) + " columns but given data has " +
          util::Format()( DATA.GetSize()) + " so cannot be inserted!"
        );
        return;
      }

      // since format of given data is ensured insert it with the correct header
      m_Data.PushBack( Pair< std::string, Row< t_DataType> >( ROW_NAME, Row< t_DataType>( m_Header, DATA)));
    }

    //! @brief insert an empty row with given ROW_NAME
    //! @param ROW_NAME name of the row to be created
    //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
    //! @return reference to inserted row that was created
    template< typename t_DataType>
    Row< t_DataType> &Table< t_DataType>::InsertRow
    (
      const std::string &ROW_NAME,
      const bool INCLUDE_SIMILAR_ROW_NAMES
    )
    {
      if( !INCLUDE_SIMILAR_ROW_NAMES)
      {
        // assert that there is no row in table with given name
        BCL_Assert( !HasRow( ROW_NAME), "This table already contains a row with given name: " + ROW_NAME);
      }

      // since format of given data is ensured insert it with the correct header
      m_Data.PushBack
      (
        Pair< std::string, Row< t_DataType> >
        (
          ROW_NAME,
          Row< t_DataType>( m_Header, Vector< t_DataType>( m_Header->GetSize()))
        )
      );

      // return reference to inserted row
      return m_Data.LastElement().Second();
    }

    //! @brief insert empty rows with given ROW_NAMES
    //! @param ROW_NAMES vector of names of the rows to be created
    template< typename t_DataType>
    void Table< t_DataType>::InsertRows( const Vector< std::string> &ROW_NAMES)
    {
      // iterate over all row names
      for
      (
        Vector< std::string>::const_iterator name_itr( ROW_NAMES.Begin()), name_itr_end( ROW_NAMES.End());
        name_itr != name_itr_end;
        ++name_itr
      )
      {
        InsertRow( *name_itr);
      }
    }

    //! @brief remove row with given ROW_NAME
    //! @param ROW_NAME name of the row to be removed
    template< typename t_DataType>
    void Table< t_DataType>::RemoveRow( const std::string &ROW_NAME)
    {
      // find the row with the given identifier
      iterator itr( std::find_if( m_Data.Begin(), m_Data.End(), PairEqualFirst< std::string>( ROW_NAME)));

      // if found
      if( itr != m_Data.End())
      {
        // remove
        m_Data.Remove( itr);
      }
      else
      {
        BCL_MessageCrt( "This table does not have a row with given name: " + ROW_NAME);
      }
    }

    //! @brief delete single element at ITR
    //! @param ITR iterator pointing to the element that will be destroyed
    //! @return iterator pointing to the element immediately following the deleted one
    template< typename t_DataType>
    typename Table< t_DataType>::iterator Table< t_DataType>::Remove( iterator ITR)
    {
      return ITR == End() ? End() : m_Data.Remove( ITR);
    }

    //! @brief delete a range of elements [FIRST, LAST)
    //! @param ITR_FIRST first element of range to be deleted
    //! @param ITR_LAST last element, in the range, will be kept
    template< typename t_DataType>
    void Table< t_DataType>::Remove( iterator ITR_FIRST, iterator ITR_LAST)
    {
      m_Data.Remove( ITR_FIRST, ITR_LAST);
    }

    //! @brief append a TABLE to the end of this Table
    //! @param TABLE Table to be appended of same format
    //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
    template< typename t_DataType>
    void Table< t_DataType>::Append( const Table< t_DataType> &TABLE, const bool INCLUDE_SIMILAR_ROW_NAMES)
    {
      // check that table headers are the same
      BCL_Assert
      (
        *TABLE.m_Header == *m_Header,
        "The headers of two tables do not match\n this table has\n" + util::Format()( *m_Header) +
          "\n the given table has\n" + util::Format()( *TABLE.m_Header)
      );

      // iterate over every row in the given table
      for( const_iterator itr( TABLE.Begin()), itr_end( TABLE.End()); itr != itr_end; ++itr)
      {
        // call the insert on this row
        InsertRow( itr->First(), itr->Second().GetData(), INCLUDE_SIMILAR_ROW_NAMES);
      }
    }

    //! @brief insert pair of row with identical header
    template< typename t_DataType>
    void Table< t_DataType>::PushBack( const Pair< std::string, Row< t_DataType> > &NAME_ROW)
    {
      BCL_Assert
      (
        NAME_ROW.Second().GetHeader() == *m_Header,
        "unable to insert row with different header"
      );

      InsertRow( NAME_ROW.First(), NAME_ROW.Second().GetData());
    }

    //! @brief sort the table by given column
    //! @param COLUMN_NAME the name of the column according to which the data will be sorted
    //! @param BINARY_PREDICATE binary predicate to be used in sorting the table using a column's values
    template< typename t_DataType>
    void Table< t_DataType>::SortByColumn
    (
      const std::string &COLUMN_NAME,
      const util::BinaryFunctionInterface< t_DataType, t_DataType, bool> &BINARY_PREDICATE
    )
    {
      m_Data.Sort
      (
        PairBinaryPredicateSecond< std::string, Row< t_DataType> >
        (
          RowComparison< t_DataType>( ( *m_Header)[ COLUMN_NAME], BINARY_PREDICATE)
        )
      );
    }

    //! @brief extract a new table composed of the columns with the names in the specified table header
    //! @param SUB_TABLE_HEADER header of the subtable to be extracted
    //! @return a new table that only contains the given column names
    template< typename t_DataType>
    Table< t_DataType> Table< t_DataType>::ExtractSubTableFromColumns( const TableHeader &SUB_TABLE_HEADER) const
    {
      // store the number of columns in the table header
      const size_t nr_columns( SUB_TABLE_HEADER.GetSize());

      // initialize vector of indices that correspond to the given column names in the original table
      Vector< size_t> column_indices;
      column_indices.AllocateMemory( nr_columns);

      // iterate over the columns in the given table header
      for
      (
        Vector< std::string>::const_iterator
          col_itr( SUB_TABLE_HEADER.Begin()), col_itr_end( SUB_TABLE_HEADER.End());
        col_itr != col_itr_end; ++col_itr
      )
      {
        // insert into indices
        column_indices.PushBack( m_Header->operator []( *col_itr));
      }

      // create a new table
      Table< t_DataType> sub_table( SUB_TABLE_HEADER);

      // iterate over this table
      for( const_iterator row_itr( Begin()), row_itr_end( End()); row_itr != row_itr_end; ++row_itr)
      {
        // create a new data vector
        Vector< t_DataType> new_row_data;
        new_row_data.AllocateMemory( nr_columns);

        // iterate over the column indices
        for
        (
          Vector< size_t>::const_iterator index_itr( column_indices.Begin()), index_itr_end( column_indices.End());
          index_itr != index_itr_end; ++index_itr
        )
        {
          // pushback the value at the given index into the new row data vector
          new_row_data.PushBack( row_itr->Second()( *index_itr));
        }

        // insert the new row
        sub_table.InsertRow( row_itr->First(), new_row_data);
      }

      // end
      return sub_table;
    }

    //! @brief extract two columns from the table as a list of pairs
    //! @param COL_NAME_A first column of interest
    //! @param COL_NAME_B second column of interest
    template< typename t_DataType>
    List< Pair< t_DataType, t_DataType> > Table< t_DataType>::ExtractDataPairs
    (
      const std::string &COL_NAME_A,
      const std::string &COL_NAME_B
    ) const
    {
      // initialize list
      List< Pair< t_DataType, t_DataType> > data_pairs;

      // initialize indices to items in row
      const size_t index_a( m_Header->operator[]( COL_NAME_A));
      const size_t index_b( m_Header->operator[]( COL_NAME_B));

      // iterate over all rows
      for( const_iterator row_itr( Begin()), row_itr_end( End()); row_itr != row_itr_end; ++row_itr)
      {
        // insert pair
        data_pairs.PushBack(
          Pair< t_DataType, t_DataType>
          (
            row_itr->Second()( index_a),
            row_itr->Second()( index_b)
          )
        );
      }

      // end
      return data_pairs;
    }

    //! @brief return a List of tables of identical sizes that are different subsets of that table
    //! @param NUMBER_TABLES number fo tables that are created
    //! @param NUMBER_ROWS number of rows per table
    template< typename t_DataType>
    List< Table< t_DataType> > Table< t_DataType>::CreateDifferentSubTables
    (
      const size_t NUMBER_TABLES, const size_t NUMBER_ROWS
    ) const
    {
      // check that there are enough rows in that table
      if( NUMBER_ROWS > GetSize())
      {
        BCL_MessageCrt
        (
          "NUMBER_ROWS is larger than the number rows within that table " +
          util::Format()( NUMBER_ROWS) + " >" + util::Format()( GetSize())
        );

        return List< Table< t_DataType> >();
      }

      // check that resulting tables are actually different -
      // the requested number of rows has to be smaller by n -1 to the number of rows in that table to get n different tables
      if( GetSize() - NUMBER_ROWS < NUMBER_TABLES - 1)
      {
        BCL_MessageCrt
        (
          "cannot create " + util::Format()( NUMBER_TABLES) +
          " different tables since there are not enough rows. Difference between requested NUMBER_ROWS and the number"
          " of rows in that table has to be larger than " + util::Format()( NUMBER_TABLES - 1) + " but is: " +
          util::Format()( GetSize() - NUMBER_ROWS)
        );

        return List< Table< t_DataType> >();
      }

      // check that there are not more rows in that table than requested in total
      if( GetSize() > NUMBER_TABLES * NUMBER_ROWS)
      {
        BCL_MessageCrt
        (
          "total number of rows requested " + util::Format()( NUMBER_TABLES * NUMBER_ROWS) +
          " is less than size of table " + util::Format()( GetSize()) + " - resulting tables will miss rows"
        );
      }

      // calculate the shift
      const size_t shift( ( GetSize() - NUMBER_ROWS) / NUMBER_TABLES);

      // initialize list
      List< Table< t_DataType> > tables( NUMBER_TABLES, Table< t_DataType>( m_Header));
      typename List< Table< t_DataType> >::iterator table_itr( tables.Begin()), table_itr_end( tables.End());

      // initialize an iterator to the end
      const_iterator row_itr_end( End());

      // iterate over individual cross validation tables
      for
      (
        size_t index_table( 0);
        index_table < NUMBER_TABLES && table_itr != table_itr_end;
        ++index_table, ++table_itr
      )
      {
        // initialize iterator to the beginning
        const_iterator row_itr( Begin());
        // move to the starting position for this table
        AdvanceIterator( row_itr, row_itr_end, index_table * shift);

        // until the requested NUMBER_ROWS entries are inserted
        for( size_t nr_inserted_rows( 0); nr_inserted_rows < NUMBER_ROWS; ++nr_inserted_rows)
        {
          // push this row
          table_itr->PushBack( *row_itr);
          // move to the next row
          ++row_itr;
          if( row_itr == row_itr_end)
          {
            row_itr = Begin();
          }
        }
      }

      // end
      return tables;
    }

    //! @brief returns a transposed table where rows and columns are swapped
    //! @return transposed table where rows and columns are swapped
    template< typename t_DataType>
    Table< t_DataType> Table< t_DataType>::GetTransposedTable() const
    {
      // initialize new column names
      Vector< std::string> new_column_names;

      // store the number of columns
      const size_t nr_rows( m_Data.GetSize());

      // iterate over rows and collect row names that are to be used as new column names
      for( const_iterator row_itr( m_Data.Begin()), row_itr_end( m_Data.End()); row_itr != row_itr_end; ++row_itr)
      {
        new_column_names.PushBack( row_itr->First());
      }

      // create a table header from the column names
      util::ShPtr< TableHeader> new_table_header( new TableHeader( new_column_names));

      // initialize storage for the new row data with empty values
      List< Pair< std::string, Vector< t_DataType> > > new_rows_list;

      // iterate over current table header and initialize new rows
      for
      (
        Vector< std::string>::const_iterator col_itr( m_Header->Begin()), col_itr_end( m_Header->End());
        col_itr != col_itr_end; ++col_itr
      )
      {
        new_rows_list.PushBack
        (
          Pair< std::string, Vector< t_DataType> >
          (
            *col_itr, Vector< t_DataType>( nr_rows)
          )
        );
      }

      // initialize new_col_ctr
      size_t new_col_ctr( 0);
      // now iterate over the rows in the current table while incrementing new column counter
      for
      (
        const_iterator row_itr( m_Data.Begin()), row_itr_end( m_Data.End());
          row_itr != row_itr_end; ++row_itr,
        ++new_col_ctr
      )
      {
        // make a reference on this row
        const Row< t_DataType> &this_row( row_itr->Second());

        // now iterate over the columns while iterating at the same time over the new rows
        typename List< Pair< std::string, Vector< t_DataType> > >::iterator
          new_row_itr( new_rows_list.Begin()), new_row_itr_end( new_rows_list.End());

        size_t new_row_ctr( 0);

        for
        (
          typename Vector< t_DataType>::const_iterator
            col_itr( this_row.GetData().Begin()), col_itr_end( this_row.GetData().End());
          col_itr != col_itr_end && new_row_itr != new_row_itr_end; ++col_itr, ++new_row_itr
        )
        {
          // insert the data from row at the correct column for the new row
          new_row_itr->Second()( new_col_ctr) = *col_itr;
          ++new_row_ctr;
        }
      }

      // create a new table
      Table< t_DataType> new_table( new_table_header);

      // now iterate over the new row data and construct the actual rows and insert them into the new table
      for
      (
        typename List< Pair< std::string, Vector< t_DataType> > >::iterator
          row_itr( new_rows_list.Begin()), row_itr_end( new_rows_list.End());
        row_itr != row_itr_end; ++row_itr
      )
      {
        // construct and insert a new row
        new_table.InsertRow( row_itr->First(), row_itr->Second());
      }

      // return table
      return new_table;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read container from io::IFStream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Table< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read in members
      ISTREAM >> m_Header
              >> m_Data;

      // return
      return ISTREAM;
    }

    //! @brief read container from io::IFStream
    //! @param ISTREAM input stream
    //! @param INCLUDE_SIMILAR_ROW_NAMES allow more than one row with the same name in the table
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Table< t_DataType>::ReadFormatted
    (
      std::istream &ISTREAM,
      const bool INCLUDE_SIMILAR_ROW_NAMES
    )
    {
      std::string temp;

      // skip first item
      ISTREAM >> temp;

      // get rest of line
      std::getline( ISTREAM, temp);
      util::TrimString( temp);

      // read col headings
      Vector< std::string> items( util::SplitString( temp));
      m_Header = util::ShPtr< TableHeader>( new TableHeader( items));

      do
      {
        // acquire the next line and trim
        std::getline( ISTREAM, temp);
        util::TrimString( temp);

        // no entries in line
        if( temp.empty())
        {
          break;
        }

        // get all items in row
        items = util::SplitString( temp);

        // first item is row name
        const std::string row_name( items.FirstElement());

        // remaining items in row
        Vector< t_DataType> row_items;

        // iterate over all items and convert them to t_DataType
        row_items.AllocateMemory( items.GetSize() - 1);
        for
        (
          Vector< std::string>::const_iterator itr( items.Begin() + 1), itr_end( items.End());
          itr != itr_end;
          ++itr
        )
        {
          row_items.PushBack( util::ConvertStringToNumericalValue< t_DataType>( *itr));
        }

        // insert row
        InsertRow( row_name, row_items, INCLUDE_SIMILAR_ROW_NAMES);

        // break if end of file is reached
      } while( ISTREAM.good() && !ISTREAM.eof());

      // return
      return ISTREAM;
    }

    //! @brief write container to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Table< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Header
              << m_Data << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write table in a formatted way to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param FORMAT_COLS the format of the columns
    //! @param TITLE An optional title to be used instead of storage::Table< t_DataType>
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &
    Table< t_DataType>::WriteFormatted
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_COLS,
      const std::string &TITLE
    ) const
    {
      // initialize the column formats for each header
      const std::vector< util::Format> column_formats( m_Header->GetColumnFormats( FORMAT_COLS));

      // initialize format for row names
      const util::Format format_row
      (
        util::Format().L().W( std::max( GetLongestRowName().size(), TITLE.size()))
      );

      // if the title is not empty
      if( !TITLE.empty())
      {
        // output empty space of specified WIDTH
        OSTREAM << format_row( TITLE) << TableHeader::s_ColSeparator;

        // write table header
        m_Header->WriteFormatted( OSTREAM, column_formats) << '\n';
      }

      // iterate over every entry in the table
      for( const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr)
      {
        // print the row name
        OSTREAM << format_row( itr->First()) << TableHeader::s_ColSeparator;

        // print the row data
        itr->Second().WriteFormatted( OSTREAM, column_formats);
      }

      // end
      return OSTREAM;
    }

    //! @brief write table in a formatted way to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param FORMAT_COLS column formating flags
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &
    Table< t_DataType>::WriteFormattedWithoutNames
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_COLS
    ) const
    {
      // iterate over every entry in the table
      for( const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr)
      {
        // print the row data
        itr->Second().WriteFormatted( OSTREAM, FORMAT_COLS);
      }

      // output endline
      OSTREAM << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write sub table
    //! @param OSTREAM the stream to write to
    //! @param COL_NAMES the column names to be considered
    //! @param ROW_NAMES the row name to be considered
    //! @param FORMAT formatting to be used when writing the values in the subtable specified
    //! @param WRITE_HEADER whether or not to write header for that sub table
    template< typename t_DataType>
    std::ostream &Table< t_DataType>::WriteSubTable
    (
      std::ostream &OSTREAM,
      const Vector< std::string> &COL_NAMES,
      const Vector< std::string> &ROW_NAMES,
      const util::Format &FORMAT,
      const bool WRITE_HEADER
    ) const
    {
      return WriteSubTable
             (
               OSTREAM,
               COL_NAMES,
               ROW_NAMES,
               std::vector< util::Format>( COL_NAMES.GetSize(), FORMAT),
               WRITE_HEADER
             );
    }

    //! @brief write sub table
    //! @param OSTREAM the stream to write to
    //! @param COL_NAMES the column names to be considered
    //! @param ROW_NAMES the row name to be considered
    //! @param FORMATS formatting for each col in the sub table
    //! @param WRITE_HEADER whether or not to write header for that sub table
    template< typename t_DataType>
    std::ostream &Table< t_DataType>::WriteSubTable
    (
      std::ostream &OSTREAM,
      const Vector< std::string> &COL_NAMES,
      const Vector< std::string> &ROW_NAMES,
      const std::vector< util::Format> &FORMATS,
      const bool WRITE_HEADER
    ) const
    {
      const std::vector< util::Format>::const_iterator format_itr_end( FORMATS.end());

      // write header
      if( WRITE_HEADER)
      {
        std::vector< util::Format>::const_iterator format_itr( FORMATS.begin());
        // iterate over col names and check if they exist
        for
        (
          Vector< std::string>::const_iterator
            col_name_itr( COL_NAMES.Begin()), col_name_itr_end( COL_NAMES.End());
          col_name_itr != col_name_itr_end;
          ++col_name_itr
        )
        {
          if( m_Header->HasColumn( *col_name_itr) || *col_name_itr == GetStaticClassName( *this))
          {
            BCL_Assert( format_itr != format_itr_end, "not enough format objects given");
            OSTREAM << format_itr->operator ()( *col_name_itr) << TableHeader::s_ColSeparator;
            ++format_itr;
          }
          // message if a col with that name does not exist
          else
          {
            BCL_MessageVrb
            (
              "given col name \"" + *col_name_itr + "\" does not exist in this table! Skipping!"
            );
          }
        }

        // line break after header
        OSTREAM << '\n';
      }

      // once header is written write all rows
      for
      (
        Vector< std::string>::const_iterator
          row_name_itr( ROW_NAMES.Begin()), row_name_itr_end( ROW_NAMES.End());
        row_name_itr != row_name_itr_end;
        ++row_name_itr
      )
      {
        // find a row with the name
        const_iterator row_itr
        (
          std::find_if
          (
            m_Data.Begin(), m_Data.End(), PairEqualFirst< std::string>( *row_name_itr)
          )
        );

        // if row of that name exists
        if( row_itr == m_Data.End())
        {
          BCL_MessageVrb
          (
            "given row name \"" + *row_name_itr + "\" does not exist in this table! Skipping!"
          );
          continue;
        }

        std::vector< util::Format>::const_iterator format_itr( FORMATS.begin());
        // iterate over all col names
        for
        (
          Vector< std::string>::const_iterator
            col_name_itr( COL_NAMES.Begin()), col_name_itr_end( COL_NAMES.End());
          col_name_itr != col_name_itr_end;
          ++col_name_itr
        )
        {
          // write value
          if( m_Header->HasColumn( *col_name_itr))
          {
            BCL_Assert( format_itr != format_itr_end, "not enough format objects given");
            OSTREAM << format_itr->operator ()( row_itr->Second()[ *col_name_itr]) << TableHeader::s_ColSeparator;
            ++format_itr;
          }
          // write the row name
          else if( *col_name_itr == GetStaticClassName( *this))
          {
            BCL_Assert( format_itr != format_itr_end, "not enough format objects given");
            OSTREAM << format_itr->operator ()( row_itr->First()) << TableHeader::s_ColSeparator;
            ++format_itr;
          }
          // message if a col with that name does not exist
          else
          {
            BCL_MessageVrb
            (
              "given col name \"" + *col_name_itr + "\" does not exist in this table! Skipping!"
            );
          }
        }

        // add line break
        OSTREAM << '\n';
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief report the widest row name ( also compares the static class name)
    //! @return the widest row name ( also compares the static class name)
    template< typename t_DataType>
    const std::string &Table< t_DataType>::GetWidestRowName() const
    {
      const std::string *widest_name( &GetStaticClassName< Table< t_DataType> >());
      size_t max_width( widest_name->length());

      // iterate over every entry in the table
      for( const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr)
      {
        // compare length
        if( itr->First().length() > max_width)
        {
          max_width = itr->First().length();
          widest_name = &itr->First();
        }
      }

      // end
      return *widest_name;
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_TABLE_HPP_
