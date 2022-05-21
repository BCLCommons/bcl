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

#ifndef BCL_STORAGE_ROW_H_
#define BCL_STORAGE_ROW_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_map.h"
#include "bcl_storage_table_header.h"
#include "math/bcl_math_comparisons.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Row
    //! @brief This is a class for storing a row of values with their associated column names
    //! @details This class can is used by Table class. In addition, it should be preferred of Vector or other sequence
    //! containers where a string ( a column name) is used to access data
    //!
    //! @tparam t_DataType indicates the type of data that the row will hold
    //!
    //! @see @link example_storage_row.cpp @endlink
    //! @author karakam
    //! @date 06/24/07
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Row :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

      friend class Table< t_DataType>;

    private:

    //////////
    // data //
    //////////

      util::ShPtr< TableHeader> m_Header; //! header vector
      Vector< t_DataType>       m_Data;   //! actual data stored in a vector

    public:

      //! typedef for iterator
      typedef typename std::vector< t_DataType>::iterator               iterator;

      //! typedef for const_iterator
      typedef typename std::vector< t_DataType>::const_iterator         const_iterator;

      //! typedef for reverse_iterator
      typedef typename std::vector< t_DataType>::reverse_iterator       reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef typename std::vector< t_DataType>::const_reverse_iterator const_reverse_iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      Row< t_DataType>() :
        m_Header(),
        m_Data()
      {
      }

    protected:

      //! @brief construct a Row from number of size and pointer to data
      //! @param DATA data to be inserted into row
      //! @param SP_HEADER ShPtr to table header
      Row< t_DataType>( const util::ShPtr< TableHeader> &SP_HEADER, const Vector< t_DataType> &DATA) :
        m_Header( SP_HEADER),
        m_Data( DATA)
      {
      }

    public:

      //! @brief construct a Row from another ROW
      //! @param ROW Row< t_DataType> to copy
      Row< t_DataType>( const Row< t_DataType> &ROW) :
        m_Header( ROW.m_Header),
        m_Data( ROW.m_Data)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to a copy of the actual object
      Row< t_DataType> *Clone() const
      {
        return new Row< t_DataType>( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns data
      //! @return data
      const Vector< t_DataType> &GetData() const
      {
        return m_Data;
      }

      //! @brief return changeable reference to element at POS
      //! @param POS the position of the returned element
      //! @return changeable reference to element at POS
      t_DataType &operator()( const size_t POS)
      {
        return m_Data( POS);
      }

      //! @brief return const_reference to element at POS
      //! @param POS the position of the returned element
      //! @return const_reference to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        return m_Data( POS);
      }

      //! @brief return changeable reference to element for given COLUMN_NAME
      //! @param COLUMN_NAME name of the column of the returned element
      //! @return changeable reference to element for given COLUMN_NAME
      t_DataType &operator[]( const std::string &COLUMN_NAME)
      {
        return m_Data( m_Header->operator[]( COLUMN_NAME));
      }

      //! @brief return const_reference to element for given COLUMN_NAME
      //! @param COLUMN_NAME name of the column of the returned element
      //! @return const_reference to element for given COLUMN_NAME
      t_DataType const &operator[]( const std::string &COLUMN_NAME) const
      {
        return m_Data( m_Header->operator[]( COLUMN_NAME));
      }

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief returns a const reference to TableHeader
      //! @return const reference to TableHeader
      const TableHeader &GetHeader() const
      {
        return *m_Header;
      }

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      const Vector< t_DataType> &InternalData() const
      {
        return m_Data;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin()
      {
        return m_Data.Begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const
      {
        return m_Data.Begin();
      }

      //! @brief return changeable iterator to last element
      //! @return changeable iterator to last element
      iterator Last()
      {
        return m_Data.Last();
      }

      //! @brief return const iterator to last element
      //! @return const iterator to last element
      const_iterator Last() const
      {
        return m_Data.Last();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End()
      {
        return m_Data.End();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const
      {
        return m_Data.End();
      }

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin()
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd()
      {
        return m_Data.ReverseEnd();
      }

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const
      {
        return m_Data.ReverseEnd();
      }

      //! @brief return const reference to first element
      //! @return const reference to first element of type t_DataType
      t_DataType const &FirstElement() const
      {
        return m_Data.FirstElement();
      }

      //! @brief return a changeable reference to first element
      t_DataType &FirstElement()
      {
        return m_Data.FirstElement();
      }

      //! @brief return const reference to last element
      //! @return const reference to last element of type t_DataType
      t_DataType const &LastElement() const
      {
        return m_Data.LastElement();
      }

      //! @brief return changeable reference to last element
      //! @return changeable reference to last element of type t_DataType
      t_DataType &LastElement()
      {
        return m_Data.LastElement();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief set all elements of the sequence to one given value X
      //! @param ELEMENT the element of type t_DataType which will be assigned to all elements, default is t_DataType()
      void SetAllElements( const t_DataType &ELEMENT = t_DataType())
      {
        m_Data.SetAllElements( ELEMENT);
      }

      //! @brief converts the row data to a map with headers as the keys pointing to values in the row
      //! @return a map that contains the information in this row with headers as the keys
      Map< std::string, t_DataType> ConvertToMap() const
      {
        // initialize the map to be returned
        Map< std::string, t_DataType> this_map;

        // initialize iterators
        typename Vector< t_DataType>::const_iterator row_itr( m_Data.Begin()), row_itr_end( m_Data.End());
        Vector< std::string>::const_iterator header_itr( m_Header->Begin()), header_itr_end( m_Header->End());

        // iterate over the row and the header at the same time
        for( ; row_itr != row_itr_end && header_itr != header_itr_end; ++row_itr, ++header_itr)
        {
          // insert the corresponding map and row iterator
          this_map[ *header_itr] = *row_itr;
        }

        // end
        return this_map;
      }

      //! @brief converts the row data to a list with headers as the first in pair pointing to values in the row
      //!        This is similar to ConvertToMap but since this uses a list it does not change order of header strings
      //! @return a list of pairs that contains the information in this row with headers as the first pair element
      List< Pair< std::string, t_DataType> > ConvertToPairList() const
      {
        // initialize the list to be returned
        List< Pair< std::string, t_DataType> > this_list;

        // initialize iterators
        typename Vector< t_DataType>::const_iterator row_itr( m_Data.Begin()), row_itr_end( m_Data.End());
        Vector< std::string>::const_iterator header_itr( m_Header->Begin()), header_itr_end( m_Header->End());

        // iterate over the row and the header at the same time
        for( ; row_itr != row_itr_end && header_itr != header_itr_end; ++row_itr, ++header_itr)
        {
          // insert the corresponding header and row iterator
          this_list.PushBack( Pair< std::string, t_DataType>( *header_itr, *row_itr));
        }

        // end
        return this_list;
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read container from io::IFStream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read in size of container
        io::Serialize::Read( m_Header, ISTREAM);
        io::Serialize::Read( m_Data, ISTREAM);

        // return
        return ISTREAM;
      }

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size of container
        io::Serialize::Write( m_Header, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    public:

      //! @brief write row in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param FORMAT_COLS the format of the columns
      //! @return output stream which was written to
      std::ostream &WriteFormatted
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_COLS = util::Format().W( 10)
      ) const
      {

        // iterate over columns
        for
        (
          const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
        )
        {
          OSTREAM << FORMAT_COLS( *itr) << TableHeader::s_ColSeparator;
        }

        // put endline
        OSTREAM << '\n';

        // end
        return OSTREAM;
      }

      //! @brief write row in a formatted way to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param COLUMN_FORMATS vector of column formats
      //! @return output stream which was written to
      std::ostream &WriteFormatted
      (
        std::ostream &OSTREAM,
        const std::vector< util::Format> &COLUMN_FORMATS
      ) const
      {
        // make sure there are enough formats
        BCL_Assert
        (
          m_Data.GetSize() <= COLUMN_FORMATS.size(), "The provided column formats does not have enough entries"
        );

        // initialize iterator on the format
        std::vector< util::Format>::const_iterator format_itr( COLUMN_FORMATS.begin());

        // iterate over columns
        for
        (
          const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr, ++format_itr
        )
        {
          OSTREAM << format_itr->operator()( *itr) << TableHeader::s_ColSeparator;
        }

        // put endline
        OSTREAM << '\n';

        // end
        return OSTREAM;
      }

    /////////////
    // sorting //
    /////////////

    }; // template class Row

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RowComparison
    //! @brief class to compare two rows according to their values for a specified column
    //!
    //! @see @link example_storage_row.cpp @endlink
    //! @author karakam
    //! @date 06/24/07
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    struct RowComparison :
      public util::BinaryFunctionInterface< Row< t_DataType>, Row< t_DataType>, bool>
    {
    private:

      //! Index of column to be compared
      const size_t m_ColumnIndex;

      //! binary predicate function
      util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> > m_BinaryPredicate;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param COLUMN_INDEX index of the column to be compared
      //! @param BINARY_PREDICATE binary predicate for comparing the values at given columns of two rows
      RowComparison< t_DataType>
      (
        const size_t COLUMN_INDEX,
        const util::BinaryFunctionInterface< t_DataType, t_DataType, bool> &BINARY_PREDICATE =
        ( **math::Comparisons< t_DataType>::GetEnums().e_Less)
      ) :
        m_ColumnIndex( COLUMN_INDEX),
        m_BinaryPredicate( BINARY_PREDICATE.Clone())
      {
      }

      //! virtual copy constructor
      RowComparison< t_DataType> *Clone() const
      {
        return new RowComparison< t_DataType>( *this);
      }

      //! @brief compares two given rows for given column index
      //! @param ROW_A first row
      //! @param ROW_B second row
      //! @return true if ROW_A comes before ROW_B when sorted by the specified column ascending
      bool operator()( const Row< t_DataType> &ROW_A, const Row< t_DataType> &ROW_B) const
      {
        return m_BinaryPredicate->operator()( ROW_A( m_ColumnIndex), ROW_B( m_ColumnIndex));
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // template class Row

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Row< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Row< t_DataType>())
    );

    //! @brief operator == checks if two Rows are the same
    //! @param ROW_A first container
    //! @param ROW_B second container
    //! @return true, if Rows are identical
    template< typename t_DataType>
    inline bool operator ==
    (
      const Row< t_DataType> &ROW_A,
      const Row< t_DataType> &ROW_B
    )
    {
      return ROW_A.InternalData() == ROW_B.InternalData();
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_ROW_H_
