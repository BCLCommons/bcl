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

#ifndef BCL_CLUSTER_INPUT_TABLE_H_
#define BCL_CLUSTER_INPUT_TABLE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_distances_stored.h"
#include "bcl_cluster_input_interface.h"
#include "io/bcl_io_ifstream.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InputTable
    //! @brief InputTable is a method for getting putting the contents of a storage::Table into the data construct
    //!         used by cluster::DistanceStored.
    //!
    //! @see @link example_cluster_input_table.cpp @endlink
    //! @author alexanns
    //! @date August 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_PrecisionType>
    class InputTable :
      public InputInterface< std::string, t_PrecisionType>
    {
    private:

    //////////
    // data //
    //////////

      //! true if the only the upper triangle needs to be read
      bool m_ReadUpperTriangle;

      //! true if the only the lower triangle needs to be read
      bool m_ReadLowerTriangle;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      InputTable();

      //! @brief constructor from a file containing the data, the objects, and a boolean
      //! @param READ_UPPER_TRIANGLE is true if the only the upper triangle of the table needs to be read
      //! @param READ_LOWER_TRIANGLE is true if the only the lower triangle of the table needs to be read
      InputTable
      (
        const bool READ_UPPER_TRIANGLE,
        const bool READ_LOWER_TRIANGLE
      );

      //! @brief copy constructor
      //! @param INPUT_TABLE the InputTable from which this InputTable will be copied
      InputTable( const InputTable &INPUT_TABLE);

      //! @brief Clone function
      //! @return pointer to new InputMatrix< std::string>
      virtual InputTable *Clone() const;

      //! @brief virtual destructor
      virtual ~InputTable();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
      //!         in m_Objects
      //! @param IFSTREAM is the stream from which the input will be read
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      virtual
      util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType
        >
      >
      HandleInput( io::IFStream &IFSTREAM);

      //! @brief StoreEntireMatrix puts the entire matrix from a Table into memory
      //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
      //! @param DATA the construct where the data from TABLE will be inserted
      //! @param COLUMN_LAST_OBJECT_MAKER_ITERATOR
      void
      StoreEntireMatrix
      (
        const storage::Table< t_PrecisionType> &TABLE,
        storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA,
        const storage::List< std::string>::const_iterator &COLUMN_LAST_OBJECT_MAKER_ITERATOR
      ) const;

      //! @brief StoreUpperTriangle puts the upper triangle half matrix from a Table into memory
      //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
      //! @param DATA the construct where the data from TABLE will be inserted
      void
      StoreUpperTriangle
      (
        const storage::Table< t_PrecisionType> &TABLE,
        storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
      ) const;

      //! @brief StoreLowerTriangle puts the upper triangle half matrix from a Table into memory
      //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
      //! @param DATA the construct where the data from TABLE will be inserted
      void StoreLowerTriangle
      (
        const storage::Table< t_PrecisionType> &TABLE,
        storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
      ) const;

    }; // template class InputTable

    // instantiate s_Instance
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> InputTable< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new InputTable())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_PrecisionType>
    InputTable< t_PrecisionType>::InputTable() :
      InputInterface< std::string, t_PrecisionType>(),
      m_ReadUpperTriangle(),
      m_ReadLowerTriangle()
    {
    }

    //! @brief constructor from a file containing the data, the objects, and a boolean
    //! @param READ_UPPER_TRIANGLE is true if the only the upper triangle of the table needs to be read
    //! @param READ_LOWER_TRIANGLE is true if the only the lower triangle of the table needs to be read
    template< typename t_PrecisionType>
    InputTable< t_PrecisionType>::InputTable
    (
      const bool READ_UPPER_TRIANGLE,
      const bool READ_LOWER_TRIANGLE
    ) :
      InputInterface< std::string, t_PrecisionType>(),
      m_ReadUpperTriangle( READ_UPPER_TRIANGLE),
      m_ReadLowerTriangle( READ_LOWER_TRIANGLE)
    {
      BCL_Assert( m_ReadUpperTriangle || m_ReadLowerTriangle, "Must need at least one of the table triangles");
    }

    //! @brief copy constructor
    //! @param INPUT_TABLE the InputTable from which this InputTable will be copied
    template< typename t_PrecisionType>
    InputTable< t_PrecisionType>::InputTable( const InputTable &INPUT_TABLE) :
      InputInterface< std::string, t_PrecisionType>( INPUT_TABLE),
      m_ReadUpperTriangle( INPUT_TABLE.m_ReadUpperTriangle),
      m_ReadLowerTriangle( INPUT_TABLE.m_ReadLowerTriangle)
    {
    }

    //! @brief Clone function
    //! @return pointer to new InputTable
    template< typename t_PrecisionType>
    InputTable< t_PrecisionType> *InputTable< t_PrecisionType>::Clone() const
    {
      return new InputTable< t_PrecisionType>( *this);
    }

    //! @brief virtual destructor
    template< typename t_PrecisionType>
    InputTable< t_PrecisionType>::~InputTable< t_PrecisionType>()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_PrecisionType>
    const std::string &InputTable< t_PrecisionType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_PrecisionType>
    std::istream &InputTable< t_PrecisionType>::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( InputInterface< std::string, t_PrecisionType>::m_Objects, ISTREAM);
      io::Serialize::Read( m_ReadUpperTriangle, ISTREAM);
      io::Serialize::Read( m_ReadLowerTriangle, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_PrecisionType>
    std::ostream &InputTable< t_PrecisionType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( InputInterface< std::string, t_PrecisionType>::m_Objects, OSTREAM, INDENT);
      io::Serialize::Write( m_ReadUpperTriangle, OSTREAM, INDENT);
      io::Serialize::Write( m_ReadLowerTriangle, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
    //!         in m_Objects
    //! @param IFSTREAM is the stream from which the input will be read
    //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
    template< typename t_PrecisionType>
    util::ShPtr
    <
      math::FunctionInterface
      <
        storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType
      >
    >
    InputTable< t_PrecisionType>::HandleInput( io::IFStream &IFSTREAM)
    {
      util::Stopwatch clustering_timer
      (
        " Input Table HandleInput ", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical
      );

      // create Table "table" and read it in from "read"
      storage::Table< t_PrecisionType> table;

      // read in the contents of "table" from ISTREAM
      {
        util::Stopwatch read_table( "Read in table", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical);
        table.ReadFormatted( IFSTREAM);
        BCL_MessageDbg( "table is " + util::Format()( table));
      }

      util::Stopwatch define_objects( "Define objects");
      // create "object_list" to contain all objects and initialize with header objects of "table"
      util::ShPtr< storage::List< std::string> > object_list
      (
        new storage::List< std::string>( table.GetHeader().Begin(), table.GetHeader().End())
      );

      // create iterator "column_last_object_marker" which marks the end of the objects in the table header
      // this is needed if the entire matrix needs to be read in, so that we can distinguish in the "object_list"
      // between the objects in the table header and the objects in the row headers
      const storage::List< std::string>::iterator column_last_object_marker( --object_list->End());

      // true if whole table is needed
      // indicating the objects of the rows are different than the objects of the columns
      if( m_ReadUpperTriangle && m_ReadLowerTriangle)
      {
        // append the row names of "table" to "object_list"
        object_list->Append( table.GetRowNames());
      }

      // set the input objects
      InputInterface< std::string, t_PrecisionType>::SetInputObjects( object_list);

      util::Stopwatch data_watch( "fill data construct");
      // create "data" which will hold the data be returned
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > data;

      // true if the entire table must be read into memory
      if( m_ReadUpperTriangle && m_ReadLowerTriangle)
      {
        StoreEntireMatrix( table, data, column_last_object_marker);
      }
      else if( m_ReadUpperTriangle && !m_ReadLowerTriangle)//< only the upper triangle of the matrix is sufficient
      {
        BCL_MessageDbg( "store upper triangle");
        BCL_MessageDbg( "table is " + util::Format()( table));
        StoreUpperTriangle( table, data);
      }
      else if( !m_ReadUpperTriangle && m_ReadLowerTriangle)//< only the lower triangle of the matrix is sufficient
      {
        StoreLowerTriangle( table, data);
      }
      else
      {
        BCL_Exit( "don't know which parts of the matrix are needed!", -1);
      }

      const util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType
        >
      > calculator( new DistancesStored< std::string, t_PrecisionType>( data));

      return calculator;
    }

    //! @brief StoreEntireMatrix puts the entire matrix from a Table into memory
    //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
    //! @param DATA the construct where the data from TABLE will be inserted
    template< typename t_PrecisionType>
    void InputTable< t_PrecisionType>::StoreEntireMatrix
    (
      const storage::Table< t_PrecisionType> &TABLE,
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA,
      const storage::List< std::string>::const_iterator &COLUMN_LAST_OBJECT_MAKER_ITERATOR
    ) const
    {
      // create size_t "current_row" to keep track of what row we are on and initialize with zero
      size_t current_row( 0);

      storage::List< std::string>::const_iterator itr_copy( COLUMN_LAST_OBJECT_MAKER_ITERATOR);

      storage::List< std::string>::const_iterator
        data_row_itr( ++itr_copy), data_row_itr_end( InputInterface< std::string, t_PrecisionType>::m_Objects->End());

      // iterate through the rows of "TABLE"
      for
      (
        typename storage::Table< t_PrecisionType>::const_iterator row_itr( TABLE.Begin()), row_itr_end( TABLE.End());
        row_itr != row_itr_end && data_row_itr != data_row_itr_end;
        ++row_itr, ++current_row, ++data_row_itr
      )
      {
        BCL_MessageStd( "Current row data is -> " + util::Format()( *data_row_itr));

        storage::List< std::string>::const_iterator data_column_itr( InputInterface< std::string, t_PrecisionType>::m_Objects->Begin());
        storage::List< std::string>::const_iterator data_column_itr_end( itr_copy);

        BCL_MessageStd( "Current column data is -> " + util::Format()( *data_column_itr));

        // iterate through the columns of "TABLE"
        for
        (
          // create size_t "start_column" and initialize with "current_row" plus one
          size_t column( 0), end_column( TABLE.GetHeader().GetSize());
          column < end_column && data_column_itr != data_column_itr_end;
          ++column, ++data_column_itr
        )
        {
          DATA[ size_t( &( *data_row_itr))][ size_t( &( *data_column_itr))] = row_itr->Second()( column);
        }
      }
    }

    //! @brief StoreUpperTriangle puts the upper triangle half matrix from a Table into memory
    //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
    //! @param DATA the construct where the data from TABLE will be inserted
    template< typename t_PrecisionType>
    void InputTable< t_PrecisionType>::StoreUpperTriangle
    (
      const storage::Table< t_PrecisionType> &TABLE,
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
    ) const
    {
      // create size_t "current_row" to keep track of what row we are on and initialize with zero
      size_t current_row( 0);

      const storage::TableHeader table_header( TABLE.GetHeader());

      storage::List< std::string>::const_iterator
        data_row_itr( InputInterface< std::string, t_PrecisionType>::m_Objects->Begin()),
        data_row_itr_end( InputInterface< std::string, t_PrecisionType>::m_Objects->End());

      // iterate through the rows of "TABLE"
      for
      (
        typename storage::Table< t_PrecisionType>::const_iterator row_itr( TABLE.Begin()), row_itr_end( TABLE.End());
        row_itr != row_itr_end && data_row_itr != data_row_itr_end;
        ++row_itr, ++current_row, ++data_row_itr
      )
      {
        storage::List< std::string>::const_iterator data_column_itr( data_row_itr);
        ++data_column_itr;
        // iterate through the columns of "TABLE"
        for
        (
          // create size_t "start_column" and initialize with "current_row" plus one
          size_t column( current_row + 1), end_column( TABLE.GetHeader().GetSize());
          column < end_column && data_column_itr != data_row_itr_end;
          ++column, ++data_column_itr
        )
        {
          DATA[ size_t( &( *data_row_itr))][ size_t( &( *data_column_itr))] = row_itr->Second()( column);
        }
      }
    }

    //! @brief StoreLowerTriangle puts the lower triangle half matrix from a Table into memory
    //! @param TABLE storage::Table< t_PrecisionType> which has the data for "m_Data"
    //! @param DATA the construct where the data from TABLE will be inserted
    template< typename t_PrecisionType>
    void InputTable< t_PrecisionType>::StoreLowerTriangle
    (
      const storage::Table< t_PrecisionType> &TABLE,
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
    ) const
    {
      // create size_t "current_row" to keep track of what row we are on and initialize with zero
      size_t current_row( 0);

      storage::List< std::string>::const_iterator
        data_row_itr( InputInterface< std::string, t_PrecisionType>::m_Objects->Begin()),
        data_row_itr_end( InputInterface< std::string, t_PrecisionType>::m_Objects->End());

      // iterate through the rows of "TABLE"
      for
      (
        typename storage::Table< t_PrecisionType>::const_iterator row_itr( TABLE.Begin()), row_itr_end( TABLE.End());
        row_itr != row_itr_end && data_row_itr != data_row_itr_end;
        ++row_itr, ++current_row, ++data_row_itr
      )
      {
        storage::List< std::string>::const_iterator data_column_itr( InputInterface< std::string, t_PrecisionType>::m_Objects->Begin());
        storage::List< std::string>::const_iterator data_column_itr_end( InputInterface< std::string, t_PrecisionType>::m_Objects->End());

        // iterate through the columns of "TABLE"
        for
        (
          // create size_t "start_column" and initialize with "current_row" plus one
          size_t column( 0), end_column( current_row);
          column < end_column && data_column_itr != data_column_itr_end;
          ++column, ++data_column_itr
        )
        {
          DATA[ size_t( &( *data_row_itr))][ size_t( &( *data_column_itr))] = row_itr->Second()( column);
        }
      }
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_INPUT_TABLE_H_
