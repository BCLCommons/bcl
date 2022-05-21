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

#ifndef BCL_CLUSTER_INPUT_PAIRWISE_LIST_H_
#define BCL_CLUSTER_INPUT_PAIRWISE_LIST_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_distances_stored.h"
#include "bcl_cluster_input_interface.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_hash_map.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InputPairwiseList
    //! @brief Reads in pairwise distances between objects from a file
    //! @details The file format allows the objects and their distances to be on any column. The columns are specified
    //!          in the first three lines of the file. The file format is :
    //!          <column for first objects>
    //!          <column for second objects>
    //!          <column for distance between objects>
    //!          Then the pairwise distances are provided one per line.
    //!          The columns start numbering with 0, so the first column is column 0, the second column is column 1...
    //!          A simple example of a file would be :
    //!          0
    //!          1
    //!          2
    //!          a b 1.0
    //!          a c 2.0
    //!          b c 0.5
    //!
    //! @see @link example_cluster_input_pairwise_list.cpp @endlink
    //! @author alexanns
    //! @date Mar 18, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_PrecisionType>
    class InputPairwiseList :
      public InputInterface< std::string, t_PrecisionType>
    {

    private:

    //////////
    // data //
    //////////

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
      InputPairwiseList();

      //! @brief Clone function
      //! @return pointer to new InputPairwiseList
      virtual InputPairwiseList *Clone() const;

      //! @brief virtual destructor
      virtual ~InputPairwiseList();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
      //!        in m_Objects
      //! @param IFSTREAM is the stream from which the input will be read
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      virtual util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType
        >
      > HandleInput( io::IFStream &IFSTREAM);

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
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class InputPairwiseList

  //////////
  // data //
  //////////

    //! instantiate s_Instance
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> InputPairwiseList< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new InputPairwiseList< t_PrecisionType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_PrecisionType> InputPairwiseList< t_PrecisionType>::InputPairwiseList() :
      InputInterface< std::string, t_PrecisionType>()
    {
    }

    //! @brief Clone function
    //! @return pointer to new InputPairwiseList
    template< typename t_PrecisionType>
    InputPairwiseList< t_PrecisionType> *InputPairwiseList< t_PrecisionType>::Clone() const
    {
      return new InputPairwiseList< t_PrecisionType>( *this);
    }

    //! @brief virtual destructor
    template< typename t_PrecisionType>
    InputPairwiseList< t_PrecisionType>::~InputPairwiseList< t_PrecisionType>()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_PrecisionType>
    const std::string &InputPairwiseList< t_PrecisionType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
    //         in m_Objects
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
    InputPairwiseList< t_PrecisionType>::HandleInput( io::IFStream &IFSTREAM)
    {
      // set m_Objects to a new storage list
      InputInterface< std::string, t_PrecisionType>::m_Objects = util::ShPtr< storage::List< std::string> >
      (
        new storage::List< std::string>()
      );

      // create data hashmap that will hold the pairwise addresses and distances of the objects
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > data;

      // create map that will hold the unique objects found in the list of objects
      storage::Map< std::string, size_t> unique_objects_address;

      // get the column where the first object in the pair can be found
      int column_object_a;
      io::Serialize::Read( column_object_a, IFSTREAM);

      // get the column where the second object in the pair can be found
      int column_object_b;
      io::Serialize::Read( column_object_b, IFSTREAM);

      // get the column where the distance between the pair of objects can be found
      int column_data;
      io::Serialize::Read( column_data, IFSTREAM);

      // message the columns where information can be found
      BCL_MessageDbg
      (
        "column_object_a " + util::Format()( column_object_a) + " column_object_b "
        + util::Format()( column_object_b) + " column_data " + util::Format()( column_data)
      );

      // true while not at the end of the pairwise distance list file
      while( !IFSTREAM.eof())
      {
        // get the current line and trim it
        std::string current_line;
        std::getline( IFSTREAM, current_line);
        current_line = util::TrimString( current_line);

        // message current line
        BCL_MessageDbg( "current line |" + current_line + "|");

        // if the current line is empty then continue to next iteration of while loop
        if( current_line.empty())
        {
          continue;
        }

        // split the current line
        const storage::Vector< std::string> split_line( util::SplitString( current_line));

        // get the two objects from the current line
        const std::string object_a( split_line( column_object_a));
        const std::string object_b( split_line( column_object_b));

        // get the distance between the two objects from the appropriate column
        const t_PrecisionType data_value( util::ConvertStringToNumericalValue< t_PrecisionType>( split_line( column_data)));

        // message the two objects and their distance
        BCL_MessageDbg
        (
          "object_a " + util::Format()( object_a) + " object_b " + util::Format()( object_b) +
          " distance " + util::Format()( data_value)
        );

        // try to find the two objects in the map of unique objects
        storage::Map< std::string, size_t>::const_iterator itr_obj_a
        (
          unique_objects_address.Find( object_a)
        );
        storage::Map< std::string, size_t>::const_iterator itr_obj_b
        (
          unique_objects_address.Find( object_b)
        );

        // size_t will hold the address of object a
        size_t address_obj_a;

        // true if object a was not found in "unique_objects_address"
        if( itr_obj_a == unique_objects_address.End())
        {
          // add object a to "m_Objects"
          InputInterface< std::string, t_PrecisionType>::m_Objects->PushBack( object_a);

          // get the address of object_a as it is in "m_Objects"
          address_obj_a = size_t( &( InputInterface< std::string, t_PrecisionType>::m_Objects->LastElement()));

          // add object_a and its address in "m_Objects" into "unique_objects_address"
          unique_objects_address.Insert( std::pair< std::string, size_t>( object_a, address_obj_a));
        }
        else //< object_a already exists in "unique_objects_address"
        {
          // assign the address of object_a to "address_obj_a"
          address_obj_a = itr_obj_a->second;
        }

        // size_t will hold the address of object b
        size_t address_obj_b;

        // true if object b was not found in "unique_objects_address"
        if( itr_obj_b == unique_objects_address.End())
        {
          // add object b to "m_Objects"
          InputInterface< std::string, t_PrecisionType>::m_Objects->PushBack( object_b);

          // get the address of object_b as it is in "m_Objects"
          address_obj_b = size_t( &( InputInterface< std::string, t_PrecisionType>::m_Objects->LastElement()));

          // add object_b and its address in "m_Objects" into "unique_objects_address"
          unique_objects_address.Insert( std::pair< std::string, size_t>( object_b, address_obj_b));
        }
        else //< object_b already exists in "unique_objects_address"
        {
          // assign the address of object_b to "address_obj_b"
          address_obj_b = itr_obj_b->second;
        }

        // add the pair of object a and b addresses and their distance into "data"
        data[ address_obj_a][ address_obj_b] = data_value;

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
    std::istream &InputPairwiseList< t_PrecisionType>::Read( std::istream &ISTREAM)
    {
      // read member variables
      InputInterface< std::string, t_PrecisionType>::Read( ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_PrecisionType>
    std::ostream &InputPairwiseList< t_PrecisionType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member variables
      InputInterface< std::string, t_PrecisionType>::Write( OSTREAM, INDENT) << '\n';

      return OSTREAM;
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_INPUT_PAIRWISE_LIST_H_
