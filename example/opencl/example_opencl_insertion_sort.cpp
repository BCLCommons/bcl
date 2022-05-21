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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "opencl/bcl_opencl_insertion_sort.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_insertion_sort.cpp
  //!
  //! @author loweew
  //! @date Mar 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclInsertionSort :
    public ExampleInterface
  {
  public:

    ExampleOpenclInsertionSort *Clone() const
    {
      return new ExampleOpenclInsertionSort( *this);
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

    int Run() const
    {
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      // creating data set
      float data[ 15] = { 5, 4, 6, 4, 5, 3, 2, 3, 1, 4, 1, 4, 2, 3, 4};
      linal::Matrix< float> data_to_sort( 3, 5, data);

      // creating command queue
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      // place data into buffer
      opencl::Matrix< float> data_buffer( data_to_sort, queue);
      opencl::Matrix< float> padded_data_buffer( data_to_sort, queue, 3, 5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor from command queue
      opencl::InsertionSort sorter_from_queue( queue);

      // clone
      util::ShPtr< opencl::InsertionSort> sp_sorter( sorter_from_queue.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( sorter_from_queue.GetClassIdentifier(), GetStaticClassName< opencl::InsertionSort>());

    ///////////////
    // operators //
    ///////////////

      // sorting data of non-padded input
      opencl::Matrix< float> sorted_buffer( sorter_from_queue( data_buffer, 1));

      // print unsorted data
      BCL_MessageStd( "unsorted data: " + util::Format()( data_to_sort));

      linal::Matrix< float> sorted_data( sorted_buffer.GetHostMatrix());

      // print k sorted data
      BCL_MessageStd( "sorted data: " + util::Format()( sorted_data));

      // checking sorted result
      BCL_ExampleCheck( sorted_data( 0, 0), data_to_sort( 2, 0));

      // create index buffer
      opencl::Matrix< int> index_buffer( sorter_from_queue.GetIndexMatrix());

      // create index matrix container for index buffer
      linal::Matrix< int> index_matrix( index_buffer.GetHostMatrix());

      // print index matrix
      BCL_MessageStd( "index matrix: " + util::Format()( index_matrix));

      // check index matrix is correct
      BCL_ExampleCheck( index_matrix( 0, 2), 2);

      // sorting padded data
      opencl::Matrix< float> padded_sorted_buffer( sorter_from_queue( padded_data_buffer, 1));

      // unsorted padded data
      linal::Matrix< float> padded_data_to_sort( padded_data_buffer.GetHostMatrix());

      // create padded matrix for sorted data
      linal::Matrix< float> padded_sorted_data( padded_sorted_buffer.GetHostMatrix());

      // print padded unsorted data
      BCL_MessageStd( "padded unsorted data: " + util::Format()( padded_data_to_sort));

      // print padded sorted data
      BCL_MessageStd( "padded sorted data: " + util::Format()( padded_sorted_data));

      // check padded sorted data is correct
      BCL_ExampleCheck( padded_sorted_data( 0, 0), data_to_sort( 2, 0));

      // get index buffer
      opencl::Matrix< int> padded_index_buffer( sorter_from_queue.GetIndexMatrix());

      // create index buffer for result of padded sort which is still just k x unpadded cols
      linal::Matrix< int> padded_index_matrix( padded_index_buffer.GetHostMatrix());

      // print index matrix for padded result which itself is not padded
      BCL_MessageStd( "padded index matrix: " + util::Format()( padded_index_matrix));

      // check index result
      BCL_ExampleCheck( padded_index_matrix( 0, 2), 2);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclInsertionSort

  const ExampleClass::EnumType ExampleOpenclInsertionSort::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclInsertionSort())
  );
  
} // namespace bcl
