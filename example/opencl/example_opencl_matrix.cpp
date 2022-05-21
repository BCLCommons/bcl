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
#include "opencl/bcl_opencl_matrix.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_matrix.cpp
  //!
  //! @author loweew
  //! @date Mar 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclMatrix :
    public ExampleInterface
  {
  public:

    ExampleOpenclMatrix *Clone() const
    {
      return new ExampleOpenclMatrix( *this);
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

      // create some data
      const size_t s_number_elements( 6);
      const size_t rows( 2);
      const size_t cols( 3);
      const size_t padding( 2);

      float data_fill[ s_number_elements] = { 1, 2, 3, 4, 5, 6};
      linal::Matrix< float> data( rows, cols, data_fill);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      // construct from size
      opencl::Matrix< float> empty_buffer_matrix( rows, cols, queue, 0, 0, 0);

      // construct from matrix
      opencl::Matrix< float> buffer_matrix_from_matrix( data, queue);

      // construct from matrix and add padding
      opencl::Matrix< float> buffer_from_mat_with_padding( data, queue, padding, padding);

      // construct from pointer
      opencl::Matrix< float> buffer_mat_from_pointer( rows, cols, data_fill, queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( empty_buffer_matrix.GetClassIdentifier(), GetStaticClassName< opencl::Matrix< float> >());

      const size_t get_row( buffer_matrix_from_matrix.GetNumberRows());
      const size_t get_col( buffer_matrix_from_matrix.GetNumberCols());
      const size_t get_row_padding( buffer_from_mat_with_padding.GetRowPadding());
      const size_t get_col_padding( buffer_from_mat_with_padding.GetColPadding());
      const size_t get_nr_elements( buffer_matrix_from_matrix.GetNumberOfElements());

      // checking get functions
      BCL_ExampleCheck( get_row, rows);
      BCL_ExampleCheck( get_col, cols);
      BCL_ExampleCheck( get_row_padding, padding);
      BCL_ExampleCheck( get_col_padding, padding);
      BCL_ExampleCheck( get_nr_elements, s_number_elements);

      linal::Matrix< float> read_back_mat( rows, cols);
      queue.enqueueReadBuffer( buffer_matrix_from_matrix.GetData(), CL_TRUE, 0, sizeof( float) * buffer_matrix_from_matrix.GetNumberOfElements(), read_back_mat.Begin());

      // checking that the same data is returned from the gpu that was sent to it
      BCL_MessageStd( "original matrix: " + util::Format()( data) + "\nafter host <-> gpu <-> host: " + util::Format()( read_back_mat));
      BCL_ExampleCheck( data( 0, 0), read_back_mat( 0, 0));
      BCL_ExampleCheck( data( 1, 1), read_back_mat( 1, 1));
      BCL_ExampleCheck( data( 1, 1), buffer_matrix_from_matrix( 1, 1));
      BCL_ExampleCheck( float( 5), buffer_matrix_from_matrix( 1, 1));
      buffer_matrix_from_matrix.SetValue( 1, 1, float( 3));
      BCL_ExampleCheck( float( 3), buffer_matrix_from_matrix( 1, 1));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclMatrix

  const ExampleClass::EnumType ExampleOpenclMatrix::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclMatrix())
  );
  
} // namespace bcl
