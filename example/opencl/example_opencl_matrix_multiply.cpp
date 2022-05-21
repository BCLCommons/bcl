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
#include "opencl/bcl_opencl_matrix_multiply.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_matrix_multiply.cpp
  //!
  //! @author loweew
  //! @date Mar 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclMatrixMultiply :
    public ExampleInterface
  {
  public:

    ExampleOpenclMatrixMultiply *Clone() const
    {
      return new ExampleOpenclMatrixMultiply( *this);
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
      const size_t rows_a    ( 2);
      const size_t cols_a    ( 3);
      const size_t rows_b    ( 3);
      const size_t cols_b    ( 2);

      float data_fill_a[ s_number_elements] = { 1, 2, 3, 4, 5, 6};
      float data_fill_b[ s_number_elements] = { 6, 5, 4, 3, 2, 1};
      linal::Matrix< float> data_a( rows_a, cols_a, data_fill_a);
      linal::Matrix< float> data_b( rows_b, cols_b, data_fill_b);

      // calculating row and col padding
      const size_t block_size( 16);
      const size_t rows_a_pad( ( block_size - ( rows_a % block_size)) % block_size);
      const size_t cols_a_pad( ( block_size - ( cols_a % block_size)) % block_size);
      const size_t rows_b_pad( ( block_size - ( rows_b % block_size)) % block_size);
      const size_t cols_b_pad( ( block_size - ( cols_b % block_size)) % block_size);

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      opencl::Matrix< float> data_a_mat( data_a, queue, rows_a_pad, cols_a_pad);
      opencl::Matrix< float> data_b_mat( data_b, queue, rows_b_pad, cols_b_pad);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from size
      opencl::MatrixMultiply< float> mmult( queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( mmult.GetClassIdentifier(), GetStaticClassName< opencl::MatrixMultiply< float> >());

    ///////////////
    // operators //
    ///////////////

      // calculating cpu result for comparison
      linal::Matrix< float> cpu_result( data_a * data_b);

      // calculating gpu result and copying back to host
      opencl::Matrix< float> opencl_result( data_a_mat.GetNumberRows() - data_a_mat.GetRowPadding(), data_b_mat.GetNumberCols() - data_b_mat.GetColPadding(), queue, data_a_mat.GetRowPadding(), data_b_mat.GetColPadding());
      mmult( data_a_mat, data_b_mat, opencl_result);
      linal::Matrix< float> read_back_mat( opencl_result.GetHostMatrix( opencl_result.GetRowPadding(), opencl_result.GetColPadding()));

      // checking that the cpu and gpu agree
      BCL_MessageStd( "cpu result matrix: " + util::Format()( cpu_result) + "\nopencl result matrix: " + util::Format()( read_back_mat));
      BCL_ExampleCheck( cpu_result( 0, 0), read_back_mat( 0, 0));
      BCL_ExampleCheck( cpu_result( 1, 1), read_back_mat( 1, 1));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclMatrixMultiply

  const ExampleClass::EnumType ExampleOpenclMatrixMultiply::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclMatrixMultiply())
  );

} // namespace bcl
