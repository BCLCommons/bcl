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
#include "opencl/bcl_opencl_matrix_add.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_matrix_add.cpp
  //!
  //! @author loweew
  //! @date Mar 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by kothiwsk on June 11, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclMatrixAdd :
    public ExampleInterface
  {
  public:

    ExampleOpenclMatrixAdd *Clone() const
    {
      return new ExampleOpenclMatrixAdd( *this);
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

      float data_fill_a[ s_number_elements] = { 1, 2, 3, 4, 5, 6};
      float data_fill_b[ s_number_elements] = { 6, 5, 4, 3, 2, 1};
      linal::Matrix< float> data_a( rows, cols, data_fill_a);
      linal::Matrix< float> data_b( rows, cols, data_fill_b);

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());
      opencl::Matrix< float> data_mat_a( data_a, queue);
      opencl::Matrix< float> data_mat_b( data_b, queue);
      opencl::Matrix< float> large_mat_a( 3859, 27, queue, 0, 0, 34);
      opencl::Matrix< float> large_mat_b( 3859, 27, queue, 0, 0, 15);

      opencl::Matrix< float> data_mat_a_padded( data_a, queue, 5, 5);
      opencl::Matrix< float> data_mat_b_padded( data_b, queue, 5, 5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from size
      opencl::MatrixAdd< float> add( queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( add.GetClassIdentifier(), GetStaticClassName< opencl::MatrixAdd< float> >());

    ///////////////
    // operators //
    ///////////////

      // calculating cpu result for comparison
      linal::Matrix< float> cpu_result( data_a + data_b);

      // calculating gpu result and copying back to host
      opencl::Matrix< float> opencl_result
      (
        data_mat_a.GetNumberRows() - data_mat_a.GetRowPadding(),
        data_mat_b.GetNumberCols() - data_mat_b.GetColPadding(),
        queue,
        data_mat_a.GetRowPadding(),
        data_mat_b.GetColPadding()
      );

      // operator
      add( data_mat_a, data_mat_b);
      opencl_result = data_mat_a;

      opencl::Matrix< float> opencl_padded_result
      (
        data_mat_a_padded.GetNumberRows() - data_mat_a_padded.GetRowPadding(),
        data_mat_b_padded.GetNumberCols() - data_mat_b_padded.GetColPadding(),
        queue,
        data_mat_a_padded.GetRowPadding(),
        data_mat_b_padded.GetColPadding()
      );

      add( data_mat_a_padded, data_mat_b_padded);
      opencl_padded_result = data_mat_a_padded;
      linal::Matrix< float> read_back_mat( opencl_result.GetHostMatrix());

      // checking that the cpu and gpu agree
      BCL_MessageStd
      (
        "cpu result matrix: " + util::Format()( cpu_result)
        + "\nopencl result matrix: " + util::Format()( read_back_mat)
      );

      BCL_MessageStd
      (
        "opencl padded result matrix: "
        + util::Format()( opencl_padded_result.GetHostMatrix())
      );

      // preallocate result
      opencl::Matrix< float> large_results_ocl
      (
        large_mat_a.GetNumberRows() - large_mat_a.GetRowPadding(),
        large_mat_b.GetNumberCols() - large_mat_b.GetColPadding(),
        queue,
        large_mat_a.GetRowPadding(),
        large_mat_b.GetColPadding()
      );

      // operator
      add( large_mat_a, large_mat_b);
      large_results_ocl = large_mat_a;
      linal::Matrix< float> large_results( large_results_ocl.GetHostMatrix());

      // example checks
      BCL_ExampleCheck( cpu_result( 0, 0), read_back_mat( 0, 0));
      BCL_ExampleCheck( cpu_result( 1, 1), read_back_mat( 1, 1));
      BCL_ExampleCheck( float( 49), large_results( 0, 0));
      BCL_ExampleCheck( float( 49), large_results( 33, 21));
      BCL_ExampleCheck( float( 49), large_results( 2338, 18));
      BCL_ExampleCheck( float( 49), large_results( 3323, 4));
      BCL_ExampleCheck( float( 49), large_results( 3858, 26));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclMatrixAdd

  const ExampleClass::EnumType ExampleOpenclMatrixAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclMatrixAdd())
  );

} // namespace bcl
