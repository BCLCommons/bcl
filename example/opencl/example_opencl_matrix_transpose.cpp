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
#include "opencl/bcl_opencl_matrix_transpose.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_matrix_transpose.cpp
  //!
  //! @author loweew
  //! @date Mar 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by kothiwsk on June 11, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclMatrixTranspose :
    public ExampleInterface
  {
  public:

    ExampleOpenclMatrixTranspose *Clone() const
    {
      return new ExampleOpenclMatrixTranspose( *this);
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
      const size_t s_number_elements( 36);
      const size_t rows( 4);
      const size_t cols( 9);

      float data_fill[ s_number_elements] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};
      linal::Matrix< float> data( rows, cols, data_fill);
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      opencl::Matrix< float> data_mat( data, queue);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from size
      opencl::MatrixTranspose< float> transpose( queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( transpose.GetClassIdentifier(), GetStaticClassName< opencl::MatrixTranspose< float> >());

    ///////////////
    // operators //
    ///////////////

      // calculating cpu result for comparison
      linal::Matrix< float> cpu_result( data.Transposed());

      // calculating gpu result and copying back to host
      opencl::Matrix< float> opencl_result
      (
        data_mat.GetNumberCols() - data_mat.GetColPadding(),
        data_mat.GetNumberRows() - data_mat.GetRowPadding(),
        queue,
        data_mat.GetColPadding(),
        data_mat.GetRowPadding()
      );

      transpose( data_mat, opencl_result);
      linal::Matrix< float> read_back_mat( opencl_result.GetHostMatrix());

      // checking that the cpu and gpu agree
      BCL_MessageStd
      (
        "original matrix: " + util::Format()( data) + "\ncpu result matrix: "
        + util::Format()( cpu_result) + "\nopencl result matrix: " + util::Format()( read_back_mat)
      );

      BCL_ExampleCheck( data_mat( 0, 1), cpu_result( 1, 0));
      BCL_ExampleCheck( data_mat( 0, 2), cpu_result( 2, 0));
      BCL_ExampleCheck( data_mat( 0, 1), read_back_mat( 1, 0));
      BCL_ExampleCheck( data_mat( 0, 2), read_back_mat( 2, 0));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclMatrixTranspose

  const ExampleClass::EnumType ExampleOpenclMatrixTranspose::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclMatrixTranspose())
  );

} // namespace bcl
