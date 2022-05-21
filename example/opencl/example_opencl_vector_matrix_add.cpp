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
#include "opencl/bcl_opencl_vector_matrix_add.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_vector_matrix_add.cpp
  //!
  //! @author loweew
  //! @date Apr 01, 2011
  //! @remarks status complete
  //! @remarks reviewed by kothiwsk on June 11, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclVectorMatrixAdd :
    public ExampleInterface
  {
  public:

    ExampleOpenclVectorMatrixAdd *Clone() const
    {
      return new ExampleOpenclVectorMatrixAdd( *this);
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
      float data_fill_b[ 3] = { 6, 5, 4};
      linal::Matrix< float> data_a( rows, cols, data_fill_a);
      linal::Vector< float> data_b( cols, data_fill_b);
      linal::Matrix< float> large_mat( 143, 3, 3);
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      opencl::Matrix< float> data_mat_a( data_a, queue);
      opencl::Vector< float> data_vec_b( data_b, queue);
      opencl::Matrix< float> large_ocl_mat( large_mat, queue);

      opencl::Matrix< float> data_mat_a_padded( data_a, queue, 5, 5);
      opencl::Vector< float> data_vec_b_padded( data_b, queue, 5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from size
      opencl::VectorMatrixAdd< float> add( queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( add.GetClassIdentifier(), GetStaticClassName< opencl::VectorMatrixAdd< float> >( queue));

    ///////////////
    // operators //
    ///////////////

      // calculating cpu result for comparison
      linal::Matrix< float> cpu_result( data_a);

      for( size_t row( 0), row_end( data_a.GetNumberRows()); row < row_end; ++row)
      {
        for( size_t col( 0), col_end( data_a.GetNumberCols()); col < col_end; ++col)
        {
          cpu_result( row, col) += data_b( col);
        }
      }

      // calculating gpu result and copying back to host
      add( data_vec_b, data_mat_a);
      add( data_vec_b_padded, data_mat_a_padded);
      linal::Matrix< float> read_back_mat( data_mat_a.GetHostMatrix());

      // checking that the cpu and gpu agree
      BCL_MessageStd
      (
        "cpu result matrix: "
        + util::Format()( cpu_result) + "\nopencl result matrix: " + util::Format()( read_back_mat)
      );
      BCL_MessageStd
      (
        "opencl padded result matrix: "
        + util::Format()( data_mat_a.GetHostMatrix())
      );

      // checking by example check
      add( data_vec_b, large_ocl_mat);
      linal::Matrix< float> large_result( large_ocl_mat.GetHostMatrix());
      BCL_ExampleCheck( cpu_result( 0, 0), read_back_mat( 0, 0));
      BCL_ExampleCheck( cpu_result( 1, 1), read_back_mat( 1, 1));
      BCL_ExampleCheck( float( 9), large_result( 29, 0));
      BCL_ExampleCheck( float( 8), large_result( 102, 1));
      BCL_ExampleCheck( float( 7), large_result( 78, 2));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclVectorMatrixAdd

  const ExampleClass::EnumType ExampleOpenclVectorMatrixAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclVectorMatrixAdd())
  );
  
} // namespace bcl
