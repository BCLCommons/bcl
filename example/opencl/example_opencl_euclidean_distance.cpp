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
#include "opencl/bcl_opencl_euclidean_distance.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_euclidean_distance.cpp
  //!
  //! @author loweew
  //! @date Mar 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclEuclideanDistance :
    public ExampleInterface
  {
  public:

    ExampleOpenclEuclideanDistance *Clone() const
    {
      return new ExampleOpenclEuclideanDistance( *this);
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
      float data_fill_a[ 15] = { 5, 4, 6, 4, 5, 3, 2, 3, 1, 4, 1, 4, 2, 3, 4};
      linal::Matrix< float> data_a( 3, 5, data_fill_a);
      linal::Matrix< float> data_b( 3, 5, 5);

      // creating command queue
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      // place data into buffer
      opencl::Matrix< float> padded_data_buffer_a( data_a, queue, 13, 11);
      opencl::Matrix< float> padded_data_buffer_b( data_b, queue, 13, 11);

      opencl::Matrix< float> big_a( 160, 160, queue, 0, 0, 5);
      opencl::Matrix< float> big_b( 160, 160, queue, 0, 0, 9);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor from command queue
      opencl::EuclideanDistance< float> euclidean_distance( queue);

      // clone
      util::ShPtr< opencl::EuclideanDistance< float> > sp_euclidean_distance( euclidean_distance.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( euclidean_distance.GetClassIdentifier(), GetStaticClassName< opencl::EuclideanDistance< float> >());

    ///////////////
    // operators //
    ///////////////

      // checking operator
      opencl::Matrix< float> padded_dist_mat_buffer( euclidean_distance( padded_data_buffer_a, padded_data_buffer_b));

      // calculating gold value on cpu
      linal::Matrix< float> cpu_euclidean_distance( data_a.GetNumberRows(), data_b.GetNumberRows());

      float tmp( 0);
      // calculating distance matrix on cpu for comparison
      for( size_t row_a( 0), nr_rows_a( data_a.GetNumberRows()); row_a < nr_rows_a; ++row_a)
      {
        for( size_t row_b( 0), nr_rows_b( data_b.GetNumberRows()); row_b < nr_rows_b; ++row_b)
        {
          tmp = 0;
          for( size_t col( 0), nr_cols( data_a.GetNumberCols()); col < nr_cols; ++col)
          {
            tmp += math::Sqr( data_a( row_a, col) - data_b( row_b, col));
          }
          cpu_euclidean_distance( row_a, row_b) = math::Sqrt( tmp);
        }
      }

      linal::Matrix< float> euc_dist_padded( padded_dist_mat_buffer.GetHostMatrix( padded_dist_mat_buffer.GetRowPadding(), padded_dist_mat_buffer.GetColPadding()));

      BCL_MessageStd
      (
        "cpu euclidean_distance calculation give: " + util::Format()( cpu_euclidean_distance)
        + "\nopencl on padded matrix gives: " + util::Format()( euc_dist_padded)
      );

      // checking non-padded euclidean_distance
      BCL_ExampleCheckWithinTolerance( double( float( euc_dist_padded( 0, 2))), double( float( cpu_euclidean_distance( 0, 2))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( float( euc_dist_padded( 1, 2))), double( float( cpu_euclidean_distance( 1, 2))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( float( euc_dist_padded( 2, 2))), double( float( cpu_euclidean_distance( 2, 2))), 0.001);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclEuclideanDistance

  const ExampleClass::EnumType ExampleOpenclEuclideanDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclEuclideanDistance())
  );

} // namespace bcl
