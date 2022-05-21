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
#include "opencl/bcl_opencl_rmsd.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_rmsd.cpp
  //!
  //! @author loweew
  //! @date Mar 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclRMSD :
    public ExampleInterface
  {
  public:

    ExampleOpenclRMSD *Clone() const
    {
      return new ExampleOpenclRMSD( *this);
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
      const size_t s_elements( 5);
      float data_a[ s_elements] = { 5, 4, 6, 4, 5};
      float data_b[ s_elements] = { 5, 8, 6, 4, 5};
      linal::Matrix< float> data_vec_a( 1, s_elements, data_a);
      linal::Matrix< float> data_vec_b( 1, s_elements, data_b);

      // creating command queue
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      // place data into buffer
      opencl::Matrix< float> data_buffer_a( data_vec_a, queue);
      opencl::Matrix< float> padded_data_buffer_a( data_vec_a, queue, 0, 5);
      opencl::Matrix< float> data_buffer_b( data_vec_b, queue);
      opencl::Matrix< float> padded_data_buffer_b( data_vec_b, queue, 0, 5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor from command queue
      opencl::RMSD rmsd( queue);

      // clone
      util::ShPtr< opencl::RMSD> sp_rmsd( rmsd.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( rmsd.GetClassIdentifier(), GetStaticClassName< opencl::RMSD>());

    ///////////////
    // operators //
    ///////////////

      // checking operator
      float rmsd_non_padded_result( rmsd( data_buffer_a, data_buffer_b));
      float rmsd_padded_result( rmsd( padded_data_buffer_a, padded_data_buffer_b));

      // calculating gold value on cpu
      float sum( 0), tmp( 0);
      for( size_t count( 0); count < s_elements; ++count)
      {
        tmp = data_vec_a( 0, count) - data_vec_b( 0, count);
        sum += tmp * tmp;
      }
      float cpu_rmsd( math::Sqrt( sum / ( s_elements)));

      BCL_MessageStd( "cpu rmsd calculation give: " + util::Format()( cpu_rmsd)
        + "\nopencl calculation gives: " + util::Format()( rmsd_non_padded_result));

      // checking non-padded rmsd
      BCL_ExampleCheck( rmsd_non_padded_result, cpu_rmsd);
      BCL_ExampleCheck( rmsd_padded_result, cpu_rmsd);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclRMSD

  const ExampleClass::EnumType ExampleOpenclRMSD::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclRMSD())
  );

} // namespace bcl
