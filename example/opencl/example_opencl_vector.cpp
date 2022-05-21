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
#include "opencl/bcl_opencl_vector.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_vector.cpp
  //!
  //! @author loweew
  //! @date Mar 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclVector :
    public ExampleInterface
  {
  public:

    ExampleOpenclVector *Clone() const
    {
      return new ExampleOpenclVector( *this);
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
      const size_t padding( 2);

      float data_fill[ s_number_elements] = { 1, 2, 3, 4, 5, 6};
      linal::Vector< float> data( s_number_elements, data_fill);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      // construct from size
      opencl::Vector< float> empty_buffer_vector( s_number_elements, queue);

      // construct from vector
      opencl::Vector< float> buffer_vector_from_vector( data, queue);

      // construct from vector and add padding
      opencl::Vector< float> buffer_from_mat_with_padding( data, queue, padding);

      // construct from pointer
      opencl::Vector< float> buffer_mat_from_pointer( s_number_elements, data_fill, queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( empty_buffer_vector.GetClassIdentifier(), GetStaticClassName< opencl::Vector< float> >());

      const size_t get_size( buffer_vector_from_vector.GetSize());
      const size_t get_padding( buffer_from_mat_with_padding.GetPadding());

      // checking get functions
      BCL_ExampleCheck( get_size, s_number_elements);
      BCL_ExampleCheck( get_padding, padding);

      // read back data from gpu to make sure it's the same
      linal::Vector< float> read_back_vec( s_number_elements);
      queue.enqueueReadBuffer( buffer_vector_from_vector.GetData(), CL_TRUE, 0, sizeof( float) * buffer_vector_from_vector.GetSize(), read_back_vec.Begin());

      // checking that the same data is returned from the gpu that was sent to it
      BCL_MessageStd( "original vector: " + util::Format()( data) + "\nafter host <-> gpu <-> host: " + util::Format()( read_back_vec));
      BCL_ExampleCheck( data( 1), read_back_vec( 1));

      BCL_ExampleCheck( data( 2), buffer_vector_from_vector( 2));
      buffer_vector_from_vector.SetValue( 2, float( 5));
      BCL_ExampleCheck( float( 5), buffer_vector_from_vector( 2));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclVector

  const ExampleClass::EnumType ExampleOpenclVector::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclVector())
  );
  
} // namespace bcl
