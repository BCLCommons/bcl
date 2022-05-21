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
#include "opencl/bcl_opencl_arg_min.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_arg_min.cpp
  //!
  //! @author loweew
  //! @date Mar 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclArgMin :
    public ExampleInterface
  {
  public:

    ExampleOpenclArgMin *Clone() const
    {
      return new ExampleOpenclArgMin( *this);
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
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

      float data_fill[ s_number_elements] = { 1, 2, 3, 4, 5, 6};
      linal::Vector< float> data( s_number_elements, data_fill);
      opencl::Vector< float> vec( data, queue);

      linal::Vector< float> big_vec( 536, float( 3));
      big_vec( 346) = 0;
      opencl::Vector< float> ocl_big_vec( big_vec, queue);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from size
      opencl::ArgMin< float> arg_min( queue);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( arg_min.GetClassIdentifier(), GetStaticClassName< opencl::ArgMin< float> >());

    ///////////////
    // operators //
    ///////////////

      storage::Pair< float, int> result( arg_min( vec));
      storage::Pair< float, int> big_result( arg_min( ocl_big_vec));

      BCL_MessageStd( "result: " + util::Format()( result));
      BCL_MessageStd( "big result: " + util::Format()( big_result));

      BCL_ExampleCheck( result.First(), data( 0));
      BCL_ExampleCheck( result.Second(), 0);
      BCL_ExampleCheck( big_result.First(), float( 0));
      BCL_ExampleCheck( big_result.Second(), 346);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclArgMin

  const ExampleClass::EnumType ExampleOpenclArgMin::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclArgMin())
  );
  
} // namespace bcl
