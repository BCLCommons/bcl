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
#include "opencl/bcl_opencl_platform.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_platform.cpp
  //! @details if no platform is available it returns prematurely without testing
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by loweew on Nov 6, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclPlatform :
    public ExampleInterface
  {
  public:

    ExampleOpenclPlatform *Clone() const
    {
      return new ExampleOpenclPlatform( *this);
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

      cl_int error( CL_SUCCESS);

      // standardized platform names
      const storage::Vector< std::string> &platform_names( opencl::Platform::QueryPlatformNamesStandardized( &error));
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "QueryPlatformNamesStandardized");

      // query all available platforms
      opencl::Platform::QueryPlatforms( &error);
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "QueryPlatforms");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::Platform platform_empty;

      // construct from name
      error = CL_SUCCESS;
      storage::Vector< opencl::Platform> opencl_platforms_from_name;
      for
      (
        storage::Vector< std::string>::const_iterator itr( platform_names.Begin()), itr_end( platform_names.End());
        itr != itr_end && error == CL_SUCCESS;
        ++itr
      )
      {
        cl_int current_error( CL_SUCCESS);
        opencl_platforms_from_name.PushBack( opencl::Platform( *itr, &current_error));
        error &= current_error;
      }
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "construct from name");
      BCL_ExampleIndirectCheck( platform_names.GetSize(), opencl_platforms_from_name.GetSize(), "construct from name");

      // query all platforms
      error = CL_SUCCESS;
      storage::Vector< opencl::Platform> opencl_platforms( opencl::Platform::QueryPlatforms( &error));
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "QueryPlatforms");
      BCL_ExampleIndirectCheck( platform_names.GetSize(), opencl_platforms.GetSize(), "QueryPlatforms");

      if( opencl_platforms.IsEmpty())
      {
        return 0;
      }

      // copy constructor
      util::ShPtr< opencl::Platform> sp_platform( opencl_platforms.FirstElement().Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( sp_platform->GetClassIdentifier(), GetStaticClassName< opencl::Platform>());

    ////////////////
    // operations //
    ////////////////

      // standardized name
      BCL_ExampleCheck( sp_platform->StandardizedName(), platform_names.FirstElement());
      BCL_MessageStd( "non standardized platform name: " + sp_platform->Name());

      // get description
      BCL_MessageStd( "description for first platform:\n" + sp_platform->Description());

      // devices for that platform
      storage::Vector< opencl::Device> platforms_devices( sp_platform->Devices( CL_DEVICE_TYPE_ALL, &error));
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "GetDevices");
      BCL_MessageStd
      (
        "platform " + sp_platform->StandardizedName() + " has " + util::Format()( platforms_devices.GetSize()) +
        " devices"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // read and write platform out to a string stream
      // machines with different gpus give different results, so do not put this file in the svn
      std::stringstream input_output;
      input_output << opencl_platforms.FirstElement();
      opencl::Platform read_platform;
      input_output >> read_platform;

      // check if platforms are the same
      BCL_ExampleIndirectCheck
      (
        read_platform(), // convert to cl_type_id through ()
        opencl_platforms.FirstElement()(), // convert to cl_type_id through ()
        "read write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // indirect test if device can be constructed from cl::Device
      opencl::Platform platform;
      storage::Vector< opencl::Device> command_line_devices;
      opencl::Platform::InitializeFromCommandLine( platform, command_line_devices, &error);

      // commandline platform
      BCL_ExampleIndirectCheck( platform.Name().empty(), false, "InitializeFromCommandLine");

      // commandline devices
      BCL_ExampleIndirectCheck
      (
        command_line_devices.GetSize() <= platform.Devices( CL_DEVICE_TYPE_ALL).GetSize(),
        true,
        "InitializeFromCommandLine"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclPlatform

  const ExampleClass::EnumType ExampleOpenclPlatform::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclPlatform())
  );

} // namespace bcl
