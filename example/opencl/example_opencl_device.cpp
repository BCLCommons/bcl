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
#include "opencl/bcl_opencl_device.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_device.cpp
  //! @details this example will exit prematurely if no platform or device is available. reading and writing cannot
  //! be supported with the current opencl standard 1.1
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by woetzen on Nov 6, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclDevice :
    public ExampleInterface
  {
  public:

    ExampleOpenclDevice *Clone() const
    {
      return new ExampleOpenclDevice( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::Device device_empty;
      device_empty.MaxComputeUnits( &error);
      BCL_ExampleIndirectCheck( error == CL_SUCCESS, false, "queries on empty constructed device");

      // reset the error
      error = CL_SUCCESS;
      // indirect test if device can be constructed from cl::Device
      opencl::Platform platform;
      storage::Vector< opencl::Device> command_line_devices;
      opencl::Platform::InitializeFromCommandLine( platform, command_line_devices, &error);
      BCL_ExampleIndirectCheck( error, CL_SUCCESS, "command line devices");

      if( command_line_devices.IsEmpty())
      {
        return 0;
      }

      // clone
      util::ShPtr< opencl::Device> sp_device( command_line_devices.LastElement().Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( sp_device->GetClassIdentifier(), GetStaticClassName< opencl::Device>());

      // cl device
      BCL_ExampleCheck( sp_device->operator()(), command_line_devices.LastElement()());

      // extension set
      BCL_MessageStd
      (
        "nr of devices extensions: " + util::Format()( sp_device->Extensions().GetSize())
      );

      // max compute units
      BCL_ExampleIndirectCheck( sp_device->MaxComputeUnits() > 0, true, "max number of compute units");

      // device type
      const cl_device_type current_device_type( sp_device->DeviceType());
      BCL_MessageStd( opencl::Device::TypeToString( current_device_type));

    ////////////////
    // operations //
    ////////////////

      // device description
      BCL_MessageStd( "device description:\n" + sp_device->GetDescription());

    //////////////////////
    // input and output //
    //////////////////////

      // no reading and writing possible yet, since the devices cannot be identified uniquely

    //////////////////////
    // helper functions //
    //////////////////////

      // type to and from string
      BCL_ExampleIndirectCheck
      (
        opencl::Device::TypeFromString( opencl::Device::TypeToString( current_device_type)),
        current_device_type,
        "device type from and to string"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclDevice

  const ExampleClass::EnumType ExampleOpenclDevice::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclDevice())
  );

} // namespace bcl
