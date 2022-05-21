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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opencl/bcl_opencl_platform.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////
  // data //
  //////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlatformFlagHelper
    //! @brief Adds Opencl flags to default flags, but only if -opencl Disable was not given on the command line
    //!        This is necessary to avoid querying the opencl platform if opencl was disabled, because even querying the
    //!        platform causes the cuda runtime to initialize, which on cluster jobs results in the processes' virtual
    //!        memory appearing to be equal to the total (phys + virtual + GPU) memory available on the given node due
    //!        to cuda's universal virtual addressing system. This is bad because it can cause pbs cluster jobs to crash
    //!        because they don't (or can't) request all virtual memory on the given node
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 09, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BCL_API PlatformFlagHelper :
      public signal::Slots
    {
      static const PlatformFlagHelper s_Instance;

      PlatformFlagHelper()
      {
        command::CommandState::GetGlobalCommandState().GetParseArgumentsSignal().Connect
        (
          this,
          &PlatformFlagHelper::AddOpenclFlagsToDefaultFlags
        );
      }

      //! @brief function to add opencl flags in a well defined order to the app default flags enum, but only if opencl
      //!        was not explicitly disabled
      void AddOpenclFlagsToDefaultFlags( const command::CommandState &STATE)
      {
        const storage::Vector< std::string> &opencl_args( STATE.GetArguments( "opencl"));
        if
        (
          ( opencl_args.GetSize() == size_t( 1) && opencl_args( 0) == "Disable")
          || ( opencl_args.IsEmpty() && !STATE.GetState().Has( "opencl"))
        )
        {
          Platform::GetIsOpenclDisabled() = true;
        }
        if( Platform::GetPlatformFlag().IsDefined())
        {
          // add opencl flags
          command::GetAppDefaultFlags().AddDefaultFlag( Platform::GetPlatformFlag(), command::e_Opencl);
          // add these flags only if opencl functionality is not prohibited due to incomplete opencl installation
          if( Platform::GetPlatformFlag()->GetParameterList().GetSize() > size_t( 1))
          {
            command::GetAppDefaultFlags().AddDefaultFlag( KernelSources::GetKernelsSourcePathFlag(), command::e_Opencl);
            command::GetAppDefaultFlags().AddDefaultFlag( KernelSources::GetKernelsBinaryPathFlag(), command::e_Opencl);

            // Ensure that the last flag has the signal to update the compiled opencl kernels and enums
            BCL_Assert
            (
              command::GetAppDefaultFlags().GetDefaultFlagsOfType( command::e_Opencl).LastElement()->GetSignal()
              == &Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag,
              "The last opencl flag must have the update signal on it"
            );
          }
        }
      }
    };

    //! @brief add the opencl flags to the default flags
    const PlatformFlagHelper PlatformFlagHelper::s_Instance;

    //! @brief access to the commandline flag
    //! @return flag to select platform and processor type
    const util::ShPtr< command::FlagInterface> &Platform::GetPlatformFlag()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "opencl",
          "choice of opencl platform and device type"
        )
      );
      static bool s_attempted_initialization( false); // keep track of whether we have attempted initialization

      if( !s_attempted_initialization)
      {
        // cast to FlagStatic for pushback
        util::ShPtr< command::FlagStatic> flag( s_flag);
        s_attempted_initialization = true;
        cl_int error( CL_SUCCESS);
        storage::Vector< std::string> names
        (
          GetIsOpenclDisabled()
          ? storage::Vector< std::string>()
          : QueryPlatformNamesStandardized( &error)
        );
        // warn if platform names could not be retrieved
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to get opencl platform names: " + Tools::ErrorString( error));
        }
        names.PushBack( "Disable");
        if( GetIsOpenclDisabled() || error != CL_SUCCESS)
        {
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; Disabled because this machine lacks libOpenCL or does not have the appropriate "
                "/etc/OpenCL/vendors .icd files, or Disable was already given",
                command::ParameterCheckAllowed( names),
                names.LastElement() // Disable
              )
            )
          );
          return s_flag;
        }
        cl_device_type type;
        Platform optimal;
        optimal = GetOptimalPlatform( type, &error);
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to get optimal platform: " + Tools::ErrorString( error));
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; Disabled because this machine lacks libOpenCL or does not have the appropriate "
                "/etc/OpenCL/vendors .icd files",
                command::ParameterCheckAllowed( names),
                names.LastElement() // Disable
              )
            )
          );
        }
        else
        {
          BCL_MessageVrb( "optimal platform: " + optimal.Name() + " " + Device::TypeToString( type));
        }

        // only insert parameters if there is at least one platform with either CPU or GPU
        if( error == CL_SUCCESS)
        {
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; select Disable on older machines with little or buggy opencl support",
                command::ParameterCheckAllowed( names),
                StandardizeName( optimal.Name())
              )
            )
          );
          flag->PushBack( GetDeviceTypeParam( &type));
          flag->PushBack( GetDeviceIDsParam());
        }
      }

      // end
      return s_flag;
    }

    //! @brief access to the command line param for processor type
    //! @param DEVICE_TYPE_PTR the device type to be default, if NULL, use the GPU or last setting
    //! @return param to select processor type
    util::ShPtr< command::ParameterInterface> &Platform::GetDeviceTypeParam( const cl_device_type *DEVICE_TYPE_PTR)
    {
      static util::ShPtr< command::ParameterInterface> s_param;

      // first time
      if( !s_param.IsDefined())
      {
        s_param = util::ShPtr< command::ParameterInterface>
        (
          new command::Parameter
          (
            "device_type",
            "choice of device type",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                Device::TypeToString( CL_DEVICE_TYPE_GPU),
                Device::TypeToString( CL_DEVICE_TYPE_CPU),
                Device::TypeToString( CL_DEVICE_TYPE_ACCELERATOR),
                Device::TypeToString( CL_DEVICE_TYPE_ALL),
                Device::TypeToString( CL_DEVICE_TYPE_DEFAULT)
              )
            ),
            Device::TypeToString( CL_DEVICE_TYPE_GPU)
          )
        );
      }

      if( DEVICE_TYPE_PTR != NULL)
      {
        s_param->SetDefaultParameter( Device::TypeToString( *DEVICE_TYPE_PTR));
      }

      // end
      return s_param;
    }

    //! @brief access to the command line param for device vendor id's
    //! @return param to select devices by vendor id's
    const util::ShPtr< command::ParameterInterface> &Platform::GetDeviceIDsParam()
    {
      static util::ShPtr< command::ParameterInterface> s_param
      (
        new command::Parameter
        (
          "device_ids",
          "comma separated list of device ids to be used for selected platform, e.g \"0,1\"",
          ""
        )
      );

      // end
      return s_param;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Platform::Platform()
    {
    }

    //! @brief construct from name
    //! @param NAME platform name e.g. "ATI", "NVIDIA", "ATI Stream", "NVIDIA CUDA"
    Platform::Platform( const std::string &NAME, cl_int *ERROR_PTR)
    {
      const std::string standardized_name( StandardizeName( NAME));
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return;
      }

      // iterate over platforms and find the one with the correct name
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        if
        (
          StandardizeName( itr->Name()) == standardized_name
        )
        {
          *this = *itr;
          return;
        }
      }

      Tools::AssignError( ERROR_PTR, CL_INVALID_VALUE);
    }

    //! @brief construct from an opencl platform
    Platform::Platform( const cl::Platform &PLATFORM, cl_int *ERROR_PTR) :
      cl::Platform( PLATFORM)
    {}

    //! @brief Clone function
    //! @return pointer to new Platform
    Platform *Platform::Clone() const
    {
      return new Platform( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Platform::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief access the name defined by CL_PLATFORM_NAME
    //! @param ERROR_PTR error will be written to this location
    //! @return the name of that platform
    std::string Platform::Name( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_NAME>( ERROR_PTR);
    }

    //! @brief get the name defined by CL_PLATFORM_NAME as standardized name
    //! @param ERROR_PTR error will be written to this location
    //! @return bcl standardized platform name (without spaces)
    std::string Platform::StandardizedName( cl_int *ERROR_PTR) const
    {
      return StandardizeName( Name( ERROR_PTR));
    }

    //! @brief get the version defined by CL_PLATFORM_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the version of the platform
    std::string Platform::Version( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_VERSION>( ERROR_PTR);
    }

    //! @brief get the vendor defined by CL_PLATFORM_VENDOR
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the platform
    std::string Platform::Vendor( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_VENDOR>( ERROR_PTR);
    }

    //! @brief get the profile defined by CL_PLATFORM_PROFILE
    //! @param ERROR_PTR error will be written to this location
    //! @return the profile of the platform
    std::string Platform::Profile( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_PROFILE>( ERROR_PTR);
    }

    //! @brief get the extensions defined by CL_PLATFORM_EXTENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return set of extensions
    storage::Set< Extension> Platform::Extensions( cl_int *ERROR_PTR) const
    {
      return GetExtensions().ExtensionsFromString( getInfo< CL_PLATFORM_EXTENSIONS>( ERROR_PTR));
    }

    //! @brief does it support ICD extension and what is the suffix
    //! @param ERROR_PTR error will be written to this location
    //! @return pair of bool, wheter is supports the icd extensions, and a string representing the icd suffix
    std::pair< bool, std::string> Platform::ICDExtension( cl_int *ERROR_PTR) const
    {
      std::pair< bool, std::string> icd_extension;
      storage::Set< Extension> extensions( Extensions( ERROR_PTR));

      icd_extension.first = extensions.Find( GetExtensions().e_khr_icd) != extensions.End();
      if( icd_extension.first)
      {
        Tools::AssignError( ERROR_PTR, getInfo( CL_PLATFORM_ICD_SUFFIX_KHR, &icd_extension.second));
      }

      // end
      return icd_extension;
    }

    //! @brief get a description string for the device
    //! @return string with a description for the device
    std::string Platform::Description( cl_int *ERROR_PTR) const
    {
      std::string description;

      description +=   "PROFILE        " + Profile( ERROR_PTR);
      description += "\nVERSION        " + Version( ERROR_PTR);
      description += "\nNAME           " + Name( ERROR_PTR);
      description += "\nVENDOR         " + Vendor( ERROR_PTR)    ;
      description += "\nEXTENSIONS     " + Extensions::ExtensionsToString( Extensions( ERROR_PTR));

      const std::pair< bool, std::string> icd_extension( ICDExtension( ERROR_PTR));
      if( icd_extension.first)
      {
        description += "\nICD_SUFFIX_KHR " + icd_extension.second;
      }

      // end
      return description;
    }

    //! @brief get devices for that platform
    //! @param DEVICE_TYPE the device type to be queried
    //! @return list of devices that match the queried types
    storage::Vector< Device> Platform::Devices( cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR) const
    {
      cl_int error( CL_SUCCESS);
      storage::Vector< Device> final_devices;
      std::vector< cl::Device> devices;

      // get devices of queried type
      error = getDevices( DEVICE_TYPE, &devices);

      // check error
      if( error == CL_DEVICE_NOT_FOUND)
      {
        BCL_MessageVrb
        (
          "cannot find device for platform and device type: " +
          Name() + ' ' + Device::TypeToString( DEVICE_TYPE)
        );

        return final_devices;
      }
      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
      }

      // iterate over devices and insert
      for( std::vector< cl::Device>::const_iterator itr( devices.begin()), itr_end( devices.end()); itr != itr_end; ++itr)
      {
        final_devices.PushBack( Device( *itr));
      }

      // end
      return final_devices;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Platform::Read( std::istream &ISTREAM)
    {
      // read the name
      std::string tmp_name;
      io::Serialize::Read( tmp_name, ISTREAM);
      if( !tmp_name.empty()) // plaform name was given - check if there is one with that name
      {
        cl_int error( CL_SUCCESS);
        Platform tmp_platform( tmp_name, &error);
        BCL_Assert( error == CL_SUCCESS, "unable to construct platform with name: " + tmp_name);
        *this = tmp_platform;
      }
      else // empty platform
      {
        *this = Platform();
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &Platform::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write the platform name
      io::Serialize::Write( Name(), OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief query all platforms available
    //! @param ERROR_PTR error will be written to this location
    //! @return vector of all platforms
    const storage::Vector< Platform> &Platform::QueryPlatforms( cl_int *ERROR_PTR)
    {
      static cl_int s_error( CL_SUCCESS);

      // final list
      static storage::Vector< Platform> s_final_platform_list;

      // if we already have a result, return it
      if( !s_final_platform_list.IsEmpty() || s_error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, s_error);
        return s_final_platform_list;
      }

      // get all platforms
      std::vector< cl::Platform> platform_list;
      s_error = cl::Platform::get( &platform_list);
      if( s_error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, s_error);
        return s_final_platform_list;
      }

      // final list
      // iterate over platforms and insert
      for( std::vector< cl::Platform>::const_iterator itr( platform_list.begin()), itr_end( platform_list.end()); itr != itr_end; ++itr)
      {
        s_final_platform_list.PushBack( Platform( *itr, ERROR_PTR));
      }

      // end
      return s_final_platform_list;
    }

    //! @brief query all platform names
    //! @param ERROR_PTR error will be written to this location
    //! @return vector of all platform names
    const storage::Vector< std::string> &Platform::QueryPlatformNamesStandardized( cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      // name
      static storage::Vector< std::string> s_names;

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return s_names;
      }
      else if( !s_names.IsEmpty() || platforms.IsEmpty())
      {
        // already initialized names
        return s_names;
      }

      // iterate over platforms
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        s_names.PushBack( StandardizeName( itr->Name()));
      }

      // end
      return s_names;
    }

    //! @brief get default device type as set in command line
    //! @brief return device type
    cl_device_type Platform::CommandLineDeviceType()
    {
      return Device::TypeFromString( GetDeviceTypeParam()->GetValue());
    }

    //! @brief get platform and devices from the command line
    //! @param PLATFORM platform, will be set to the command line platform
    //! @param DEVICES vector of devices to initialize
    //! @param ERROR_PTR error will be written to this location
    void Platform::InitializeFromCommandLine
    (
      Platform &PLATFORM,
      storage::Vector< Device> &DEVICES,
      cl_int *ERROR_PTR
    )
    {
      if( !GetPlatformFlag().IsDefined())
      {
        Tools::AssignError( ERROR_PTR, CL_PLATFORM_NOT_FOUND_KHR);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      const std::string standardized_name( StandardizeName( GetPlatformFlag()->GetFirstParameter()->GetValue()));
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        PLATFORM = Platform();
        DEVICES.Reset();
        Tools::AssignError( ERROR_PTR, error);
        return;
      }

      // get the command line device type
      cl_device_type device_type( CommandLineDeviceType());

      // keep track of whether a platform with the given name was found
      bool found_platform_with_name( false);

      // if a platform with the given name was found, assign error accordingly
      cl_int devices_err( CL_SUCCESS);

      // iterate over platforms and find the one with the correct name
      for
      (
        storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End());
        itr != itr_end;
        ++itr
      )
      {
        if( StandardizeName( itr->Name()) != standardized_name)
        {
          continue;
        }
        found_platform_with_name = true;
        cl_int local_devices_err( CL_SUCCESS);
        DEVICES = itr->Devices( device_type, &local_devices_err);

        // if no devices are available for that platform, issue a warning message that the default type has been updated
        if( DEVICES.IsEmpty() && !GetDeviceTypeParam()->GetWasSetInCommandLine())
        {
          DEVICES = itr->Devices( CL_DEVICE_TYPE_ALL, &local_devices_err);
          if( local_devices_err != CL_SUCCESS)
          {
            BCL_MessageStd
            (
              "Warning: Found platform with no devices (bad OpenCL driver setup?)"
            );
          }
          else if( !DEVICES.IsEmpty())
          {
            BCL_MessageStd
            (
              "Note: no devices with the default type = " + GetDeviceTypeParam()->GetDefaultValue()
              + " was found; instead taking the first device of the given platform"
            );
            GetDeviceTypeParam()->SetDefaultParameter( Device::TypeToString( CL_DEVICE_TYPE_ALL));
          }
        }
        if( local_devices_err != CL_SUCCESS)
        {
          devices_err = local_devices_err;
          continue;
        }
        if( !DEVICES.IsEmpty())
        {
          PLATFORM = *itr;
          break;
        }
      }

      // platform not found
      if( !found_platform_with_name)
      {
        Tools::AssignError( ERROR_PTR, CL_INVALID_VALUE);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      // no devices found for the platform
      if( devices_err != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, devices_err);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      // filter devices, if device_vendor_ids parameter was given
      if( GetDeviceIDsParam()->GetWasSetInCommandLine())
      {
        // vendor ids from commandline
        const storage::Vector< size_t> ids_vector( util::SplitStringToNumerical< size_t>( GetDeviceIDsParam()->GetValue(), ","));

        // create set to condense to unique ids
        const storage::Set< size_t> ids( ids_vector.Begin(), ids_vector.End());

        storage::Vector< Device> new_devices;

        // iterate over ids
        for( storage::Set< size_t>::const_iterator id_itr( ids.Begin()), id_itr_end( ids.End()); id_itr != id_itr_end; ++id_itr)
        {
          if( *id_itr >= DEVICES.GetSize())
          {
            BCL_MessageCrt
            (
              "there is no device with id: " + util::Format()( *id_itr) + " for platform: " +
              PLATFORM.Name() + " => ignoring this id"
            );
          }
          else
          {
            new_devices.PushBack( DEVICES( *id_itr));
          }
        }

        // end
        DEVICES = new_devices;
      }
    }

    //! @brief first platform with at least one device of given type
    //! @param DEVICE_TYPE the device type to be queried
    //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
    //! @return platform with gpu
    Platform Platform::FirstPlatformWithDeviceType( const cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      // get all platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return Platform();
      }

      const Platform *optimal_platform( NULL);
      // iterate over platforms and find first with given device type
      // iterate over platforms
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        storage::Vector< Device> devices( itr->Devices( DEVICE_TYPE));
        if( !devices.IsEmpty())
        {
          Tools::AssignError( ERROR_PTR, CL_SUCCESS);
          optimal_platform = &( *itr);

          // prefer Intel over AMD, since AMD writes all binaries to /tmp/OCL*.so, which causes the tmp folder to overflow
          if( optimal_platform->Name() == "Intel(R) OpenCL")
          {
            return *optimal_platform;
          }
        }
      }

      // no optimal platform found
      if( optimal_platform != NULL)
      {
        Tools::AssignError( ERROR_PTR, CL_SUCCESS);
        return *optimal_platform;
      }

      Tools::AssignError( ERROR_PTR, CL_INVALID_PLATFORM);
      // end - no platform found
      return Platform();
    }

    //! @brief identify optimal platform
    //! searches first for plaform with gpu, if there is non, the first one with cpu
    //! @param DEVICE_TYPE location where the device type of the optimal platform is stored
    //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
    //! @return optimal Platform
    Platform Platform::GetOptimalPlatform( cl_device_type &DEVICE_TYPE, cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      Platform optimal;

      // try gpu first
      optimal = FirstPlatformWithDeviceType( CL_DEVICE_TYPE_GPU, &error);
      if( error == CL_SUCCESS)
      {
        DEVICE_TYPE = CL_DEVICE_TYPE_GPU;
        return optimal;
      }

      // try cpu
      error = CL_SUCCESS;
      optimal = FirstPlatformWithDeviceType( CL_DEVICE_TYPE_CPU, &error);
      if( error == CL_SUCCESS)
      {
        DEVICE_TYPE = CL_DEVICE_TYPE_CPU;
        return optimal;
      }

      // no such platform
      Tools::AssignError( ERROR_PTR, CL_INVALID_PLATFORM);
      return optimal;
    }

    //! @brief standardize platform names, by replacing spaces
    //! @param PLATFORM_NAME name of platform that contains spaces like "ATI Stream"
    //! @return NAME with '_' instead of ' ' like "ATI_Stream"
    std::string Platform::StandardizeName( const std::string &PLATFORM_NAME)
    {
      std::string new_name( PLATFORM_NAME);
      // iterate over string
      for( std::string::iterator itr( new_name.begin()), itr_end( new_name.end()); itr != itr_end; ++itr)
      {
        if( *itr == ' ')
        {
          *itr = '_';
        }
      }

      // end
      return new_name;
    }

  } // namespace opencl
} // namespace bcl
