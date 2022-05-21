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
#include "opencl/bcl_opencl_tools.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_operations.h"
#include "opencl/bcl_opencl_kernel_source_file.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Tools::Tools()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Tools::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the command line platform
    //! @return Platform reference to the platform selected on the command line
    const Platform &Tools::GetPlatform()
    {
      return m_Platform;
    }

    //! @brief access the default context with the commandline platform and devices
    //! @return Context reference to the context created from the command line platform and devices
    const Context &Tools::GetContext()
    {
      return m_Context;
    }

    //! @brief access to the command line devices
    //! @return storage::Vector< Device> vector of devices selected on the command line
    const storage::Vector< Device> &Tools::GetDevices()
    {
      return m_Devices;
    }

    //! @brief access to the command line queues
    //! @return storage::Vector< CommandQueue> vector of command queues selected on the command line
    const storage::Vector< CommandQueue> &Tools::GetCommandQueues()
    {
      return m_CommandQueues;
    }

    //! @brief ceck if any opencl commandqueue is available
    //! @return true if the is at least one command queue available, false if not e.g. not platform or no device on the platform was found
    bool Tools::HasCommandQueues()
    {
      return !m_CommandQueues.IsEmpty();
    }

    //! @brief convenience function to access the first command queue
    //! @return CommandQueue the first command queue
    const CommandQueue &Tools::GetFirstCommandQueue()
    {
      return m_CommandQueues.FirstElement();
    }

    //! @brief return a buffer program for that datatype
    //! @param DATA_TYPE cpp data type
    //! @return program associated with first command queues device
    const cl::Program &Tools::GetBufferProgram( const util::CPPDataTypes::Types &DATA_TYPE, const CommandQueue &QUEUE)
    {
      std::pair< util::CPPDataTypes::Types, util::SiPtr< const CommandQueue> >
        search_pair( DATA_TYPE, util::SiPtr< const CommandQueue>( QUEUE));
      std::map
      <
        std::pair< util::CPPDataTypes::Types, util::SiPtr< const CommandQueue> >,
        cl::Program
      >::const_iterator itr( m_BufferPrograms.find( search_pair));

      if( itr != m_BufferPrograms.end())
      {
        return itr->second;
      }

      // try to compile program
      cl_int error_number( CL_SUCCESS);
      itr =
        m_BufferPrograms.insert
        (
          std::make_pair
          (
            search_pair,
            KernelSources::Compile( GetKernelSources().e_Buffer, DATA_TYPE, QUEUE, std::string(), &error_number)
          )
        ).first;
      BCL_Assert( error_number == CL_SUCCESS, "cannot compile buffer program for: " + util::Format()( DATA_TYPE));

      return itr->second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Initialize the platform, devices and the commandqueue from the command line flag
    //! This function should be called by all member functions that access queues, platform, or context
    void Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag()
    {
      GetTools().UpdateCurrentPlatformDevicesQueuesFromCommandLineFlagImpl();
    }

    //! @brief Initialize the platform, devices and the commandqueue from the command line flag
    //! This function should be called by all member functions that access queues, platform, or context
    void Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlagImpl()
    {
      static storage::Vector< std::string> s_platform;

      // if no platform is available, then just return
      if( !Platform::GetPlatformFlag().IsDefined())
      {
        return;
      }

      // if opencl was disabled, then just return
      if( Platform::GetPlatformFlag()->GetFirstParameter()->GetValue() == "Disable")
      {
        return;
      }

      // first, check whether anything has changed
      storage::Vector< std::string> platform_info( Platform::GetPlatformFlag()->GetStringList());
      platform_info.Append( KernelSources::GetKernelsBinaryPathFlag()->GetStringList());
      platform_info.Append( KernelSources::GetKernelsSourcePathFlag()->GetStringList());
      if( platform_info == s_platform)
      {
        return;
      }

      s_platform = platform_info;

      BCL_MessageVrb( "Updating current platform device queues for OpenCL");

      m_BufferPrograms.clear();
      cl_int error( CL_SUCCESS);

      Platform::InitializeFromCommandLine( m_Platform, m_Devices, &error);
      if( error != CL_SUCCESS)
      {
        BCL_MessageCrt( "unable to initialize platform for devices with error: " + ErrorString( error));
        return;
      }

      for( size_t gpus( 0); gpus < m_Devices.GetSize(); ++gpus)
      {
        BCL_MessageDbg
        (
          "UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag() gives this as device( " + util::Format()( gpus) + "): \n\n"
          + util::Format()( m_Devices( gpus).GetDescription())
        );
      }

      BCL_MessageVrb( "platform created: " + m_Platform.Name());

      // update the context
      m_Context = Context( m_Devices, &error);
      if( error != CL_SUCCESS)
      {
        BCL_MessageCrt( "unable to initialize context from devices with error: " + ErrorString( error));
        m_Platform = Platform();
        m_Devices.Reset();
        return;
      }

      // update s_platform, since the device may have been updated
      s_platform = Platform::GetPlatformFlag()->GetStringList();
      s_platform.Append( KernelSources::GetKernelsBinaryPathFlag()->GetStringList());
      s_platform.Append( KernelSources::GetKernelsSourcePathFlag()->GetStringList());

      BCL_MessageDbg( "Initialized opencl context");

      // reset the command queues
      m_CommandQueues.Reset();

      // iterate over devices and create a command queue each
      m_CommandQueues.AllocateMemory( m_Devices.GetSize());
      for
      (
        storage::Vector< Device>::const_iterator dev_itr( m_Devices.Begin()), dev_itr_end( m_Devices.End());
        dev_itr != dev_itr_end;
        ++dev_itr
      )
      {
        m_CommandQueues.PushBack( CommandQueue( m_Context, *dev_itr, &error));
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to initialize command queue from context and device with error: " + ErrorString( error));
          m_CommandQueues.Reset();
          m_Context = Context();
          m_Devices.Reset();
          m_Platform = Platform();

          return;
        }
      }

      // notify listeners that command queue has changed
      if( Platform::GetPlatformFlag().IsDefined())
      {
        m_QueueUpdateSignal.Emit( *this);
      }
      if( GetPlatform().GetPlatformFlag()->GetFlag())
      {
        linal::GetOperationsNonConst< float>().SetDefaultOperationsType( linal::Operations< float>::Operator( "OpenCL"));
        if( GetFirstCommandQueue().GetDevice( NULL).Extensions( NULL).Contains( GetExtensions().e_khr_fp64))
        {
          linal::GetOperationsNonConst< double>().SetDefaultOperationsType( linal::Operations< double>::Operator( "OpenCL"));
        }
      }
    }

    //! list all platforms and devices
    //! @param MESSAGE_LEVEL the message level, to use
    void Tools::ListPlatformsWithDevices( const util::Message::MessageLevel &MESSAGE_LEVEL)
    {
      cl_int error( CL_SUCCESS);

      // get all platforms
      const storage::Vector< Platform> &platform_list( Platform::QueryPlatforms( &error));
      BCL_Assert( error == CL_SUCCESS, "unable to acquire opencl platforms: " + ErrorString( error));

      // iterate over all platforms and display info
      for( storage::Vector< Platform>::const_iterator itr( platform_list.Begin()), itr_end( platform_list.End()); itr != itr_end; ++itr)
      {
        BCL_Message( MESSAGE_LEVEL, "OpenCL platform:\n" + itr->Description());

        // list all devices for that platform
        const storage::Vector< Device> devices( itr->Devices( CL_DEVICE_TYPE_ALL, &error));

        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to query devices");
          continue;
        }
        // iterate over devices and display properties
        for( storage::Vector< Device>::const_iterator dev_itr( devices.Begin()), dev_itr_end( devices.End()); dev_itr != dev_itr_end; ++dev_itr)
        {
          BCL_Message( MESSAGE_LEVEL, "device:\n" + dev_itr->GetDescription());
        }
      }
    }

    //! Get and log the binary (PTX) from the OpenCL compiler for the requested program & device
    //! @param cpProgram                   OpenCL program
    //! @param cPtxFileName   optional PTX file name
    void Tools::LogPtx( const cl::Program &cpProgram, const std::string &cPtxFileName)
    {
      cl_int error_number; // Error code var

      std::vector< char *> binaries( cpProgram.getInfo< CL_PROGRAM_BINARIES>( &error_number));
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "could not acquire the binaries error: " + ErrorString( error_number));
        return;
      }

      // successfully acquired the binary information
      BCL_MessageDbg( "writing " + util::Format()( binaries.size()) + " binaries for opencl operations to file: " + cPtxFileName);

      // write binaries to file
      io::OFStream write;
      io::File::MustOpenOFStream( write, cPtxFileName);
      // iterate over all binaries and write
      for( std::vector< char *>::iterator itr( binaries.begin()), itr_end( binaries.end()); itr != itr_end; ++itr)
      {
        if( *itr != NULL)
        {
          write.write( *itr, std::strlen( *itr));
          delete[] ( *itr);
        }
      }
      io::File::CloseClearFStream( write);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Tools::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Tools::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief assign error to a pointer
    //! @param ERROR_PTR pointer to error storage - can be NULL
    //! @param ACTUAL_ERROR the actual error to be assigned
    void Tools::AssignError( cl_int *ERROR_PTR, const cl_int ACTUAL_ERROR)
    {
      if( ERROR_PTR != NULL)
      {
        *ERROR_PTR = ACTUAL_ERROR;
      }
    }

    //! @brief calculate smallest multiple of block size larger equal than global_size needed
    //! @param GROUP_SIZE the number of threads in that block dimension
    //! @param GLOBAL_SIZE the total number of threads needed in that dimension
    //! @return the number of threads needed
    int Tools::RoundUp( const int GROUP_SIZE, const int GLOBAL_SIZE)
    {
      return ( ( ( GLOBAL_SIZE - 1) / GROUP_SIZE) + 1) * GROUP_SIZE;
    }

    //! @brief calculate the padding to add
    //! @param BLOCK_SIZE
    //! @param CURRENT_DIMENSION
    //! @return the amount of padding to add
    size_t Tools::CalcPadding( const size_t BLOCK_SIZE, const size_t CURRENT_DIMENSION)
    {
      return ( BLOCK_SIZE - ( CURRENT_DIMENSION % BLOCK_SIZE)) % BLOCK_SIZE;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @param QUEUE the command queue
    //! @param PROGRAM the cl program
    //! @param KERNEL_FILENAME the filename of the kernel
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int Tools::CompilePrograms
    (
      const util::CPPDataTypes::Types &PRECISION,
      const CommandQueue &QUEUE,
      cl::Program &PROGRAM,
      const std::string &KERNEL_FILENAME
    )
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;
      KernelSourceFile source_file( KERNEL_FILENAME);

      // construct opencl device
      const Device device( QUEUE.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      BCL_MessageDbg( "device description: " + util::Format()( device.GetDescription()));

      // construct kernel source strings
      const std::string kernels_source( source_file.GetSource( PRECISION, device.Extensions()));

      // check if kernel strings are empty
      if( kernels_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }

      // pushback strings to program sources vector
      source.push_back( std::make_pair( kernels_source.c_str(), kernels_source.length()));

      // create the program
      cl::Program &current_program( PROGRAM);
      current_program = cl::Program( QUEUE.GetContext(), source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
        return error_number;
      }

      // build the program
      error_number = current_program.build( std::vector< cl::Device>( 1, device), s_CLCompilerOptions);
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "build program error: " + Tools::ErrorString( error_number));
        std::string build_info;
        error_number = current_program.getBuildInfo( device, CL_PROGRAM_BUILD_LOG, &build_info);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "get build info error: " + Tools::ErrorString( error_number));
        }
        else
        {
          BCL_MessageCrt( "build log: " + build_info);
        }
        return error_number;
      }

      // end
      return error_number;
    }

    //! @brief common compiler options
    const char *Tools::s_CLCompilerOptions( "-cl-mad-enable -cl-fast-relaxed-math");

    // Helper function to get error string
    const std::string &Tools::ErrorString( cl_int error)
    {
      static const std::string s_error_string[] = {
        "CL_SUCCESS",                         //   0
        "CL_DEVICE_NOT_FOUND",                //  -1
        "CL_DEVICE_NOT_AVAILABLE",            //  -2
        "CL_COMPILER_NOT_AVAILABLE",          //  -3
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",   //  -4
        "CL_OUT_OF_RESOURCES",                //  -5
        "CL_OUT_OF_HOST_MEMORY",              //  -6
        "CL_PROFILING_INFO_NOT_AVAILABLE",    //  -7
        "CL_MEM_COPY_OVERLAP",                //  -8
        "CL_IMAGE_FORMAT_MISMATCH",           //  -9
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",      // -10
        "CL_BUILD_PROGRAM_FAILURE",           // -11
        "CL_MAP_FAILURE",                     // -12
        "",                                   // -13
        "",                                   // -14
        "",                                   // -15
        "",                                   // -16
        "",                                   // -17
        "",                                   // -18
        "",                                   // -19
        "",                                   // -20
        "",                                   // -21
        "",                                   // -22
        "",                                   // -23
        "",                                   // -24
        "",                                   // -25
        "",                                   // -26
        "",                                   // -27
        "",                                   // -28
        "",                                   // -29
        "CL_INVALID_VALUE",                   // -30
        "CL_INVALID_DEVICE_TYPE",             // -31
        "CL_INVALID_PLATFORM",                // -32
        "CL_INVALID_DEVICE",                  // -33
        "CL_INVALID_CONTEXT",                 // -34
        "CL_INVALID_QUEUE_PROPERTIES",        // -35
        "CL_INVALID_COMMAND_QUEUE",           // -36
        "CL_INVALID_HOST_PTR",                // -37
        "CL_INVALID_MEM_OBJECT",              // -38
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR", // -39
        "CL_INVALID_IMAGE_SIZE",              // -40
        "CL_INVALID_SAMPLER",                 // -41
        "CL_INVALID_BINARY",                  // -42
        "CL_INVALID_BUILD_OPTIONS",           // -43
        "CL_INVALID_PROGRAM",                 // -44
        "CL_INVALID_PROGRAM_EXECUTABLE",      // -45
        "CL_INVALID_KERNEL_NAME",             // -46
        "CL_INVALID_KERNEL_DEFINITION",       // -47
        "CL_INVALID_KERNEL",                  // -48
        "CL_INVALID_ARG_INDEX",               // -49
        "CL_INVALID_ARG_VALUE",               // -50
        "CL_INVALID_ARG_SIZE",                // -51
        "CL_INVALID_KERNEL_ARGS",             // -52
        "CL_INVALID_WORK_DIMENSION",          // -53
        "CL_INVALID_WORK_GROUP_SIZE",         // -54
        "CL_INVALID_WORK_ITEM_SIZE",          // -55
        "CL_INVALID_GLOBAL_OFFSET",           // -56
        "CL_INVALID_EVENT_WAIT_LIST",         // -57
        "CL_INVALID_EVENT",                   // -58
        "CL_INVALID_OPERATION",               // -59
        "CL_INVALID_GL_OBJECT",               // -60
        "CL_INVALID_BUFFER_SIZE",             // -61
        "CL_INVALID_MIP_LEVEL",               // -62
        "CL_INVALID_GLOBAL_WORK_SIZE"         // -63
      };

      static const int s_error_count( sizeof( s_error_string) / sizeof( s_error_string[ 0]));

      const int index( -error);

      static const std::string s_unknwon( "ERROR_OUT_OF_KNOWN_RANGE");

      // valid error
      if( index >= 0 && index < s_error_count)
      {
        return s_error_string[ index];
      }
      else
      {
        if( error == CL_PLATFORM_NOT_FOUND_KHR)
        {
          static const std::string s_CL_PLATFORM_NOT_FOUND_KHR_string( "CL_PLATFORM_NOT_FOUND_KHR");
          return s_CL_PLATFORM_NOT_FOUND_KHR_string;
        }
        BCL_MessageCrt( "given opencl error out of known range: " + util::Format()( error));
        return s_unknwon;
      }
    }

    //! @brief access to the tools singleton
    //! @return reference to the only instance of Tools
    Tools &GetTools()
    {
      static Tools s_tools;
      return s_tools;
    }

  } // namespace opencl
} // namespace bcl
