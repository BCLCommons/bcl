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

#ifndef BCL_OPENCL_TOOLS_H_
#define BCL_OPENCL_TOOLS_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_context.h"
#include "bcl_opencl_device.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_platform.h"
#include "signal/bcl_signal_signal.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <map>
#if defined (__APPLE__) || defined(MACOSX)
  #include <OpenCL/opencl.h>
#else
  #include <CL/cl.hpp>
  #include <CL/opencl.h>
#endif

namespace bcl
{
  namespace opencl
  {

    //! @brief access to the tools singleton
    //! @return reference to the only instance of Tools
    BCL_API Tools &GetTools();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Tools
    //! @brief the Tools class contains various opencl helper functions. It also contains the default platform, context,
    //!        command queue and devices as defined on the commandline
    //!
    //! @see @link example_opencl_tools.cpp @endlink
    //! @author woetzen
    //! @date Jun 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Tools
    {
    private:

    /////////////
    // friends //
    /////////////

      friend Tools &GetTools();

    //////////
    // data //
    //////////

      //! platform from the commandline
      Platform m_Platform;

      //! devices from the commandline
      storage::Vector< Device> m_Devices;

      //! default context
      Context m_Context;

      //! command queues for the devices
      storage::Vector< CommandQueue> m_CommandQueues;

      //! @brief buffer programs for the devices
      mutable std::map< std::pair< util::CPPDataTypes::Types, util::SiPtr< const CommandQueue> >, cl::Program>
        m_BufferPrograms;

      //! signal handler for updates on the default platform, command queue and context
      mutable signal::Signal1< Tools &> m_QueueUpdateSignal;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Tools();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the command line platform
      //! @return Platform reference to the platform selected on the command line
      const Platform &GetPlatform();

      //! @brief access the default context with the commandline platform and devices
      //! @return Context reference to the context created from the command line platform and devices
      const Context &GetContext();

      //! @brief access to the command line devices
      //! @return storage::Vector< Device> vector of devices selected on the command line
      const storage::Vector< Device> &GetDevices();

      //! @brief access to the command line queues
      //! @return storage::Vector< CommandQueue> vector of command queues selected on the command line
      const storage::Vector< CommandQueue> &GetCommandQueues();

      //! @brief ceck if any opencl commandqueue is available
      //! @return true if the is at least one command queue available, false if not e.g. not platform or no device on the platform was found
      bool HasCommandQueues();

      //! @brief convenience function to access the first command queue
      //! @return CommandQueue the first command queue
      const CommandQueue &GetFirstCommandQueue();

      //! @brief return a buffer program for that datatype
      //! @param DATA_TYPE cpp data type
      //! @return program associated with first command queues device
      const cl::Program &GetBufferProgram( const util::CPPDataTypes::Types &DATA_TYPE, const CommandQueue &QUEUE);

      //! @brief get the update signal handler, that emits if the default commandqueue and dependent objects are updated
      //! @return signal::Signal1< Tools&> the handler that emit a signal on updating to queue
      signal::Signal1< Tools &> &GetQueueUpdateSignal() const
      {
        return m_QueueUpdateSignal;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Initialize the platform, devices and the commandqueue from the command line flag
      //! This function should be called by all member functions that access queues, platform, or context
      static void UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag();

      //! list all platforms and devices
      //! @param MESSAGE_LEVEL the message level, to use
      static void ListPlatformsWithDevices( const util::Message::MessageLevel &MESSAGE_LEVEL);

      //! Get and log the binary (PTX) from the OpenCL compiler for the requested program & device
      //! @param cpProgram      OpenCL program
      //! @param cPtxFileName   optional PTX file name
      static void LogPtx( const cl::Program &cpProgram, const std::string &cPtxFileName);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief Initialize the platform, devices and the commandqueue from the command line flag
      //! This function should be called by all member functions that access queues, platform, or context
      void UpdateCurrentPlatformDevicesQueuesFromCommandLineFlagImpl();

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief assign error to a pointer
      //! @param ERROR_PTR pointer to error storage - can be NULL
      //! @param ACTUAL_ERROR the actual error to be assigned
      static void AssignError( cl_int *ERROR_PTR, const cl_int ACTUAL_ERROR);

      //! @brief calculate smallest multiple of block size larger equal than global_size needed
      //! @param GROUP_SIZE the number of threads in that block dimension
      //! @param GLOBAL_SIZE the total number of threads needed in that dimension
      //! @return the number of threads needed
      static int RoundUp( const int GROUP_SIZE, const int GLOBAL_SIZE);

      //! @brief calculate the padding to add
      //! @param BLOCK_SIZE
      //! @param CURRENT_DIMENSION
      //! @return the amount of padding to add
      static size_t CalcPadding( const size_t BLOCK_SIZE, const size_t CURRENT_DIMENSION);

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @param QUEUE the command queue
      //! @param PROGRAM the cl program
      //! @param KERNEL_FILENAME the filename of the kernel
      //! @return ERROR error that occured, CL_SUCCESS if no error
      static cl_int CompilePrograms
      (
        const util::CPPDataTypes::Types &PRECISION,
        const CommandQueue &QUEUE,
        cl::Program &PROGRAM,
        const std::string &KERNEL_FILENAME
      );

      //! @brief common compiler options
      static const char *s_CLCompilerOptions;

      //! @brief Helper function to get error string
      static const std::string &ErrorString( cl_int error);

    }; // class Tools

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_TOOLS_H_
