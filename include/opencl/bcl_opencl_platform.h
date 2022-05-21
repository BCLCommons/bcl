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

#ifndef BCL_OPENCL_PLATFORM_H_
#define BCL_OPENCL_PLATFORM_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_extensions.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <CL/cl.hpp>

namespace bcl
{
  namespace opencl
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Platform
    //! @brief this class provides easy access to the opencl platform and derives from cl::Platform
    //! it provides convenience functions for querying platform information and implements the util::ObjectInterface
    //!
    //! @see @link example_opencl_platform.cpp @endlink
    //! @author woetzen
    //! @date Aug 15, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Platform :
      public cl::Platform,
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief access to the command line param for device type
      //! @param DEVICE_TYPE_PTR the device type to be default, if NULL, use the GPU or last setting
      //! @return param to select processor type
      static util::ShPtr< command::ParameterInterface> &GetDeviceTypeParam( const cl_device_type *DEVICE_TYPE_PTR = NULL);

      //! @brief access to the command line param for device vendor id's
      //! @return param to select devices by vendor id's
      static const util::ShPtr< command::ParameterInterface> &GetDeviceIDsParam();

    public:

    //////////
    // data //
    //////////

      //! @brief access to the command line flag for platform
      //! @return flag to select platform and processor type
      static const util::ShPtr< command::FlagInterface> &GetPlatformFlag();

      //! @brief access whether OpenCL is explicitly disabled by flag
      //! @return true if opencl was explicitly disabled over the command line
      static bool &GetIsOpenclDisabled()
      {
        static bool s_is_disabled( false);
        return s_is_disabled;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Platform();

      //! @brief construct from name
      //! @param NAME platform name e.g. "ATI Stream", "NVIDIA CUDA", "ATI_Stream", "NVIDIA_CUDA"
      //! @param ERROR_PTR error will be written to this location
      Platform( const std::string &NAME, cl_int *ERROR_PTR = NULL);

      //! @brief construct from an opencl platform
      //! @param PLATFORM
      //! @param ERROR_PTR error will be written to this location
      Platform( const cl::Platform &PLATFORM, cl_int *ERROR_PTR = NULL);

      //! @brief Clone function
      //! @return pointer to new Platform
      Platform *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief access the name defined by CL_PLATFORM_NAME
      //! @param ERROR_PTR error will be written to this location
      //! @return the name of that platform
      std::string Name( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the name defined by CL_PLATFORM_NAME as standardized name
      //! @param ERROR_PTR error will be written to this location
      //! @return bcl standardized platform name (without spaces)
      std::string StandardizedName( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the version defined by CL_PLATFORM_VERSION
      //! @param ERROR_PTR error will be written to this location
      //! @return the version of the platform
      std::string Version( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the vendor defined by CL_PLATFORM_VENDOR
      //! @param ERROR_PTR error will be written to this location
      //! @return the vendor of the platform
      std::string Vendor( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the profile defined by CL_PLATFORM_PROFILE
      //! @param ERROR_PTR error will be written to this location
      //! @return the profile of the platform
      std::string Profile( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the extensions defined by CL_PLATFORM_EXTENSIONS
      //! @param ERROR_PTR error will be written to this location
      //! @return set of extensions
      storage::Set< Extension> Extensions( cl_int *ERROR_PTR = NULL) const;

      //! @brief does it support ICD extension and what is the suffix
      //! @param ERROR_PTR error will be written to this location
      //! @return pair of bool, wheter is supports the icd extensions, and a string representing the icd suffix
      std::pair< bool, std::string> ICDExtension( cl_int *ERROR_PTR = NULL) const;

      //! @brief get a description string for the platform containing all information
      //! @param ERROR_PTR error will be written to this location
      //! @return string with a description for the device
      std::string Description( cl_int *ERROR_PTR = NULL) const;

      //! @brief get devices for that platform
      //! @param DEVICE_TYPE the device type to be queried
      //! @param ERROR_PTR error will be written to this location
      //! @return list of devices that match the queried types
      storage::Vector< Device> Devices( cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR = NULL) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief query all platforms available
      //! @param ERROR_PTR error will be written to this location
      //! @return vector of all platforms
      static const storage::Vector< Platform> &QueryPlatforms( cl_int *ERROR_PTR = NULL);

      //! @brief query all platform names (standardized names)
      //! @param ERROR_PTR error will be written to this location
      //! @return vector of all platform names
      static const storage::Vector< std::string> &QueryPlatformNamesStandardized( cl_int *ERROR_PTR = NULL);

      //! @brief get default device type as set in command line
      //! @brief return device type
      static cl_device_type CommandLineDeviceType();

      //! @brief get platform and devices from the command line
      //! @param PLATFORM platform, will be set to the command line platform
      //! @param DEVICES vector of devices to initialize
      //! @param ERROR_PTR error will be written to this location
      static void InitializeFromCommandLine
      (
        Platform &PLATFORM,
        storage::Vector< Device> &DEVICES,
        cl_int *ERROR_PTR = NULL
      );

      //! @brief first platform with at least one device of given type
      //! @param DEVICE_TYPE the device type to be queried
      //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
      //! @return platform with gpu
      static Platform FirstPlatformWithDeviceType( const cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR = NULL);

    private:

      //! @brief identify optimal platform
      //! searches first for plaform with gpu, if there is non, the first one with cpu
      //! @param DEVICE_TYPE location where the device type of the optimal platform is stored
      //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
      //! @return optimal Platform
      static Platform GetOptimalPlatform( cl_device_type &DEVICE_TYPE, cl_int *ERROR_PTR = NULL);

      //! @brief standardize platform names, by replacing spaces
      //! @param PLATFORM_NAME name of platform that contains spaces like "ATI Stream"
      //! @return NAME with '_' instead of ' ' like "ATI_Stream"
      static std::string StandardizeName( const std::string &PLATFORM_NAME);

    }; // class Platform

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_PLATFORM_H_
