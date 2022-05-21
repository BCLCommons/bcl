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
#include "opencl/bcl_opencl_device.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief DataType as string
    //! @param DATA_TYPE the data type
    //! @return the DataType as string
    const std::string &Device::GetDataTypeString( const DataType &DATA_TYPE)
    {
      static const std::string s_data_type_strings[] =
      {
        "CHAR", "SHORT", "INT", "LONG", "FLOAT", "DOUBLE", "HALF",
        GetStaticClassName< DataType>()
      };

      return s_data_type_strings[ DATA_TYPE];
    }

    //! @brief convert DataType to cl device info preferred vector width
    //! @param DATA_TYPE the data type
    //! @return cl_device_info preferred vector width for data type
    cl_device_info Device::DataTypeToDeviceInfoPreferredVectorWidth( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Char:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR;
        case e_Short:  return CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT;
        case e_Int:    return CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT;
        case e_Long:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG;
        case e_Float:  return CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT;
        case e_Double: return CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE;
        case e_Half:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF;
        default:       return 0;
      }
    }

    //! @brief convert DataType to cl device info native vector width
    //! @param DATA_TYPE the data type
    //! @return cl_device_info native vector width for data type
    cl_device_info Device::DataTypeToDeviceInfoNativeVectorWidth( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Char:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR;
        case e_Short:  return CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT;
        case e_Int:    return CL_DEVICE_NATIVE_VECTOR_WIDTH_INT;
        case e_Long:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG;
        case e_Float:  return CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT;
        case e_Double: return CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE;
        case e_Half:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF;
        default:       return 0;
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Device::Device() :
      cl::Device()
    {
    }

    //! @brief construct from cl::Device
    //! @param DEVICE the cl::Device
    Device::Device( const cl::Device &DEVICE) :
      cl::Device( DEVICE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Device
    Device *Device::Clone() const
    {
      return new Device( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Device::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief access the name defined by CL_DEVICE_NAME
    //! @param ERROR_PTR error will be written to this location
    //! @return name of device
    std::string Device::Name( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_NAME>( ERROR_PTR);
    }

    //! @brief get the vendor defined by CL_DEVICE_VENDOR
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the device
    std::string Device::Vendor( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR>( ERROR_PTR);
    }

    //! @brief get the version defined by CL_DEVICE_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the version of the device
    std::string Device::Version( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VERSION>( ERROR_PTR);
    }

    //! @brief get the driver version defined by CL_DRIVER_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the device
    std::string Device::DriverVersion( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DRIVER_VERSION>( ERROR_PTR);
    }

    //! @brief returns device type defined by CL_DEVICE_TYPE
    //! @param ERROR_PTR error will be written to this location
    //! @return device type
    cl_device_type Device::DeviceType( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_TYPE>( ERROR_PTR);
    }

    //! @brief vendor id as defined by CL_DEVICE_VENDOR_ID
    //! @param ERROR_PTR error will be written to this location
    //! @return vendor id
    cl_uint Device::VendorID( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR_ID>( ERROR_PTR);
    }

    //! @brief the platform this device is associated with
    //! @param ERROR_PTR error will be written to this location
    //! @return the platform
    Platform Device::GetPlatform( cl_int *ERROR_PTR) const
    {
      return Platform( getInfo< CL_DEVICE_PLATFORM>( ERROR_PTR));
    }

    //! @brief return number of compute units on device
    //! @param ERROR_PTR error will be written to this location
    //! @return number of compute units
    cl_uint Device::MaxComputeUnits( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR_ID>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return max work item dimension
    cl_uint Device::MaxWorkItemDimension( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_SIZES
    //! @param ERROR_PTR error will be written to this location
    //! @return the max work item size for all dimensions
    std::vector< size_t> Device::MaxWorkItemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_ITEM_SIZES>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_GROUP_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return the max group size - which means the max sum of work item sizes
    size_t Device::MaxWorkGroupSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_GROUP_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_CLOCK_FREQUENCY in MHz
    //! @param ERROR_PTR error will be written to this location
    //! @return max clock for the device
    cl_uint Device::MaxClockFrequency( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_CLOCK_FREQUENCY>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_ADDRESS_BITS
    //! @param ERROR_PTR error will be written to this location
    //! @return a bitfiled for the address bits
    cl_bitfield Device::AddressBits( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ADDRESS_BITS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_MEM_ALLOC_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return maximal size for allocated memory
    cl_ulong Device::MaxMemAllocSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_MEM_ALLOC_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_GLOBAL_MEM_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return total memory size for global memory
    cl_ulong Device::GlobalMemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_GLOBAL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_ERROR_CORRECTION_SUPPORT
    //! @param ERROR_PTR error will be written to this location
    //! @return is error correcte memory supported (ECC)
    cl_bool Device::ErrorCorrection( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ERROR_CORRECTION_SUPPORT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_LOCAL_MEM_TYPE
    //! @param ERROR_PTR error will be written to this location
    //! @return type of local memory
    cl_device_local_mem_type Device::LocalMemType( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_LOCAL_MEM_TYPE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_LOCAL_MEM_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return size of local memory
    cl_ulong Device::LocalMemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_LOCAL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return size of constant memory buffer
    cl_ulong Device::MaxConstantBufferSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MEM_BASE_ADDR_ALIGN
    //! @param ERROR_PTR error will be written to this location
    //! @return address alignment
    cl_uint Device::MemBaseAddrAlign( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MEM_BASE_ADDR_ALIGN>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return
    cl_uint Device::MinDataTypeAlignSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return does device support out of order execution with multiple command queues
    cl_bool Device::QueueOutOfOrderExecution( cl_int *ERROR_PTR) const
    {
      const cl_command_queue_properties device_command_queue_properties( getInfo< CL_DEVICE_QUEUE_PROPERTIES>( ERROR_PTR));
      return device_command_queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
    }

    //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_PROFILING_ENABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is queue profiling available
    cl_bool Device::QueueProfiling( cl_int *ERROR_PTR) const
    {
      const cl_command_queue_properties device_command_queue_properties( getInfo< CL_DEVICE_QUEUE_PROPERTIES>( ERROR_PTR));
      return device_command_queue_properties & CL_QUEUE_PROFILING_ENABLE;
    }

    //! @brief defined by CL_DEVICE_PROFILING_TIMER_RESOLUTION
    //! @param ERROR_PTR error will be written to this location
    //! @return timer resolution for queue profiling
    size_t Device::ProfilingTimerResolution( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_PROFILING_TIMER_RESOLUTION>();
    }

    //! @brief defined by CL_DEVICE_ENDIAN_LITTLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is device little endian (false = big endian)
    cl_bool Device::EndianLittle( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ENDIAN_LITTLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_AVAILABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is device available
    cl_bool Device::Available( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_AVAILABLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_COMPILER_AVAILABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is compiler for device available
    cl_bool Device::CompilerAvailable( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_AVAILABLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_KERNEL
    //! @param ERROR_PTR error will be written to this location
    //! @return can kernel be executed
    cl_bool Device::ExecKernel( cl_int *ERROR_PTR) const
    {
      const cl_device_exec_capabilities device_exec_capabilities( getInfo< CL_DEVICE_EXECUTION_CAPABILITIES>( ERROR_PTR));
      return device_exec_capabilities & CL_EXEC_KERNEL;
    }

    //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_NATIVE_KERNEL
    //! @param ERROR_PTR error will be written to this location
    //! @return does device execute native kernels
    cl_bool Device::ExecNativeKernel( cl_int *ERROR_PTR) const
    {
      const cl_device_exec_capabilities device_exec_capabilities( getInfo< CL_DEVICE_EXECUTION_CAPABILITIES>( ERROR_PTR));
      return device_exec_capabilities & CL_EXEC_NATIVE_KERNEL;
    }

    //! @brief defined by CL_DEVICE_IMAGE_SUPPORT
    //! @param ERROR_PTR error will be written to this location
    //! @return does device support iamge data
    cl_bool Device::ImageSupport( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE_SUPPORT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_READ_IMAGE_ARGS
    //! @param ERROR_PTR error will be written to this location
    //! @return max number of image arguments
    cl_uint Device::MaxReadImageArgs( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_READ_IMAGE_ARGS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WRITE_IMAGE_ARGS
    //! @param ERROR_PTR error will be written to this location
    //! @return maximal number of image write arguments
    cl_uint Device::MaxWriteImageArgs( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WRITE_IMAGE_ARGS>( ERROR_PTR);
    }

    //! @brief defined by CL_FP_DENORM
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the denorm support for this data type (Float, Double, Half)
    //! @return floating point denorms supported
    cl_bool Device::FPDenorms( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_DENORM;
    }

    //! @brief defined by CL_FP_INF_NAN
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the inf nan support for this data type (Float, Double, Half)
    //! @return floating point inf and nan supported
    cl_bool Device::FPInfNan( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_INF_NAN;
    }

    //! @brief defined by CL_FP_ROUND_TO_NEAREST
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to nearest support for this data type (Float, Double, Half)
    //! @return floating point round to nearest supported
    cl_bool Device::FPRoundToNearest( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_NEAREST;
    }

    //! @brief defined by CL_FP_ROUND_TO_ZERO
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to zero support for this data type (Float, Double, Half)
    //! @return floating point round to sero supported
    cl_bool Device::FPRoundToZero( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_ZERO;
    }

    //! @brief defined by CL_FP_ROUND_TO_INF
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to inf support for this data type (Float, Double, Half)
    //! @return floating point round to +inf or -inf supported
    cl_bool Device::FPRoundToInf( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_INF;
    }

    //! @brief defined by CL_FP_FMA
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the fma support for this data type (Float, Double, Half)
    //! @return floating mointed fused multiply add supported
    cl_bool Device::FPFusedMultiplyAdd( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_FMA;
    }

    //! @brief defined by CL_DEVICE_IMAGE2D_MAX_WIDTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max width for 2d image
    size_t Device::Image2DMaxWidth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE2D_MAX_WIDTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE2D_MAX_HEIGHT
    //! @param ERROR_PTR error will be written to this location
    //! @return max height for 2d image
    size_t Device::Image2DMaxHeight( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE2D_MAX_HEIGHT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_WIDTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max width for 3d image
    size_t Device::Image3DMaxWidth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_WIDTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_HEIGHT
    //! @param ERROR_PTR error will be written to this location
    //! @return max height for 3d image
    size_t Device::Image3DMaxHeight( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_HEIGHT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_DEPTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max depth for 3d image
    size_t Device::Image3DMaxDepth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_DEPTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_SAMPLERS
    //! @param ERROR_PTR error will be written to this location
    //! @return max samlers for image
    cl_uint Device::MaxSamplers( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_SAMPLERS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_PARAMETER_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return max kernel parameter size
    size_t Device::MaxParameterSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_PARAMETER_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_EXTENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return all extensions the device supports
    storage::Set< Extension> Device::Extensions( cl_int *ERROR_PTR) const
    {
      return GetExtensions().ExtensionsFromString( getInfo< CL_DEVICE_EXTENSIONS>( ERROR_PTR));
    }

    //! @brief defined if extensions contain "cl_nv_device_attribute_query"
    //! @param ERROR_PTR error will be written to this location
    //! @return is nvidia device
    cl_bool Device::NVDevice( cl_int *ERROR_PTR) const
    {
      const storage::Set< Extension> extensions( Extensions( ERROR_PTR));
      return extensions.Find( GetExtensions().e_nv_device_attribute_query) != extensions.End();
    }

    //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return major compute capability for nvidia device
    cl_uint Device::ComputeCapabilityMajorNV( cl_int *ERROR_PTR) const
    {
      cl_uint cc_major( 999);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, &cc_major));
      return cc_major;
    }

    //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return minor compute capability for nvidia device
    cl_uint Device::ComputeCapabilityMinorNV( cl_int *ERROR_PTR) const
    {
      cl_uint cc_minor( 999);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, &cc_minor));
      return cc_minor;
    }

    //! @brief defined by CL_DEVICE_REGISTERS_PER_BLOCK_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return registers per block for nvidia device
    cl_uint Device::RegistersPerBlockNV( cl_int *ERROR_PTR) const
    {
      cl_uint registers( 0);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_REGISTERS_PER_BLOCK_NV, &registers));
      return registers;
    }

    //! @brief defined by CL_DEVICE_WARP_SIZE_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return warp size for nvidia device
    cl_uint Device::WarpSizeNV( cl_int *ERROR_PTR) const
    {
      cl_uint warp_size( 0);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_WARP_SIZE_NV, &warp_size));
      return warp_size;
    }

    //! @brief defined by CL_DEVICE_GPU_OVERLAP_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return overlapping gpu for nvidia device
    cl_bool Device::GPUOverlapNV( cl_int *ERROR_PTR) const
    {
      cl_bool gpu_overlap( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_GPU_OVERLAP_NV, &gpu_overlap));
      return gpu_overlap;
    }

    //! @brief defined by CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return kernel execution timeout for nvidia device
    cl_bool Device::KernelExecTimeoutNV( cl_int *ERROR_PTR) const
    {
      cl_bool exec_timeout( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV, &exec_timeout));
      return exec_timeout;
    }

    //! @brief defined by CL_DEVICE_INTEGRATED_MEMORY_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return integrated memory for nvidia device
    cl_bool Device::IntegratedMemoryNV( cl_int *ERROR_PTR) const
    {
      cl_bool integrated_memory( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_INTEGRATED_MEMORY_NV, &integrated_memory));
      return integrated_memory;
    }

    //! @brief preferred and native vector width - native and DataType half only available from opencl 1.1
    //! @param DATA_TYPE for which datatype
    //! @param ERROR_PTR error will be written to this location
    //! @return pair of preferred and native vector width, undefined values for opencl < 1.1
    std::pair< cl_uint, cl_uint> Device::PreferredNativeVectorWidth( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      std::pair< cl_uint, cl_uint> preferred_native( util::GetUndefined< cl_uint>(), util::GetUndefined< cl_uint>());

      Tools::AssignError( ERROR_PTR, getInfo( DataTypeToDeviceInfoPreferredVectorWidth( DATA_TYPE), &preferred_native.first));
      const cl_bool opencl11( OpenclCVersion( ERROR_PTR).find( "1.1") != std::string::npos);

      // TODO make this if statement test whether half type is available
      if( opencl11)
      {
        Tools::AssignError( ERROR_PTR, getInfo( DataTypeToDeviceInfoNativeVectorWidth( DATA_TYPE), &preferred_native.second));
      }

      return preferred_native;
    }

  ////////////////
  // opencl 1.1 //
  ////////////////

    //! defined by CL_DEVICE_OPENCL_C_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return OPENCLCVERSION - empty for opencl < 1.1
    std::string Device::OpenclCVersion( cl_int *ERROR_PTR) const
    {
      std::string opencl_version;
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_OPENCL_C_VERSION, &opencl_version));

      return opencl_version;
    }

    //! defined by CL_DEVICE_HOST_UNIFIED_MEMORY
    //! @param ERROR_PTR error will be written to this location
    //! @return
    cl_bool Device::HostUnifiedMemory( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_HOST_UNIFIED_MEMORY>( ERROR_PTR);
    }

    //! @brief get a description string for the device
    //! @param ERROR_PTR error will be written to this location
    //! @return string with a description for the device
    std::string Device::GetDescription( cl_int *ERROR_PTR) const
    {
      std::string description;

      description +=   "NAME                       " + Name( ERROR_PTR);
      description += "\nVENDOR                     " + Vendor( ERROR_PTR);
      description += "\nDRIVER_VERSION             " + DriverVersion( ERROR_PTR);

      description += "\nTYPE                       ";
      const cl_device_type type( DeviceType( ERROR_PTR));
      if( type & CL_DEVICE_TYPE_CPU)         description += "TYPE_CPU ";
      if( type & CL_DEVICE_TYPE_GPU)         description += "TYPE_GPU ";
      if( type & CL_DEVICE_TYPE_ACCELERATOR) description += "TYPE_ACCELERATOR ";
      if( type & CL_DEVICE_TYPE_DEFAULT)     description += "TYPE_DEFAULT ";

      description += "\nVERSION                    " + Version( ERROR_PTR);
      description += "\nVENDOR_ID                  " + util::Format()( VendorID( ERROR_PTR));
      description += "\nMAX_COMPUTE_UNITS          " + util::Format()( MaxComputeUnits( ERROR_PTR));
      description += "\nMAX_WORK_ITEM_DIMENSIONS   " + util::Format()( MaxWorkItemDimension( ERROR_PTR));

      description += "\nMAX_WORK_ITEM_SIZES        ";
      const std::vector< size_t> max_workitem_size( MaxWorkItemSize( ERROR_PTR));
      for( std::vector< size_t>::const_iterator itr( max_workitem_size.begin()), itr_end( max_workitem_size.end()); itr != itr_end; ++itr)
      {
        description += util::Format()( *itr) + " ";
      }

      description += "\nMAX_WORK_GROUP_SIZE        " + util::Format()( MaxWorkGroupSize( ERROR_PTR));
      description += "\nMAX_CLOCK_FREQUENCY        " + util::Format()( MaxClockFrequency( ERROR_PTR)) + " MHz";
      description += "\nADDRESS_BITS               " + util::Format()( AddressBits( ERROR_PTR));
      description += "\nMAX_MEM_ALLOC_SIZE         " + util::Format()( MaxMemAllocSize( ERROR_PTR) / ( 1024 * 1024)) + " MByte";
      description += "\nGLOBAL_MEM_SIZE            " + util::Format()( GlobalMemSize( ERROR_PTR) / ( 1024 * 1024)) + " MByte";
      description += "\nERROR_CORRECTION_SUPPORT   " + std::string(    ErrorCorrection( ERROR_PTR)   ? "yes" : "no");

      description += "\nLOCAL_MEM_TYPE             ";
      switch( LocalMemType( ERROR_PTR))
      {
        case CL_LOCAL:
          description += "local";
          break;
        case CL_GLOBAL:
          description += "global";
          break;
        default:
          description += "unknown";
          break;
      }

      description += "\nLOCAL_MEM_SIZE             " + util::Format()( LocalMemSize( ERROR_PTR) / 1024) + " kByte";
      description += "\nMAX_CONSTANT_BUFFER_SIZE   " + util::Format()( MaxConstantBufferSize( ERROR_PTR) / 1024) + " kByte";
      description += "\nMEM_BASE_ADDR_ALIGN        " + util::Format()( MemBaseAddrAlign( ERROR_PTR));
      description += "\nMIN_DATA_TYPE_ALIGN_SIZE   " + util::Format()( MinDataTypeAlignSize( ERROR_PTR));
      description += "\nQUEUE_OUT_OF_ORDER_EXEC    " + std::string(    QueueOutOfOrderExecution( ERROR_PTR) ? "yes" : "no");
      description += "\nQUEUE_PROFILING            " + std::string(    QueueProfiling( ERROR_PTR) ? "yes" : "no");
      description += "\nPROFILING_TIMER_RESOLUTION " + util::Format()( ProfilingTimerResolution( ERROR_PTR)) + " ns";

      description += "\nENDIAN_LITTLE              " + std::string( EndianLittle( ERROR_PTR)      ? "yes" : "no");
      description += "\nAVAILABLE                  " + std::string( Available( ERROR_PTR)         ? "yes" : "no");
      description += "\nCOMPILER_AVAILABLE         " + std::string( CompilerAvailable( ERROR_PTR) ? "yes" : "no");

      description += "\nEXEC_KERNEL                " + std::string( ExecKernel( ERROR_PTR)       ? "yes" : "no");
      description += "\nEXEC_NATIVE_KERNEL         " + std::string( ExecNativeKernel( ERROR_PTR) ? "yes" : "no");

      description += "\nIMAGE_SUPPORT              " + std::string(    QueueProfiling( ERROR_PTR) ? "yes" : "no");
      description += "\nMAX_READ_IMAGE_ARGS        " + util::Format()( MaxReadImageArgs( ERROR_PTR));
      description += "\nMAX_WRITE_IMAGE_ARGS       " + util::Format()( MaxWriteImageArgs( ERROR_PTR));

      for( size_t i( e_Float); i <= e_Half; ++i)
      {
        description += '\n' + GetDataTypeString( DataType( i));
        description += "\n  FP_DENORM                " + std::string( FPDenorms(          DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_INF_NAN               " + std::string( FPInfNan(           DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_NEAREST      " + std::string( FPRoundToNearest(   DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_ZERO         " + std::string( FPRoundToZero(      DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_INF          " + std::string( FPRoundToInf(       DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_FMA                   " + std::string( FPFusedMultiplyAdd( DataType( i), ERROR_PTR) ? "yes" : "no");
      }

      description += "\nIMAGE2D_MAX_WIDTH          " + util::Format()( Image2DMaxWidth( ERROR_PTR) );
      description += "\nIMAGE2D_MAX_HEIGHT         " + util::Format()( Image2DMaxHeight( ERROR_PTR));
      description += "\nIMAGE3D_MAX_WIDTH          " + util::Format()( Image3DMaxWidth( ERROR_PTR) );
      description += "\nIMAGE3D_MAX_HEIGHT         " + util::Format()( Image3DMaxHeight( ERROR_PTR));
      description += "\nIMAGE3D_MAX_DEPTH          " + util::Format()( Image3DMaxDepth( ERROR_PTR) );

      description += "\nMAX_SAMPLERS               " + util::Format()( MaxSamplers( ERROR_PTR));

      description += "\nMAX_PARAMETER_SIZE         " + util::Format()( MaxParameterSize( ERROR_PTR));

      description += "\nEXTENSIONS                 " + Extensions::ExtensionsToString( Extensions( ERROR_PTR));

      // for nvidia device
      if( NVDevice( ERROR_PTR))
      {
        description += "\nREGISTERS_PER_BLOCK_NV     " + util::Format()( RegistersPerBlockNV( ERROR_PTR));
        description += "\nWARP_SIZE_NV               " + util::Format()( WarpSizeNV( ERROR_PTR));
        description += "\nGPU_OVERLAP_NV             " + std::string(    GPUOverlapNV( ERROR_PTR)        ? "yes" : "no");
        description += "\nKERNEL_EXEC_TIMEOUT_NV     " + std::string(    KernelExecTimeoutNV( ERROR_PTR) ? "yes" : "no");
        description += "\nINTEGRATED_MEMORY_NV       " + std::string(    IntegratedMemoryNV( ERROR_PTR)  ? "yes" : "no");
      }

      description += "\nPREFERRED_VECTOR_WIDTH_CHAR   " + util::Format()( PreferredNativeVectorWidth( e_Char  , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_SHORT  " + util::Format()( PreferredNativeVectorWidth( e_Short , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_INT    " + util::Format()( PreferredNativeVectorWidth( e_Int   , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_LONG   " + util::Format()( PreferredNativeVectorWidth( e_Long  , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_FLOAT  " + util::Format()( PreferredNativeVectorWidth( e_Float , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_DOUBLE " + util::Format()( PreferredNativeVectorWidth( e_Double, ERROR_PTR).first);

      if( OpenclCVersion( ERROR_PTR).find( "1.1") == std::string::npos)
      {
        return description;
      }

      // opencl 1.1
      description += "\nPREFFERED_VECTOR_WIDTH_HALF   " + util::Format()( PreferredNativeVectorWidth( e_Half  , ERROR_PTR).first);
      description += "\nNATIVE_VECTOR_WIDTH_CHAR      " + util::Format()( PreferredNativeVectorWidth( e_Char  , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_SHORT     " + util::Format()( PreferredNativeVectorWidth( e_Short , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_INT       " + util::Format()( PreferredNativeVectorWidth( e_Int   , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_LONG      " + util::Format()( PreferredNativeVectorWidth( e_Long  , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_FLOAT     " + util::Format()( PreferredNativeVectorWidth( e_Float , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_DOUBLE    " + util::Format()( PreferredNativeVectorWidth( e_Double, ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_HALF      " + util::Format()( PreferredNativeVectorWidth( e_Half  , ERROR_PTR).second);
      description += "\nOPENCL_C_VERSION              " + OpenclCVersion( ERROR_PTR);
      description += "\nHOST_UNIFIED_MEMORY           " + std::string(    HostUnifiedMemory( ERROR_PTR) ? "yes" : "no");

      // end
      return description;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Device::Read( std::istream &ISTREAM)
    {
      // write

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &Device::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief the config flag to get the info for fp config
    //! @param DATA_TYPE float, Double or Half
    //! @return  CL_DEVICE_SINGLE_FP_CONFIG, CL_DEVICE_DOUBLE_FP_CONFIG, CL_DEVICE_HALF_FP_CONFIG, 0 for non available config
    cl_int Device::FpConfigFromDataType( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Float : return CL_DEVICE_SINGLE_FP_CONFIG;
        case e_Double: return CL_DEVICE_DOUBLE_FP_CONFIG;
        case e_Half  : return CL_DEVICE_HALF_FP_CONFIG;
        default      : return 0;
      }

      // end
      return 0;
    }

    //! @brief convert device type from string
    //! @param TYPE device type
    //! @return string for that device type
    const std::string &Device::TypeToString( const cl_device_type TYPE)
    {
      switch( TYPE)
      {
        case CL_DEVICE_TYPE_CPU:
        {
          static const std::string s_type_string_cpu( "TYPE_CPU");
          return s_type_string_cpu;
        }
        case CL_DEVICE_TYPE_GPU:
        {
          static const std::string s_type_string_gpu( "TYPE_GPU");
          return s_type_string_gpu;
        }
        case CL_DEVICE_TYPE_ACCELERATOR:
        {
          static const std::string s_type_string_accelerator( "TYPE_ACCELERATOR");
          return s_type_string_accelerator;
        }
        case CL_DEVICE_TYPE_ALL:
        {
          static const std::string s_type_string_all( "TYPE_ALL");
          return s_type_string_all;
        }
        case CL_DEVICE_TYPE_DEFAULT:
        default:
        {
          static const std::string s_type_string_default( "TYPE_DEFAULT");
          return s_type_string_default;
        }
      }
    }

    //! @brief convert string to device type
    //! @param TYPE_STRING string for device type
    //! @return string for that device type
    cl_device_type Device::TypeFromString( const std::string &TYPE_STRING)
    {
      if( TYPE_STRING == "TYPE_CPU")
      {
        return CL_DEVICE_TYPE_CPU;
      }
      else if( TYPE_STRING == "TYPE_GPU")
      {
        return CL_DEVICE_TYPE_GPU;
      }
      else if( TYPE_STRING == "TYPE_ACCELERATOR")
      {
        return CL_DEVICE_TYPE_ACCELERATOR;
      }
      else if( TYPE_STRING == "TYPE_ALL")
      {
        return CL_DEVICE_TYPE_ALL;
      }

      // end
      return CL_DEVICE_TYPE_DEFAULT;
    }

  } // namespace opencl
} // namespace bcl
