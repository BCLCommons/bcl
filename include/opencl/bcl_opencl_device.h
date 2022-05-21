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

#ifndef BCL_OPENCL_DEVICE_H_
#define BCL_OPENCL_DEVICE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

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
    //! @class Device
    //! @brief Device is wrapped around cl::Device and caches properties that would have to be queried otherwise
    //!
    //! @see @link example_opencl_device.cpp @endlink
    //! @author woetzen
    //! @date Aug 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Device :
      public cl::Device,
      public util::ObjectInterface
    {
    public:

      //! @enum DataType
      //! @brief datatypes that can be used on device
      enum DataType
      {
        e_Char,
        e_Short,
        e_Int,
        e_Long,
        e_Float,
        e_Double,
        e_Half,   //!< half type - only from opencl 1.1
        s_NumberDataTypes
      };

      //! @brief DataType as string
      //! @param DATA_TYPE the data type
      //! @return the DataType as string
      static const std::string &GetDataTypeString( const DataType &DATA_TYPE);

      //! @brief convert DataType to cl device info preferred vector width
      //! @param DATA_TYPE the data type
      //! @return cl_device_info preferred vector width for data type
      static cl_device_info DataTypeToDeviceInfoPreferredVectorWidth( const DataType &DATA_TYPE);

      //! @brief convert DataType to cl device info native vector width
      //! @param DATA_TYPE the data type
      //! @return cl_device_info native vector width for data type
      static cl_device_info DataTypeToDeviceInfoNativeVectorWidth( const DataType &DATA_TYPE);

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Device();

      //! @brief construct from cl::Device
      //! @param DEVICE the cl::Device
      Device( const cl::Device &DEVICE);

      //! @brief Clone function
      //! @return pointer to new Device
      Device *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // opencl 1.0 //
    ////////////////

      //! @brief access the name defined by CL_DEVICE_NAME
      //! @param ERROR_PTR error will be written to this location
      //! @return name of device
      std::string Name( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the vendor defined by CL_DEVICE_VENDOR
      //! @param ERROR_PTR error will be written to this location
      //! @return the vendor of the device
      std::string Vendor( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the version defined by CL_DEVICE_VERSION
      //! @param ERROR_PTR error will be written to this location
      //! @return the version of the device
      std::string Version( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the driver version defined by CL_DRIVER_VERSION
      //! @param ERROR_PTR error will be written to this location
      //! @return the vendor of the device
      std::string DriverVersion( cl_int *ERROR_PTR = NULL) const;

      //! @brief returns device type defined by CL_DEVICE_TYPE
      //! @param ERROR_PTR error will be written to this location
      //! @return device type
      cl_device_type DeviceType( cl_int *ERROR_PTR = NULL) const;

      //! @brief vendor id as defined by CL_DEVICE_VENDOR_ID
      //! @param ERROR_PTR error will be written to this location
      //! @return vendor id
      cl_uint VendorID( cl_int *ERROR_PTR = NULL) const;

      //! @brief the platform this device is associated with
      //! @param ERROR_PTR error will be written to this location
      //! @return the platform
      Platform GetPlatform( cl_int *ERROR_PTR = NULL) const;

      //! @brief return number of compute units on device
      //! @param ERROR_PTR error will be written to this location
      //! @return number of compute units
      cl_uint MaxComputeUnits( cl_int *ERROR_PTR = NULL) const;

      //! @brief return all extensions for that device defined by CL_DEVICE_EXTENSIONS
      //! @param ERROR_PTR error will be written to this location
      //! @return set of extensions that this device supports
      storage::Set< Extension> Extensions( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
      //! @param ERROR_PTR error will be written to this location
      //! @return max work item dimension
      cl_uint MaxWorkItemDimension( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_SIZES
      //! @param ERROR_PTR error will be written to this location
      //! @return the max work item size for all dimensions
      std::vector< size_t> MaxWorkItemSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_WORK_GROUP_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return the max group size - which means the max sum of work item sizes
      size_t MaxWorkGroupSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_CLOCK_FREQUENCY in MHz
      //! @param ERROR_PTR error will be written to this location
      //! @return max clock for the device
      cl_uint MaxClockFrequency( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_ADDRESS_BITS
      //! @param ERROR_PTR error will be written to this location
      //! @return a bitfiled for the address bits
      cl_bitfield AddressBits( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_MEM_ALLOC_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return maximal size for allocated memory
      cl_ulong MaxMemAllocSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_GLOBAL_MEM_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return total memory size for global memory
      cl_ulong GlobalMemSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_ERROR_CORRECTION_SUPPORT
      //! @param ERROR_PTR error will be written to this location
      //! @return is error correcte memory supported (ECC)
      cl_bool ErrorCorrection( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_LOCAL_MEM_TYPE
      //! @param ERROR_PTR error will be written to this location
      //! @return type of local memory
      cl_device_local_mem_type LocalMemType( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_LOCAL_MEM_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return size of local memory
      cl_ulong LocalMemSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return size of constant memory buffer
      cl_ulong MaxConstantBufferSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MEM_BASE_ADDR_ALIGN
      //! @param ERROR_PTR error will be written to this location
      //! @return address alignment
      cl_uint MemBaseAddrAlign( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return
      cl_uint MinDataTypeAlignSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
      //! @param ERROR_PTR error will be written to this location
      //! @return does device support out of order execution with multiple command queues
      cl_bool QueueOutOfOrderExecution( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_PROFILING_ENABLE
      //! @param ERROR_PTR error will be written to this location
      //! @return is queue profiling available
      cl_bool QueueProfiling( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_PROFILING_TIMER_RESOLUTION
      //! @param ERROR_PTR error will be written to this location
      //! @return timer resolution for queue profiling
      size_t ProfilingTimerResolution( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_ENDIAN_LITTLE
      //! @param ERROR_PTR error will be written to this location
      //! @return is device little endian (false = big endian)
      cl_bool EndianLittle( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_AVAILABLE
      //! @param ERROR_PTR error will be written to this location
      //! @return is device available
      cl_bool Available( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_COMPILER_AVAILABLE
      //! @param ERROR_PTR error will be written to this location
      //! @return is compiler for device available
      cl_bool CompilerAvailable( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_KERNEL
      //! @param ERROR_PTR error will be written to this location
      //! @return can kernel be executed
      cl_bool ExecKernel( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_NATIVE_KERNEL
      //! @param ERROR_PTR error will be written to this location
      //! @return does device execute native kernels
      cl_bool ExecNativeKernel( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE_SUPPORT
      //! @param ERROR_PTR error will be written to this location
      //! @return does device support iamge data
      cl_bool ImageSupport( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_READ_IMAGE_ARGS
      //! @param ERROR_PTR error will be written to this location
      //! @return max number of image arguments
      cl_uint MaxReadImageArgs( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_WRITE_IMAGE_ARGS
      //! @param ERROR_PTR error will be written to this location
      //! @return maximal number of image write arguments
      cl_uint MaxWriteImageArgs( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_DENORM
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the denorm support for this data type (Float, Double, Half)
      //! @return floating point denorms supported
      cl_bool FPDenorms( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_INF_NAN
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the inf nan support for this data type (Float, Double, Half)
      //! @return floating point inf and nan supported
      cl_bool FPInfNan( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_ROUND_TO_NEAREST
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the round to nearest support for this data type (Float, Double, Half)
      //! @return floating point round to nearest supported
      cl_bool FPRoundToNearest( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_ROUND_TO_ZERO
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the round to zero support for this data type (Float, Double, Half)
      //! @return floating point round to sero supported
      cl_bool FPRoundToZero( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_ROUND_TO_INF
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the round to inf support for this data type (Float, Double, Half)
      //! @return floating point round to +inf or -inf supported
      cl_bool FPRoundToInf( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_FP_FMA
      //! @param ERROR_PTR error will be written to this location
      //! @param DATA_TYPE the fma support for this data type (Float, Double, Half)
      //! @return floating mointed fused multiply add supported
      cl_bool FPFusedMultiplyAdd( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE2D_MAX_WIDTH
      //! @param ERROR_PTR error will be written to this location
      //! @return max width for 2d image
      size_t Image2DMaxWidth( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE2D_MAX_HEIGHT
      //! @param ERROR_PTR error will be written to this location
      //! @return max height for 2d image
      size_t Image2DMaxHeight( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE3D_MAX_WIDTH
      //! @param ERROR_PTR error will be written to this location
      //! @return max width for 3d image
      size_t Image3DMaxWidth( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE3D_MAX_HEIGHT
      //! @param ERROR_PTR error will be written to this location
      //! @return max height for 3d image
      size_t Image3DMaxHeight( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_IMAGE3D_MAX_DEPTH
      //! @param ERROR_PTR error will be written to this location
      //! @return max depth for 3d image
      size_t Image3DMaxDepth( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_SAMPLERS
      //! @param ERROR_PTR error will be written to this location
      //! @return max samlers for image
      cl_uint MaxSamplers( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_MAX_PARAMETER_SIZE
      //! @param ERROR_PTR error will be written to this location
      //! @return max kernel parameter size
      size_t MaxParameterSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined if extensions contain "cl_nv_device_attribute_query"
      //! @param ERROR_PTR error will be written to this location
      //! @return is nvidia device
      cl_bool NVDevice( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return major compute capability for nvidia device
      cl_uint ComputeCapabilityMajorNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return minor compute capability for nvidia device
      cl_uint ComputeCapabilityMinorNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_REGISTERS_PER_BLOCK_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return registers per block for nvidia device
      cl_uint RegistersPerBlockNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_WARP_SIZE_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return warp size for nvidia device
      cl_uint WarpSizeNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_GPU_OVERLAP_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return overlapping gpu for nvidia device
      cl_bool GPUOverlapNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return kernel execution timeout for nvidia device
      cl_bool KernelExecTimeoutNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief defined by CL_DEVICE_INTEGRATED_MEMORY_NV
      //! @param ERROR_PTR error will be written to this location
      //! @return integrated memory for nvidia device
      cl_bool IntegratedMemoryNV( cl_int *ERROR_PTR = NULL) const;

      //! @brief preferred and native vector width - native and DataType half only available from opencl 1.1
      //! @param DATA_TYPE for which datatype
      //! @param ERROR_PTR error will be written to this location
      //! @return pair of preferred and native vector width, undefined values for opencl < 1.1
      std::pair< cl_uint, cl_uint> PreferredNativeVectorWidth( const DataType DATA_TYPE, cl_int *ERROR_PTR = NULL) const;

    ////////////////
    // opencl 1.1 //
    ////////////////

      //! defined by CL_DEVICE_OPENCL_C_VERSION
      //! @param ERROR_PTR error will be written to this location
      //! @return OPENCLCVERSION - empty for opencl < 1.1
      std::string OpenclCVersion( cl_int *ERROR_PTR = NULL) const;

      //! defined by CL_DEVICE_HOST_UNIFIED_MEMORY
      //! @param ERROR_PTR error will be written to this location
      //! @return
      cl_bool HostUnifiedMemory( cl_int *ERROR_PTR = NULL) const;

      //! @brief get a description string for the device
      //! @param ERROR_PTR error will be written to this location
      //! @return string with a description for the device
      std::string GetDescription( cl_int *ERROR_PTR = NULL) const;

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

    private:

      //! @brief the config flag to get the info for fp config
      //! @param DATA_TYPE float, Double or Half
      //! @return  CL_DEVICE_SINGLE_FP_CONFIG, CL_DEVICE_DOUBLE_FP_CONFIG, CL_DEVICE_HALF_FP_CONFIG, 0 for non available config
      static cl_int FpConfigFromDataType( const DataType &DATA_TYPE);

    public:

      //! @brief convert device type from string
      //! @param TYPE device type
      //! @return string for that device type
      static const std::string &TypeToString( const cl_device_type TYPE);

      //! @brief convert string to device type
      //! @param TYPE_STRING string for device type
      //! @return string for that device type
      static cl_device_type TypeFromString( const std::string &TYPE_STRING);

    }; // class Device

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_DEVICE_H_ 
