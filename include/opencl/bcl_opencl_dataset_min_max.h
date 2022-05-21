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

#ifndef BCL_OPENCL_DATASET_MIN_MAX_H_
#define BCL_OPENCL_DATASET_MIN_MAX_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_device.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_platform.h"
#include "bcl_opencl_tools.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetMinMax
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_opencl_dataset_min_max.cpp @endlink
    //! @author woetzen
    //! @date Nov 16, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class DataSetMinMax :
      public math::FunctionInterfaceSerializable< linal::MatrixConstInterface< t_DataType>, linal::Vector< t_DataType> >
    {

    private:

    //////////
    // data //
    //////////

      CommandQueue m_CommandQueue; //!< command queue
      cl::Program  m_Program;      //!< OpenCL program that contains kernels for each precision

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DataSetMinMax() :
        m_CommandQueue()
      {
      }

      //! @brief construct from command queue
      DataSetMinMax( const CommandQueue &COMMAND_QUEUE) :
        m_CommandQueue()
      {
        Initialize( COMMAND_QUEUE);
      }

      //! @brief Clone function
      //! @return pointer to new DataSetMinMax
      DataSetMinMax< t_DataType> *Clone() const
      {
        return new DataSetMinMax< t_DataType>( *this);
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

    ////////////////
    // operations //
    ////////////////

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const
      {
        return true;
      }

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE)
      {
        m_CommandQueue = COMMAND_QUEUE;

        cl_int error_number( CL_SUCCESS);
        error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>());

        if( error_number != CL_SUCCESS)
        {
          BCL_MessageDbg( "error compiling programs: " + Tools::ErrorString( error_number));
          m_CommandQueue = CommandQueue();
          return false;
        }

        return true;
      }

      linal::Vector< t_DataType> Max( const Matrix< t_DataType> &DEVICE_MATRIX) const
      {
        cl_int error_number( CL_SUCCESS);

        // default block size
        const size_t block_size( 256);
        const size_t vector_size( DEVICE_MATRIX.GetNumberCols());

        // number of rows in matrix
        const size_t number_rows( DEVICE_MATRIX.GetNumberRows());

        const size_t num_groups( ( number_rows - 1) / block_size + 1);

        // create the result buffer
        Matrix< t_DataType> device_result( vector_size, num_groups, m_CommandQueue, 0, 0, -std::numeric_limits< t_DataType>::max());

        // Create the kernel
        cl::Kernel kernel( m_Program, "DataSetStatisticMax", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const cl::NDRange block_dimensions( block_size, 1);
        const cl::NDRange offset;
        const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, number_rows), vector_size);

        error_number = kernel.setArg( 0, DEVICE_MATRIX.GetData());
        error_number |= kernel.setArg( 1, cl_uint( number_rows));
        error_number |= kernel.setArg( 2, cl_uint( vector_size));
        error_number |= kernel.setArg( 3, device_result.GetData());
        error_number |= kernel.setArg( 4, block_size * sizeof( t_DataType), 0);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        linal::Matrix< t_DataType> result_max( device_result.GetHostMatrix());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // return result without padding
        linal::Vector< t_DataType> result_vector( vector_size);
        for( size_t i( 0); i < vector_size; ++i)
        {
          result_vector( i) = *std::max_element( result_max.Begin() + i * num_groups, result_max.Begin() + ( i + 1) * num_groups);
        }

        // end
        return result_vector;
      }

      linal::Vector< t_DataType> Min( const Matrix< t_DataType> &DEVICE_MATRIX) const
      {
        cl_int error_number( CL_SUCCESS);

        // default block size
        const size_t block_size( 256);
        const size_t vector_size( DEVICE_MATRIX.GetNumberCols());

        // number of rows in matrix
        const size_t number_rows( DEVICE_MATRIX.GetNumberRows());

        const size_t num_groups( ( number_rows - 1) / block_size + 1);

        // create the result buffer
        Matrix< t_DataType> device_result( vector_size, num_groups, m_CommandQueue, 0, 0, std::numeric_limits< t_DataType>::max());

        // Create the kernel
        cl::Kernel kernel( m_Program, "DataSetStatisticMin", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const cl::NDRange block_dimensions( block_size, 1);
        const cl::NDRange offset;
        const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, number_rows), vector_size);

        error_number = kernel.setArg( 0, DEVICE_MATRIX.GetData());
        error_number |= kernel.setArg( 1, cl_uint( number_rows));
        error_number |= kernel.setArg( 2, cl_uint( vector_size));
        error_number |= kernel.setArg( 3, device_result.GetData());
        error_number |= kernel.setArg( 4, block_size * sizeof( t_DataType), 0);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        linal::Matrix< t_DataType> result_min( device_result.GetHostMatrix());
        // return result without padding
        linal::Vector< t_DataType> result_vector( vector_size);
        for( size_t i( 0); i < vector_size; ++i)
        {
          result_vector( i) = *std::min_element( result_min.Begin() + i * num_groups, result_min.Begin() + ( i + 1) * num_groups);
        }

        // end
        return result_vector;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief determine the min and max value of each column in given matrix
      //! @param MATRIX the matrix of interest
      linal::Vector< t_DataType> operator()( const linal::MatrixConstInterface< t_DataType> &MATRIX) const
      {
        Matrix< t_DataType> device_matrix( MATRIX, m_CommandQueue);

        return Max( device_matrix);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occurred, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
      {
        cl_int error_number( CL_SUCCESS);

        // compile the program
        cl::Program::Sources source;
        KernelSourceFile max_source_file( "dataset_statistic_max.cl");
        KernelSourceFile min_source_file( "dataset_statistic_min.cl");

        const Device device( m_CommandQueue.GetDevice( &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "cannot get device for commandline");

        const std::string reduction_max_source( max_source_file.GetSource( PRECISION, device.Extensions()));
        const std::string reduction_min_source( min_source_file.GetSource( PRECISION, device.Extensions()));
        if( reduction_max_source.empty() || reduction_min_source.empty())
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        source.push_back( std::make_pair( reduction_max_source.c_str(), reduction_max_source.length()));
        source.push_back( std::make_pair( reduction_min_source.c_str(), reduction_min_source.length()));

        const Context context( m_CommandQueue.GetContext( &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "cannot get context from command queue");

        // create the program
        cl::Program &current_program( m_Program);
        current_program = cl::Program( context, source, &error_number);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
          return error_number;
        }
        // build the program
        error_number = current_program.build( std::vector< cl::Device>( 1, device));
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
        }

        // end
        return error_number;
      }

    }; // template class DataSetMinMax

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_DATASET_MIN_MAX_H_
