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

#ifndef BCL_OPENCL_TRANSFER_FUNCTION_GAUSSIAN_H_
#define BCL_OPENCL_TRANSFER_FUNCTION_GAUSSIAN_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_device.h"
#include "bcl_opencl_kernel_source_alternative.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_kernel_source_string.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_platform.h"
#include "bcl_opencl_tools.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_const_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransferFunctionGaussian
    //! @brief this class performs a gaussian transfer function using opencl
    //! @tparam data type of float, t_DataType, int, complex, etc...
    //!
    //! @see @link example_opencl_transfer_function_gaussian.cpp @endlink
    //! @author loweew
    //! @date 09/16/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class TransferFunctionGaussian :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! queue
      CommandQueue m_Queue;

      //! opencl program
      cl::Program m_Program;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      static const char* s_CLGaussianF;

      static const char* s_CLGaussiandF;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TransferFunctionGaussian()
      {
        if( GetTools().HasCommandQueues())
        {
          UpdateQueue( GetTools().GetFirstCommandQueue());
        }
      }

      //! @brief default constructor
      TransferFunctionGaussian( const CommandQueue &QUEUE) :
        m_Queue()
      {
        UpdateQueue( QUEUE);
      }

      //! @brief virtual copy constructor
      //! @return pointer to a new TransferFunctionGaussian copied from this one
      TransferFunctionGaussian *Clone() const
      {
        return new TransferFunctionGaussian( *this);
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

      //! @brief get the working x range of that function
      //! @return x math::Range this function works on
      const math::Range< t_DataType> &GetOutputRange() const
      {
        static const math::Range< t_DataType> s_OutputRange( 0, 1);
        return s_OutputRange;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief function for one argument - depends on x ( ARGUMENT( 0))
      //! @param ELEMENT element on which to perform transfer function
      //! @return the transformed element
      t_DataType F( const t_DataType &ELEMENT) const
      {
        return exp( -math::Sqr( ELEMENT) / t_DataType( 2));
      }

      //! @brief derivative for one argument - depends on x ( ARGUMENT( 0)) and F( x) ( ARGUMENT( 1))
      //! @param ELEMENT_X x
      //! @param ELEMENT_Y F( x)
      //! @return the derivative of the transfer function
      t_DataType dF( const t_DataType &ELEMENT_X, const t_DataType &ELEMENT_Y) const
      {
        return ELEMENT_X * ELEMENT_Y;
      }

      //! @brief overloaded F for vectors, calls F(vector) of the base class
      //! @param VECTOR the vector of elements to perform transfer function on
      //! @return the resulting vector of values
      linal::Vector< t_DataType> F( const linal::Vector< t_DataType> &VECTOR) const
      {
        return GaussianF( VECTOR);
      }

      //! @brief overloaded dF for Vectors, calls dF(vector, vector) of the base class
      //! @param VECTOR_X the vector of elements x
      //! @param VECTOR_Y the vector of elements F( x)
      //! @return the resulting vector of elements from the transfer function derivative
      linal::Vector< t_DataType> dF( const linal::Vector< t_DataType> &VECTOR_X, const linal::Vector< t_DataType> &VECTOR_Y) const
      {
        return GaussiandF( VECTOR_X, VECTOR_Y);
      }

      //! @brief overloaded F for matrices
      //! @param MATRIX the matrix of elements to perform transfer function on
      //! @return the resulting matrix of values
      linal::Matrix< t_DataType> F( const linal::MatrixConstInterface< t_DataType> &MATRIX) const
      {
        return GaussianF( MATRIX);
      }

      //! @brief overloaded dF for Matrices, calls dF(matrix, matrix) of the base class
      //! @param MATRIX_X the matrix of elements
      //! @param MATRIX_Y the matrix of elements
      //! @return the resulting matrix of elements from the transfer function derivative
      linal::Matrix< t_DataType> dF( const linal::MatrixConstInterface< t_DataType> &MATRIX_X, const linal::MatrixConstInterface< t_DataType> &MATRIX_Y) const
      {
        return GaussiandF( MATRIX_X, MATRIX_Y);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return the kernel source for the F( x)
      //! @return kernel source of F( x)
      const KernelSource &GetKernelF()
      {
        static const KernelSource e_gaussian_kernel
        (
          GetKernelSources().AddEnum
          (
            "TransferFunctionGaussianF",
            util::ShPtr< KernelSourceInterface>
            (
              new KernelSourceAlternative
              (
                KernelSourceFile( "transfer_function_gaussian_F.cl"),
                KernelSourceString( s_CLGaussianF)
              )
            )
          )
        );
        return e_gaussian_kernel;
      }

      //! @brief return the kernel source for dF( x)
      //! @return kernel source of dF( x)
      const KernelSource &GetKerneldF()
      {
        static const KernelSource e_gaussian_derivative_kernel( GetKernelSources().AddEnum( "TransferFunctionGaussiandF", util::ShPtr< KernelSourceInterface>( new KernelSourceAlternative( KernelSourceFile( "transfer_function_gaussian_dF.cl"), KernelSourceString( s_CLGaussiandF)))));
        return e_gaussian_derivative_kernel;
      }

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occured, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
      {
        cl_int error_number( CL_SUCCESS);

        const Device device( m_Queue.GetDevice( &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

        // compile the program
        cl::Program::Sources source;
        const std::string gaussian_source_F( ( *GetKernelF())->GetSource( PRECISION, device.Extensions()));
        if
        (
          gaussian_source_F.empty()
        )
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        const std::string gaussian_source_dF( ( *GetKerneldF())->GetSource( PRECISION, device.Extensions()));
        if
        (
          gaussian_source_dF.empty()
        )
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        source.push_back( std::make_pair( gaussian_source_F.c_str()      , gaussian_source_F.length()));
        source.push_back( std::make_pair( gaussian_source_dF.c_str()      , gaussian_source_dF.length()));

        // create the program
        cl::Program &current_program( m_Program);
        current_program = cl::Program( m_Queue.GetContext(), source, &error_number);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
          return error_number;
        }
        // build the program
//        error_number = current_program.build( std::vector< cl::Device>( 1, device), KernelSourceInterface::GetPrecisionCompilerOptions( PRECISION).c_str());
        error_number = current_program.build( std::vector< cl::Device>( 1, device), KernelSourceInterface::GetAdditionalCompilerOptions().c_str());
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

      //! @brief perform gaussian transfer function on a vector of elements
      //! @param VECTOR the vector of elements
      //! @return vector of gaussian results
      linal::Vector< t_DataType> GaussianF
      (
        const linal::Vector< t_DataType> &VECTOR
      ) const
      {
        const int number_elements( VECTOR.GetSize());
        const int block_size( 256);
        cl_int error_number = CL_SUCCESS;

        Vector< t_DataType> device_input( VECTOR, m_Queue);
        Vector< t_DataType> device_output( number_elements, m_Queue);

        cl::Kernel kernel( m_Program, "GaussianKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::NDRange local_worksize( block_size); //< all thread blocks have same dimensions
        cl::NDRange global_worksize( Tools::RoundUp( block_size, number_elements));
        cl::NDRange offset;
        error_number  = kernel.setArg( 0,  device_output.GetData());
        error_number |= kernel.setArg( 1,  device_input.GetData());
        error_number |= kernel.setArg( 2,  number_elements);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

        linal::Vector< t_DataType> output_vector( device_output.GetHostVector());

        // wait for all commands to finish in queue
        m_Queue.finish();

        return output_vector;
      }

      //! @brief perform gaussian derivative
      //! @param VECTOR_A the vector of elements - x
      //! @param VECTOR_B the vector of elements - F( x)
      //! @return vector of derivative results
      linal::Vector< t_DataType> GaussiandF
      (
        const linal::Vector< t_DataType> &VECTOR_A,
        const linal::Vector< t_DataType> &VECTOR_B
      ) const
      {
        BCL_Assert( VECTOR_A.GetSize() == VECTOR_B.GetSize(), "vector sizes do not match!");

        const int number_elements( VECTOR_A.GetSize());
        const int block_size( 256);

        cl_int error_number = CL_SUCCESS;

        Vector< t_DataType> device_input_a( VECTOR_A, m_Queue);
        Vector< t_DataType> device_input_b( VECTOR_B, m_Queue);
        Vector< t_DataType> device_output( number_elements, m_Queue);

        // Core sequence... copy input data to GPU, compute, copy results back
        cl::Kernel kernel( m_Program, "GaussianDerivativeKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::NDRange local_worksize( block_size); //< all thread blocks have same dimensions
        cl::NDRange global_worksize( Tools::RoundUp( block_size, number_elements));
        cl::NDRange offset;
        error_number  = kernel.setArg( 0,  device_output.GetData());
        error_number |= kernel.setArg( 1,  device_input_a.GetData());
        error_number |= kernel.setArg( 2,  device_input_b.GetData());
        error_number |= kernel.setArg( 3,  number_elements);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

        linal::Vector< t_DataType> output_vector( device_output.GetHostVector());

        // wait for all commands to finish in queue
        m_Queue.finish();

        return output_vector;
      }

      //! @brief perform gaussian transfer function on a vector of elements
      //! @param MATRIX the matrix of elements
      //! @return matrix of gaussian results
      linal::Matrix< t_DataType> GaussianF
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX
      ) const
      {
        const int number_elements( MATRIX.GetNumberOfElements());
        const int block_size( 256);
        cl_int error_number = CL_SUCCESS;

        Matrix< t_DataType> device_input( MATRIX, m_Queue);
        Matrix< t_DataType> device_output( MATRIX.GetNumberRows(), MATRIX.GetNumberCols(), m_Queue);

        // Core sequence... copy input data to GPU, compute, copy results back
        cl::Kernel kernel( m_Program, "GaussianKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::NDRange local_worksize( block_size); //< all thread blocks have same dimensions
        cl::NDRange global_worksize( Tools::RoundUp( block_size, number_elements));
        cl::NDRange offset;
        error_number  = kernel.setArg( 0,  device_output.GetData());
        error_number |= kernel.setArg( 1,  device_input.GetData());
        error_number |= kernel.setArg( 2,  number_elements);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

        linal::Matrix< t_DataType> output_matrix( device_output.GetHostMatrix());

        // wait for all commands to finish in queue
        m_Queue.finish();

        return output_matrix;
      }

      //! @brief perform gaussian derivative
      //! @param MATRIX_A the matrix of elements - x
      //! @param MATRIX_B the matrix of elements - F( x)
      //! @return matrix of derivative results
      linal::Matrix< t_DataType> GaussiandF
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX_A,
        const linal::MatrixConstInterface< t_DataType> &MATRIX_B
      ) const
      {
        BCL_Assert( MATRIX_A.GetNumberOfElements() == MATRIX_B.GetNumberOfElements(), "matrix sizes do not match!");

        const int number_elements( MATRIX_A.GetNumberOfElements());
        const int block_size( 256);
        cl_int error_number = CL_SUCCESS;

        Matrix< t_DataType> device_input_a( MATRIX_A, m_Queue);
        Matrix< t_DataType> device_input_b( MATRIX_B, m_Queue);
        Matrix< t_DataType> device_output( MATRIX_A.GetNumberRows(), MATRIX_A.GetNumberCols(), m_Queue);

        // Core sequence... copy input data to GPU, compute, copy results back
        cl::Kernel kernel( m_Program, "GaussianDerivativeKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::NDRange local_worksize( block_size); //< all thread blocks have same dimensions
        cl::NDRange global_worksize( Tools::RoundUp( block_size, number_elements));
        cl::NDRange offset;
        error_number  = kernel.setArg( 0,  device_output.GetData());
        error_number |= kernel.setArg( 1,  device_input_a.GetData());
        error_number |= kernel.setArg( 2,  device_input_b.GetData());
        error_number |= kernel.setArg( 3,  number_elements);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

        linal::Matrix< t_DataType> output_matrix( device_output.GetHostMatrix());

        // wait for all commands to finish in queue
        m_Queue.finish();

        return output_matrix;
      }

      //! @brief updates the queue
      //! @param QUEUE the new queue
      void UpdateQueue( const CommandQueue &QUEUE)
      {
        m_Queue = QUEUE;

        cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>());
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

    }; // class TransferFunctionGaussian

    template< typename t_DataType>
    const char *TransferFunctionGaussian< t_DataType>::s_CLGaussianF =
      "                                                                                         \n"
      "/////////////////////////////////////////////////////////////////////                    \n"
      "//! Gaussian kernel                                                                      \n"
      "/////////////////////////////////////////////////////////////////////                    \n"
      "__kernel void GaussianKernel                                                             \n"
      "(                                                                                        \n"
      "  __global PRECISION *OUTPUT, __global PRECISION *INPUT, int NUMBER_ELEMENTS             \n"
      ")                                                                                        \n"
      "{                                                                                        \n"
      "  int index = get_global_id( 0);                                                         \n"
      "  if( index >= NUMBER_ELEMENTS)                                                          \n"
      "  {                                                                                      \n"
      "    return;                                                                              \n"
      "  }                                                                                      \n"
      "  PRECISION temp;                                                                        \n"
      "  temp = INPUT[ index];                                                                  \n"
      "  OUTPUT[ index] = exp( -( temp * temp) / 2.0);                                          \n"
      "}                                                                                        \n"
      ; // gaussian kernel

    template< typename t_DataType>
    const char *TransferFunctionGaussian< t_DataType>::s_CLGaussiandF =
      "                                                                                                             \n"
      "/////////////////////////////////////////////////////////////////////                                        \n"
      "//! Gaussian derivative kernel                                                                               \n"
      "/////////////////////////////////////////////////////////////////////                                        \n"
      "__kernel void GaussianDerivativeKernel                                                                       \n"
      "(                                                                                                            \n"
      "  __global PRECISION *OUTPUT, __global PRECISION *INPUT_A, __global PRECISION *INPUT_B, int NUMBER_ELEMENTS  \n"
      ")                                                                                                            \n"
      "{                                                                                                            \n"
      "  int index = get_global_id( 0);                                                                             \n"
      "  if( index >= NUMBER_ELEMENTS)                                                                              \n"
      "  {                                                                                                          \n"
      "    return;                                                                                                  \n"
      "  }                                                                                                          \n"
      "  PRECISION temp_a, temp_b;                                                                                  \n"
      "  temp_a = INPUT_A[ index];                                                                                  \n"
      "  temp_b = INPUT_B[ index];                                                                                  \n"
      "  OUTPUT[ index] = temp_a * temp_b;                                                                          \n"
      "}                                                                                                            \n"
      ; // gaussian kernel

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_TRANSFER_FUNCTION_GAUSSIAN_H_
