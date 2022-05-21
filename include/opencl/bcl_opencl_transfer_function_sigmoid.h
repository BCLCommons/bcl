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

#ifndef BCL_OPENCL_TRANSFER_FUNCTION_SIGMOID_H_
#define BCL_OPENCL_TRANSFER_FUNCTION_SIGMOID_H_

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
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_const_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransferFunctionSigmoid
    //! @brief this class performs a sigmoid transfer function using opencl
    //! @tparam data type of float, t_DataType, int, complex, etc...
    //!
    //! @see @link example_opencl_transfer_function_sigmoid.cpp @endlink
    //! @author loweew
    //! @date Apr 4, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class TransferFunctionSigmoid :
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
      static const char* s_CLSigmoidF;
      static const char* s_CLSigmoiddF;
      static const cl_uint s_BlockSize = 16;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param QUEUE the command queue
      TransferFunctionSigmoid()
      {
      }

      //! @brief default constructor
      //! @param QUEUE the command queue
      TransferFunctionSigmoid( const CommandQueue &QUEUE) :
        m_Queue( QUEUE),
        m_Program()
      {
        cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>());
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

      //! @brief virtual copy constructor
      //! @return pointer to a new TransferFunctionSigmoid copied from this one
      TransferFunctionSigmoid *Clone() const
      {
        return new TransferFunctionSigmoid( *this);
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

      //! @brief performs sigmoid of all elements in matrix
      //! @param INPUT matrix for which to perform sigmoid
      //! @return the resulting matrix
      void F( const Matrix< t_DataType> &INPUT, Matrix< t_DataType> &OUTPUT) const
      {
        SigmoidF( INPUT, OUTPUT);
      }

      //! @brief derivative of sigmoid
      //! @param INPUT the input matrix for which to perform the derivative of sigmoid
      //! @return the resulting matrix
      void dF( const Matrix< t_DataType> &INPUT, Matrix< t_DataType> &OUTPUT) const
      {
        SigmoiddF( INPUT, OUTPUT);
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
        static const KernelSource e_sigmoid_kernel
        (
          GetKernelSources().AddEnum
          (
            "TransferFunctionSigmoidF",
            util::ShPtr< KernelSourceInterface>
            (
              new KernelSourceString( s_CLSigmoidF)
            )
          )
        );
        return e_sigmoid_kernel;
      }

      //! @brief return the kernel source for dF( x)
      //! @return kernel source of dF( x)
      const KernelSource &GetKerneldF()
      {
        static const KernelSource e_sigmoid_derivative_kernel
        (
          GetKernelSources().AddEnum
          (
            "TransferFunctionSigmoiddF",
            util::ShPtr< KernelSourceInterface>
            (
              new KernelSourceString( s_CLSigmoiddF)
            )
          )
        );
        return e_sigmoid_derivative_kernel;
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
        const std::string sigmoid_source_F( ( *GetKernelF())->GetSource( PRECISION, device.Extensions()));
        if
        (
          sigmoid_source_F.empty()
        )
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        const std::string sigmoid_source_dF( ( *GetKerneldF())->GetSource( PRECISION, device.Extensions()));
        if
        (
          sigmoid_source_dF.empty()
        )
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        source.push_back( std::make_pair( sigmoid_source_F.c_str()      , sigmoid_source_F.length()));
        source.push_back( std::make_pair( sigmoid_source_dF.c_str()      , sigmoid_source_dF.length()));

        // create the program
        cl::Program &current_program( m_Program);
        current_program = cl::Program( m_Queue.GetContext(), source, &error_number);
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
          return error_number;
        }

        // end
        return error_number;
      }

      //! @brief perform sigmoid transfer function on a matrix
      //! @param INPUT the matrix of elements
      //! @param OUTPUT matrix in which to store the sigmoid results
      void SigmoidF
      (
        const Matrix< t_DataType> &INPUT,
              Matrix< t_DataType> &OUTPUT
      ) const
      {
        const cl_uint rows( INPUT.GetNumberRows());
        const cl_uint cols( INPUT.GetNumberCols());
        const cl_uint col_pad( INPUT.GetColPadding());
        const cl_uint row_pad( INPUT.GetRowPadding());

        cl_int error_number = CL_SUCCESS;

        cl::Kernel kernel( m_Program, "SigmoidKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        cl::NDRange block_dims( s_BlockSize, s_BlockSize); //< all thread blocks have same dimensions
        cl::NDRange worksize( Tools::RoundUp( s_BlockSize, cols), Tools::RoundUp( s_BlockSize, rows));
        cl::NDRange offset;

        error_number  = kernel.setArg( 0,  OUTPUT.GetData());
        error_number |= kernel.setArg( 1,  INPUT.GetData());
        error_number |= kernel.setArg( 2,  rows);
        error_number |= kernel.setArg( 3,  cols);
        error_number |= kernel.setArg( 4,  col_pad);
        error_number |= kernel.setArg( 5,  row_pad);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));
      }

      //! @brief perform sigmoid derivative
      //! @param INPUT the matrix of elements - F( x)
      //! @param OUTPUT matrix in which to store the sigmoid results
      //! @param OUTPUT matrix in which to store the sigmoid derivative results
      void SigmoiddF
      (
        const Matrix< t_DataType> &INPUT,
              Matrix< t_DataType> &OUTPUT
      ) const
      {
        const cl_uint rows( INPUT.GetNumberRows());
        const cl_uint cols( INPUT.GetNumberCols());
        const cl_uint col_pad( INPUT.GetColPadding());
        const cl_uint row_pad( INPUT.GetRowPadding());

        cl_int error_number = CL_SUCCESS;

        cl::Kernel kernel( m_Program, "SigmoidDerivativeKernel", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        cl::NDRange block_dims( s_BlockSize, s_BlockSize); //< all thread blocks have same dimensions
        cl::NDRange worksize( Tools::RoundUp( s_BlockSize, cols), Tools::RoundUp( s_BlockSize, rows));
        cl::NDRange offset;

        error_number  = kernel.setArg( 0, OUTPUT.GetData());
        error_number |= kernel.setArg( 1, INPUT.GetData());
        error_number |= kernel.setArg( 2, rows);
        error_number |= kernel.setArg( 3, cols);
        error_number |= kernel.setArg( 4, col_pad);
        error_number |= kernel.setArg( 5, row_pad);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        // launching kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));
      }

    }; //class TransferFunctionSigmoid

      template< typename t_DataType>
      const char *TransferFunctionSigmoid< t_DataType>::s_CLSigmoidF =
        "//! @brief calculates sigmoid of all elements in matrix except for those added as padding        \n"
        "//! @param  OUTPUT sigmoid output matrix                                                         \n"
        "//! @param  INPUT  input matrix                                                                  \n"
        "//! @param  ROWS   rows in input matrix                                                          \n"
        "//! @param  COLS   cols in input matrix                                                          \n"
        "//! @param  XPAD   padding added to cols                                                         \n"
        "//! @param  YPAD   padding added to rows                                                         \n"
        "__kernel void SigmoidKernel                                                                      \n"
        "(                                                                                                \n"
        "          __global PRECISION *OUTPUT,                                                            \n"
        "  __const __global PRECISION *INPUT,                                                             \n"
        "  __const          uint        ROWS,                                                             \n"
        "  __const          uint        COLS,                                                             \n"
        "  __const          uint        XPAD,                                                             \n"
        "  __const          uint        YPAD                                                              \n"
        ")                                                                                                \n"
        "{                                                                                                \n"
        "  const uint block_size = get_local_size( 0);                                                    \n"
        "  uint xindex = get_global_id( 0);                                                               \n"
        "  uint yindex = get_global_id( 1);                                                               \n"
        "  uint size_x = COLS - XPAD;                                                                     \n"
        "  uint size_y = ROWS - YPAD;                                                                     \n"
        "  uint index = yindex * COLS + xindex;                                                           \n"
        "  if( xindex < size_x && yindex < size_y)                                                        \n"
        "  {                                                                                              \n"
        "    PRECISION temp;                                                                              \n"
        "    temp = INPUT[ index];                                                                        \n"
        "    OUTPUT[ index] = native_recip( 1 + native_exp( -temp));                                      \n"
        "  }                                                                                              \n"
        "}"
        ; // sigmoid kernel

      template< typename t_DataType>
      const char *TransferFunctionSigmoid< t_DataType>::s_CLSigmoiddF =
        "//! @brief calculates sigmoid derivative of all elements in matrix except for padding            \n"
        "//! @param  OUTPUT sigmoid output matrix                                                         \n"
        "//! @param  INPUT  input matrix                                                                  \n"
        "//! @param  ROWS   rows in input matrix                                                          \n"
        "//! @param  COLS   cols in input matrix                                                          \n"
        "//! @param  XPAD   padding added to cols                                                         \n"
        "//! @param  YPAD   padding added to rows                                                         \n"
        "__kernel void SigmoidDerivativeKernel                                                            \n"
        "(                                                                                                \n"
        "          __global PRECISION *OUTPUT,                                                            \n"
        "  __const __global PRECISION *INPUT,                                                             \n"
        "  __const          uint        ROWS,                                                             \n"
        "  __const          uint        COLS,                                                             \n"
        "  __const          uint        XPAD,                                                             \n"
        "  __const          uint        YPAD                                                              \n"
        ")                                                                                                \n"
        "{                                                                                                \n"
        "  const uint block_size = get_local_size( 0);                                                    \n"
        "  uint xindex = get_global_id( 0);                                                               \n"
        "  uint yindex = get_global_id( 1);                                                               \n"
        "  uint size_x = COLS - XPAD;                                                                     \n"
        "  uint size_y = ROWS - YPAD;                                                                     \n"
        "  uint index = yindex * COLS + xindex;                                                           \n"
        "  if( xindex < size_x && yindex < size_y)                                                        \n"
        "  {                                                                                              \n"
        "    PRECISION temp;                                                                              \n"
        "    temp = INPUT[ index];                                                                        \n"
        "    OUTPUT[ index] = temp * ( 1 - temp);                                                         \n"
        "  }                                                                                              \n"
        "}"
        ; // sigmoid kernel

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_TRANSFER_FUNCTION_SIGMOID_H_
