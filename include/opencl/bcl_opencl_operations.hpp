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

#ifndef BCL_OPENCL_OPERATIONS_HPP_
#define BCL_OPENCL_OPERATIONS_HPP_

// include header of this class
#include "bcl_opencl_operations.h"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_matrix_transpose.h"
#include "bcl_opencl_tools.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_operations.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_const_interface.h"
#include "math/bcl_math.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OperationsEnumHandler
    //! @brief handler class for adding the operations enum handler
    //! @author woetzen, loweew
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class OperationsEnumHandler :
      public signal::Slots
    {

    public:

    //////////
    // data //
    //////////

      //! ShPtr to the instance of the Operations implementation
      util::ShPtr< Operations< t_DataType> > m_Instance;

      //! the enum in the linal::OperationsInterface
      typename linal::Operations< t_DataType>::EnumType e_Operations;

      //! the only instance of this class
      static const OperationsEnumHandler< t_DataType> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      OperationsEnumHandler();

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS);

    }; // template class OperationsEnumHandler

    //! instance of OperationsEnumHandler
    template< typename t_DataType>
    const OperationsEnumHandler< t_DataType> OperationsEnumHandler< t_DataType>::s_Instance = OperationsEnumHandler< t_DataType>();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    OperationsEnumHandler< t_DataType>::OperationsEnumHandler() :
      m_Instance(),
      e_Operations( linal::Operations< t_DataType>::GetEnums().e_Undefined)
    {
      GetTools().GetQueueUpdateSignal().Connect( this, &OperationsEnumHandler< t_DataType>::UpdateEnum);
    }

    //! @brief update the enum with the command queue from the Tools
    //! @param TOOLS the tolls to get the commandqueue from
    template< typename t_DataType>
    void OperationsEnumHandler< t_DataType>::UpdateEnum( Tools &TOOLS)
    {
      util::ShPtr< Operations< t_DataType> > sp_operations( new Operations< t_DataType>());
      if( !TOOLS.HasCommandQueues())
      {
        m_Instance = util::ShPtr< Operations< t_DataType> >();
        *e_Operations = m_Instance;
        return;
      }

      // try to initialize
      if( sp_operations->Initialize( TOOLS.GetFirstCommandQueue()))
      {
        if( e_Operations.IsDefined())
        {
          // just update the existing one with the new one
          m_Instance = sp_operations;
          *e_Operations = m_Instance;
        }
        else
        {
          m_Instance = sp_operations;
          e_Operations = linal::Operations< t_DataType>::GetEnums().AddEnum( "OpenCL", m_Instance);
        }
      }
      else
      {
        BCL_MessageVrb( "unable to initialize enum: OpenCL");
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Operations< t_DataType>::Operations()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Operations
    template< typename t_DataType>
    Operations< t_DataType> *Operations< t_DataType>::Clone() const
    {
      return new Operations< t_DataType>( *this);
    }

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    template< typename t_DataType>
    bool Operations< t_DataType>::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      cl_int error_number( CL_SUCCESS);
      const Device device( COMMAND_QUEUE.GetDevice( &error_number));

      // can get device
      if( error_number != CL_SUCCESS)
      {
        return false;
      }

      const storage::Set< Extension> extensions( device.Extensions( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "unable to get extensions from device");

      return KernelSourceInterface::PrecisionCompatibleWithExtensions
             (
               util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(),
               extensions
             );
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    template< typename t_DataType>
    bool Operations< t_DataType>::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      // check if this is a compatible command queue
      if( !IsCompatible( COMMAND_QUEUE))
      {
        return false;
      }

      cl_int error_number( CL_SUCCESS);

      // update the command queue
      m_CommandQueue = COMMAND_QUEUE;

      // update the context
      m_Context = m_CommandQueue.GetContext( &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "cannot get context for command queue");

      // for precision type
      error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>());
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageDbg( "error compiling programs:\n" + Tools::ErrorString( error_number));
        return false;
      }

      return true;
    }

    //! @brief access to instance
    //! @return reference to instance of this class
    template< typename t_DataType>
    const Operations< t_DataType> &Operations< t_DataType>::GetInstance()
    {
      return *OperationsEnumHandler< t_DataType>::s_Instance.m_Instance;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Operations< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the program
    //! @return reference to the program
    template< typename t_DataType>
    const cl::Program &Operations< t_DataType>::GetProgram() const
    {
      return m_Program;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief dot product of two vectors
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A * VECTOR_B
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::DotProduct
    (
      const linal::VectorConstInterface< t_DataType> &VECTOR_A,
      const linal::VectorConstInterface< t_DataType> &VECTOR_B
    ) const
    {
      BCL_Assert( VECTOR_A.GetSize() == VECTOR_B.GetSize(), "non-matching dimensions!");

      // OpenCL Vars
      // set and log Global and Local work size dimensions
      const size_t szLocalWorkSize( 256); // # of work items in the 1D work group
      const size_t szGlobalWorkSize( Tools::RoundUp( szLocalWorkSize, VECTOR_A.GetSize()));      // Total # of work items in the 1D range
      const size_t nr_groups( szGlobalWorkSize / szLocalWorkSize);
      cl_int error_number( CL_SUCCESS);              // Error code var

      // Allocate the OpenCL buffer memory objects for source and result on the device GMEM
      Buffer cmDevSrcA( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * szGlobalWorkSize, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      Buffer cmDevSrcB( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * szGlobalWorkSize, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      Buffer cmDevDst( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * nr_groups, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "VectorDotProduct", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // initialize result on host
      linal::Vector< t_DataType> result( nr_groups, t_DataType( 0));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( cmDevSrcA, CL_FALSE, 0, sizeof( t_DataType) * szGlobalWorkSize, ( void*)VECTOR_A.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmDevSrcB, CL_FALSE, 0, sizeof( t_DataType) * szGlobalWorkSize, ( void*)VECTOR_B.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Launch kernel
      error_number  = kernel.setArg( 0, cmDevSrcA);
      error_number |= kernel.setArg( 1, cmDevSrcB);
      error_number |= kernel.setArg( 2, cl_uint( VECTOR_A.GetSize()));
      error_number |= kernel.setArg( 3, cmDevDst);
      error_number |= kernel.setArg( 4, szLocalWorkSize * sizeof( t_DataType), NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, cl::NDRange(), cl::NDRange( Tools::RoundUp( szLocalWorkSize, VECTOR_A.GetSize())), cl::NDRange( szLocalWorkSize), NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Read back results and check error
      error_number = m_CommandQueue.enqueueReadBuffer( cmDevDst, CL_TRUE, 0, sizeof( t_DataType) * result.GetSize(), ( void *)result.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return result.Sum();
    }

    //! @brief norm of a vector
    //! @param VECTOR vector
    //! @return the euclidean norm of the VECTOR
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Norm( const linal::VectorConstInterface< t_DataType> &VECTOR) const
    {
      // OpenCL Vars
      const size_t szLocalWorkSize( 256); // # of work items in the 1D work group
      const size_t szGlobalWorkSize( Tools::RoundUp( szLocalWorkSize, VECTOR.GetSize()));
      const size_t nr_groups( szGlobalWorkSize / szLocalWorkSize);
      cl_int error_number( CL_SUCCESS);          // Error code var

      // Allocate the OpenCL buffer memory objects for source and result on the device GMEM
      Buffer cmDevSrc( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * szGlobalWorkSize, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      Buffer cmDevDst( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * nr_groups, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "VectorInnerProduct", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // initialize result on host
      linal::Vector< t_DataType> result( nr_groups, t_DataType( 0));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( cmDevSrc, CL_FALSE, 0, sizeof( t_DataType) * szGlobalWorkSize, ( void *)VECTOR.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Launch kernel
      error_number  = kernel.setArg( 0, cmDevSrc);
      error_number |= kernel.setArg( 1, cl_uint( VECTOR.GetSize()));
      error_number |= kernel.setArg( 2, cmDevDst);
      error_number |= kernel.setArg( 3, szLocalWorkSize * sizeof( t_DataType), NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, cl::NDRange(), cl::NDRange( szGlobalWorkSize), cl::NDRange( szLocalWorkSize));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Read back results and check error
      error_number = m_CommandQueue.enqueueReadBuffer( cmDevDst, CL_TRUE, 0, sizeof( t_DataType) * result.GetSize(), ( void *)result.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return math::Sqrt( result.Sum());
    }

    //! @brief matrix-vector multiplication
    //! @param MATRIX matrix to be multiplied
    //! @param VECTOR vector to be multiplied
    //! @return resulting linal::Vector< t_DataType>
    template< typename t_DataType>
    linal::Vector< t_DataType> Operations< t_DataType>::Multiply
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const linal::VectorConstInterface< t_DataType> &VECTOR
    ) const
    {
      BCL_Assert( MATRIX.GetNumberCols() == VECTOR.GetSize(), "non-matching dimensions!");

      cl_int error_number( CL_SUCCESS);              // Error code var

      // create buffer for matrix
      Buffer cmM( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * MATRIX.GetNumberOfElements(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      // creat buffer for vector
      Buffer cmV( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * VECTOR.GetSize(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Vector< t_DataType> result( MATRIX.GetNumberRows(), t_DataType( 0.0));
      Buffer cmDevDst( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * result.GetSize(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "MatrixVectorMultiplication", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( cmM, CL_FALSE, 0, sizeof( t_DataType) * MATRIX.GetNumberOfElements(), ( void *)MATRIX.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmV, CL_TRUE, 0, sizeof( t_DataType) * VECTOR.GetSize(), ( void *)VECTOR.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      int block_size( 256);
      cl::NDRange block_dimensions( block_size); // # of work items in the 1D work group
      cl::NDRange kernel_worksize( Tools::RoundUp( block_size, int( MATRIX.GetNumberRows()))); // Total # of work items in the 1D range
      cl::NDRange offset;

      // Launch kernel
      error_number =  kernel.setArg( 0, cmM);
      error_number |= kernel.setArg( 1, cmV);
      error_number |= kernel.setArg( 2, cl_uint( MATRIX.GetNumberCols()));
      error_number |= kernel.setArg( 3, cl_uint( MATRIX.GetNumberRows()));
      error_number |= kernel.setArg( 4, cmDevDst);
      error_number |= kernel.setArg( 5, sizeof( t_DataType) * block_size, 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Read back results and check error
      error_number = m_CommandQueue.enqueueReadBuffer( cmDevDst, CL_TRUE, 0, sizeof( t_DataType) * int( result.GetSize()), ( void *)result.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return result;
    }

    //! @brief performs matrix-matrix multiplication using level 3 cublas
    //! @param MATRIX_A the first matrix in the multiplication
    //! @param MATRIX_B the second matrix in the multiplication
    //! @return the resulting linal::Matrix< t_DataType>
    template< typename t_DataType>
    linal::Matrix< t_DataType> Operations< t_DataType>::Multiply
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX_A,
      const linal::MatrixConstInterface< t_DataType> &MATRIX_B
    ) const
    {
      const size_t block_size( 16);
      const size_t pad_a_row( ( block_size - ( MATRIX_A.GetNumberRows() % block_size)) % block_size);
      const size_t pad_a_col( ( block_size - ( MATRIX_A.GetNumberCols() % block_size)) % block_size);
      const size_t pad_b_row( ( block_size - ( MATRIX_B.GetNumberRows() % block_size)) % block_size);
      const size_t pad_b_col( ( block_size - ( MATRIX_B.GetNumberCols() % block_size)) % block_size);

      // create padded matrix a and b
      linal::Matrix< t_DataType> matrix_pad_a( MATRIX_A, pad_a_row, pad_a_col);
      linal::Matrix< t_DataType> matrix_pad_b( MATRIX_B, pad_b_row, pad_b_col);

      cl_int error_number( CL_SUCCESS);

      // Allocate the OpenCL buffer memory objects for source and result on the device GMEM and copy data
      Buffer cmDevSrcA( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * matrix_pad_a.GetNumberOfElements(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error creating buffer a: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmDevSrcA, CL_FALSE, 0, sizeof( t_DataType) * matrix_pad_a.GetNumberOfElements(), matrix_pad_a.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error writing buffer a: " + Tools::ErrorString( error_number));

      Buffer cmDevSrcB( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * matrix_pad_b.GetNumberOfElements(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error creating buffer b: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmDevSrcB, CL_FALSE, 0, sizeof( t_DataType) * matrix_pad_b.GetNumberOfElements(), matrix_pad_b.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error writing buffer b: " + Tools::ErrorString( error_number));

      // create kernel
      cl::Kernel kernel( m_Program, "MatrixMultiplication", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Matrix for the result
      linal::Matrix< t_DataType> result( matrix_pad_a.GetNumberRows(), matrix_pad_b.GetNumberCols(), 0.0);

      const cl_uint number_cols_a( matrix_pad_a.GetNumberCols());
      const cl_uint number_cols_b( matrix_pad_b.GetNumberCols());

      // Output buffer
      Buffer buffer_c( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, result.GetNumberOfElements() * sizeof( t_DataType), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, cmDevSrcA);
      error_number |= kernel.setArg( 1, cmDevSrcB);
      error_number |= kernel.setArg( 2, number_cols_a);
      error_number |= kernel.setArg( 3, number_cols_b);
      error_number |= kernel.setArg( 4, buffer_c);
      error_number |= kernel.setArg( 5, block_size * block_size * sizeof( t_DataType), NULL); //shared memory
      error_number |= kernel.setArg( 6, block_size * block_size * sizeof( t_DataType), NULL); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange localWorkSize( block_size, block_size);
      const cl::NDRange globalOffset;
      const cl::NDRange globalWorkSize( Tools::RoundUp( block_size, result.GetNumberCols()), Tools::RoundUp( block_size, result.GetNumberRows()));

      // Multiplication - non-blocking execution
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, globalOffset, globalWorkSize, localWorkSize, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // blocking copy of result from device to host
      error_number = m_CommandQueue.enqueueReadBuffer
      (
        buffer_c, CL_TRUE, 0, result.GetNumberOfElements() * sizeof( t_DataType),
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // return result without padding
      return result.CreateSubMatrix( MATRIX_A.GetNumberRows(), MATRIX_B.GetNumberCols());
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS list of vectors
    //! @return triangular matrix of values of the distances
    template< typename t_DataType>
    linal::Matrix< t_DataType> Operations< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS
    ) const
    {
      const int num_vectors( LIST_VECTORS.GetSize());
      const int vector_size( LIST_VECTORS.FirstElement()->GetSize());
      const int block_size( 16);

      linal::Matrix< t_DataType> input_matrix( num_vectors, vector_size, t_DataType( 0.0));
      t_DataType *itr( input_matrix.Begin());

      for( int ctr( 0); ctr < num_vectors; ++ctr)
      {
        const t_DataType *vec_itr( LIST_VECTORS( ctr)->Begin());
        for( int vec( 0); vec < vector_size; ++vec, ++itr, ++vec_itr)
        {
          ( *itr) = ( *vec_itr);
        }
      }

      const size_t pad_a_row( ( block_size - ( input_matrix.GetNumberRows() % block_size)) % block_size);
      const size_t pad_a_col( ( block_size - ( input_matrix.GetNumberCols() % block_size)) % block_size);

      // create padded matrix
      linal::Matrix< t_DataType> matrix_padded( input_matrix, pad_a_row, pad_a_col);
      int num_rows_padded( matrix_padded.GetNumberRows());
      int num_cols_padded( matrix_padded.GetNumberCols());

      cl_int error_number( CL_SUCCESS);              // Error code var
      // create buffer for matrix
      Buffer device_input_matrix( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * num_rows_padded * num_cols_padded, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Matrix< t_DataType> result( num_rows_padded, num_rows_padded, t_DataType( 0.0));
      Buffer device_result( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * num_rows_padded * num_rows_padded, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "EuclideanDistance", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( device_input_matrix, CL_FALSE, 0, sizeof( t_DataType) * num_rows_padded * num_cols_padded, matrix_padded.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size, block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, num_rows_padded), Tools::RoundUp( block_size, num_rows_padded));

      error_number =  kernel.setArg( 0, device_result);
      error_number |= kernel.setArg( 1, device_input_matrix);
      error_number |= kernel.setArg( 2, device_input_matrix);
      error_number |= kernel.setArg( 3, sizeof( t_DataType) * block_size * block_size, 0);
      error_number |= kernel.setArg( 4, sizeof( t_DataType) * block_size * block_size, 0);
      error_number |= kernel.setArg( 5, num_rows_padded);
      error_number |= kernel.setArg( 6, num_cols_padded);
      error_number |= kernel.setArg( 7, num_rows_padded);
      error_number |= kernel.setArg( 8, num_cols_padded);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueReadBuffer
      (
        device_result, CL_TRUE, 0, sizeof( t_DataType) * num_rows_padded * num_rows_padded,
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      return result.CreateSubMatrix( num_vectors, num_vectors);
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS_A list of vectors
    //! @param LIST_VECTORS_B list of vectors
    //! @return matrix of values of the distances
    template< typename t_DataType>
    linal::Matrix< t_DataType> Operations< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
      const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS_B
    ) const
    {
      const int num_vectors_a( LIST_VECTORS_A.GetSize());
      const int vector_size_a( LIST_VECTORS_A.FirstElement()->GetSize());
      const int num_vectors_b( LIST_VECTORS_B.GetSize());
      const int vector_size_b( LIST_VECTORS_B.FirstElement()->GetSize());

      linal::Matrix< t_DataType> input_matrix_a( num_vectors_a, vector_size_a, t_DataType( 0.0));
      t_DataType *itr_a( input_matrix_a.Begin());

      for( int ctr( 0); ctr < num_vectors_a; ++ctr)
      {
        const t_DataType *vec_itr( LIST_VECTORS_A( ctr)->Begin());
        for( int vec( 0); vec < vector_size_a; ++vec, ++itr_a, ++vec_itr)
        {
          ( *itr_a) = ( *vec_itr);
        }
      }

      linal::Matrix< t_DataType> input_matrix_b( num_vectors_b, vector_size_b, t_DataType( 0.0));
      t_DataType *itr_b( input_matrix_b.Begin());

      for( int ctr( 0); ctr < num_vectors_b; ++ctr)
      {
        const t_DataType *vec_itr( LIST_VECTORS_B( ctr)->Begin());
        for( int vec( 0); vec < vector_size_b; ++vec, ++itr_b, ++vec_itr)
        {
          ( *itr_b) = ( *vec_itr);
        }
      }

      const size_t block_size( 16);
      const size_t pad_a_row( ( block_size - ( input_matrix_a.GetNumberRows() % block_size)) % block_size);
      const size_t pad_a_col( ( block_size - ( input_matrix_a.GetNumberCols() % block_size)) % block_size);
      const size_t pad_b_row( ( block_size - ( input_matrix_b.GetNumberRows() % block_size)) % block_size);
      const size_t pad_b_col( ( block_size - ( input_matrix_b.GetNumberCols() % block_size)) % block_size);

      // create padded matrix a and b
      linal::Matrix< t_DataType> matrix_pad_a( input_matrix_a, pad_a_row, pad_a_col);
      linal::Matrix< t_DataType> matrix_pad_b( input_matrix_b, pad_b_row, pad_b_col);
      int num_rows_a( matrix_pad_a.GetNumberRows());
      int num_cols_a( matrix_pad_a.GetNumberCols());
      int num_rows_b( matrix_pad_b.GetNumberRows());
      int num_cols_b( matrix_pad_b.GetNumberCols());

      cl_int error_number( CL_SUCCESS);              // Error code var
      // create buffer for matrix a
      Buffer device_input_matrix_a( cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof( t_DataType) * num_rows_a * num_cols_a, matrix_pad_a.Begin(), &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create buffer for matrix b
      Buffer device_input_matrix_b( cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof( t_DataType) * num_rows_b * num_cols_b, matrix_pad_b.Begin(), &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Matrix< t_DataType> result( num_rows_a, num_rows_b, t_DataType( 0.0));
      Buffer device_result( cl::Buffer( m_Context, CL_MEM_READ_WRITE, sizeof( t_DataType) * num_rows_a * num_rows_b, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "EuclideanDistance", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size, block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, num_rows_a), Tools::RoundUp( block_size, num_rows_b));

      error_number =  kernel.setArg( 0, device_result);
      error_number |= kernel.setArg( 1, device_input_matrix_a);
      error_number |= kernel.setArg( 2, device_input_matrix_b);
      error_number |= kernel.setArg( 3, sizeof( t_DataType) * block_size * block_size, 0);
      error_number |= kernel.setArg( 4, sizeof( t_DataType) * block_size * block_size, 0);
      error_number |= kernel.setArg( 5, num_rows_a);
      error_number |= kernel.setArg( 6, num_cols_a);
      error_number |= kernel.setArg( 7, num_rows_b);
      error_number |= kernel.setArg( 8, num_cols_b);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueReadBuffer
      (
        device_result, CL_TRUE, 0, sizeof( t_DataType) * num_rows_a * num_rows_b,
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      return result.CreateSubMatrix( num_vectors_a, num_vectors_b);
    }

    //! @brief reduction sum kernel
    //! @param VECTOR the vector to get the sum of
    //! @return the reduced result sum
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Sum( const linal::VectorConstInterface< t_DataType> &VECTOR) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      // create buffer for matrix
      const size_t block_size( 256);
      const size_t padding( ( block_size - ( VECTOR.GetSize() % block_size)) % block_size);
      linal::Vector< t_DataType> vector_padded( VECTOR, padding);
      const size_t vector_size( vector_padded.GetSize());
      const size_t num_groups( vector_size / block_size);
      Buffer device_vector( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * vector_size, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Vector< t_DataType> result( num_groups, t_DataType( 0.0));
      Buffer device_result( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * num_groups, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "ReductionSum", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( device_vector, CL_TRUE, 0, sizeof( t_DataType) * vector_size, ( void *)vector_padded.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, vector_size));

      error_number  = kernel.setArg( 0, device_vector);
      error_number |= kernel.setArg( 1, cl_uint( vector_size));
      error_number |= kernel.setArg( 2, device_result);
      error_number |= kernel.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueReadBuffer
      (
        device_result, CL_TRUE, 0, sizeof( t_DataType) * num_groups,
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // return result without padding
      return result.Sum();
    }

    //! @brief reduction sum kernel
    //! @param VECTOR the vector to get the sum of
    //! @return the reduced result sum
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Sum( const Vector< t_DataType> &VECTOR) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var

      // create buffer for result
      const size_t block_size( 256);
      const size_t vector_size( VECTOR.GetSize());
      const size_t num_groups( Tools::RoundUp( block_size, vector_size) / block_size);

      // create the result buffer
      Vector< t_DataType> result( num_groups, m_CommandQueue);

      // Create the kernel
      cl::Kernel kernel( m_Program, "ReductionSum", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, vector_size));

      error_number  = kernel.setArg( 0, VECTOR.GetData());
      error_number |= kernel.setArg( 1, vector_size);
      error_number |= kernel.setArg( 2, result.GetData());
      error_number |= kernel.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // return result as sum of the group results
      return result.GetHostVector().Sum();
    }

    //! @brief reduction min kernel
    //! @param VECTOR the vector to get the min of
    //! @return the reduced resulting min
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Min( const linal::VectorConstInterface< t_DataType> &VECTOR) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      // create buffer for matrix
      int block_size( 256);
      const size_t padding( ( block_size - ( VECTOR.GetSize() % block_size)) % block_size);
      linal::Vector< t_DataType> vector_padded( VECTOR, padding);
      const int vector_size( VECTOR.GetSize()); //( vector_padded.GetSize());
      const int num_groups( ( vector_size % block_size == 0 ? 0 : 1) + vector_size / block_size);
      Buffer device_vector( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * vector_size, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Vector< t_DataType> result( num_groups, t_DataType( 0.0));
      Buffer device_result( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * num_groups, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "ReductionMin", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( device_vector, CL_TRUE, 0, sizeof( t_DataType) * vector_size, ( void *)VECTOR.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, vector_size));

      error_number =  kernel.setArg( 0, device_vector);
      error_number |= kernel.setArg( 1, device_result);
      error_number |= kernel.setArg( 2, vector_size);
      error_number |= kernel.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueReadBuffer
      (
        device_result, CL_TRUE, 0, sizeof( t_DataType) * num_groups,
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // return result without padding
      return *std::min_element( result.Begin(), result.End());
    }

    //! @brief reduction max kernel
    //! @param VECTOR the vector to get the max of
    //! @return the reduced resulting max
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Max( const linal::VectorConstInterface< t_DataType> &VECTOR) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      // create buffer for matrix
      int block_size( 256);
      const size_t padding( ( block_size - ( VECTOR.GetSize() % block_size)) % block_size);
      linal::Vector< t_DataType> vector_padded( VECTOR, padding);
      const int vector_size( vector_padded.GetSize());
      const int num_groups( vector_size / block_size);
      Buffer device_vector( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * vector_size, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // create the result buffer
      linal::Vector< t_DataType> result( num_groups, t_DataType( 0.0));
      Buffer device_result( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * num_groups, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "ReductionMax", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( device_vector, CL_TRUE, 0, sizeof( t_DataType) * vector_size, ( void *)vector_padded.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, vector_size));

      error_number =  kernel.setArg( 0, device_vector);
      error_number |= kernel.setArg( 1, device_result);
      error_number |= kernel.setArg( 2, vector_size);
      error_number |= kernel.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueReadBuffer
      (
        device_result, CL_TRUE, 0, sizeof( t_DataType) * num_groups,
        result.Begin(), NULL, NULL
      );
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // return result without padding
      return *std::max_element( result.Begin(), result.End());
    }

    //! @brief gets the min and max for each column in a matrix
    //! @param MATRIX the matrix input
    //! @return a storage vector of math ranges with min and max for each column in the matrix
    template< typename t_DataType>
    storage::Vector< math::Range< t_DataType> > Operations< t_DataType>::MinMax( const linal::MatrixConstInterface< t_DataType> &MATRIX) const
    {
      if( MATRIX.GetNumberOfElements() == 0)
      {
        return storage::Vector< math::Range< t_DataType> >();
      }

      const size_t number_data_pts( MATRIX.GetNumberRows());
      const size_t nr_cols( MATRIX.GetNumberCols());
      storage::Vector< math::Range< t_DataType> > ranges;

      ranges.AllocateMemory( nr_cols);
      cl_int error_number( CL_SUCCESS);              // Error code var

      // create buffer for matrix
      const size_t block_size( 256);
      const size_t num_groups( ( number_data_pts / block_size) + 1);

      linal::Matrix< t_DataType> trans( MATRIX.Transposed());

      Vector< t_DataType> dev_matrix( trans.GetNumberCols(), m_CommandQueue);
      // create the result buffer
      linal::Vector< t_DataType> result( num_groups, t_DataType( 0.0));
      Vector< t_DataType> dev_result_min( result, m_CommandQueue);
      Vector< t_DataType> dev_result_max( result, m_CommandQueue);

      // Create the kernel
      cl::Kernel kernel_max( m_Program, "ReductionMax", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      cl::Kernel kernel_min( m_Program, "ReductionMin", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      const cl::NDRange block_dimensions( block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, number_data_pts));

      error_number |= kernel_max.setArg( 1, dev_result_max.GetData());
      error_number |= kernel_max.setArg( 2, cl_uint( number_data_pts));
      error_number |= kernel_max.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error max kernel args: " + Tools::ErrorString( error_number));

      error_number |= kernel_min.setArg( 1, dev_result_min.GetData());
      error_number |= kernel_min.setArg( 2, cl_uint( number_data_pts));
      error_number |= kernel_min.setArg( 3, block_size * sizeof( t_DataType), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error min kernel args: " + Tools::ErrorString( error_number));

      for( size_t col( 0), number_features( nr_cols); col < number_features; ++col)
      {
        dev_matrix = Vector< t_DataType>( number_data_pts, trans[ col], m_CommandQueue);

        error_number = kernel_max.setArg( 0, dev_matrix.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = kernel_min.setArg( 0, dev_matrix.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_CommandQueue.enqueueNDRangeKernel( kernel_min, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_CommandQueue.enqueueNDRangeKernel( kernel_max, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        ranges.PushBack( math::Range< t_DataType>( dev_result_min.GetHostVector().Min(), dev_result_max.GetHostVector().Max()));
      }

      return ranges;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void Operations< t_DataType>::VectorEqualsVectorTimesMatrix
    (
      linal::VectorInterface< t_DataType> &STORAGE,
      const linal::VectorConstInterface< t_DataType> &FEATURE,
      const linal::MatrixConstInterface< t_DataType> &MATRIX
    ) const
    {
      BCL_Assert( MATRIX.GetNumberRows() == FEATURE.GetSize() && MATRIX.GetNumberCols() == STORAGE.GetSize(), "non-matching dimensions!");

      cl_int error_number( CL_SUCCESS);              // Error code var

      Matrix< t_DataType> cmM( MATRIX, m_CommandQueue);
      Vector< t_DataType> cmV( FEATURE, m_CommandQueue);
      Vector< t_DataType> cmDevDst( STORAGE, m_CommandQueue);

      // Create the kernel
      cl::Kernel kernel( m_Program, "MatrixVectorMultiplication", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      MatrixTranspose< t_DataType> trans( m_CommandQueue);
      cmM = trans( cmM);

      int block_size( 256);
      cl::NDRange block_dimensions( block_size); // # of work items in the 1D work group
      cl::NDRange kernel_worksize( Tools::RoundUp( block_size, int( cmM.GetNumberRows()))); // Total # of work items in the 1D range
      cl::NDRange offset;

      // Launch kernel
      error_number =  kernel.setArg( 0, cmM.GetData());
      error_number |= kernel.setArg( 1, cmV.GetData());
      error_number |= kernel.setArg( 2, cl_uint( cmM.GetNumberCols()));
      error_number |= kernel.setArg( 3, cl_uint( cmM.GetNumberRows()));
      error_number |= kernel.setArg( 4, cmDevDst.GetData());
      error_number |= kernel.setArg( 5, sizeof( t_DataType) * block_size, 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Read back results and check error
      error_number = m_CommandQueue.enqueueReadBuffer( cmDevDst.GetData(), CL_TRUE, 0, sizeof( t_DataType) * int( STORAGE.GetSize()), ( void *)STORAGE.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void Operations< t_DataType>::VectorPlusEqualsMatrixTimesVector
    (
      linal::VectorInterface< t_DataType> &STORAGE,
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const linal::VectorConstInterface< t_DataType> &FEATURE
    ) const
    {
      BCL_Assert( MATRIX.GetNumberCols() == FEATURE.GetSize() && MATRIX.GetNumberRows() == STORAGE.GetSize(), "non-matching dimensions!");

      cl_int error_number( CL_SUCCESS);              // Error code var

      // create buffer for matrix
      Buffer cmM( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * MATRIX.GetNumberOfElements(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      // creat buffer for vector
      Buffer cmV( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * FEATURE.GetSize(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      // creat buffer for vector to add
      Buffer cmV_A( cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof( t_DataType) * STORAGE.GetSize(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      Buffer cmDevDst( cl::Buffer( m_Context, CL_MEM_WRITE_ONLY, sizeof( t_DataType) * STORAGE.GetSize(), NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Create the kernel
      cl::Kernel kernel( m_Program, "MatrixVectorMultiplicationPlusVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      linal::Vector< t_DataType> tmp( STORAGE.GetSize());

      // Asynchronous write of data to GPU device
      error_number = m_CommandQueue.enqueueWriteBuffer( cmM, CL_FALSE, 0, sizeof( t_DataType) * MATRIX.GetNumberOfElements(), ( void *)MATRIX.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmV, CL_FALSE, 0, sizeof( t_DataType) * FEATURE.GetSize(), ( void *)FEATURE.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      error_number = m_CommandQueue.enqueueWriteBuffer( cmV_A, CL_FALSE, 0, sizeof( t_DataType) * STORAGE.GetSize(), ( void *)STORAGE.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      int block_size( 256);
      cl::NDRange block_dimensions( block_size); // # of work items in the 1D work group
      cl::NDRange kernel_worksize( Tools::RoundUp( block_size, int( MATRIX.GetNumberRows()))); // Total # of work items in the 1D range
      cl::NDRange offset;

      // Launch kernel
      error_number =  kernel.setArg( 0, cmM);
      error_number |= kernel.setArg( 1, cmV);
      error_number |= kernel.setArg( 2, cmV_A);
      error_number |= kernel.setArg( 3, cl_uint( MATRIX.GetNumberCols()));
      error_number |= kernel.setArg( 4, cl_uint( MATRIX.GetNumberRows()));
      error_number |= kernel.setArg( 5, cmDevDst);
      error_number |= kernel.setArg( 6, sizeof( t_DataType) * block_size, 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Read back results and check error
      error_number = m_CommandQueue.enqueueReadBuffer( cmDevDst, CL_TRUE, 0, sizeof( t_DataType) * int( tmp.GetSize()), ( void *)tmp.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      for( typename linal::VectorInterface< t_DataType>::iterator itr( STORAGE.Begin()), itr_end( STORAGE.End()), itr_add( tmp.Begin()); itr != itr_end; ++itr, ++itr_add)
      {
        *itr += *itr_add;
      }

      return;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Operations< t_DataType>::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &Operations< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @return ERROR error that occured, CL_SUCCESS if no error
    template< typename t_DataType>
    cl_int Operations< t_DataType>::CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
    {
      cl_int error_number( CL_SUCCESS);

      const Device device( m_CommandQueue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      // determine the extensions supported by the device
      const storage::Set< Extension> s_device_extensions( device.Extensions( &error_number));
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt
        (
          "cannot get device extensions, error: " + Tools::ErrorString( error_number)
        );
        return error_number;
      }

      // compile the program
      cl::Program::Sources source;
      const std::string linal( ( *GetKernelSources().e_Linal)->GetSource( PRECISION, s_device_extensions));
      const std::string euclidean_source( ( *GetKernelSources().e_EuclideanDistance)->GetSource( PRECISION, s_device_extensions));
      if( linal.empty() || euclidean_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }
      source.push_back( std::make_pair( linal.c_str()      , linal.length()));
      source.push_back( std::make_pair( euclidean_source.c_str(), euclidean_source.length()));

      // create the program
      cl::Program &current_program( m_Program);
      current_program = cl::Program( m_Context, source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
        return error_number;
      }
      // build the program
      Device queue_device( m_CommandQueue.GetDevice( &error_number));
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "create program error getting device: " + Tools::ErrorString( error_number));
        return error_number;
      }
      error_number = current_program.build( std::vector< cl::Device>( 1, queue_device), KernelSourceInterface::GetPrecisionCompilerOptions( PRECISION).c_str());
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

      // for debugging, write the binary
//        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
//        {
//          Tools::LogPtx( current_program, "OpenCLProgram.ptx");
//        }

      // end
      return error_number;
    }

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_OPERATIONS_HPP_
