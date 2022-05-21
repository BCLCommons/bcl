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

#ifndef BCL_OPENCL_MATRIX_MULTIPLY_H_
#define BCL_OPENCL_MATRIX_MULTIPLY_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixMultiply
    //! @brief for matrix-matrix multiplication optimized for the gpu
    //!
    //! @tparam t_DataType can be float, double, complex, int, etc...
    //!
    //! @see @link example_opencl_matrix_multiply.cpp @endlink
    //! @author loweew
    //! @date Mar 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixMultiply :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! command queue
      CommandQueue m_Queue;

      //! program
      cl::Program m_Program;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MatrixMultiply()
      {
      }

      //! @brief constructor
      //! @param QUEUE the command queue
      MatrixMultiply( const CommandQueue &QUEUE) :
        m_Queue( QUEUE),
        m_Program()
      {
        cl_int error_number( CL_SUCCESS);
        m_Program = KernelSources::Compile( GetKernelSources().e_Linal, util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue, std::string(), &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

      //! @brief Clone function
      //! @return pointer to new MatrixMultiply
      MatrixMultiply< t_DataType> *Clone() const
      {
        return new MatrixMultiply< t_DataType>( *this);
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

      //! @brief return the programs
      //! @return reference to the programs
      const cl::Program &GetPrograms() const
      {
        return m_Program;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator for performing matrix-matrix multiplication
      //! @param INPUT_A, INPUT_B, OUTPUT the two matrices to multiply and the results
      void operator()( const Matrix< t_DataType> &INPUT_A, const Matrix< t_DataType> &INPUT_B, Matrix< t_DataType> &OUTPUT) const
      {
        // error catching
        cl_int error_number = CL_SUCCESS;
        const cl_uint block_size( 16);

        BCL_Assert
        (
          INPUT_A.GetNumberRows() % block_size == 0 &&
          INPUT_A.GetNumberCols() % block_size == 0 &&
          INPUT_B.GetNumberRows() % block_size == 0 &&
          INPUT_B.GetNumberCols() % block_size == 0 &&
          OUTPUT.GetNumberRows() == INPUT_A.GetNumberRows() &&
          OUTPUT.GetNumberCols() == INPUT_B.GetNumberCols(),
          "input buffers are not padded correctly!"
        );

        cl::Kernel kernel( m_Program, "MatrixMultiplication", &error_number);

        // set the args values
        error_number  = kernel.setArg( 0, INPUT_A.GetData());
        error_number |= kernel.setArg( 1, INPUT_B.GetData());
        error_number |= kernel.setArg( 2, cl_uint( INPUT_A.GetNumberCols()));
        error_number |= kernel.setArg( 3, cl_uint( INPUT_B.GetNumberCols()));
        error_number |= kernel.setArg( 4, OUTPUT.GetData());
        error_number |= kernel.setArg( 5, block_size * block_size * sizeof( t_DataType), 0); //shared memory
        error_number |= kernel.setArg( 6, block_size * block_size * sizeof( t_DataType), 0); //shared memory
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute Multiplication
        const cl::NDRange block_dims( block_size, block_size);
        const cl::NDRange offset;
        const cl::NDRange worksize( OUTPUT.GetNumberCols(), OUTPUT.GetNumberRows());

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      //! @brief operator for performing matrix-matrix multiplication
      //! @param INPUT_A, INPUT_B the two matrices to multiply
      //! @return the resulting matrix
      Matrix< t_DataType> operator()( const Matrix< t_DataType> &INPUT_A, const Matrix< t_DataType> &INPUT_B) const
      {
        // error catching
        cl_int error_number = CL_SUCCESS;
        const cl_uint block_size( 16);

        BCL_Assert
        (
          INPUT_A.GetNumberRows() % block_size == 0 &&
          INPUT_A.GetNumberCols() % block_size == 0 &&
          INPUT_B.GetNumberRows() % block_size == 0 &&
          INPUT_B.GetNumberCols() % block_size == 0 &&
          INPUT_A.GetNumberCols() == INPUT_B.GetNumberRows(),
          "input buffers are not padded correctly!"
        );

        cl::Kernel kernel( m_Program, "MatrixMultiplication", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        Matrix< t_DataType> output
        (
          INPUT_A.GetNumberRows() - INPUT_A.GetRowPadding(),
          INPUT_B.GetNumberCols() - INPUT_B.GetColPadding(),
          m_Queue,
          INPUT_A.GetRowPadding(),
          INPUT_B.GetColPadding()
        );

        // set the args values
        error_number  = kernel.setArg( 0, INPUT_A.GetData());
        error_number |= kernel.setArg( 1, INPUT_B.GetData());
        error_number |= kernel.setArg( 2, cl_uint( INPUT_A.GetNumberCols()));
        error_number |= kernel.setArg( 3, cl_uint( INPUT_B.GetNumberCols()));
        error_number |= kernel.setArg( 4, output.GetData());
        error_number |= kernel.setArg( 5, block_size * block_size * sizeof( t_DataType), 0); //shared memory
        error_number |= kernel.setArg( 6, block_size * block_size * sizeof( t_DataType), 0); //shared memory
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute Multiplication
        const cl::NDRange block_dims( block_size, block_size);
        const cl::NDRange offset;
        const cl::NDRange worksize( output.GetNumberCols(), output.GetNumberRows());

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        m_Queue.finish();

        return output;
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

    }; // class MatrixMultiply

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_MATRIX_MULTIPLY_H_ 
