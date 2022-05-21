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

#ifndef BCL_OPENCL_MATRIX3X3_HPP_
#define BCL_OPENCL_MATRIX3X3_HPP_

// include the header of this class
#include "bcl_opencl_matrix3x3.h"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_operations.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix3x3.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3()
    {
    }

    //! @brief construct from filler
    //! @param QUEUE command queue
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const CommandQueue &QUEUE, const t_DataType &FILL_VALUE) :
      m_Queue( QUEUE),
      m_Data( Buffer::CreateBufferFromMatrix< t_DataType>( linal::Matrix< t_DataType>( s_BlockSize, s_BlockSize, FILL_VALUE), m_Queue))
    {
    }

    //! @brief copy constructor
    //! @param MATRIX the matrix to copy
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const Matrix3x3< t_DataType> &MATRIX) :
      m_Queue( MATRIX.m_Queue),
      m_Data( MATRIX.m_Data)
    {
    }

    //! @brief constructor from linal matrix
    //! @param MATRIX matrix to create buffer from
    //! @param QUEUE command queue
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const CommandQueue &QUEUE
    ) :
      m_Queue( QUEUE),
      m_Data
      (
        Buffer::CreateBufferFromMatrix< t_DataType>
        (
          linal::Matrix< t_DataType>( MATRIX, s_Padding, s_Padding), m_Queue
        )
      )
    {
      BCL_Assert
      (
        MATRIX.GetNumberCols() == s_Dimension && MATRIX.GetNumberRows() == s_Dimension,
        "construct matrix from matrix of incompatible size: " + util::Format()( MATRIX)
      );
    }

    //! @brief Clone function
    //! @return pointer to new Matrix< t_DataType>
    template< typename t_DataType>
    Matrix3x3< t_DataType> *Matrix3x3< t_DataType>::Clone() const
    {
      return new Matrix3x3< t_DataType>( *this);
    }

    //! @brief hard copy
    //! @return matrix with new copied buffer
    template< typename t_DataType>
    Matrix3x3< t_DataType> Matrix3x3< t_DataType>::HardCopy() const
    {
      Matrix3x3< t_DataType> copy;
      copy.m_Queue = m_Queue;
      copy.m_Data = Buffer( cl::Buffer( copy.m_Queue.GetContext(), CL_MEM_READ_WRITE, s_BlockSize * s_BlockSize * sizeof( t_DataType), NULL));
      m_Queue.enqueueCopyBuffer( m_Data, copy.m_Data, 0, 0, s_BlockSize * s_BlockSize * sizeof( t_DataType));

      return copy;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Matrix3x3< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get number of rows
    //! @return number of rows
    template< typename t_DataType>
    size_t Matrix3x3< t_DataType>::GetNumberRows() const
    {
      return s_Dimension;
    }

    //! @brief get number of columns
    //! @return number of columns
    template< typename t_DataType>
    size_t Matrix3x3< t_DataType>::GetNumberCols() const
    {
      return s_Dimension;
    }

    //! @brief gets buffer object
    //! @return buffer
    template< typename t_DataType>
    const Buffer &Matrix3x3< t_DataType>::GetData() const
    {
      return m_Data;
    }

    //! @brief gets command queue associated with this buffer
    //! @return queue
    template< typename t_DataType>
    const CommandQueue &Matrix3x3< t_DataType>::GetQueue() const
    {
      return m_Queue;
    }

    //! @brief returns linal::Matrix3x3
    //! @return the linal::Matrix3x3
    template< typename t_DataType>
    linal::Matrix3x3< t_DataType> Matrix3x3< t_DataType>::GetHostMatrix() const
    {
      linal::Matrix< t_DataType> tmp( s_BlockSize, s_BlockSize);
      cl_int error_number = m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * s_BlockSize * s_BlockSize, tmp.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));
      return linal::Matrix3x3< t_DataType>( tmp.CreateSubMatrix( s_Dimension, s_Dimension));
    }

    //! @brief returns linal::Matrix with padding
    //! @return the linal::Matrix with padding
    template< typename t_DataType>
    linal::Matrix< t_DataType> Matrix3x3< t_DataType>::GetHostMatrixWithPadding() const
    {
      linal::Matrix< t_DataType> tmp( s_BlockSize, s_BlockSize);
      cl_int error_number = m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * s_BlockSize * s_BlockSize, tmp.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));
      return tmp;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is matrix a square matrix
    //! @return true if number of cols and rows are identical
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsSquare() const
    {
      return true;
    }

    //! @brief swap elements between two rows
    //! @param ROW_A the first row
    //! @param ROW_B the second row
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::SwapRows( const size_t ROW_A, const size_t ROW_B)
    {
      BCL_Exit( "not implemented yet", -1);
    }

    //! @brief replace row
    //! @param ROW the row to replace
    //! @param VECTOR the vector to replace the current row with
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::ReplaceRow( const size_t ROW, const Vector< t_DataType> &VECTOR)
    {
      BCL_Exit( "not implemented yet", -1);
    }

    //! @brief sort rows and given vectors (less than)
    //! @param VECTOR vector with 3 elements
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::SortRowsAndVector( Vector< t_DataType> &VECTOR)
    {
      // Create the kernel
      cl_int error_number( CL_SUCCESS);              // Error code var
      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "Matrix3x3SortRowsAndVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, VECTOR.GetData());
      error_number |= kernel.setArg( 2, s_BlockSize * s_BlockSize * sizeof( t_DataType), 0); //shared memory
      error_number |= kernel.setArg( 3, s_BlockSize * sizeof( t_DataType), 0); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
    }

    //! @brief orthogonalize row with cross product of the other two
    //! @param ROW the row to replace
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::Orthogonalize( const size_t ROW)
    {
      // Create the kernel
      cl_int error_number( CL_SUCCESS);              // Error code var
      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "Matrix3x3OrthogonalizeRow", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, cl_uint( ROW));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
    }

    //! @brief normalize rows by square of elements in vector
    //! @param VECTOR vector with normalization values
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::NormalizeRows( const Vector< t_DataType> &VECTOR)
    {
      // Create the kernel
      cl_int error_number( CL_SUCCESS);              // Error code var
      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "Matrix3x3NormalizeRows", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, VECTOR.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
    }

    //! @brief determinant of this matrix
    //! @return the determinant of the matrix
    template< typename t_DataType>
    Vector< t_DataType> Matrix3x3< t_DataType>::Determinant() const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      Vector< t_DataType> determinat( 1, m_Queue);

      // Create the kernel
      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "Matrix3x3Determinant", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, determinat.GetData());
      error_number |= kernel.setArg( 2, s_BlockSize * s_BlockSize * sizeof( t_DataType), 0); //shared memory
      error_number |= kernel.setArg( 3, s_BlockSize * sizeof( t_DataType), 0); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return determinat;
    }

    //! @brief transpose the matrix
    //! @return reference to the transposed matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::Transpose()
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "MatrixTranspose", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, cl_uint( s_BlockSize));
      error_number |= kernel.setArg( 2, cl_uint( s_BlockSize));
      error_number |= kernel.setArg( 3, m_Data);
      error_number |= kernel.setArg( 4, s_BlockSize * ( s_BlockSize + 1) * sizeof( t_DataType), NULL); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return *this;
    }

    //! @brief multiply matrix with its transposed
    //! @return symmetrized matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::MultiplyWithTransposed()
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "Matrix3x3MultiplyWithTransposed", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, s_BlockSize * s_BlockSize * sizeof( t_DataType), NULL); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return *this;
    }

    //! @brief eigenvalues
    //! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
    //! @return the eigenvalues
    template< typename t_DataType>
    Vector< t_DataType> Matrix3x3< t_DataType>::EigenValues( const bool SQRT) const
    {
      return EigenValuesSymmetric( SQRT);
    }

    //! @brief eigenvalues for a symmetric matrix
    //! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
    //! @return the eigenvalues
    template< typename t_DataType>
    Vector< t_DataType> Matrix3x3< t_DataType>::EigenValuesSymmetric( const bool SQRT) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      Vector< t_DataType> eigenvalues( s_Dimension, m_Queue, s_Padding);

      // Create the kernel
      cl::Kernel kernel( Operations< t_DataType>().GetInstance().GetProgram(), "Matrix3x3SymmetricEigenvalues", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, cl_uint( SQRT));
      error_number |= kernel.setArg( 2, eigenvalues.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      error_number = m_Queue.enqueueTask( kernel);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return eigenvalues;
    }

    //! @brief eigenvectors of symmetric matrix
    //! @param EIGEN_VECTORS reference to the matrix that will hold the eigenvectors
    //! @param EIGEN_VALUES reference to the vector that will hold the eigenvalues
    //! @return true is successful
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::EigenVectorsSymmetric( Matrix3x3< t_DataType> &EIGEN_VECTORS, Vector< t_DataType> &EIGEN_VALUES) const
    {
      cl_int error_number( CL_SUCCESS);              // Error code var
      {
        // Create the kernel
        cl::Kernel kernel( Operations< t_DataType>().GetInstance().GetProgram(), "Matrix3x3SymmetricEigenvalues", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // set the args values
        error_number  = kernel.setArg( 0, m_Data);
        error_number |= kernel.setArg( 1, cl_uint( 0));
        error_number |= kernel.setArg( 2, EIGEN_VALUES.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute Multiplication
        error_number = m_Queue.enqueueTask( kernel);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // Create the kernel
      cl::Kernel kernel( Operations< t_DataType>().GetInstance().GetProgram(), "Matrix3x3SymmetricEigenvectors", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, EIGEN_VECTORS.GetData());
      error_number |= kernel.setArg( 2, EIGEN_VALUES.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute
      error_number = m_Queue.enqueueTask( kernel);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment
    //! @param MATRIX the matrix to copy
    //! @return reference to this matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator =( const Matrix3x3< t_DataType> &MATRIX)
    {
      m_Queue = MATRIX.m_Queue;
      m_Data  = MATRIX.m_Data;

      // end
      return *this;
    }

    //! @brief *= multiply this with given matrix
    //! @param MATRIX the matrix to multiply onto
    //! @return reference to this matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator *=( const Matrix3x3< t_DataType> &MATRIX)
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      cl::Kernel kernel( Operations< t_DataType>::GetInstance().GetProgram(), "MatrixMultiplication", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // set the args values
      error_number  = kernel.setArg( 0, m_Data);
      error_number |= kernel.setArg( 1, MATRIX.m_Data);
      error_number |= kernel.setArg( 2, cl_uint( s_BlockSize));
      error_number |= kernel.setArg( 3, cl_uint( s_BlockSize));
      error_number |= kernel.setArg( 4, m_Data);
      error_number |= kernel.setArg( 5, s_BlockSize * s_BlockSize * sizeof( t_DataType), 0); //shared memory
      error_number |= kernel.setArg( 6, s_BlockSize * s_BlockSize * sizeof( t_DataType), 0); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute Multiplication
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Matrix3x3< t_DataType>::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Matrix3x3< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_MATRIX3X3_HPP_
