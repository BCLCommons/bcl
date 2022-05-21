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

#ifndef BCL_OPENCL_MATRIX3X3_H_
#define BCL_OPENCL_MATRIX3X3_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Matrix3x3
    //! @brief is an implementation of a matrix with 3x3 elements
    //! @details
    //!
    //! @see @link example_opencl_matrix3x3.cpp @endlink
    //! @author woetzen
    //! @date Mar 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Matrix3x3 :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      static const size_t s_Dimension      =  3; //!< dimension of this matrix
      static const size_t s_Padding        = 13; //!< padding of that matrix
      static const size_t s_BlockSize      = s_Dimension + s_Padding;

      //! command queue
      CommandQueue m_Queue;

      //! opencl buffer
      Buffer       m_Data;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Matrix3x3();

      //! @brief construct from filler
      //! @param QUEUE command queue
      //! @param FILL_VALUE assign every element to that value
      Matrix3x3( const CommandQueue &QUEUE, const t_DataType &FILL_VALUE);

      //! @brief copy constructor
      //! @param MATRIX the matrix to copy
      Matrix3x3( const Matrix3x3< t_DataType> &MATRIX);

      //! @brief constructor from linal matrix
      //! @param MATRIX matrix to create buffer from
      //! @param QUEUE command queue
      Matrix3x3
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const CommandQueue &QUEUE
      );

      //! @brief Clone function
      //! @return pointer to new Matrix< t_DataType>
      Matrix3x3< t_DataType> *Clone() const;

      //! @brief hard copy
      //! @return matrix with new copied buffer
      Matrix3x3< t_DataType> HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get number of rows
      //! @return number of rows
      size_t GetNumberRows() const;

      //! @brief get number of columns
      //! @return number of columns
      size_t GetNumberCols() const;

      //! @brief gets buffer object
      //! @return buffer
      const Buffer &GetData() const;

      //! @brief gets command queue associated with this buffer
      //! @return queue
      const CommandQueue &GetQueue() const;

      //! @brief returns linal::Matrix3x3
      //! @return the linal::Matrix3x3
      linal::Matrix3x3< t_DataType> GetHostMatrix() const;

      //! @brief returns linal::Matrix with padding
      //! @return the linal::Matrix with padding
      linal::Matrix< t_DataType> GetHostMatrixWithPadding() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief is matrix a square matrix
      //! @return true if number of cols and rows are identical
      bool IsSquare() const;

      //! @brief swap elements between two rows
      //! @param ROW_A the first row
      //! @param ROW_B the second row
      void SwapRows( const size_t ROW_A, const size_t ROW_B);

      //! @brief replace row
      //! @param ROW the row to replace
      //! @param VECTOR the vector to replace the current row with
      void ReplaceRow( const size_t ROW, const Vector< t_DataType> &VECTOR);

      //! @brief sort rows and given vectors (less than)
      //! @param VECTOR vector with 3 elements
      void SortRowsAndVector( Vector< t_DataType> &VECTOR);

      //! @brief orthogonalize row with cross product of the other two
      //! @param ROW the row to replace
      void Orthogonalize( const size_t ROW);

      //! @brief normalize rows by square of elements in vector
      //! @param VECTOR vector with normalization values
      void NormalizeRows( const Vector< t_DataType> &VECTOR);

      //! @brief determinant of this matrix
      //! @return the determinant of the matrix
      Vector< t_DataType> Determinant() const;

      //! @brief transpose the matrix
      //! @return reference to the transposed matrix
      Matrix3x3< t_DataType> &Transpose();

      //! @brief multiply matrix with its transposed
      //! @return symmetrized matrix
      Matrix3x3< t_DataType> &MultiplyWithTransposed();

      //! @brief eigenvalues
      //! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
      //! @return the eigenvalues
      Vector< t_DataType> EigenValues( const bool SQRT = false) const;

    private:

      //! @brief eigenvalues for a symmetric matrix
      //! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
      //! @return the eigenvalues
      Vector< t_DataType> EigenValuesSymmetric( const bool SQRT) const;

    public:

      //! @brief eigenvectors of symmetric matrix
      //! @param EIGEN_VECTORS reference to the matrix that will hold the eigenvectors
      //! @param EIGEN_VALUES reference to the vector that will hold the eigenvalues
      //! @return true is successful
      bool EigenVectorsSymmetric( Matrix3x3< t_DataType> &EIGEN_VECTORS, Vector< t_DataType> &EIGEN_VALUES) const;

    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment
      //! @param MATRIX the matrix to copy
      //! @return reference to this matrix
      Matrix3x3< t_DataType> &operator =( const Matrix3x3< t_DataType> &MATRIX);

      //! @brief *= multiply this with given matrix
      //! @param MATRIX the matrix to multiply onto
      //! @return reference to this matrix
      Matrix3x3< t_DataType> &operator *=( const Matrix3x3< t_DataType> &MATRIX);

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
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

    }; // template class Matrix

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix3x3< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix3x3< double>;

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_MATRIX3X3_H_
