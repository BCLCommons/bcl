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

#ifndef BCL_LINAL_MATRIX_OPERATIONS_H_
#define BCL_LINAL_MATRIX_OPERATIONS_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_vector.h"
#include "bcl_linal_vector_const_reference.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically
#include "bcl_linal_operations.h"
#include "bcl_linal_operations_interface.h"
#include <algorithm>
#include <numeric>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_linal_matrix_operations.h
  //! @brief mathematical operators/operations for linal::Matrix
  //!
  //! @see @link example_linal_matrix_operations.cpp @endlink
  //! @author woetzen, loweew, mendenjl
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace linal
  {

  //////////////////////
  // Matrix functions //
  //////////////////////

    //! @brief check dimension agreement of two Matrices
    //! @param MATRIX_LHS rhs matrix
    //! @param MATRIX_RHS lhs matrix
    //! @return true if number rows and cols are the same between both Matrices
    template< typename t_DataType>
    inline
    bool
    SameDimensions
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      return    MATRIX_LHS.GetNumberRows() == MATRIX_RHS.GetNumberRows()
             && MATRIX_LHS.GetNumberCols() == MATRIX_RHS.GetNumberCols();
    }

    //! @brief check inverse dimension agreement of two Matrices
    //! comapre number ros of lhs with number cols of rhs and number cols of lhs with number rows of rhs
    //! @param MATRIX_LHS rhs matrix
    //! @param MATRIX_RHS lhs matrix
    //! @return true if number rows with cols and cols with rows agree between both Matrices
    template< typename t_DataType>
    inline
    bool
    InverseDimensions
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      return    MATRIX_LHS.GetNumberRows() == MATRIX_RHS.GetNumberCols()
             && MATRIX_LHS.GetNumberCols() == MATRIX_RHS.GetNumberRows();
    }

    //! @brief check dimensions for multiplication A*B
    //! compare number cols of lhs with number rows of rhs
    //! @param MATRIX_LHS rhs matrix
    //! @param MATRIX_RHS lhs matrix
    //! @return true if number cols rhs and number cols lhs agree
    template< typename t_DataType>
    inline
    bool
    MultiplicationDimension
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      return MATRIX_LHS.GetNumberCols() == MATRIX_RHS.GetNumberRows();
    }

  //////////////////////
  // binary operators //
  //////////////////////

    //! @brief add one matrix to another
    //! @param MATRIX_LHS matrix to add to
    //! @param MATRIX_RHS matrix to add
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator +=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      BCL_Assert
      (
        SameDimensions( MATRIX_LHS, MATRIX_RHS),
        "provided Matrices have different dimensions; rhs != lhs : " + util::Format()( MATRIX_LHS) + " != " +
        util::Format()( MATRIX_RHS)
      );

      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_RHS.Begin(),
        MATRIX_LHS.Begin(),
        std::plus< t_DataType>()
      );

      return MATRIX_LHS;
    }

    //! @brief subtract one matrix from another
    //! @param MATRIX_LHS matrix to subtract from
    //! @param MATRIX_RHS matrix to subtract
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator -=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      BCL_Assert
      (
        SameDimensions< t_DataType>( MATRIX_LHS, MATRIX_RHS),
        "provided Matrices have different dimensions; rhs != lhs : " + util::Format()( MATRIX_LHS) + " != " +
        util::Format()( MATRIX_RHS)
      );

      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_RHS.Begin(),
        MATRIX_LHS.Begin(),
        std::minus< t_DataType>()
      );

      return MATRIX_LHS;
    }

    //! @brief add scalar to matrix
    //! @param MATRIX_LHS matrix to add to
    //! @param VALUE scalar to be added
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator +=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE
    )
    {
      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_LHS.Begin(),
        std::bind2nd( std::plus< t_DataType>(), VALUE)
      );

      return MATRIX_LHS;
    }

    //! @brief subtract scalar from matrix
    //! @param MATRIX_LHS matrix to subtract from
    //! @param VALUE scalar to be added
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator -=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE
    )
    {
      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_LHS.Begin(),
        std::bind2nd( std::minus< t_DataType>(), VALUE)
      );

      return MATRIX_LHS;
    }

    //! @brief multiply matrix with scalar
    //! @param MATRIX_LHS matrix to multiply to
    //! @param SCALAR scalar to be multiplied
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator *=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &SCALAR
    )
    {
      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_LHS.Begin(),
        std::bind2nd( std::multiplies< t_DataType>(), SCALAR)
      );

      return MATRIX_LHS;
    }

    //! @brief divide matrix by scalar
    //! @param MATRIX_LHS matrix to divide
    //! @param SCALAR scalar to divide by
    //! @return the changed lhs matrix
    template< typename t_DataType>
    inline
    MatrixInterface< t_DataType> &
    operator /=
    (
      MatrixInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &SCALAR
    )
    {
      std::transform
      (
        MATRIX_LHS.Begin(), MATRIX_LHS.End(),
        MATRIX_LHS.Begin(),
        std::bind2nd( std::divides< t_DataType>(), SCALAR)
      );

      return MATRIX_LHS;
    }

  //////////////////////////////
  // binary logical operators //
  //////////////////////////////

    //! @brief compare if all items in matrix are equal to a given VALUE
    //! @param MATRIX_LHS matrix with values
    //! @param VALUE_RHS value that is compared against
    //! @return true if matrix is empty are all elements in matrix are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator ==
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      return std::find_if
             (
               MATRIX_LHS.Begin(), MATRIX_LHS.End(),
               std::bind2nd( std::not_equal_to< t_DataType>(), VALUE_RHS)
             ) == MATRIX_LHS.End();
    }

    //! @brief compare if all items in matrix are equal to a given VALUE
    //! @param VALUE_LHS value that is compared against
    //! @param MATRIX_RHS matrix with values
    //! @return true if matrix is empty are all elements in matrix are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator ==
    (
      const t_DataType &VALUE_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      return ( MATRIX_RHS == VALUE_LHS);
    }

    //! @brief compare if all items in matrix are not equal to a given VALUE
    //! @param MATRIX_LHS matrix with values
    //! @param VALUE_RHS value that is compared against
    //! @return false if matrix is empty are all elements in matrix are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator !=
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      return std::find_if
             (
               MATRIX_LHS.Begin(), MATRIX_LHS.End(),
               std::bind2nd( std::not_equal_to< t_DataType>(), VALUE_RHS)
             ) != MATRIX_LHS.End();
    }

    //! @brief compare if all items in matrix are not equal to a given VALUE
    //! @param VALUE_LHS value that is compared against
    //! @param MATRIX_RHS matrix with values
    //! @return false if matrix is empty are all elements in matrix are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator !=
    (
      const t_DataType &VALUE_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      return ( MATRIX_RHS != VALUE_LHS);
    }

  //////////////////////
  // binary operators //
  //////////////////////

    //! @brief sum two matrixs of equal size
    //! @param MATRIX_LHS lhs matrix
    //! @param MATRIX_RHS rhs matrix
    //! @return matrix with all individual summed elements of lhs and rhs matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator +
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      BCL_Assert
      (
        SameDimensions( MATRIX_LHS, MATRIX_RHS),
        "matrices of different sizes supplied: " +
        util::Format()( MATRIX_LHS) + " != " + util::Format()( MATRIX_RHS)
      );

      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix += MATRIX_RHS;
    }

    //! @brief subtract two matrixs of equal size
    //! @param MATRIX_LHS lhs matrix
    //! @param MATRIX_RHS rhs matrix
    //! @return matrix with all individual subtracted elements of rhs from lhs matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator -
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      BCL_Assert
      (
        SameDimensions< t_DataType>( MATRIX_LHS, MATRIX_RHS),
        "matrices of different sizes supplied: " +
        util::Format()( MATRIX_LHS) + " != " + util::Format()( MATRIX_RHS)
      );

      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix -= MATRIX_RHS;
    }

    namespace
    {
      // Anonymous namespace for optimized matrix multiplication functions for various types of matrices
      // all functions in this namespace assume that the matrices have multipliable dimensions

      //! @brief outer product of two vectors
      //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
      //! @param VECTOR_U left hand side vector with m elements
      //! @param VECTOR_V right hand side vector with n elements
      //! @return Matrix with m*n elements (outer product of vector u and v)
      template< typename t_DataType>
      Matrix< t_DataType> OuterProduct
      (
        const MatrixConstInterface< t_DataType> &VECTOR_U,
        const MatrixConstInterface< t_DataType> &VECTOR_V
      )
      {
        const size_t v_size( VECTOR_V.GetNumberCols());
        Matrix< t_DataType> outer_product( VECTOR_U.GetNumberRows(), VECTOR_V.GetNumberCols());

        // ptr to first element of first row in matrix
        t_DataType *ptr_matrix( outer_product.Begin());

        // each result row is the product of the current VECTOR_U element multiplied with VECTOR_V
        for( const t_DataType *ptr_u( VECTOR_U.Begin()), *ptr_u_end( VECTOR_U.End()); ptr_u != ptr_u_end; ++ptr_u)
        {
          std::transform
          (
            VECTOR_V.Begin(), VECTOR_V.End(),
            ptr_matrix,
            std::bind2nd( std::multiplies< t_DataType>(), *ptr_u)
          );
          // ptr_matrix points to first element in next row
          ptr_matrix += v_size;
        }

        // end
        return outer_product;
      }

      //! @brief determine whether a matrix multiplication would likely be faster to compute using transposition than
      //! @param ROWS number of rows in the resulting matrix (== number rows in the LHS matrix)
      //! @param DIM number of columns of the LHS matrix == number of rows of the RHS matrix
      //! @param COLS number of columns in the resulting matrix (== number cols in RHS matrix)
      bool MatrixMultiplyPreferTransposition( const size_t &ROWS, const size_t &DIM, const size_t &COLS)
      {
        // this heuristic was determined on a benchmark of around 100 data points with various numbers for each of these
        // parameters ranging from 2 to 36.
        // A linear curve was fit to the ratio of transposed matrix multiply speed vs. various metrics, then simplified
        // to minimize the cost of the heuristic itself.  In practice, the only cases the heuristic is wrong on are
        // where there is minimal difference between using the transposed vs. the non-transposed version.
        return size_t( 3 * std::min( ROWS, DIM) + 2 * std::min( DIM, COLS) + ROWS) >= size_t( 25);
      }

    }

    //! @brief dense matrix multiplications; no optimizations used; appropriate for very small matrices
    //! @param STORAGE the matrix where the result will be stored
    //! @param MATRIX_A lhs matrix
    //! @param MATRIX_B rhs matrix
    //! @return matrix computed using rows of A times columns of B, no transposition
    template< typename t_DataType>
    void MultiplySmallDenseMatrices
    (
      MatrixInterface< t_DataType> &STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX_A,
      const MatrixConstInterface< t_DataType> &MATRIX_B
    )
    {
      t_DataType *mat_ptr( STORAGE.Begin());
      const size_t n_rows( STORAGE.GetNumberRows()), n_cols( STORAGE.GetNumberCols()), dim( MATRIX_A.GetNumberCols());
      // MATRIX B is small and likely will all be cached at once; so doing the computation via B transpose is likely
      // insignificant or even more costly in terms of performance
      for( size_t i( 0); i < n_rows; ++i)
      {
        const t_DataType *const row_a( MATRIX_A[ i]);
        for( size_t j( 0); j < n_cols; ++j, ++mat_ptr)
        {
          t_DataType sum( 0);
          for( size_t k( 0); k < dim; ++k)
          {
            sum += row_a[ k] * MATRIX_B( k, j);
          }
          *mat_ptr = sum;
        }
      }
    }

    //! @brief compute a matrix times its transpose, without the need to compute the transpose, e.g. M * M^T
    //! @param A the matrix for which to compute A * A^T
    //! @return M * M^T
    template< typename t_DataType>
    Matrix< t_DataType> MatrixTimesItselfTransposed( const MatrixConstInterface< t_DataType> &A)
    {
      const size_t n_rows( A.GetNumberRows());
      const size_t dim( A.GetNumberCols());
      Matrix< t_DataType> mat( n_rows, n_rows);
      // ptr to the current element of the results matrix (mat) being edited
      t_DataType *mat_ptr( mat.Begin());
      for( size_t i( 0); i < n_rows; ++i)
      {
        const t_DataType *const a_row_i( A[ i]), *const a_row_i_end( a_row_i + dim);
        size_t j( 0);
        // copy the lower diagonal from the already-computed upper diagonal
        for( ; j < i; ++j, ++mat_ptr)
        {
          *mat_ptr = mat( j, i);
        }
        // compute upper diagonal as A[row i] * A[row j]
        for( ; j < n_rows; ++j, ++mat_ptr)
        {
          // compute
          *mat_ptr = std::inner_product( a_row_i, a_row_i_end, A[ j], *mat_ptr);
        }
      }
      return mat;
    }

    //! @brief A^T * A, without the need to compute A^T
    //! @param A the matrix for which to compute A^T * A
    //! @return A^T * A
    template< typename t_DataType>
    Matrix< t_DataType> MatrixTransposeTimesMatrix( const MatrixConstInterface< t_DataType> &A)
    {
      const size_t n_rows( A.GetNumberRows()), n_cols( A.GetNumberCols());
      Matrix< t_DataType> mat( n_cols, n_cols, t_DataType( 0));
      // ptr to the current element of the results matrix (mat) being edited
      for( size_t i( 0); i < n_rows; ++i)
      {
        const t_DataType *const a_row_i( A[ i]);
        // compute upper diagonal as A[row i] * A[row j]
        for( size_t j( 0); j < n_cols; ++j)
        {
          t_DataType *itr_mat_row_j( mat[ j]);
          const t_DataType *itr_a_row_i( a_row_i);
          const t_DataType a_j( a_row_i[ j]);
          for( size_t k( 0); k <= j; ++k, ++itr_mat_row_j, ++itr_a_row_i)
          {
            *itr_mat_row_j += a_j * *itr_a_row_i;
          }
        }
      }
      for( size_t i( 0); i < n_cols; ++i)
      {
        // copy lower triangular to upper triangular
        for( size_t j( 0); j < i; ++j)
        {
          mat( j, i) = mat( i, j);
        }
      }
      return mat;
    }

    //! @brief A^T * B, without the need to compute A^T
    //! @param A the matrix for which to compute A^T * B
    //! @param B the rhs of the * in A^T * B
    //! @return A^T * B
    template< typename t_DataType>
    Matrix< t_DataType> MatrixTransposeTimesMatrix
    (
      const MatrixConstInterface< t_DataType> &A,
      const MatrixConstInterface< t_DataType> &B
    )
    {
      // A is an N X M matrix
      // B is an N X K matrix
      const size_t n_rows( A.GetNumberRows()), n_cols_a( A.GetNumberCols()), n_cols_b( B.GetNumberCols());
      Matrix< t_DataType> mat( n_cols_a, n_cols_b, t_DataType( 0));
      // ptr to the current element of the results matrix (mat) being edited
      for( size_t i( 0); i < n_rows; ++i)
      {
        const t_DataType *const a_row_i( A[ i]);
        const t_DataType *const b_row_i( B[ i]), *const b_row_i_end( b_row_i + n_cols_b);
        t_DataType *itr_mat( mat.Begin());
        // compute upper diagonal as A[row i] * A[row j]
        for( size_t j( 0); j < n_cols_a; ++j)
        {
          const t_DataType a_j( a_row_i[ j]);
          for( const t_DataType *itr_b_row_i( b_row_i); itr_b_row_i != b_row_i_end; ++itr_mat, ++itr_b_row_i)
          {
            *itr_mat += a_j * *itr_b_row_i;
          }
        }
      }
      return mat;
    }

    //! @brief A * B^T, without the need to compute B^T
    //! @param A the matrix for which to compute A * B^T
    //! @param B the rhs of the * in A * B^T
    //! @return A * B^T
    template< typename t_DataType>
    Matrix< t_DataType> MatrixTimesMatrixTranspose
    (
      const MatrixConstInterface< t_DataType> &MATRIX_A,
      const MatrixConstInterface< t_DataType> &MATRIX_B
    )
    {
      const size_t n_rows( MATRIX_A.GetNumberRows()), n_cols( MATRIX_B.GetNumberRows());
      const size_t dim( MATRIX_A.GetNumberCols());
      BCL_Assert
      (
        MATRIX_A.GetNumberCols() == MATRIX_B.GetNumberCols(),
        "matrices have non-matching dimension! [" + util::Format()( n_rows) + " X " + util::Format()( dim) + " ] * [ "
        + util::Format()( MATRIX_B.GetNumberCols()) + " X " + util::Format()( n_cols) + " ]"
      );
      Matrix< t_DataType> mat( n_rows, n_cols);

      // ptr to the current element of the results matrix (mat) being edited
      t_DataType *mat_ptr( mat.Begin());
      for( size_t i( 0); i < n_rows; ++i)
      {
        for( size_t j( 0); j < n_cols; ++j, ++mat_ptr)
        {
          // compute MATRIX_A[row i] * MATRIX_B[ col j]
          t_DataType sum( 0);
          // get a pointer to row i of MATRIX_A, row j of MATRIX_B'
          const t_DataType *mat_b_ptr( MATRIX_B[ j]);
          for
          (
            const t_DataType *mat_a_ptr( MATRIX_A[ i]), *mat_a_ptr_end( mat_a_ptr + dim);
            mat_a_ptr != mat_a_ptr_end;
            ++mat_a_ptr, ++mat_b_ptr
          )
          {
            sum += *mat_a_ptr * ( *mat_b_ptr);
          }
          // set the new matrix value up correctly
          *mat_ptr = sum;
        }
      }
      return mat;
    }

    //! @brief multiply two matrixs of equal size by building the inner product yielding the scalar product
    //! @param MATRIX_A lhs matrix
    //! @param MATRIX_B rhs matrix
    //! @return scalar representing root of inner product of the two ranges
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator *
    (
      const MatrixConstInterface< t_DataType> &MATRIX_A,
      const MatrixConstInterface< t_DataType> &MATRIX_B
    )
    {
      const size_t n_rows( MATRIX_A.GetNumberRows()), n_cols( MATRIX_B.GetNumberCols());
      const size_t dim( MATRIX_A.GetNumberCols());
      BCL_Assert
      (
        MultiplicationDimension< t_DataType>( MATRIX_A, MATRIX_B),
        "matrices have non-matching dimension! [" + util::Format()( n_rows) + " X " + util::Format()( dim) + " ] * [ "
        + util::Format()( MATRIX_B.GetNumberRows()) + " X " + util::Format()( n_cols) + " ]"
      );

      // result will be empty matrix
      if( dim == 0 || n_cols == 0)
      {
        return Matrix< t_DataType>( 0, 0);
      }
      else if( dim == 1)
      {
        // outer product
        return OuterProduct( MATRIX_A, MATRIX_B);
      }
      else if( n_rows == 1 && n_cols == 1)
      {
        // scalar product
        return
          Matrix< t_DataType>
          (
            1,
            1,
            std::inner_product( MATRIX_A.Begin(), MATRIX_A.End(), MATRIX_B.Begin(), t_DataType( 0))
          );
      }
      else if( n_rows <= 2 && n_cols <= 2 && dim <= 2)
      {
        // tiny (usually 2x2) matrix multiply.  Even checking IsDiagonal in these cases can be slower than just running
        // them through dense matrix multiplication
        Matrix< t_DataType> result( n_rows, n_cols);
        MultiplySmallDenseMatrices( result, MATRIX_A, MATRIX_B);
        return result;
      }

      // handle diagonal cases:
      // 1. Diagonal X Diagonal
      // 2. Diagonal X Dense
      // 3. Dense X Diagonal
      if( MATRIX_A.IsDiagonal())
      {
        if( MATRIX_B.IsDiagonal())
        {
          // multiplying two diagonal matrices, trivial case
          Matrix< t_DataType> square_matrix( MATRIX_A);
          for( size_t i( 0); i < dim; ++i)
          {
            square_matrix( i, i) *= MATRIX_B( i, i);
          }
          return square_matrix;
        }
        // A-diag * B -> Yields matrix B where every value B(i,j) is multiplied by A(i,i)
        Matrix< t_DataType> mat( MATRIX_B);
        t_DataType *itr_mat( mat.Begin());
        for( size_t i( 0); i < dim; ++i)
        {
          const t_DataType a_value( MATRIX_A( i, i));
          for( t_DataType *itr_mat_end( itr_mat + n_cols); itr_mat != itr_mat_end; ++itr_mat)
          {
            *itr_mat *= a_value;
          }
        }
        return mat;
      }
      else if( MATRIX_B.IsDiagonal())
      {
        // A * B-diag -> Yields matrix A where every value A(i,j) is multiplied by B(j,j)
        // N X M * M X M
        Matrix< t_DataType> mat( MATRIX_A);
        Vector< t_DataType> b_diag( MATRIX_B.GetDiagonal());
        const t_DataType *itr_b_diag_end( b_diag.End());
        t_DataType *itr_mat( mat.Begin());
        for( size_t i( 0); i < n_rows; ++i)
        {
          for( const t_DataType *itr_b( b_diag.Begin()); itr_b != itr_b_diag_end; ++itr_mat, ++itr_b)
          {
            *itr_mat *= *itr_b;
          }
        }
        return mat;
      }

      // handle the cases where MATRIX_B and/or MATRIX_A are too small for tranposition to improve performance
      if( !MatrixMultiplyPreferTransposition( n_rows, dim, n_cols))
      {
        Matrix< t_DataType> mat( n_rows, n_cols);
        MultiplySmallDenseMatrices( mat, MATRIX_A, MATRIX_B);
        return mat;
      }

      // Memory non-locality (cache misses) tends to be the dominant speed bottleneck in matrix multiplication
      // for non-trivially size matrices.  Computing A x B using Transpose(B) is faster, because memory access is
      // contiguous.
      // Computing Tranpose(B) is not always desired though since it may consume extra memory, and is not always
      // necessary, so don't do it for the following special cases:
      // 1. Matrix B is symmetric (-> so MatrixB == its own transpose)
      // 2. Matrix A is Matrix B transposed (e.g. when computing SVD / tridiagonalization / eigensystems)
      // Otherwise, construct Transpose(B)

      // create an object to hold the transpose of matrix B
      Matrix< t_DataType> matrix_b_transposed_storage;
      util::SiPtr< const MatrixConstInterface< t_DataType> > si_ptr_matrix_b_transposed( matrix_b_transposed_storage);
      if( MATRIX_B.IsSymmetric())
      {
        // special case where there is no need to tranpose B
        si_ptr_matrix_b_transposed = util::SiPtr< const MatrixConstInterface< t_DataType> >( MATRIX_B);
      }
      else if( dim == MATRIX_B.GetNumberRows() && n_rows == n_cols)
      {
        // check whether B == A transpose
        bool b_is_a_tranpose( true);
        for( size_t i( 0); i < n_rows && b_is_a_tranpose; ++i)
        {
          for( size_t j( 0); j < dim; ++j)
          {
            if( MATRIX_A( i, j) != MATRIX_B( j, i))
            {
              b_is_a_tranpose = false;
              break;
            }
          }
        }
        if( b_is_a_tranpose)
        {
          // return the matrix times its transpose using this specialized algorithm
          BCL_MessageVrb
          (
            "Matrix * Matrix Transpose detected; prefer use of linal::MatrixTimesItsTranspose"
            "or linal::MatrixTransposeTimesMatrix to save computation of the transpose!"
          );
          return MatrixTimesItselfTransposed( MATRIX_A);
        }
        else
        {
          // asymmetric B matrix, transpose 1st
          matrix_b_transposed_storage = MATRIX_B.Transposed();
        }
      }
      else
      {
        // asymmetric B matrix, transpose 1st
        matrix_b_transposed_storage = MATRIX_B.Transposed();
      }

      Matrix< t_DataType> mat( n_rows, n_cols);

      // create a reference on the matrix
      const MatrixConstInterface< t_DataType> &matrix_b_transposed( *si_ptr_matrix_b_transposed);

      // ptr to the current element of the results matrix (mat) being edited
      t_DataType *mat_ptr( mat.Begin());
      for( size_t i( 0); i < n_rows; ++i)
      {
        for( size_t j( 0); j < n_cols; ++j, ++mat_ptr)
        {
          // compute MATRIX_A[row i] * MATRIX_B[ col j]
          t_DataType sum( 0);
          // get a pointer to row i of MATRIX_A, row j of MATRIX_B'
          const t_DataType *mat_b_ptr( matrix_b_transposed[ j]);
          for
          (
            const t_DataType *mat_a_ptr( MATRIX_A[ i]), *mat_a_ptr_end( mat_a_ptr + dim);
            mat_a_ptr != mat_a_ptr_end;
            ++mat_a_ptr, ++mat_b_ptr
          )
          {
            sum += *mat_a_ptr * ( *mat_b_ptr);
          }
          // set the new matrix value up correctly
          *mat_ptr = sum;
        }
      }

      return mat;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void VectorPlusEqualsMatrixTimesVector
    (
      VectorInterface< t_DataType> &STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX,
      const VectorConstInterface< t_DataType> &FEATURE
    )
    {
      const size_t num_rows( MATRIX.GetNumberRows());
      const size_t num_cols( MATRIX.GetNumberCols());
      BCL_Assert
      (
        num_cols == FEATURE.GetSize() && num_rows == STORAGE.GetSize(),
        "non-matching dimensions! " + util::Format()( num_rows) + "X" + util::Format()( num_cols)
        + " * " + util::Format()( FEATURE.GetSize()) + " * " + util::Format()( STORAGE.GetSize())
      );

      Vector< t_DataType> local_feature;
      const bool feature_overwrites_storage( size_t( &FEATURE) == size_t( &STORAGE));
      if( feature_overwrites_storage)
      {
        local_feature = FEATURE; // STORAGE is the same as FEATURE, so a copy of the feature vector is necessary
      }

      const t_DataType *vector_data( feature_overwrites_storage ? local_feature.Begin() : FEATURE.Begin());

      for( size_t i( 0); i < num_rows; ++i)
      {
        t_DataType &sum_of_products( STORAGE( i));
        const t_DataType *row( MATRIX[ i]);
        for( size_t j( 0); j < num_cols; ++j)
        {
          sum_of_products += row[ j] * vector_data[ j];
        }
      }
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void VectorEqualsVectorTimesMatrix
    (
      VectorInterface< t_DataType> &STORAGE,
      const VectorConstInterface< t_DataType> &FEATURE,
      const MatrixConstInterface< t_DataType> &MATRIX
    )
    {
      const size_t num_rows( MATRIX.GetNumberRows());
      const size_t num_cols( MATRIX.GetNumberCols());
      BCL_Assert
      (
        num_cols == STORAGE.GetSize() && num_rows == FEATURE.GetSize(),
        "non-matching dimensions! " + util::Format()( num_rows) + "X" + util::Format()( num_cols)
        + " * " + util::Format()( FEATURE.GetSize()) + " * " + util::Format()( STORAGE.GetSize())
      );

      Vector< t_DataType> local_feature;
      const bool feature_overwrites_storage( FEATURE.Begin() == STORAGE.Begin());
      if( feature_overwrites_storage)
      {
        local_feature = FEATURE; // STORAGE is the same vector as FEATURE, so a copy of the feature vector is necessary
      }

      // initialize the storage
      STORAGE = 0;

      // get a pointer to the beginning of the feature array, accounting for whether a copy of the feature
      // had to be made
      const t_DataType *feature_ref( feature_overwrites_storage ? local_feature.Begin() : FEATURE.Begin());

      // STORAGE(i) = ( sum of column i of MATRIX) * FEATURE( i)
      // start by setting STORAGE(i) = ( sum of column i of MATRIX)
      for( size_t i( 0); i < num_rows; ++i)
      {
        const t_DataType *itr_row( MATRIX[ i]);
        const t_DataType vector_value( feature_ref[ i]);
        for
        (
          t_DataType *itr_result( STORAGE.Begin()), *itr_result_end( STORAGE.End());
          itr_result != itr_result_end;
          ++itr_result, ++itr_row
        )
        {
          *itr_result += vector_value * ( *itr_row);
        }
      }
    }

    //! operator Matrix * Vector
    template< typename t_DataType>
    inline Vector< t_DataType> operator *( const MatrixConstInterface< t_DataType> &MATRIX, const VectorConstInterface< t_DataType> &VECTOR)
    {
      // While valid, the following line of code will produce a huge waste of memory, as it copies what may be a very
      // large matrix, and creates a vector of the resulting size at least twice.
      // return Matrix< t_DataType> ( MATRIX).operator *=( VECTOR);
      Vector< t_DataType> storage( MATRIX.GetNumberRows(), 0.0);
      VectorPlusEqualsMatrixTimesVector( storage, MATRIX, VECTOR);
      return storage;
    }

    //! operator Vector * Matrix
    template< typename t_DataType>
    inline Vector< t_DataType> operator *( const VectorConstInterface< t_DataType> &VECTOR, const MatrixConstInterface< t_DataType> &MATRIX)
    {
      Vector< t_DataType> result( MATRIX.GetNumberCols());
      VectorEqualsVectorTimesMatrix( result, VECTOR, MATRIX);
      return result;
    }

    //! @brief add value to matrix
    //! @param MATRIX_LHS lhs matrix
    //! @param VALUE_RHS rhs value to be added
    //! @return matrix that has the value added to each value of the lhs given matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator +
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix += VALUE_RHS;
    }

    //! @brief add matrix to value
    //! @param VALUE_LHS lhs value to be added
    //! @param MATRIX_RHS rhs matrix
    //! @return matrix that has the value added to each value of the lhs given matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator +
    (
      const t_DataType &VALUE_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_RHS);
      return new_matrix += VALUE_LHS;
    }

    //! @brief subtract value from matrix
    //! @param MATRIX_LHS lhs matrix
    //! @param VALUE_RHS rhs value to be subtracted
    //! @return matrix that has the value subtracted from each value of the lhs given matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator -
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix -= VALUE_RHS;
    }

    //! @brief subtract matrix from value
    //! @param VALUE_LHS rhs value to be subtracted
    //! @param MATRIX_RHS lhs matrix
    //! @return matrix that has the values in the matrix subtracted from the value
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator -
    (
      const t_DataType &VALUE_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_RHS.GetSize(), VALUE_LHS);
      return new_matrix -= MATRIX_RHS;
    }

    //! @brief multiply scalar with matrix
    //! @param SCALAR_LHS lhs value to be multiplied
    //! @param MATRIX_RHS rhs matrix
    //! @return matrix that has the values multiplied with the scalar
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator *
    (
      const t_DataType &SCALAR_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_RHS);
      return new_matrix *= SCALAR_LHS;
    }

    //! @brief multiply matrix with scalar
    //! @param MATRIX_LHS lhs matrix
    //! @param SCALAR_RHS rhs value to be multiplied
    //! @return matrix that has the values multiplied with the scalar
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator *
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &SCALAR_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix *= SCALAR_RHS;
    }

    //! @brief divide matrix with scalar
    //! @param MATRIX_LHS lhs matrix
    //! @param SCALAR_RHS rhs value to be divided by
    //! @return matrix that has the values divided by the scalar
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator /
    (
      const MatrixConstInterface< t_DataType> &MATRIX_LHS,
      const t_DataType &SCALAR_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_LHS);
      return new_matrix /= SCALAR_RHS;
    }

    //! @brief divide scalar by matrix
    //! @param SCALAR_LHS lhs value to be divided
    //! @param MATRIX_RHS rhs matrix to be used to divide the scalar
    //! @return matrix that has the values of scalar divided by each according value of the matrix
    template< typename t_DataType>
    inline
    Matrix< t_DataType>
    operator /
    (
      const t_DataType &SCALAR_LHS,
      const MatrixConstInterface< t_DataType> &MATRIX_RHS
    )
    {
      Matrix< t_DataType> new_matrix( MATRIX_RHS.GetNumberRows(), MATRIX_RHS.GetNumberCols(), SCALAR_LHS);
      return new_matrix /= MATRIX_RHS;
    }

  //////////////////////
  // Matrix functions //
  //////////////////////

    //! length of one matrix
    template< typename t_DataType>
    inline
    t_DataType
    Norm( const MatrixConstInterface< t_DataType> &MATRIX)
    {
      return MATRIX.Norm();
    }

    //! copies elements of argument VECTOR into diagonal of this object
    template< typename t_DataType>
    void ReplaceDiagonal
    (
      Matrix< t_DataType> &MATRIX,
      const Vector< t_DataType> &VECTOR
    )
    {
      //check that matrix is a square matrix
      BCL_Assert( MATRIX.IsSquare(), "applied on non-square matrix!");
      //check that argument vector has the same size as the width and height of matrix
      BCL_Assert( VECTOR.GetSize() == MATRIX.GetNumberRows(), "non-matching dimensions!");

      //pointer on first element of argument vector
      const t_DataType *dat( VECTOR.Begin()), *dat_end( VECTOR.End());
      const size_t increment( MATRIX.GetNumberCols() + 1);

      //iterate over diagonal elements
      for
      (
        t_DataType *ptr( MATRIX.Begin()), *ptr_end( MATRIX.End());
        ptr < ptr_end && dat != dat_end;
        ptr += increment, ++dat
      )
      {
        ( *ptr) = ( *dat);
      }
    }

    //! @brief find rows with undefined values in a given matrix
    //! @return a vector containing all row ids that have undefined values
    template< typename t_DataType>
    storage::Vector< size_t> FindUndefinedRows( const MatrixConstInterface< t_DataType> &MATRIX)
    {
      storage::Vector< size_t> rows_undefined;
      const size_t n_cols( MATRIX.GetNumberCols());
      for( size_t row( 0), n_rows( MATRIX.GetNumberRows()); row < n_rows; ++row)
      {
        const VectorConstReference< t_DataType> row_reference( n_cols, MATRIX[ row]);
        if( !row_reference.IsDefined())
        {
          rows_undefined.PushBack( row);
        }
      }
      return rows_undefined;
    }

    //! @brief find rows with undefined values in a given matrix
    //! @return a vector containing all row ids that have undefined values
    template< typename t_DataType>
    storage::Vector< size_t> FindCompletelyUndefinedRows( const MatrixConstInterface< t_DataType> &MATRIX)
    {
      storage::Vector< size_t> rows_undefined;
      const size_t n_cols( MATRIX.GetNumberCols());
      if( !n_cols)
      {
        return rows_undefined;
      }
      for( size_t row( 0), n_rows( MATRIX.GetNumberRows()); row < n_rows; ++row)
      {
        const VectorConstReference< t_DataType> row_reference( n_cols, MATRIX[ row]);
        bool has_defined( false);
        for( size_t col( 0); col < n_cols; ++col)
        {
          if( util::IsDefined( row_reference( col)))
          {
            has_defined = true;
            break;
          }
        }
        if( !has_defined)
        {
          rows_undefined.PushBack( row);
        }
      }
      return rows_undefined;
    }

    //! @brief outer product of two vectors given a matrix to store the result into and an assignment operation
    //! A = u x v, A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
    //! @param A the matrix to use as storage.  Must be of dimension m*n
    //! @param OPERATION e.g. Equals< t_DataType> gives A = u x v. PlusEquals< t_DataType> gives A += u x v
    //! @param VECTOR_U left hand side vector with m elements
    //! @param VECTOR_V right hand side vector with n elements
    //! @return Matrix with m*n elements (A OPERATION (outer product of vector u and v))
    template< typename t_DataType>
    void AddOuterProductToMatrix
    (
      MatrixInterface< t_DataType> &A,
      const VectorConstInterface< t_DataType> &VECTOR_U,
      const VectorConstInterface< t_DataType> &VECTOR_V
    )
    {
      const size_t num_rows( VECTOR_U.GetSize());
      const size_t num_cols( VECTOR_V.GetSize());

      BCL_Assert
      (
        num_rows == A.GetNumberRows() && num_cols == A.GetNumberCols(),
        "non-matching dimensions passed to StoreOuterProduct"
      );

      // ptr to first element of each vertex
      const t_DataType *v_row( VECTOR_V.Begin());

      // each storage_row is the product of the current VECTOR_U element multiplied with VECTOR_V
      for( size_t i( 0); i < num_rows; ++i)
      {
        // get a pointer to the first position of this row in STORAGE
        t_DataType *storage_row( A[ i]);

        // the value for VECTOR_U in this row will not change, so store it now
        const t_DataType u_value( VECTOR_U( i));

        // make the storage row by multiplying  VECTOR_U(i) by each element in VECTOR_V
        for( size_t j( 0); j < num_cols; ++j)
        {
          storage_row[ j] += u_value * v_row[ j];
        }
      }
    }

    //! @brief subtracts outer product of two vectors given a matrix to store the result into and an assignment operation
    //! A = u x v, A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
    //! @param A the matrix to use as storage.  Must be of dimension m*n
    //! @param VECTOR_U left hand side vector with m elements
    //! @param VECTOR_V right hand side vector with n elements
    //! @return Matrix with m*n elements (A OPERATION (outer product of vector u and v))
    template< typename t_DataType>
    void SubtractOuterProductFromMatrix
    (
      MatrixInterface< t_DataType> &A,
      const VectorConstInterface< t_DataType> &VECTOR_U,
      const VectorConstInterface< t_DataType> &VECTOR_V
    )
    {
      const size_t num_rows( VECTOR_U.GetSize());
      const size_t num_cols( VECTOR_V.GetSize());

      BCL_Assert
      (
        num_rows == A.GetNumberRows() && num_cols == A.GetNumberCols(),
        "non-matching dimensions passed to StoreOuterProduct"
      );

      // ptr to first element of each vertex
      const t_DataType *v_row( VECTOR_V.Begin());

      // each storage_row is the product of the current VECTOR_U element multiplied with VECTOR_V
      for( size_t i( 0); i < num_rows; ++i)
      {
        // get a pointer to the first position of this row in STORAGE
        t_DataType *storage_row( A[ i]);

        // the value for VECTOR_U in this row will not change, so store it now
        const t_DataType u_value( VECTOR_U( i));

        // make the storage row by multiplying  VECTOR_U(i) by each element in VECTOR_V
        for( size_t j( 0); j < num_cols; ++j)
        {
          storage_row[ j] -= u_value * v_row[ j];
        }
      }
    }

    //! @brief Test matrices for equality, ignoring nans (which would ordinarily cause the comparison to fail)
    //! @param LHS, RHS the matrices to test for equality
    //! @return true if the matrices are equal
    template< typename t_DataType>
    bool EqualIgnoringNans
    (
      const MatrixConstInterface< t_DataType> &A,
      const MatrixConstInterface< t_DataType> &B
    )
    {
      if( A.GetNumberRows() != B.GetNumberRows() || A.GetNumberCols() != B.GetNumberCols())
      {
        return false;
      }
      for
      (
        const t_DataType *itr_a( A.Begin()), *a_end( A.End()), *itr_b( B.Begin());
        itr_a != a_end;
        ++itr_a, ++itr_b
      )
      {
        if( *itr_a != *itr_b)
        {
          if( !util::IsDefined( *itr_a) && !util::IsDefined( *itr_b))
          {
            continue;
          }
          return false;
        }
      }
      return true;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_OPERATIONS_H_
