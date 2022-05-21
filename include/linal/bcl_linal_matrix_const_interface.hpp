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

#ifndef BCL_LINAL_MATRIX_CONST_INTERFACE_HPP_
#define BCL_LINAL_MATRIX_CONST_INTERFACE_HPP_

// include the header of this class
#include "bcl_linal_matrix_const_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_inversion_cholesky.h"
#include "bcl_linal_matrix_inversion_gauss_jordan.h"
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_vector.h"
#include "bcl_linal_vector_const_reference.h"
#include "math/bcl_math.h"
#include "math/bcl_math_statistics.h"
#include "type/bcl_type_enable_if.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////
  // Data Access //
  /////////////////

    //! @brief is the matrix empty
    //! @return true if number rows or number cols is zero
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsEmpty() const
    {
      return !GetNumberOfElements();
    }

    //! @brief is matrix a square matrix
    //! @return true if number of cols and rows are identical
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsSquare() const
    {
      return GetNumberRows() && GetNumberRows() == GetNumberCols();
    }

    //! check whether matrix is symmetric
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsSymmetric() const
    {
      if( !IsSquare())
      {
        return false;
      }

      // get the square matrix dimension
      const size_t dim( GetNumberRows());
      for( size_t i( 0); i < dim; ++i)
      {
        const t_DataType *const row( operator[]( i));
        for( size_t j( i + 1); j < dim; ++j)
        {
          if( row[ j] != operator()( j, i))
          {
            return false;
          }
        }
      }

      return true;
    }

    //! @brief is matrix a diagonal matrix
    //! @return true if all but the elements in the diagonal are 0
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsDiagonal() const
    {
      // if matrix is not square or empty
      if( !IsSquare() || GetNumberOfElements() == 0)
      {
        return false;
      }

      // create a zero object of the type
      const t_DataType zero( 0);
      const size_t dim( GetNumberRows());

      // check that all but the elements in the diagonal are 0
      for( size_t i( 0); i < dim; ++i)
      {
        const t_DataType *const row( operator[]( i));
        for( size_t j( 0); j < dim; ++j)
        {
          if( !math::EqualWithinMachineTolerance( row[ j], zero) && j != i)
          {
            return false;
          }
        }
      }

      // return true if all elements but diagonal are 0
      return true;
    }

    //! @brief is matrix a tridiagonal matrix
    //! @return if diagonal and adjacent diagonals are filled and the rest is 0
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsTriDiagonal() const
    {
      // if matrix is not square and does not have at least 2 elements in each dimension
      const size_t dim( GetNumberRows());
      if( !IsSquare() || dim < 2)
      {
        return false;
      }
      else if( dim == 2)
      {
        // 2d matrices are always tri-diagonal
        return true;
      }

      // create a zero object of the type
      const t_DataType zero( 0);

      // check that all but the inner three diagonal elements are 0
      for( size_t i( 0); i < dim; ++i)
      {
        const t_DataType *const row_start( operator[]( i));
        const t_DataType *const first_diagonal( row_start + i - 1);
        // check elements before the first of the three diagonals
        for( const t_DataType *col( row_start); col < first_diagonal; ++col)
        {
          if( !math::EqualWithinMachineTolerance( *col, zero))
          {
            return false;
          }
        }
        // check elements after the last of the three diagonals
        const t_DataType *const row_end( row_start + dim);
        const t_DataType *const last_diagonal( row_start + i + 1);
        for( const t_DataType *col( last_diagonal + 1); col < row_end; ++col)
        {
          if( !math::EqualWithinMachineTolerance( *col, zero))
          {
            return false;
          }
        }
      }

      // return true if all elements but the inner three diagonals are 0
      return true;
    }

    //! @brief checks if matrix is defined
    //! @return bool - true if matrix is defined
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::IsDefined() const
    {
      return math::Statistics::IsDefined( Begin(), End());
    }

    //! @brief get a vector reference on a particular row
    //! @param ROW the row of interest
    template< typename t_DataType>
    VectorConstReference< t_DataType> MatrixConstInterface< t_DataType>::GetRow( const size_t &ROW) const
    {
      return
        VectorConstReference< t_DataType>
        (
          this->GetNumberCols(),
          this->operator[]( ROW)
        );
    }

    //! @brief copy a column into a vector
    //! @param COL the column of interest
    template< typename t_DataType>
    Vector< t_DataType> MatrixConstInterface< t_DataType>::GetCol( const size_t &COL) const
    {
      const size_t number_cols( GetNumberCols());
      Vector< t_DataType> col( GetNumberRows());
      const_iterator itr_this( Begin() + COL);
      for
      (
        typename Vector< t_DataType>::iterator itr_col( col.Begin()), itr_col_end( col.End());
        itr_col != itr_col_end;
        ++itr_col, itr_this += number_cols
      )
      {
        *itr_col = *itr_this;
      }
      return col;
    }

    //! @brief get a vector containing the diagonal of this matrix
    template< typename t_DataType>
    Vector< t_DataType> MatrixConstInterface< t_DataType>::GetDiagonal() const
    {
      BCL_Assert( GetNumberRows() == GetNumberCols(), "Diagonal of a non-square matrix!");

      //copy elements
      const size_t number_rows( GetNumberRows());
      Vector< t_DataType> diagonal( number_rows);

      for( size_t i( 0); i < number_rows; ++i)
      {
        diagonal( i) = operator()( i, i);
      }

      //return
      return diagonal;
    }

    //! @brief get the sum of the matrix
    template< typename t_DataType>
    t_DataType MatrixConstInterface< t_DataType>::Sum() const
    {
      return math::Statistics::Sum( Begin(), End());
    }

    //! compute trace of matrix
    template< typename t_DataType>
    t_DataType MatrixConstInterface< t_DataType>::Trace() const
    {
      BCL_Assert( GetNumberRows() == GetNumberCols(), "Trace of a non-square matrix!");
      t_DataType sum( 0);
      const size_t number_rows( GetNumberRows());
      for( size_t i( 0); i < number_rows; ++i)
      {
        sum += operator()( i, i);
      }
      return sum;
    }

    //! @brief Get the matrix as a vector
    template< typename t_DataType>
    VectorConstReference< t_DataType> MatrixConstInterface< t_DataType>::AsVector() const
    {
      return VectorConstReference< t_DataType>( GetNumberOfElements(), Begin());
    }

    //! @brief Return a transposed copy of the matrix
    //! @return A transposed copy of the matrix
    template< typename t_DataType>
    Matrix< t_DataType> MatrixConstInterface< t_DataType>::Transposed() const
    {
      const size_t rows( GetNumberRows()), cols( GetNumberCols());
      Matrix< t_DataType> transposed( cols, rows);
      for( size_t i( 0); i < cols; ++i)
      {
        t_DataType *transposed_row_i( transposed[ i]);
        for( size_t j( 0); j < rows; ++j)
        {
          transposed_row_i[ j] = operator()( j, i);
        }
      }

      //end
      return transposed;
    }

    namespace
    {
      // overloads for Solve and Inverse

      //! @brief Get the inverse of the matrix; uses the fastest method that appears to work on the matrix
      //! @param A the matrix to compute the inverse of
      //! @return floating point matrix with the inverse
      template< typename t_DataType>
      typename type::EnableIf
      <
        !std::numeric_limits< t_DataType>::is_integer,
        Matrix< typename MatrixConstInterface< t_DataType>::FloatType>
      >::Type
      InverseImpl( const MatrixConstInterface< t_DataType> &A)
      {
        const size_t rows( A.GetNumberRows()), cols( A.GetNumberCols());
        if( A.IsEmpty())
        {
          return A;
        }
        else if( rows != cols)
        {
          // handle rectangular matrices
          Matrix< t_DataType> newmatrix;
          const bool transposed( rows < cols);
          newmatrix = transposed ? MatrixTimesItselfTransposed( A) : MatrixTransposeTimesMatrix( A);
          return transposed ? ( newmatrix.Inverse() * A).Transposed() : newmatrix.Inverse() * A.Transposed();
        }
        // handle diagonal/tridiagonal cases; for which tri-diagonal will always return true
        if( A.IsTriDiagonal())
        {
          Matrix< t_DataType> inverse( A);
          if( A.IsDiagonal())
          {
            BCL_Assert
            (
              MatrixInversionInterface< t_DataType>::TryInvertDiagonalMatrix( inverse, A),
              "Could not invert diagonal matrix: " + util::Format()( A)
            );
            return inverse;
          }
          BCL_Assert
          (
            MatrixInversionInterface< t_DataType>::TryInvertTridiagonalMatrix( inverse, A),
            "Could not invert tridiagonal matrix: " + util::Format()( A)
          );
          return inverse;
        }

        // ultimately, everything below this line should default to running LU-decomposition, with QR as backup
        // if the LU-decomposition fails.
        MatrixInversionGaussJordan< t_DataType> inv( A, true);
        BCL_Assert( inv.IsDefined(), "Non-invertible matrix!");
        return inv.ComputeInverse();
      }

      //! @brief Get the inverse of the matrix; uses the fastest method that appears to work on the matrix
      //! @param A the matrix to compute the inverse of
      //! @return floating point matrix with the inverse
      template< typename t_DataType>
      typename type::EnableIf
      <
        std::numeric_limits< t_DataType>::is_integer,
        Matrix< typename MatrixConstInterface< t_DataType>::FloatType>
      >::Type
      InverseImpl( const MatrixConstInterface< t_DataType> &A)
      {
        // integer types; convert to floating point first
        return InverseImpl( Matrix< typename MatrixConstInterface< t_DataType>::FloatType>( A));
      }

      //! @brief solve implementation for integer types; first converts the matrix to double, then calls solve
      //! @param A the matrix to solve in Ax=B
      //! @param B the vector to solve for in Ax=B
      //! @return vector X
      template< typename t_DataType>
      typename type::EnableIf
      <
        !std::numeric_limits< t_DataType>::is_integer,
        Vector< typename MatrixConstInterface< t_DataType>::FloatType>
      >::Type
      SolveImpl
      (
        const MatrixConstInterface< t_DataType> &A,
        const VectorConstInterface< typename MatrixConstInterface< t_DataType>::FloatType> &B
      )
      {
        // overload for floating point types
        const size_t rows( A.GetNumberRows()), cols( A.GetNumberCols()), b_dim( B.GetSize());
        if( A.IsEmpty())
        {
          return Vector< t_DataType>();
        }
        else if( rows == cols)
        {
          BCL_Assert( rows == b_dim, "Incompatible dimensions for matrix to solve!");
          // direct solve
          // handle diagonal/tridiagonal cases; for which tri-diagonal will always return true
          if( A.IsTriDiagonal())
          {
            if( A.IsDiagonal())
            {
              return MatrixInversionInterface< t_DataType>::SolveDiagonalMatrix( A, B);
            }
            return MatrixInversionInterface< t_DataType>::SolveTridiagonalMatrix( A, B);
          }
          if( A.IsSymmetric())
          {
            // Try cholesky decomposition first, since it is nearly an order of magnitude
            // faster than gauss jordan as it is currently implemented; also more stable
            MatrixInversionCholesky< t_DataType> inv_chol( A);
            if( inv_chol.IsDefined())
            {
              return inv_chol.Solve( B);
            }
          }
          // if all else fails, use gauss jordan
          return MatrixInversionGaussJordan< t_DataType>( A, true).Solve( B);
        }
        // least squares case

        // Need to test whether A or A^T was provided; for convenience we allow either to be provided, since it is
        // clear which was desired by the dimensions of B.
        // M X N matrix times N vector = N vector
        // Compute Transpose(A) * A (N X N) and solve with cholesky
        return SolveImpl( b_dim == rows ? A * A.Transposed() : A.Transposed() * A, B);
      }

      //! @brief solve implementation for integer types; first converts the matrix to double, then calls solve
      //! @param A the matrix to solve in Ax=B
      //! @param B the vector to solve for in Ax=B
      //! @return vector X
      template< typename t_DataType>
      typename type::EnableIf
      <
        std::numeric_limits< t_DataType>::is_integer,
        Vector< typename MatrixConstInterface< t_DataType>::FloatType>
      >::Type
      SolveImpl
      (
        const MatrixConstInterface< t_DataType> &A,
        const VectorConstInterface< typename MatrixConstInterface< t_DataType>::FloatType> &B
      )
      {
        // overload for integer types
        return SolveImpl( Matrix< typename MatrixConstInterface< t_DataType>::FloatType>( A), B);
      }

    }

    //! @brief Get the inverse of the matrix; uses the fastest method that appears to work on the matrix
    template< typename t_DataType>
    Matrix< typename MatrixConstInterface< t_DataType>::FloatType> MatrixConstInterface< t_DataType>::Inverse() const
    {
      return InverseImpl( *this);
    }

    //! @brief Solve the normal matrix equation Ax=B using the fastest available method
    //! @param B the rhs of the matrix equation
    //! @return x, the solution vector
    //! If this matrix is rectangular and over-determined, computes the least squares fit of A onto X
    template< typename t_DataType>
    Vector< typename MatrixConstInterface< t_DataType>::FloatType> MatrixConstInterface< t_DataType>::Solve
    (
      const VectorConstInterface< FloatType> &B
    ) const
    {
      return SolveImpl( *this, B);
    }

    //! @brief compare two matrices for equality
    //! @param MATRIX_RHS rhs matrix to compare to
    //! @return true if dimensions and content match
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::operator ==( const MatrixConstInterface< t_DataType> &MATRIX_RHS) const
    {
      if( GetNumberRows() != MATRIX_RHS.GetNumberRows() || GetNumberCols() != MATRIX_RHS.GetNumberCols())
      {
        return false;
      }

      return std::equal( Begin(), End(), MATRIX_RHS.Begin());
    }

    //! @brief compare two matrices for inequality
    //! @param MATRIX_RHS rhs matrix to compare to
    //! @return true if dimensions and content match
    template< typename t_DataType>
    bool MatrixConstInterface< t_DataType>::operator !=( const MatrixConstInterface< t_DataType> &MATRIX_RHS) const
    {
      return !operator ==( MATRIX_RHS);
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_CONST_INTERFACE_HPP_
