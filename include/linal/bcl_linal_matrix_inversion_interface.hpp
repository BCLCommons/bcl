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

#ifndef BCL_LINAL_MATRIX_INVERSION_INTERFACE_HPP_
#define BCL_LINAL_MATRIX_INVERSION_INTERFACE_HPP_

// include the header of this class
#include "bcl_linal_matrix_inversion_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_vector_const_reference.h"
#include "bcl_linal_vector_operations.h"
#include "bcl_linal_vector_reference.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //! @brief Invert diagonal matrix
    //! @param INVERSE_STORAGE storage for the inverse matrix; must have the same size as MATRIX!
    //! @param MATRIX the matrix to try and invert as a diagonal matrix
    //! @return true if an inverse of the diagonal was computed successfully, false if there are zeros on the diagonal
    template< typename t_DataType>
    bool MatrixInversionInterface< t_DataType>::TryInvertDiagonalMatrix
    (
      MatrixInterface< t_DataType> &INVERSE_STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX
    )
    {
      const t_DataType zero( 0), one( 1);
      bool zero_on_diagonal( false);
      INVERSE_STORAGE = zero;
      for( size_t i( 0), n( MATRIX.GetNumberRows()); i < n; ++i)
      {
        if( MATRIX( i, i) != zero)
        {
          INVERSE_STORAGE( i, i) = one / MATRIX( i, i);
        }
        else
        {
          zero_on_diagonal = true;
        }
      }
      return !zero_on_diagonal;
    }

    //! @brief Invert a tridiagonal matrices
    //! @param INVERSE_STORAGE storage for the inverse matrix; must have the same size as MATRIX!
    //! @param MATRIX the matrix to try and invert as a tridiagonal matrix
    //! @return true if an inverse of the tridiagonal was computed successfully; false otherwise
    template< typename t_DataType>
    bool MatrixInversionInterface< t_DataType>::TryInvertTridiagonalMatrix
    (
      MatrixInterface< t_DataType> &INVERSE_STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX
    )
    {
      if( MATRIX.IsEmpty())
      {
        // empty matrix, nothing to do!
        return true;
      }
      // handle the tridiagonal case
      // fastest way to compute is to first extract the diagonals; label them d_lo(lower diagonal),d_on(diagonal),d_up(upper diagonal)
      // see http://en.wikipedia.org/wiki/Tridiagonal_matrix
      const t_DataType zero( 0), one( 1);
      const size_t n( MATRIX.GetNumberRows()), n_off( n - 1);
      Vector< t_DataType> d_lo( n_off), d_on( MATRIX.GetDiagonal()), d_up( n_off);
      for( size_t j( 1), j_prev( 0); j < n; ++j, ++j_prev)
      {
        d_lo( j_prev) = MATRIX( j, j_prev);
        d_up( j_prev) = MATRIX( j_prev, j);
      }

      // now compute theta and phi via recurrence relationships given in
      // Usmani, Inversion of a tridiagonal Jacobi Matrix
      Vector< t_DataType> theta( n + 1, t_DataType( 0)), phi( n + 2, t_DataType( 0));
      theta( 0) = 1;
      theta( 1) = d_on( 0);
      for( size_t j( 2); j <= n; ++j)
      {
        theta( j) = d_on( j - 1) * theta( j - 1) - d_lo( j - 2) * d_up( j - 2) * theta( j - 2);
      }
      if( theta( n) == zero)
      {
        return false;
      }
      phi( n) = 1;
      phi( n - 1) = d_on( n - 1);
      for( size_t j( n - 2); j < n; --j)
      {
        phi( j) = d_on( j) * phi( j + 1) - d_lo( j) * d_up( j) * phi( j + 2);
      }

      // divide all the thetas by the last theta value, which is the determinant
      theta *= one / theta( n);

      // compute the inverse of the tri-diagonal matrix in O(N^2)
      // Pass one computes the minors product, e.g.
      // (-1)^(i+j)*product(d_up[x] for x where i<=x<j) when i < j
      // (-1)^(i+j)*product(d_lo[x] for x where j<x<=i) when i > j
      // TODO: Implement the more numerically stable variant of this algorithm using the formula in
      // Da Fonseca, C. M. (2007). "On the eigenvalues of some tridiagonal matrices".
      // Journal of Computational and Applied Mathematics 200: 283–286
      for( size_t i( 0); i < n; ++i)
      {
        t_DataType *inv_row_i( INVERSE_STORAGE[ i]);
        if( i)
        {
          // handle i > j (lower-triangle)
          const t_DataType *inv_row_i_prev( INVERSE_STORAGE[ i - 1]);
          const t_DataType d_lo_val( d_lo( i - 1));
          for( size_t j( 0); j < i; ++j)
          {
            inv_row_i[ j] = -d_lo_val * inv_row_i_prev[ j];
          }
        }
        // handle diagonal
        inv_row_i[ i] = one;
        if( i + 1 < n)
        {
          // handle i < j (upper triangle
          const t_DataType *itr_d_up( d_up[ i]);
          for( size_t j( i + 1); j < n; ++j, ++itr_d_up)
          {
            inv_row_i[ j] = -inv_row_i[ j - 1] * ( *itr_d_up);
          }
        }
      }

      // 2nd pass multiplies each element in the matrix by its theta and phis.
      // Note that we have already normalized the theta's by theta[n]
      for( size_t i( 0); i < n; ++i)
      {
        t_DataType *inv_row_i( INVERSE_STORAGE[ i]);
        const t_DataType phi_i( phi( i + 1));
        for( size_t j( 0); j < i; ++j)
        {
          inv_row_i[ j] *= theta( j) * phi_i;
        }
        const t_DataType theta_i( theta( i));
        for( size_t j( i); j < n; ++j)
        {
          inv_row_i[ j] *= theta_i * phi( j + 1);
        }
      }
      return true;
    }

    namespace
    {
      //! @brief divide with check for 0/0; in that case, return 0
      //! @param LHS lhs of the / sign
      //! @param RHS rhs of the / sign
      template< typename t_DataType>
      inline t_DataType SafeDiv( const t_DataType &LHS, const t_DataType &RHS)
      {
        return RHS || LHS ? LHS / RHS : t_DataType( 0.0);
      }
    }

    //! @brief Solve a diagonal matrix equation Ax=B
    //! @param A should be the original matrix
    //! @param B the desired output
    //! @note  There is no error checking for this function to keep it O(N), so only use it with known Diagonal matrices!
    template< typename t_DataType>
    Vector< t_DataType> MatrixInversionInterface< t_DataType>::SolveDiagonalMatrix
    (
      const MatrixConstInterface< t_DataType> &A,
      const VectorConstInterface< t_DataType> &B
    )
    {
      Vector< t_DataType> solution( B);
      for( size_t i( 0), n( A.GetNumberRows()); i < n; ++i)
      {
        solution( i) = SafeDiv( B( i), A( i, i));
      }
      return solution;
    }

    //! @brief Solve a tridiagonal matrix; this should be preferred to taking the inverse whenever possible because
    //!        the algorithm is more numerically stable and faster (O(N) vs. O(N^2)).  There is no error checking for
    //!        this function to keep it O(N), so only use it with known TriDiagonal matrices!
    //! @param A should be the original matrix
    //! @param B rhs of the equation
    //! @return Solution vector
    template< typename t_DataType>
    Vector< t_DataType> MatrixInversionInterface< t_DataType>::SolveTridiagonalMatrix
    (
      const MatrixConstInterface< t_DataType> &A,
      const VectorConstInterface< t_DataType> &B
    )
    {
      const size_t n( B.GetSize()), n_last( n - 1);
      Vector< t_DataType> c_prime( n, t_DataType( 0)), solution( B);
      if( B.IsEmpty())
      {
        return solution;
      }

      c_prime( 0) = SafeDiv( A( 0, 1), A( 0, 0));
      solution( 0) = SafeDiv( B( 0), A( 0, 0));
      for( size_t i( 1), i_prev( 0); i < n_last; ++i, ++i_prev)
      {
        const t_DataType *row_i( A[ i]);
        const t_DataType &a( row_i[ i_prev]);
        t_DataType denom( row_i[ i] - a * c_prime( i_prev));
        c_prime( i) = SafeDiv( row_i[ i + 1], denom);
        solution( i) -= a * solution( i_prev);
        solution( i) = SafeDiv( solution( i), denom);
      }
      {
        // handle last element; no need to for c_prime, which remains 0
        const t_DataType *row_i( A[ n_last]);
        const size_t i_prev( n_last - 1);
        const t_DataType &a( row_i[ i_prev]);
        solution( n_last) -= a * solution( i_prev);
        solution( n_last) = SafeDiv( solution( n_last), row_i[ n_last] - a * c_prime( i_prev));
      }

      // back substitution
      for( size_t i( n - 2); i < n; --i)
      {
        solution( i) -= c_prime( i) * solution( i + 1);
      }
      return solution;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_INTERFACE_HPP_
