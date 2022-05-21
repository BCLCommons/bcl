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

#ifndef BCL_LINAL_MATRIX_INVERSION_CHOLESKY_HPP_
#define BCL_LINAL_MATRIX_INVERSION_CHOLESKY_HPP_

// include the header of this class
#include "bcl_linal_matrix_inversion_cholesky.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_vector_const_reference.h"
#include "bcl_linal_vector_operations.h"
#include "bcl_linal_vector_reference.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
  //////////////////
  // Construction //
  //////////////////

    //! @brief Default constructor
    //! @param SMOOTHING amount of smoothing to apply; higher values result in smaller L2 norm of the solutions vector,
    //!        and more robust solutions to numerical issues and ill-conditioning, but at the cost of slightly lowered
    //!        accuracy
    template< typename t_DataType>
    MatrixInversionCholesky< t_DataType>::MatrixInversionCholesky( const t_DataType &SMOOTHING) :
      m_DecompositionIsDefined( true),
      m_Smoothing( SMOOTHING)
    {
    }

    //! @brief constructor from a matrix
    //! @param MATRIX matrix for which to compute an inverse
    template< typename t_DataType>
    MatrixInversionCholesky< t_DataType>::MatrixInversionCholesky
    (
      const MatrixConstInterface< t_DataType> &MATRIX
    ) :
      m_Decomposition(),
      m_DecompositionIsDefined( true),
      m_Smoothing( 0.000)
    {
      SetMatrix( MATRIX);
    }

    //! @brief Clone function
    //! @return pointer to new MatrixInterface
    template< typename t_DataType>
    MatrixInversionCholesky< t_DataType> *MatrixInversionCholesky< t_DataType>::Clone() const
    {
      return new MatrixInversionCholesky< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &MatrixInversionCholesky< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    template< typename t_DataType>
    const std::string &MatrixInversionCholesky< t_DataType>::GetAlias() const
    {
      static const std::string s_alias( "Cholesky");
      return s_alias;
    }

    //! @brief determine whether the solution is well-defined for all input vectors
    //! @return true if the inverse is well-defined
    template< typename t_DataType>
    bool MatrixInversionCholesky< t_DataType>::IsDefined() const
    {
      return m_DecompositionIsDefined;
    }

    //! @brief Set the source matrix.  Derived classes should store the source matrix, or a decomposition thereof
    //! @param SOURCE the matrix for which inversions and solutions will be required
    //! @return bool indicating whether matrix decomposition and storage was successful
    //! Implementations may choose to compute the matrix in this function, store the matrix directly, or compute a
    //! decomposition of the matrix
    template< typename t_DataType>
    bool MatrixInversionCholesky< t_DataType>::SetMatrix( const MatrixConstInterface< t_DataType> &A)
    {
      // reset
      m_Decomposition = A;
      // multiply the diagonal by 1%; this makes the solutions just a little off, but at the benefit of vastly improved
      // numerical stability.
      ReplaceDiagonal( m_Decomposition, A.GetDiagonal() * t_DataType( 1.0 + m_Smoothing));
      m_DecompositionIsDefined = true;

      // a small number; if diagonal is between -s_small and s_small, it will be set to s_small
      const t_DataType s_small( A.GetNumberOfElements() * std::numeric_limits< t_DataType>::min());
      if( !A.IsSymmetric())
      {
        // Non-symmetric matrix detected; inverse is not defined by this method
        m_DecompositionIsDefined = false;
        return false;
      }
      const size_t n( A.GetNumberRows());
      // for each row; perform the decomposition in-place using the Cholesky-Banachiewicz algorithm
      for( size_t j( 0); j < n; ++j)
      {
        // compute the diagonal element, C(j,j) = Sqrt(A(j,j)-Sum(C(j,k)^2 | j < k))
        t_DataType *a_row_j( m_Decomposition[ j]);
        t_DataType diagonal( a_row_j[ j]);
        for( size_t k( 0); k < j; ++k)
        {
          diagonal -= math::Sqr( a_row_j[ k]);
        }
        // test whether the diagonal is <= 0 -> no defined inverse
        if( diagonal < s_small)
        {
          if( diagonal >= -s_small)
          {
            // set values that are very near zero just slightly above zero
            diagonal = s_small;
          }
          else
          {
            BCL_MessageStd
            (
              "Matrix is not a positive definite matrix starting at row " + util::Format()( j)
              + "; cholesky diagonal entry was: Sqrt(" + util::Format()( diagonal) + "). "
              "Recommend increasing the smoothing parameter if this is being used for linear least squares"
            );
            m_DecompositionIsDefined = false;
            break;
          }
        }
        else
        {
          // take the sqrt -> not nan because diagonal > 0
          a_row_j[ j] = diagonal = math::Sqrt( diagonal);

          // fill in column J of the remaining rows
          const t_DataType inverse_diagonal( t_DataType( 1.0) / diagonal);
          for( size_t i( j + 1); i < n; ++i)
          {
            t_DataType *a_row_i( m_Decomposition[ i]);
            t_DataType a_row_i_col_j( a_row_i[ j]);
            for( size_t k( 0); k < j; ++k)
            {
              a_row_i_col_j -= a_row_i[ k] * a_row_j[ k];
            }
            a_row_i[ j] = a_row_i_col_j * inverse_diagonal;
          }
        }
      }

      // set the upper triangle of m_Decomposition = lower triangle
      for( size_t j( 0); j < n; ++j)
      {
        t_DataType *a_row_j( m_Decomposition[ j]);
        for( size_t k( j + 1); k < n; ++k)
        {
          a_row_j[ k] = m_Decomposition( k, j);
        }
      }
      return m_DecompositionIsDefined;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get/Compute the inverse matrix, which shall be empty if matrix inversion was unsuccessful
    //! Implementations should not perform this operation repeatedly on the same matrix; just cache whether the
    //! matrix in SetMatrix was inverted
    //! @note For some implementations (Cholesky, LU decomposition, etc.), matrix inverse should not be used to solve
    //!       Ax=B for just 1-2 B vectors because it is usually slower, less stable, and requires more memory than
    //!       just solving the equations directly; for those cases, use Solve below
    //! @return Inverse matrix
    template< typename t_DataType>
    Matrix< t_DataType> MatrixInversionCholesky< t_DataType>::ComputeInverse()
    {
      BCL_Assert( m_DecompositionIsDefined, "Cannot compute inverse of non-invertible matrix");
      // This is not the most efficient way to invert a cholesky decomposition, but Cholesky is primarily intended
      // for use with solving single equations anyway, so the performance hit of roughly 3x slower than the optimal
      // method (at http://arxiv.org/pdf/1111.4144.pdf) should be acceptable
      const size_t n( m_Decomposition.GetNumberRows());
      Matrix< t_DataType> inverse( IdentityMatrix< t_DataType>( n));
      for( size_t i( 0); i < n; ++i)
      {
        inverse.ReplaceRow( i, Solve( inverse.GetRow( i)));
      }
      return inverse;
    }

    //! @brief Solve Ax=B, where A is the matrix already given, x is the vector to determine, and B is the result vec
    //! @param B the result vector
    //! @return Solution vector and bool indicating whether the solution is unique.
    template< typename t_DataType>
    Vector< t_DataType> MatrixInversionCholesky< t_DataType>::Solve
    (
      const VectorConstInterface< t_DataType> &B
    ) const
    {
      const size_t n( m_Decomposition.GetNumberRows());
      BCL_Assert( m_DecompositionIsDefined, "Cannot solve on a non-invertible matrix");

      Vector< t_DataType> solution( B);

      t_DataType *b( solution.Begin());
      // Solve Lx=B substitution
      for( size_t i( 0); i < n; ++i)
      {
        t_DataType sum( b[ i]);
        const t_DataType *decomp_row_i( m_Decomposition[ i]);
        for( size_t j( 0); j < i; ++j)
        {
          sum -= decomp_row_i[ j] * b[ j];
        }
        if( m_Decomposition( i, i) < 1.0e-20)
        {
          b[ i] = 0.0;
        }
        else
        {
          b[ i] = sum / m_Decomposition( i, i);
        }
      }

      // Solve L*y=x
      for( size_t i( n - 1); i < n; --i)
      {
        t_DataType sum( b[ i]);
        const t_DataType *decomp_row_i( m_Decomposition[ i]);
        for( size_t j( i + 1); j < n; j++)
        {
          sum -= decomp_row_i[ j] * b[ j];
        }
        if( m_Decomposition( i, i) < 1.0e-20)
        {
          b[ i] = 0.0;
        }
        else
        {
          b[ i] = sum / m_Decomposition( i, i);
        }
      }
      return solution;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer MatrixInversionCholesky< t_DataType>::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Computes matrix inverse using Cholesky decomposition"
        "This is usually preferred for least squares and similar symmetric problems. "
      );
      serializer.AddInitializer
      (
        "smoothing",
        "A small value added as an offset to the adjoint matrix as (1.0 + smoothing) * diagonal.  "
        "Higher values result in more numerically stable solutions and decrease the L2 norm of the solutions. "
        "Zero should only be used for very small and well-conditioned matrices",
        io::Serialization::GetAgent( &m_Smoothing),
        "0.00"
      );
      return serializer;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_CHOLESKY_HPP_
