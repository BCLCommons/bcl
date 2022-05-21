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

#ifndef BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_HPP_
#define BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_HPP_

// include the header of this class
#include "bcl_linal_matrix_inversion_gauss_jordan.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_vector_const_reference.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////
  // Data Access //
  /////////////////

    //! @brief Constructor from whether to perform pivoting
    //! @param PARTIAL_PIVOTING whether to enable pivoting; should be true except when the matrix is guaranteed to be
    //!        well conditioned (e.g. splines).  In those cases, there is a performance benefit to setting this
    //!        parameter to false
    template< typename t_DataType>
    MatrixInversionGaussJordan< t_DataType>::MatrixInversionGaussJordan( const bool &PARTIAL_PIVOTING) :
      m_InverseIsDefined( true),
      m_PartialPivoting( PARTIAL_PIVOTING)
    {
    }

    //! @brief constructor from a matrix
    //! @param MATRIX matrix for which to compute an inverse
    //! @param PARTIAL_PIVOTING whether to enable pivoting; should be true except when the matrix is guaranteed to be
    //!        well conditioned (e.g. splines).  In those cases, there is a performance benefit to setting this
    //!        parameter to false
    template< typename t_DataType>
    MatrixInversionGaussJordan< t_DataType>::MatrixInversionGaussJordan
    (
      const MatrixConstInterface< t_DataType> &MATRIX,
      const bool &PARTIAL_PIVOTING
    ) :
      m_InverseIsDefined( true),
      m_PartialPivoting( PARTIAL_PIVOTING)
    {
      SetMatrix( MATRIX);
    }

    //! @brief Clone function
    //! @return pointer to new MatrixInterface
    template< typename t_DataType>
    MatrixInversionGaussJordan< t_DataType> *MatrixInversionGaussJordan< t_DataType>::Clone() const
    {
      return new MatrixInversionGaussJordan< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &MatrixInversionGaussJordan< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    template< typename t_DataType>
    const std::string &MatrixInversionGaussJordan< t_DataType>::GetAlias() const
    {
      static const std::string s_alias( "GaussJordan"), s_alias_no_pivot( "GaussJordanWithoutPivoting");
      return m_PartialPivoting ? s_alias : s_alias_no_pivot;
    }

    //! @brief determine whether the solution is well-defined for all input vectors
    //! @return true if the inverse is well-defined
    template< typename t_DataType>
    bool MatrixInversionGaussJordan< t_DataType>::IsDefined() const
    {
      return m_InverseIsDefined;
    }

    //! @brief Set the source matrix.  Derived classes should store the source matrix, or a decomposition thereof
    //! @param SOURCE the matrix for which inversions and solutions will be required
    //! @return bool indicating whether matrix decomposition and storage was successful
    //! Implementations may choose to compute the matrix in this function, store the matrix directly, or compute a
    //! decomposition of the matrix
    template< typename t_DataType>
    bool MatrixInversionGaussJordan< t_DataType>::SetMatrix( const MatrixConstInterface< t_DataType> &A)
    {
      const size_t nrows( A.GetNumberRows());
      m_Inverse = Matrix< t_DataType>( nrows, A.GetNumberCols());
      if( !A.IsSquare())
      {
        // Non-square matrix detected; inverse is not defined by this method
        m_InverseIsDefined = false;
        return false;
      }

      // set unitary matrix
      m_Inverse.SetZero();
      for( size_t k( 0); k < nrows; ++k)
      {
        m_Inverse( k, k) = t_DataType( 1);
      }

      // make a copy of A that will contain the current state of A
      Matrix< t_DataType> a_temp( A);

      // gauss jordan elimination with partial (row) pivoting
      for( size_t k( 0); k < nrows; ++k)
      {
        // find the highest-magnitude remaining pivot in A for column k
        t_DataType amax( math::Absolute( a_temp( k, k)));

        if( m_PartialPivoting)
        {
          size_t best_pivot( k);
          for( size_t pivot( k + 1); pivot < nrows; ++pivot)
          {
            const t_DataType absolute_value( math::Absolute( a_temp( pivot, k)));
            if( absolute_value > amax)
            {
              amax = absolute_value;
              best_pivot = pivot;
            }
          }

          // swap the rows
          if( best_pivot != k)
          {
            a_temp.SwapRows( best_pivot, k);
            m_Inverse.SwapRows( best_pivot, k);
          }
        }

        if( a_temp( k, k) == t_DataType( 0))
        {
          m_InverseIsDefined = false;
          return false;
        }

        t_DataType a1 = a_temp( k, k);
        for( size_t j( 0); j < nrows; ++j)
        {
          a_temp( k, j) /= a1;
          m_Inverse( k, j) /= a1;
        }
        for( size_t i( 0); i < nrows; ++i)
        {
          if( i == k)
          {
            continue;
          }

          const t_DataType a2( a_temp( i, k));
          for( size_t j( 0); j < nrows; ++j)
          {
            a_temp( i, j)    -= a2 * a_temp( k, j);
            m_Inverse( i, j) -= a2 * m_Inverse( k, j);
          }
        }
      }
      return true;
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
    Matrix< t_DataType> MatrixInversionGaussJordan< t_DataType>::ComputeInverse()
    {
      return m_Inverse;
    }

    //! @brief Solve Ax=B, where A is the matrix already given, x is the vector to determine, and B is the result vec
    //! @param B the result vector
    //! @return Solution vector and bool indicating whether the solution is unique.
    template< typename t_DataType>
    Vector< t_DataType> MatrixInversionGaussJordan< t_DataType>::Solve
    (
      const VectorConstInterface< t_DataType> &B
    ) const
    {
      return m_Inverse * B;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer MatrixInversionGaussJordan< t_DataType>::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        std::string( "Computes matrix inverse using standard Gauss-Jordan elimination ")
        + ( m_PartialPivoting ? "with row pivoting (appropriate for most problems)" : "with no pivoting")
      );
      return serializer;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_HPP_
