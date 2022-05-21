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

#ifndef BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_HPP_
#define BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_HPP_

// include the header of this class
#include "bcl_linal_matrix_inversion_moore_penrose.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_symmetric_eigensolver.h"
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

    //! @brief Default constructor
    template< typename t_DataType>
    MatrixInversionMoorePenrose< t_DataType>::MatrixInversionMoorePenrose()
    {
    }

    //! @brief constructor from a matrix
    //! @param MATRIX matrix for which to compute an inverse
    template< typename t_DataType>
    MatrixInversionMoorePenrose< t_DataType>::MatrixInversionMoorePenrose
    (
      const MatrixConstInterface< t_DataType> &MATRIX
    )
    {
      SetMatrix( MATRIX);
    }

    //! @brief Clone function
    //! @return pointer to new MatrixInterface
    template< typename t_DataType>
    MatrixInversionMoorePenrose< t_DataType> *MatrixInversionMoorePenrose< t_DataType>::Clone() const
    {
      return new MatrixInversionMoorePenrose< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &MatrixInversionMoorePenrose< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    template< typename t_DataType>
    const std::string &MatrixInversionMoorePenrose< t_DataType>::GetAlias() const
    {
      static const std::string s_alias( "MoorePenrose");
      return s_alias;
    }

    //! @brief determine whether the solution is well-defined for all input vectors
    //! @return true if the inverse is well-defined
    template< typename t_DataType>
    bool MatrixInversionMoorePenrose< t_DataType>::IsDefined() const
    {
      return true;
    }

    //! @brief Set the source matrix.  Derived classes should store the source matrix, or a decomposition thereof
    //! @param SOURCE the matrix for which inversions and solutions will be required
    //! @return bool indicating whether matrix decomposition and storage was successful
    //! Implementations may choose to compute the matrix in this function, store the matrix directly, or compute a
    //! decomposition of the matrix
    template< typename t_DataType>
    bool MatrixInversionMoorePenrose< t_DataType>::SetMatrix( const MatrixConstInterface< t_DataType> &A)
    {
      const size_t nrows( A.GetNumberRows()), ncols( A.GetNumberCols());
      Matrix< t_DataType> eigenvalues, eigenvectors_v( nrows, ncols), eigenvectors_u( ncols, ncols);

      // perform SVD
      eigenvalues = SingularValueDecomposition( A, eigenvectors_v, eigenvectors_u);

      // compute eigenvectors_v
      m_Eigenvalues = eigenvalues.GetDiagonal();
      Matrix< t_DataType> inv_eigenvalues( eigenvalues);
      for( size_t i( 0), nr( eigenvalues.GetNumberRows()); i < nr; ++i)
      {
        inv_eigenvalues( i, i) = inv_eigenvalues( i, i) > std::numeric_limits< t_DataType>::min() ? 1.0 / inv_eigenvalues( i, i) : 0.0;
      }
      m_Inverse = ( eigenvectors_v.Transposed()) * inv_eigenvalues * ( eigenvectors_u.Transposed());

      // return
      return true;
    }

    //! @brief get the eigenvalues for the matrix
    //! @return the eigenvalues of the square-symmetric matrix
    template< typename t_DataType>
    const Vector< t_DataType> &MatrixInversionMoorePenrose< t_DataType>::GetEigenvalues() const
    {
      return m_Eigenvalues;
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
    Matrix< t_DataType> MatrixInversionMoorePenrose< t_DataType>::ComputeInverse()
    {
      return m_Inverse;
    }

    //! @brief Solve Ax=B, where A is the matrix already given, x is the vector to determine, and B is the result vec
    //! @param B the result vector
    //! @return Solution vector and bool indicating whether the solution is unique.
    template< typename t_DataType>
    Vector< t_DataType> MatrixInversionMoorePenrose< t_DataType>::Solve
    (
      const VectorConstInterface< t_DataType> &B
    ) const
    {
      return m_Inverse * B;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer MatrixInversionMoorePenrose< t_DataType>::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Computes matrix inverse using moore-penrose pseudo-inverse.  This is the most stable, albiet the most "
        " computationally costly of the matrix inversion methods; and is best reserved for rectangular matrices or "
        "very ill-conditioned matrices"
      );
      return serializer;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_HPP_
