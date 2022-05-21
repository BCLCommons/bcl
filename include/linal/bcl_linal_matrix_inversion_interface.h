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

#ifndef BCL_LINAL_MATRIX_INVERSION_INTERFACE_H_
#define BCL_LINAL_MATRIX_INVERSION_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixInversionInterface
    //! @brief template
    //! @details The MatrixInversionInterface provides the Inverse and Solve functions for double, float, and complex
    //!          Integral and other data types do not generally make sense because the inverse of an integer-matrix is
    //!          integral only under special circumstances that are rare in our use (determinant == +-1)
    //!
    //! @see @link example_linal_matrix_inversion_interface.cpp @endlink
    //! @author mendenjl
    //! @date Jun 02, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixInversionInterface :
      public util::SerializableInterface
    {

    public:

      //! @brief Clone function
      //! @return pointer to new MatrixInterface
      virtual MatrixInversionInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief determine whether the solution is well-defined for all input vectors
      //! @return true if the inverse is well-defined
      virtual bool IsDefined() const = 0;

      //! @brief Set the source matrix.  Derived classes should store the source matrix, or a decomposition thereof
      //! @param SOURCE the matrix for which inversions and solutions will be required
      //! @return bool indicating whether matrix decomposition and storage was successful
      //! Implementations are free to store a reference to the matrix
      virtual bool SetMatrix( const MatrixConstInterface< t_DataType> &A) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief Get/Compute the inverse matrix, which shall be empty if matrix inversion was unsuccessful
      //! Implementations should not perform this operation repeatedly on the same matrix; just cache whether the
      //! matrix in SetMatrix was inverted
      //! @note Matrix inverse should not be used to solve Ax=B for any small number ofB vectors because it is slower,
      //!       less stable, and requires more memory than just solving the equations directly; for those cases, use
      //!       Solve below.
      //! @return Inverse matrix
      virtual Matrix< t_DataType> ComputeInverse() = 0;

      //! @brief Solve Ax=B, where A is the matrix already given, x is the vector to determine, and B is the result vec
      //! @param B the result vector
      //! @return Solution vector
      virtual Vector< t_DataType> Solve( const VectorConstInterface< t_DataType> &B) const = 0;

      //! @brief Invert diagonal matrix
      //! @param INVERSE_STORAGE storage for the inverse matrix; must have the same size as MATRIX!
      //! @param MATRIX the matrix to try and invert as a diagonal matrix
      //! @return true if an inverse of the diagonal was computed successfully, false if there are zeros on the diagonal
      static bool TryInvertDiagonalMatrix
      (
        MatrixInterface< t_DataType> &INVERSE_STORAGE,
        const MatrixConstInterface< t_DataType> &MATRIX
      );

      //! @brief Invert a tridiagonal matrices
      //! @param INVERSE_STORAGE storage for the inverse matrix; must have the same size as MATRIX!
      //! @param MATRIX the matrix to try and invert as a tridiagonal matrix
      //! @return true if an inverse of the tridiagonal was computed successfully; false otherwise
      static bool TryInvertTridiagonalMatrix
      (
        MatrixInterface< t_DataType> &INVERSE_STORAGE,
        const MatrixConstInterface< t_DataType> &MATRIX
      );

      //! @brief Solve a diagonal matrix equation Ax=B
      //! @param A should be the original matrix
      //! @param B the desired output
      //! @note  There is no error checking for this function to keep it O(N), so only use it with known Diagonal matrices!
      static Vector< t_DataType> SolveDiagonalMatrix
      (
        const MatrixConstInterface< t_DataType> &A,
        const VectorConstInterface< t_DataType> &B
      );

      //! @brief Solve a tridiagonal matrix; this should be preferred to taking the inverse whenever possible because
      //!        the algorithm is more numerically stable and faster (O(N) vs. O(N^2)).  There is no error checking for
      //!        this function to keep it O(N), so only use it with known TriDiagonal matrices!
      //! @param A should be the original matrix
      //! @param B rhs of the equation
      //! @return Solution vector
      static Vector< t_DataType> SolveTridiagonalMatrix
      (
        const MatrixConstInterface< t_DataType> &A,
        const VectorConstInterface< t_DataType> &B
      );

    }; // template class MatrixInversionInterface

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionInterface< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionInterface< float>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_INTERFACE_H_
