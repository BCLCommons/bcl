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

#ifndef BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_H_
#define BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_inversion_interface.h"
#include "bcl_linal_vector.h"
#include "io/bcl_io_serializer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixInversionMoorePenrose
    //! @brief Matrix pseudo-inverse; appropriate for rectangular and non-symmetric square matrices.
    //! @details Solves the eigen system of m * Transpose( m)
    //!
    //! @tparam t_DataType can be double, float, int, complex, etc.
    //!
    //! @see @link example_linal_matrix_inversion_moore_penrose.cpp @endlink
    //! @author mendenjl
    //! @date Jun 05, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixInversionMoorePenrose :
      public MatrixInversionInterface< t_DataType>
    {

    private:

      Matrix< t_DataType> m_Inverse;          //!< The computed inverse
      bool                m_InverseIsDefined; //!< True if the inverse is defined
      Vector< t_DataType> m_Eigenvalues;      //!< Eigenvalues computed as a side-effect of the computation

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

      //! @brief Default constructor
      MatrixInversionMoorePenrose();

      //! @brief constructor from a matrix
      //! @param MATRIX matrix for which to compute an inverse
      MatrixInversionMoorePenrose( const MatrixConstInterface< t_DataType> &MATRIX);

      //! @brief Clone function
      //! @return pointer to new MatrixInterface
      MatrixInversionMoorePenrose *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief determine whether the solution is well-defined for all input vectors
      //! @return true if the inverse is well-defined
      bool IsDefined() const;

      //! @brief Set the source matrix.  Derived classes should store the source matrix, or a decomposition thereof
      //! @param SOURCE the matrix for which inversions and solutions will be required
      //! @return bool indicating whether matrix decomposition and storage was successful
      //! Implementations may choose to compute the matrix in this function, store the matrix directly, or compute a
      //! decomposition of the matrix
      bool SetMatrix( const MatrixConstInterface< t_DataType> &A);

      //! @brief get the eigenvalues for the matrix
      //! @return the eigenvalues of the square-symmetric matrix
      const Vector< t_DataType> &GetEigenvalues() const;

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
      Matrix< t_DataType> ComputeInverse();

      //! @brief Solve Ax=B, where A is the matrix already given, x is the vector to determine, and B is the result vec
      //! @param B the result vector
      //! @return Solution vector and bool indicating whether the solution is unique.
      Vector< t_DataType> Solve( const VectorConstInterface< t_DataType> &B) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // template class MatrixInversionMoorePenrose

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionMoorePenrose< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionMoorePenrose< float>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_MOORE_PENROSE_H_
