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

#ifndef BCL_LINAL_MATRIX_INVERSION_CHOLESKY_H_
#define BCL_LINAL_MATRIX_INVERSION_CHOLESKY_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_inversion_interface.h"
#include "io/bcl_io_serializer.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixInversionCholesky
    //! @brief Matrix inverse using the cholesky decomposition.  Fastest numerically stable (even without pivoting)
    //!        method for linear least squares problem and other cases where the matrix is guaranteed to be symmetric
    //!        and positive definite
    //!        This class computes the inverse only if the user requests it, since back-substituting with
    //!        the cholesky decomposition is as efficient as matrix multiplication.
    //!        For very ill-conditioned matrices, cholesky decomposition may provide lower-quality solutions than using
    //!        QR with householder reductions; however, for such rank deficient matrices, the decomposition itself will
    //!        usually fail, allowing the calling class to perform the more expensive QR operation.
    //!
    //! @tparam t_DataType can be double, float, int, complex, etc.
    //!
    //! @see @link example_linal_matrix_inversion_cholesky.cpp @endlink
    //! @see @link http://en.wikipedia.org/wiki/Cholesky_decomposition @endlink
    //! @author mendenjl
    //! @date Jun 03, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixInversionCholesky :
      public MatrixInversionInterface< t_DataType>
    {

    private:

      Matrix< t_DataType> m_Decomposition;          //!< The computed decomposition; is empty if the inverse was computed
      bool                m_DecompositionIsDefined; //!< True if the decomposition is defined
      t_DataType          m_Smoothing;              //!< multiplicative factor for diagonal to avoid ill-conditioning

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_TolerantInstance;

    public:

    //////////////////
    // construction //
    //////////////////

      //! @brief Default constructor
      //! @param SMOOTHING amount of smoothing to apply; higher values result in smaller L2 norm of the solutions vector,
      //!        and more robust solutions to numerical issues and ill-conditioning, but at the cost of slightly lowered
      //!        accuracy
      MatrixInversionCholesky( const t_DataType &SMOOTHING = 0.000);

      //! @brief constructor from a matrix
      //! @param MATRIX matrix for which to compute an inverse
      MatrixInversionCholesky( const MatrixConstInterface< t_DataType> &MATRIX);

      //! @brief Clone function
      //! @return pointer to new MatrixInterface
      MatrixInversionCholesky *Clone() const;

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

    }; // template class MatrixInversionCholesky

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionCholesky< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionCholesky< float>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_CHOLESKY_H_
