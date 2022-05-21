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

#ifndef BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_H_
#define BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_inversion_interface.h"
#include "io/bcl_io_serializer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixInversionGaussJordan
    //! @brief Matrix inverse using gauss-jordan elimination
    //!
    //! @tparam t_DataType can be double, float, int, complex, etc.
    //!
    //! @see @link example_linal_matrix_inversion_gauss_jordan.cpp @endlink
    //! @author mendenjl
    //! @date Jun 02, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixInversionGaussJordan :
      public MatrixInversionInterface< t_DataType>
    {

    private:

      Matrix< t_DataType> m_Inverse;          //!< The computed inverse
      bool                m_InverseIsDefined; //!< True if the inverse is defined
      bool                m_PartialPivoting;  //!< True if row pivoting can be used

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_PivotInstance;
      static const util::SiPtr< const util::ObjectInterface> s_NoPivotInstance;

    public:

      //! @brief Constructor from whether to perform pivoting
      //! @param PARTIAL_PIVOTING whether to enable pivoting; should be true except when the matrix is guaranteed to be
      //!        well conditioned (e.g. splines).  In those cases, there is a performance benefit to setting this
      //!        parameter to false
      MatrixInversionGaussJordan( const bool &PARTIAL_PIVOTING);

      //! @brief constructor from a matrix
      //! @param MATRIX matrix for which to compute an inverse
      //! @param PARTIAL_PIVOTING whether to enable pivoting; should be true except when the matrix is guaranteed to be
      //!        well conditioned (e.g. splines).  In those cases, there is a performance benefit to setting this
      //!        parameter to false
      MatrixInversionGaussJordan( const MatrixConstInterface< t_DataType> &MATRIX, const bool &PARTIAL_PIVOTING = true);

      //! @brief Clone function
      //! @return pointer to new MatrixInterface
      MatrixInversionGaussJordan *Clone() const;

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

    }; // template class MatrixInversionGaussJordan

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionGaussJordan< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInversionGaussJordan< float>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INVERSION_GAUSS_JORDAN_H_
