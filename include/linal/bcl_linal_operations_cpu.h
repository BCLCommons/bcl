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

#ifndef BCL_LINAL_OPERATIONS_CPU_H_
#define BCL_LINAL_OPERATIONS_CPU_H_

// include header of this class
#include "bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_vector_const_interface.h"
#include "bcl_linal_vector_operations.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OperationsCPU
    //! @brief implements linear algebra operations on the cpu
    //!
    //! TODO: add an general comment to this class
    //!
    //! @tparam t_DataType TODO description of template parameter
    //!
    //! @remarks example unnecessary
    //!
    //! @author woetzen, loweew, vuot2
    //!
    //! @date Mar 14, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class OperationsCPU :
      public OperationsInterface< t_DataType>
    {
    private:

    //////////
    // data //
    //////////

      //! @brief instance of this operations class in Operations enumerator
      static const typename Operations< t_DataType>::EnumType e_CPU;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new OperationsCPU
      OperationsCPU< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief dot product of two vectors
      //! @param VECTOR_A first vector
      //! @param VECTOR_B second vector
      //! @return dot product of VECTOR_A * VECTOR_B
      t_DataType DotProduct
      (
        const VectorConstInterface< t_DataType> &VECTOR_A,
        const VectorConstInterface< t_DataType> &VECTOR_B
      ) const;

      //! @brief outer product of two vectors
      //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
      //! @param VECTOR_U left hand side vector with m elements
      //! @param VECTOR_V right hand side vector with n elements
      //! @return Matrix with m*n elements (outer product of vector u and v)
      Matrix< t_DataType> OuterProduct
      (
        const VectorConstInterface< t_DataType> &VECTOR_U,
        const VectorConstInterface< t_DataType> &VECTOR_V
      );

      //! @brief norm of a vector
      //! @param VECTOR vector
      //! @return the euclidean norm of the VECTOR
      t_DataType Norm( const VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief compute a matrix times its transpose, without the need to compute the transpose, e.g. M * M^T
      //! @param A the matrix for which to compute A * A^T
      //! @return M * M^T
      //template< typename t_DataType>
      //Matrix< t_DataType> MatrixTimesItselfTransposed( const MatrixConstInterface< t_DataType> &A);

      //! @brief A^T * A, without the need to compute A^T
      //! @param A the matrix for which to compute A^T * A
      //! @return A^T * A
      //template< typename t_DataType>
      //Matrix< t_DataType> MatrixTransposeTimesMatrix( const MatrixConstInterface< t_DataType> &A);

      //! @brief A^T * B, without the need to compute A^T
      //! @param A the matrix for which to compute A^T * B
      //! @param B the rhs of the * in A^T * B
      //! @return A^T * B
      /*template< typename t_DataType>
      Matrix< t_DataType> MatrixTransposeTimesMatrix
      (
        const MatrixConstInterface< t_DataType> &A,
        const MatrixConstInterface< t_DataType> &B
      ); */

      //! @brief matrix-matrix multiplication
      //! @param MATRIX_A matrix to be multiplied
      //! @param MATRIX_B matrix to be multiplied
      //! @return resulting linal::Matrix< t_DataType>
      Matrix< t_DataType> Multiply
      (
        const MatrixConstInterface< t_DataType> &MATRIX_A,
        const MatrixConstInterface< t_DataType> &MATRIX_B
      ) const;

      //! @brief matrix-vector multiplication
      //! @param MATRIX matrix to be multiplied
      //! @param VECTOR vector to be multiplied
      //! @return resulting linal::Vector< t_DataType>
      Vector< t_DataType> Multiply
      (
        const MatrixConstInterface< t_DataType> &MATRIX,
        const VectorConstInterface< t_DataType> &VECTOR
      ) const;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS list of vectors
      //! @return triangular matrix of values of the distances
      Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
      ) const;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS_A list of vectors
      //! @param LIST_VECTORS_B list of vectors
      //! @return matrix of values of the distances
      Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_B
      ) const;

      //! @brief returns sum of all elements in vector
      //! @param VECTOR the vector from which to get the sum
      //! @return the sum
      t_DataType Sum( const VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief reduction kernel
      //! @param VECTOR the vector to reduce
      //! @return the reduced result
      t_DataType Min( const VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief gets max
      //! @param VECTOR the vector in which to find the max
      //! @return the max value
      t_DataType Max( const VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief gets the min and max for each column in a matrix
      //! @param MATRIX the matrix input
      //! @return a storage vector of math ranges with min and max for each column in the matrix
      storage::Vector< math::Range< t_DataType> > MinMax( const MatrixConstInterface< t_DataType> &MATRIX) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorEqualsVectorTimesMatrix
      (
        VectorInterface< t_DataType> &STORAGE,
        const VectorConstInterface< t_DataType> &FEATURE,
        const MatrixConstInterface< t_DataType> &MATRIX
      ) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorPlusEqualsMatrixTimesVector
      (
        VectorInterface< t_DataType> &STORAGE,
        const MatrixConstInterface< t_DataType> &MATRIX,
        const VectorConstInterface< t_DataType> &FEATURE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class OperationsCPU

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API OperationsInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API OperationsInterface< double>;

  } // namespace linal
} // namespace bcl
#endif // BCL_LINAL_OPERATIONS_CPU_H_
