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

#ifndef BCL_LINAL_OPERATIONS_INTERFACE_H_
#define BCL_LINAL_OPERATIONS_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_interface.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OperationsInterface
    //! @brief is an interface, that defines linear algebra operations that can be found in blas
    //! @details These operations include but are not limited to
    //! @li matrix matrix operations (*, +, -)
    //! @li vector vector operations (dot, *, +, -)
    //! @li matrix vector operations (*)
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class OperationsInterface :
      public util::ObjectInterface
    {

    public:
    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new OperationsInterface
      virtual OperationsInterface< t_DataType> *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief dot product of two vectors
      //! @param VECTOR_A first vector
      //! @param VECTOR_B second vector
      //! @return dot product of VECTOR_A * VECTOR_B
      virtual t_DataType DotProduct
      (
        const VectorConstInterface< t_DataType> &VECTOR_A,
        const VectorConstInterface< t_DataType> &VECTOR_B
      ) const = 0;

      //! @brief norm of a vector
      //! @param VECTOR vector
      //! @return the Euclidean norm of the VECTOR
      virtual t_DataType Norm( const VectorConstInterface< t_DataType> &VECTOR) const = 0;

      //! @brief matrix * vector
      //! @param MATRIX matrix
      //! @param VECTOR vector
      //! @return product of MATRIX * VECTOR
      virtual Vector< t_DataType> Multiply
      (
        const MatrixConstInterface< t_DataType> &MATRIX,
        const VectorConstInterface< t_DataType> &VECTOR
      ) const = 0;

      //! @brief matrix-matrix multiplication
      //! @param MATRIX_A first matrix
      //! @param MATRIX_B second matrix
      //! @return product of MATRIX_A * MATRIX_B
      virtual Matrix< t_DataType> Multiply
      (
        const MatrixConstInterface< t_DataType> &MATRIX_A,
        const MatrixConstInterface< t_DataType> &MATRIX_B
      ) const = 0;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS list of vectors
      //! @return triangular matrix of values of the distances
      virtual Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
      ) const = 0;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS_A list of vectors
      //! @param LIST_VECTORS_B list of vectors
      //! @return matrix of values of the distances
      virtual Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_B
      ) const = 0;

      //! @brief reduction kernel
      //! @param VECTOR the vector to reduce
      //! @return the reduced result
      virtual t_DataType Sum( const VectorConstInterface< t_DataType> &VECTOR) const = 0;

      //! @brief reduction kernel
      //! @param VECTOR the vector to reduce
      //! @return the reduced result
      virtual t_DataType Min( const VectorConstInterface< t_DataType> &VECTOR) const = 0;

      //! @brief gets max
      //! @param VECTOR the vector in which to find the max
      //! @return the max value
      virtual t_DataType Max( const VectorConstInterface< t_DataType> &VECTOR) const = 0;

      //! @brief gets the min and max for each column in a matrix
      //! @param MATRIX the matrix input
      //! @return a storage vector of math ranges with min and max for each column in the matrix
      virtual storage::Vector< math::Range< t_DataType> > MinMax( const MatrixConstInterface< t_DataType> &MATRIX) const = 0;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      virtual void VectorEqualsVectorTimesMatrix
      (
        VectorInterface< t_DataType> &STORAGE,
        const VectorConstInterface< t_DataType> &FEATURE,
        const MatrixConstInterface< t_DataType> &MATRIX
      ) const = 0;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      virtual void VectorPlusEqualsMatrixTimesVector
      (
        VectorInterface< t_DataType> &STORAGE,
        const MatrixConstInterface< t_DataType> &MATRIX,
        const VectorConstInterface< t_DataType> &FEATURE
      ) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief creates a matrix from a SiPtrVector of VectorInterfaces
      //! @param LIST_VECTORS the SiPtrVector of VectorInterfaces to be transferred into a matrix
      //! @return the matrix created from the SiPtrList of VectorInterfaces
      static Matrix< t_DataType> CopyListToMatrix
      (
        const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
      );

    }; // class OperationsInterface

    BCL_EXPIMP_TEMPLATE template class BCL_API OperationsInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API OperationsInterface< double>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_OPERATIONS_INTERFACE_H_
