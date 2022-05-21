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

#ifndef BCL_OPENBLAS_OPERATIONS_H_
#define BCL_OPENBLAS_OPERATIONS_H_

// include the namespace header
#include "bcl_openblas.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_operations_cpu.h"
#include "linal/bcl_linal_operations_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_operations.h"

namespace bcl
{
  namespace openblas
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Operations
    //! @brief is an interface, that defines linear algebra operations that can be found in blas
    //! @details These operations include but are not limited to
    //! @li matrix matrix operations (*, +, -)
    //! @li vector vector operations (dot, *, +, -)
    //! @li matrix vector operations (*)
    //!
    //! @see @link example_openblas_operations.cpp @endlink
    //! @author vuot2
    //! @date Mar 13, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class Operations :
      public linal::OperationsCPU< t_DataType>
    {
    private:

    //////////
    // data //
    //////////

      //! instantiate the openblas flag
      static bool s_IsFlagInstantiated;

    public:

    //////////
    // data //
    //////////

      //! instance of this operations class in Operations enumerator
      static const typename linal::Operations< t_DataType>::EnumType e_OpenBlas;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new OperationsInterface
      Operations *Clone() const;

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
        const linal::VectorInterface< t_DataType> &VECTOR_A,
        const linal::VectorInterface< t_DataType> &VECTOR_B
      ) const;

      //! @brief outer product of two vectors
      //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
      //! @param VECTOR_U left hand side vector with m elements
      //! @param VECTOR_V right hand side vector with n elements
      //! @return Matrix with m*n elements (outer product of vector u and v)
      linal::Matrix< t_DataType> OuterProduct
      (
        const linal::VectorConstInterface< t_DataType> &VECTOR_U,
        const linal::VectorConstInterface< t_DataType> &VECTOR_V
      );

      //! @brief norm of a vector
      //! @param VECTOR vector
      //! @return the Euclidean norm of the VECTOR
      t_DataType Norm( const linal::VectorInterface< t_DataType> &VECTOR) const;

      //! @brief matrix * vector
      //! @param MATRIX matrix
      //! @param VECTOR vector
      //! @return product of MATRIX * VECTOR
      linal::Vector< t_DataType> Multiply
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const linal::VectorInterface< t_DataType> &VECTOR
      ) const;

      //! @brief matrix-matrix multiplication
      //! @param MATRIX_A first matrix
      //! @param MATRIX_B second matrix
      //! @return product of MATRIX_A * MATRIX_B
      linal::Matrix< t_DataType> Multiply
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX_A,
        const linal::MatrixConstInterface< t_DataType> &MATRIX_B
      ) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorEqualsVectorTimesMatrix
      (
        linal::VectorInterface< t_DataType> &STORAGE,
        const linal::VectorConstInterface< t_DataType> &FEATURE,
        const linal::MatrixConstInterface< t_DataType> &MATRIX
      ) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorPlusEqualsMatrixTimesVector
      (
        linal::VectorInterface< t_DataType> &STORAGE,
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const linal::VectorConstInterface< t_DataType> &FEATURE
      ) const;

    }; // class BlasLinalOperationsDoudle

    //! @brief command line flag to be used to set Logger over the command line
    //! @return ShPtr to a FlagInterface which is used to set Logger
    ShPtr< command::FlagInterface> &GetFlagOpenblas();

    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< double>;

  } // namespace openblas
} // namespace bcl

#endif // BCL_OPENBLAS_OPERATIONS_H_
