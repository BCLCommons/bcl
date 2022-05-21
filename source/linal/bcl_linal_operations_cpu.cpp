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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "linal/bcl_linal_operations_cpu.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_const_interface.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new OperationsCPU
    template< typename t_DataType>
    OperationsCPU< t_DataType>  *OperationsCPU< t_DataType>::Clone() const
    {
      return new OperationsCPU( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &OperationsCPU< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief dot product of two vectors
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A * VECTOR_B
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::DotProduct
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    ) const
    {
      // check for identical size
      if( VECTOR_A.GetSize() != VECTOR_B.GetSize())
      {
        return util::GetUndefined< t_DataType>();
      }

      // inner product
      return std::inner_product( VECTOR_A.Begin(), VECTOR_A.End(), VECTOR_B.Begin(), t_DataType( 0));
    }

    //! @brief outer product of two vectors
    //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
    //! @param VECTOR_U left hand side vector with m elements
    //! @param VECTOR_V right hand side vector with n elements
    //! @return Matrix with m*n elements (outer product of vector u and v)
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::OuterProduct
    (
      const VectorConstInterface< t_DataType> &VECTOR_U,
      const VectorConstInterface< t_DataType> &VECTOR_V
    )
    {
      Matrix< t_DataType> outer_product( VECTOR_U.GetSize(), VECTOR_V.GetSize());

      // ptr to first element of first row in matrix
      t_DataType *ptr_matrix( outer_product.Begin());

      // each result row is the product of the current VECTOR_U element multiplied with VECTOR_V
      for( const t_DataType *ptr_u( VECTOR_U.Begin()), *ptr_u_end( VECTOR_U.End()); ptr_u != ptr_u_end; ++ptr_u)
      {
        std::transform
        (
          VECTOR_V.Begin(), VECTOR_V.End(),
          ptr_matrix,
          std::bind2nd( std::multiplies< t_DataType>(), *ptr_u)
        );
        // ptr_matrix points to first element in next row
        ptr_matrix += VECTOR_V.GetSize();
      }

      // end
      return outer_product;
    }

    //! @brief norm of a vector
    //! @param VECTOR vector
    //! @return the euclidean norm of the VECTOR
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Norm( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Norm();
    }

    //! @brief matrix-matrix multiplication
    //! @param MATRIX_A matrix to be multiplied
    //! @param MATRIX_B matrix to be multiplied
    //! @return resulting linal::Matrix< t_DataType>
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::Multiply
    (
      const MatrixConstInterface< t_DataType> &MATRIX_A,
      const MatrixConstInterface< t_DataType> &MATRIX_B
    ) const
    {
      return MATRIX_A * MATRIX_B;
    }

    //! @brief matrix-vector multiplication
    //! @param MATRIX matrix to be multiplied
    //! @param VECTOR vector to be multiplied
    //! @return resulting linal::Vector< t_DataType>
    template< typename t_DataType>
    Vector< t_DataType> OperationsCPU< t_DataType>::Multiply
    (
      const MatrixConstInterface< t_DataType> &MATRIX,
      const VectorConstInterface< t_DataType> &VECTOR
    ) const
    {
      return MATRIX * VECTOR;
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS list of vectors
    //! @return triangular matrix of values of the distances
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
    ) const
    {
      const size_t num_elements( LIST_VECTORS.GetSize());
      Matrix< t_DataType> results( num_elements, num_elements, util::GetUndefined< t_DataType>());
      // iterate over rows
      for( size_t row_num( 0); row_num < num_elements; ++row_num)
      {
        const VectorConstInterface< t_DataType> &vector_a( *LIST_VECTORS( row_num));
        // iterate over cols
        for( size_t col_num( row_num); col_num < num_elements; ++col_num)
        {
          // assign to matrix position
          results( row_num, col_num) = Distance( vector_a, *LIST_VECTORS( col_num));
        }
      }
      return results;
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS_A list of vectors
    //! @param LIST_VECTORS_B list of vectors
    //! @return matrix of values of the distances
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_B
    ) const
    {
      const size_t num_rows( LIST_VECTORS_A.GetSize());
      const size_t num_cols( LIST_VECTORS_B.GetSize());
      Matrix< t_DataType> results( num_rows, num_cols, util::GetUndefined< t_DataType>());
      // iterate over rows
      for( size_t row_num( 0); row_num < num_rows; ++row_num)
      {
        const VectorConstInterface< t_DataType> &vector_a( *LIST_VECTORS_A( row_num));
        // iterate over cols
        for( size_t col_num( 0); col_num < num_cols; ++col_num)
        {
          // assign to matrix position
          results( row_num, col_num) = Distance( vector_a, *LIST_VECTORS_B( col_num));
        }
      }
      return results;
    }

    //! @brief returns sum of all elements in vector
    //! @param VECTOR the vector from which to get the sum
    //! @return the sum
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Sum( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Sum();
    }

    //! @brief reduction kernel
    //! @param VECTOR the vector to reduce
    //! @return the reduced result
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Min( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Min();
    }

    //! @brief gets max
    //! @param VECTOR the vector in which to find the max
    //! @return the max value
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Max( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Max();
    }

    //! @brief gets the min and max for each column in a matrix
    //! @param MATRIX the matrix input
    //! @return a storage vector of math ranges with min and max for each column in the matrix
    template< typename t_DataType>
    storage::Vector< math::Range< t_DataType> > OperationsCPU< t_DataType>::MinMax( const MatrixConstInterface< t_DataType> &MATRIX) const
    {
      storage::Vector< math::Range< t_DataType> > ranges;
      ranges.AllocateMemory( MATRIX.GetNumberCols());
      const size_t number_data_pts( MATRIX.GetNumberRows());

      // find ranges of each column
      t_DataType tmp( 0);
      for( size_t col( 0), number_features( MATRIX.GetNumberCols()); col < number_features; ++col)
      {
        // set initial min and max
        t_DataType min( std::numeric_limits< t_DataType>::infinity());
        t_DataType max( -std::numeric_limits< t_DataType>::infinity());

        // iterate through all values in that feature col
        for( size_t row( 0); row < number_data_pts; ++row)
        {
          // set as min or max if qualifies
          tmp = MATRIX( row, col);
          if( tmp < min)
          {
            min = tmp;
          }
          else if( tmp > max)
          {
            max = tmp;
          }
        }
        // add to range vector
        ranges.PushBack( math::Range< t_DataType>( min, max));
      }
      return ranges;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void OperationsCPU< t_DataType>::VectorEqualsVectorTimesMatrix
    (
      VectorInterface< t_DataType> &STORAGE,
      const VectorConstInterface< t_DataType> &FEATURE,
      const MatrixConstInterface< t_DataType> &MATRIX
    ) const
    {
      linal::VectorEqualsVectorTimesMatrix( STORAGE, FEATURE, MATRIX);
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void OperationsCPU< t_DataType>::VectorPlusEqualsMatrixTimesVector
    (
      VectorInterface< t_DataType> &STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX,
      const VectorConstInterface< t_DataType> &FEATURE
    ) const
    {
      linal::VectorPlusEqualsMatrixTimesVector( STORAGE, MATRIX, FEATURE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &OperationsCPU< t_DataType>::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &OperationsCPU< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief instance of this operations class in Operations enumerator
    template< typename t_DataType>
    const typename Operations< t_DataType>::EnumType OperationsCPU< t_DataType>::e_CPU
    (
      Operations< t_DataType>::GetEnums().SetDefaultOperationsType
      (
        Operations< t_DataType>::GetEnums().AddEnum
        (
          "CPU", util::ShPtr< OperationsInterface< t_DataType> >( new OperationsCPU< t_DataType>())
        )
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API OperationsCPU< float>;
    template class BCL_API OperationsCPU< double>;

  } // namespace linal
} // namespace bcl
