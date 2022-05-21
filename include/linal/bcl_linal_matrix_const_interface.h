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

#ifndef BCL_LINAL_MATRIX_CONST_INTERFACE_H_
#define BCL_LINAL_MATRIX_CONST_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_reference.h"
#include "type/bcl_type_chooser.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixConstInterface
    //! @brief template
    //! @details The MatrixConstInterface provides general const functions for matrices
    //!
    //! @tparam t_DataType can be double, float, int, complex, etc.
    //!
    //! @remarks example unnecessary
    //! @author woetzen, loweew, mendenjl
    //! @date September 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixConstInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      typedef const t_DataType* const_iterator;
      typedef const t_DataType& const_reference;

      //! Data type for matrices and vectors returned by Inverse and Solve, respectively
      typedef typename type::Chooser< std::numeric_limits< t_DataType>::is_integer, double, t_DataType>::Type FloatType;

      //! @brief get number of rows
      //! @return number of rows
      virtual size_t GetNumberRows() const = 0;

      //! @brief get number of columns
      //! @return number of columns
      virtual size_t GetNumberCols() const = 0;

      //! @brief number of elements
      //! @return total number of elements in matrix
      virtual size_t GetNumberOfElements() const = 0;

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Matrix
      virtual const_iterator Begin() const = 0;

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      virtual const_iterator End() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief is the matrix empty
      //! @return true if number rows or number cols is zero
      virtual bool IsEmpty() const;

      //! @brief is matrix a square matrix
      //! @return true if number of cols and rows are identical
      virtual bool IsSquare() const;

      //! check whether matrix is symmetric
      virtual bool IsSymmetric() const;

      //! @brief is matrix a diagonal matrix
      //! @return true if all but the elements in the diagonal are 0
      virtual bool IsDiagonal() const;

      //! @brief is matrix a tridiagonal matrix
      //! @return if diagonal and adjacent diagonals are filled and the rest is 0
      virtual bool IsTriDiagonal() const;

      //! @brief checks if matrix is defined
      //! @return bool - true if matrix is defined
      virtual bool IsDefined() const;

      //! @brief get a vector reference on a particular row
      //! @param ROW the row of interest
      virtual VectorConstReference< t_DataType> GetRow( const size_t &ROW) const;

      //! @brief copy a column into a vector
      //! @param COL the column of interest
      virtual Vector< t_DataType> GetCol( const size_t &COL) const;

      //! @brief get a vector containing the diagonal of this matrix
      virtual Vector< t_DataType> GetDiagonal() const;

      //! @brief get the sum of the matrix
      virtual t_DataType Sum() const;

      //! compute trace of matrix
      virtual t_DataType Trace() const;

      //! @brief Get the matrix as a vector
      virtual VectorConstReference< t_DataType> AsVector() const;

      //! @brief Return a transposed copy of the matrix
      //! @return A transposed copy of the matrix
      virtual Matrix< t_DataType> Transposed() const;

      //! @brief Get the inverse of the matrix; uses the fastest method appropriate for the matrix
      //! For integral matrices, first casts them to double and returned matrix is also a double
      virtual Matrix< FloatType> Inverse() const;

      //! @brief Solve the normal matrix equation Ax=B using the fastest available method
      //! @param B the rhs of the matrix equation
      //! @return x, the solution vector
      //! If this matrix is rectangular and over-determined, computes the least squares fit of A onto X
      virtual Vector< FloatType> Solve( const VectorConstInterface< FloatType> &B) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return reference to const  element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return const element defined bey ROW and COL number
      virtual const t_DataType &operator()( const size_t ROW, const size_t COL) const = 0;

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a constant pointer to the first member of that row
      virtual const_iterator operator[]( const size_t ROW) const = 0;

      //! @brief compare two matrices for equality
      //! @param MATRIX_RHS rhs matrix to compare to
      //! @return true if dimensions and content match
      bool operator ==( const MatrixConstInterface< t_DataType> &MATRIX_RHS) const;

      //! @brief compare two matrices for inequality
      //! @param MATRIX_RHS rhs matrix to compare to
      //! @return true if dimensions and content match
      bool operator !=( const MatrixConstInterface< t_DataType> &MATRIX_RHS) const;

    }; // template class MatrixConstInterface

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixConstInterface< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_CONST_INTERFACE_H_
