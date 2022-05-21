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

#ifndef BCL_LINAL_MATRIX_INTERFACE_H_
#define BCL_LINAL_MATRIX_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_const_interface.h"
#include "bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixInterface
    //! @brief MatrixInterface extends MatrixConstInterface by adding non const operations on matrices
    //!
    //! @tparam t_DataType can be double, float, int, complex, etc
    //!
    //! @see @link example_linal_matrix_interface.cpp @endlink
    //! @author woetzen, loweew, mendenjl
    //! @date Aug 18, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixInterface :
      public MatrixConstInterface< t_DataType>
    {

    public:

    /////////////////
    // data access //
    /////////////////

      typedef t_DataType* iterator;
      typedef t_DataType& reference;

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Matrix
      virtual iterator Begin() = 0;

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Matrix
      virtual typename MatrixConstInterface< t_DataType>::const_iterator Begin() const = 0;

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Matrix
      virtual iterator End() = 0;

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Matrix
      virtual typename MatrixConstInterface< t_DataType>::const_iterator End() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief swap elements between two rows
      //! @param ROW_A the first row
      //! @param ROW_B the second row
      virtual void SwapRows( const size_t ROW_A, const size_t ROW_B);

      //! @brief swap elements between two columns
      //! @param COL_A the first column
      //! @param COL_B the second column
      virtual void SwapCols( const size_t COL_A, const size_t COL_B);

      //! @brief get a vector reference on a particular row
      //! @param ROW the row of interest
      virtual VectorConstReference< t_DataType> GetRow( const size_t &ROW) const;

      //! @brief get a vector reference on a particular row
      //! @param ROW the row of interest
      virtual VectorReference< t_DataType> GetRow( const size_t &ROW);

      //! @brief Get the matrix as a vector
      virtual VectorConstReference< t_DataType> AsVector() const;

      //! @brief Get the matrix as a vector
      virtual VectorReference< t_DataType> AsVector();

      //! copies elements of argument VECTOR into this object at position (ROW)
      virtual void ReplaceRow( const size_t ROW, const VectorConstInterface< t_DataType> &VECTOR);

      //! copies elements of argument VECTOR into this object at position (COL)
      virtual void ReplaceCol( const size_t COL, const VectorConstInterface< t_DataType> &VECTOR);

      //! @brief Set all elements in a column to a single value
      //! @param COL the column to fill
      //! @param VALUE the value to fill the column with
      virtual void SetCol( const size_t COL, const t_DataType &VALUE);

    ///////////////
    // operators //
    ///////////////

      //! @brief return reference to changeable element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return changable reference to the element defined bey ROW and COL number
      virtual reference operator()( const size_t ROW, const size_t COL) = 0;

      //! @brief mutable access to a particular row
      //! @param ROW the row to access
      //! @return a pointer to the first member of that row
      virtual iterator operator[]( const size_t ROW) = 0;

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a constant pointer to the first member of that row
      virtual typename MatrixConstInterface< t_DataType>::const_iterator operator[]( const size_t ROW) const = 0;

      //! @brief Fill a matrix with a scalar value
      //! @param SCALAR the scalar value to fill the matrix with
      virtual MatrixInterface< t_DataType> &operator =( const t_DataType &SCALAR);

      //! @brief Transpose the matrix; asserts that matrix is square
      virtual void Transpose();

    private:

      //! @brief undefined assigment operator.  Matrix interfaces should not be copied with a concrete instance of the
      //! class, to prevent data slicing.
      MatrixInterface< t_DataType> &operator =( const MatrixInterface< t_DataType> &MATRIX);

    }; // template class MatrixInterface

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MatrixInterface< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INTERFACE_H_
