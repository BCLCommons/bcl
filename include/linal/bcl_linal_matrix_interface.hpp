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

#ifndef BCL_LINAL_MATRIX_INTERFACE_HPP_
#define BCL_LINAL_MATRIX_INTERFACE_HPP_

// include the header of this class
#include "bcl_linal_matrix_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_reference.h"
#include "bcl_linal_vector_reference.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief swap elements between two rows
    //! @param ROW_A the first row
    //! @param ROW_B the second row
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::SwapRows( const size_t ROW_A, const size_t ROW_B)
    {
      const size_t n_rows( this->GetNumberRows());
      BCL_Assert
      (
        ROW_A < n_rows && ROW_B < n_rows,
        "Tried to swap row #" + util::Format()( std::max( ROW_A, ROW_B)) + " on matrix with only "
        + util::Format()( n_rows) + " rows"
      );

      t_DataType *row_a = operator[]( ROW_A);

      std::swap_ranges( row_a, row_a + this->GetNumberCols(), operator[]( ROW_B));
    }

    //! @brief swap elements between two columns
    //! @param COL_A the first column
    //! @param COL_B the second column
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::SwapCols( const size_t COL_A, const size_t COL_B)
    {
      const size_t n_cols( this->GetNumberCols());

      BCL_Assert
      (
        COL_A < n_cols && COL_B < n_cols,
        "Tried to swap column #" + util::Format()( std::max( COL_A, COL_B)) + " on matrix with only "
        + util::Format()( n_cols) + " columns"
      );

      //swap each pair in cols
      for
      (
        iterator col_a( Begin() + COL_A), col_b( Begin() + COL_B), end( End());
        col_a < end;
        col_a += n_cols, col_b += n_cols
      )
      {
        std::swap( *col_a, *col_b);
      }
    }

    //! @brief Set all elements in a column to a single value
    //! @param COL the column to fill
    //! @param VALUE the value to fill the column with
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::SetCol( const size_t COL, const t_DataType &VALUE)
    {
      const size_t n_cols( this->GetNumberCols());
      for( iterator itr( Begin() + COL), end( End()); itr < end; itr += n_cols)
      {
        *itr = VALUE;
      }
    }

    //! @brief get a vector reference on a particular row
    //! @param ROW the row of interest
    template< typename t_DataType>
    VectorConstReference< t_DataType> MatrixInterface< t_DataType>::GetRow( const size_t &ROW) const
    {
      return
        VectorConstReference< t_DataType>
        (
          this->GetNumberCols(),
          this->operator[]( ROW)
        );
    }

    //! @brief get a vector reference on a particular row
    //! @param ROW the row of interest
    template< typename t_DataType>
    VectorReference< t_DataType> MatrixInterface< t_DataType>::GetRow( const size_t &ROW)
    {
      return VectorReference< t_DataType>( this->GetNumberCols(), operator[]( ROW));
    }

    //! copies elements of argument VECTOR into this object at position (ROW)
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::ReplaceRow
    (
      const size_t ROW,
      const VectorConstInterface< t_DataType> &VECTOR
    )
    {
      if( VECTOR.IsEmpty())
      {
        return;
      }
      BCL_Assert
      (
        ROW < this->GetNumberRows(),
        "ReplaceRow given row # " + util::Format()( ROW)
        + " on matrix with only " + util::Format()( this->GetNumberRows()) + " rows"
      );
      BCL_Assert
      (
        VECTOR.GetSize() <= this->GetNumberCols(),
        "ReplaceRow with vector of " + util::Format()( VECTOR.GetSize())
        + " elements on matrix with only " + util::Format()( this->GetNumberCols()) + " columns"
      );
      std::copy( VECTOR.Begin(), VECTOR.End(), operator[]( ROW));
    }

    //! copies elements of argument VECTOR into this object at position (COL)
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::ReplaceCol
    (
      const size_t COL,
      const VectorConstInterface< t_DataType> &VECTOR
    )
    {
      if( VECTOR.IsEmpty())
      {
        return;
      }
      const size_t number_rows( this->GetNumberRows());
      const size_t number_cols( this->GetNumberCols());
      BCL_Assert
      (
        COL < number_cols,
        "ReplaceCol given col # " + util::Format()( COL)
        + " on matrix with only " + util::Format()( number_cols) + " cols"
      );

      BCL_Assert
      (
        VECTOR.GetSize() <= number_rows,
        "ReplaceCol with vector of " + util::Format()( VECTOR.GetSize())
        + " elements on matrix with only " + util::Format()( number_rows) + " rows"
      );

      const t_DataType *dat( VECTOR.Begin()), *dat_end( VECTOR.End());
      // copy elements
      for( size_t i( 0); dat != dat_end; ++i, ++dat)
      {
        operator()( i, COL) = ( *dat);
      }
    }

    //! @brief Get the matrix as a vector
    template< typename t_DataType>
    VectorConstReference< t_DataType> MatrixInterface< t_DataType>::AsVector() const
    {
      return VectorConstReference< t_DataType>( this->GetNumberOfElements(), this->Begin());
    }

    //! @brief Get the matrix as a vector
    template< typename t_DataType>
    VectorReference< t_DataType> MatrixInterface< t_DataType>::AsVector()
    {
      return VectorReference< t_DataType>( this->GetNumberOfElements(), this->Begin());
    }

    //! @brief Fill a matrix with a scalar value
    //! @param SCALAR the scalar value to fill the matrix with
    template< typename t_DataType>
    MatrixInterface< t_DataType> &MatrixInterface< t_DataType>::operator =( const t_DataType &SCALAR)
    {
      std::fill( Begin(), End(), SCALAR);
      return *this;
    }

    //! @brief Transpose the matrix; asserts that matrix is square
    template< typename t_DataType>
    void MatrixInterface< t_DataType>::Transpose()
    {
      BCL_Assert
      (
        this->IsSquare(),
        this->GetClassIdentifier() + " needs to override Transpose to handle non-square matrices!"
      );
      const size_t rows( this->GetNumberRows());
      for( size_t i( 0); i < rows; ++i)
      {
        t_DataType *transposed_row_i( operator[]( i));
        for( size_t j( 0); j < i; ++j)
        {
          std::swap( transposed_row_i[ j], operator()( j, i));
        }
      }
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_INTERFACE_HPP_
