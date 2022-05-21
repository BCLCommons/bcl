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

#ifndef BCL_LINAL_MATRIX_HPP_
#define BCL_LINAL_MATRIX_HPP_

// include the header for this class
#include "bcl_linal_matrix.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector.h"
#include "bcl_linal_vector_const_interface.h"
#include "bcl_linal_vector_reference.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

    //! single instance of the Matrix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Matrix< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Matrix< t_DataType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix() :
      m_NumberRows( 0),
      m_NumberCols( 0),
      m_Data( NULL)
    {
    }

    //! @brief construct from dimension and possible filler
    //! @param NUMBER_ROWS number of rows in matrix
    //! @param NUMBER_COLS number of cols in matrix
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix
    (
      const size_t NUMBER_ROWS,
      const size_t NUMBER_COLS,
      const t_DataType &FILL_VALUE
    ) :
      m_NumberRows( NUMBER_ROWS),
      m_NumberCols( NUMBER_COLS),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // set all values to FILL_VALUE
      std::fill( m_Data, m_Data + m_NumberRows * m_NumberCols, FILL_VALUE);
    }

    //! @brief construct from dimension and pointer to data
    //! @param NUMBER_ROWS number of rows in matrix
    //! @param NUMBER_COLS number of cols in matrix
    //! @param DATA pointer to field of data
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix
    (
      const size_t NUMBER_ROWS,
      const size_t NUMBER_COLS,
      const t_DataType *DATA
    ) :
      m_NumberRows( NUMBER_ROWS),
      m_NumberCols( NUMBER_COLS),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // copy data
      std::copy( DATA, DATA + m_NumberRows * m_NumberCols, m_Data);
    }

    //! construct Matrix from number of rows, number of columns, and storage::Vector< t_DataType>
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix
    (
      const size_t &ROWS,
      const size_t &COLS,
      const storage::Vector< t_DataType> &STORAGEVECTOR
    ) :
      m_NumberRows( ROWS),
      m_NumberCols( COLS),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // copy data
      std::copy( STORAGEVECTOR.Begin(), STORAGEVECTOR.End(), m_Data);
    }

    //! @brief copy constructor from Matrix
    //! @param MATRIX matrix to be copied from
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix( const Matrix< t_DataType> &MATRIX) :
      m_NumberRows( MATRIX.m_NumberRows),
      m_NumberCols( MATRIX.m_NumberCols),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( MATRIX.m_Data, MATRIX.m_Data + m_NumberRows * m_NumberCols, m_Data);
    }

    //! @brief move constructor from Matrix
    //! @param MATRIX matrix to be copied from
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix( Matrix< t_DataType> && MATRIX) :
      m_NumberRows( MATRIX.m_NumberRows),
      m_NumberCols( MATRIX.m_NumberCols),
      m_Data( MATRIX.m_Data)
    {
      MATRIX.m_Data = nullptr;
      MATRIX.m_NumberCols = MATRIX.m_NumberRows = 0;
    }

    //! @brief copy constructor from MatrixInterface
    //! @param MATRIX_INTERFACE matrix interface to be copied from
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE) :
      m_NumberRows( MATRIX_INTERFACE.GetNumberRows()),
      m_NumberCols( MATRIX_INTERFACE.GetNumberCols()),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), m_Data);
    }

    //! @brief construct from matrix interface and add padding
    //! @param MATRIX_INTERFACE matrix interface to be copied from
    //! @param ROW_PADDING the number of additional hidden rows to pad the matrix
    //! @param COL_PADDING the number of additional hidden cols to pad the matrix
    //! @param FILL_VALUE assign every pad element to that value
    template< typename t_DataType>
    Matrix< t_DataType>::Matrix
    (
      const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE,
      const size_t ROW_PADDING,
      const size_t COL_PADDING,
      const t_DataType &FILL_VALUE
    ) :
      m_NumberRows( MATRIX_INTERFACE.GetNumberRows() + ROW_PADDING),
      m_NumberCols( MATRIX_INTERFACE.GetNumberCols() + COL_PADDING),
      m_Data( new t_DataType[ m_NumberRows * m_NumberCols])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // set all values to FILL_VALUE
      std::fill( m_Data, m_Data + m_NumberRows * m_NumberCols, FILL_VALUE);

      // copy into submatrix
      SetSubMatrix( MATRIX_INTERFACE, 0, 0);
    }

    //! @brief Clone function
    //! @return pointer to new Matrix< t_DataType>
    template< typename t_DataType>
    Matrix< t_DataType> *Matrix< t_DataType>::Clone() const
    {
      return new Matrix< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Matrix< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a specified number of rows to the matrix
    //! @param N_ROWS the number of rows to add
    //! @param FILL_VALUE value to initialize the additional rows with
    template< typename t_DataType>
    void Matrix< t_DataType>::AddNRows( const size_t &N_ROWS, const t_DataType &FILL_VALUE)
    {
      if( m_NumberCols == size_t( 0))
      {
        // no columns, just update the number of rows
        m_NumberRows += N_ROWS;
        return;
      }

      // allocate memory for new expanded matrix
      Matrix< t_DataType> addendum( m_NumberRows + N_ROWS, m_NumberCols, FILL_VALUE);

      // copy this matrix into the addendum, followed by MATRIX
      std::copy( m_Data, End(), addendum.Begin());

      // swap pointer of extended matrix with current matrix
      std::swap( m_Data, addendum.m_Data);

      // update number rows
      std::swap( m_NumberRows, addendum.m_NumberRows);
    }

    //! @brief append a matrix to the bottom of the current matrix
    //! @param MATRIX matrix to add after the current matrix
    template< typename t_DataType>
    void Matrix< t_DataType>::Append( const MatrixConstInterface< t_DataType> &MATRIX)
    {
      // If this matrix is empty, then copy the matrix directly
      if( GetNumberOfElements() == size_t( 0))
      {
        operator=( MATRIX);
        return;
      }
      else if( MATRIX.GetNumberOfElements() == size_t( 0))
      {
        return;
      }

      // determine total number rows
      BCL_Assert( MATRIX.GetNumberCols() == m_NumberCols, "The appended matrix has a different column number!");

      // allocate memory for new expanded matrix
      Matrix< t_DataType> addendum( m_NumberRows + MATRIX.GetNumberRows(), m_NumberCols);

      // copy this matrix into the addendum, followed by MATRIX
      t_DataType *old_end( std::copy( m_Data, End(), addendum.Begin()));
      std::copy( MATRIX.Begin(), MATRIX.End(), old_end);

      // swap pointer of extended matrix with current matrix
      std::swap( m_Data, addendum.m_Data);

      // update number rows
      std::swap( m_NumberRows, addendum.m_NumberRows);
    }

    //! @brief append a matrix to the bottom of the current matrix
    //! @param MATRICES container with matrices that are added onto the last row of the current matrix
    //! @return matrix with addendum
    template< typename t_DataType>
    void Matrix< t_DataType>::Append( const storage::Vector< Matrix< t_DataType> > &MATRICES)
    {
      // determine total number rows of all involved matrices
      size_t total_length( m_NumberRows);

      // if there are no matrices to append then return
      if( MATRICES.IsEmpty())
      {
        return;
      }

      // check whether initial matrix has zero cols and if so allow adding of a new matrix
      if( GetNumberOfElements() == size_t( 0))
      {
        m_NumberCols = MATRICES.FirstElement().GetNumberCols();
      }

      // iterate over all matrices and determine total number rows
      for
      (
        typename storage::Vector< Matrix< t_DataType> >::const_iterator
          itr( MATRICES.Begin()), itr_end( MATRICES.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Assert
        (
          !itr->GetNumberRows() || itr->GetNumberCols() == m_NumberCols,
          "The appended matrix has a different column number!"
        );

        // add number of rows
        total_length += itr->GetNumberRows();
      }

      // allocate memory for new expanded matrix
      Matrix< t_DataType> addendum( total_length, m_NumberCols);

      // ommit copying when m_Data has no content
      if( m_Data != NULL)
      {
        // add current matrix to expanded matrix
        std::copy( m_Data, m_Data + m_NumberRows * m_NumberCols, addendum.Begin());
      }

      // start counting rows
      size_t end_position( m_NumberRows * m_NumberCols);

      // iterate over all matrices and determine total number rows
      for
      (
        typename storage::Vector< Matrix< t_DataType> >::const_iterator
          itr( MATRICES.Begin()), itr_end( MATRICES.End());
        itr != itr_end;
        ++itr
      )
      {
        // add current matrix to expanded matrix
        std::copy( itr->Begin(), itr->End(), addendum.Begin() + end_position);

        // update end position
        end_position += itr->GetNumberRows() * itr->GetNumberCols();
      }

      // swap pointer of extended matrix with current matrix
      std::swap( m_Data, addendum.m_Data);

      // update number rows
      m_NumberRows = total_length;
    }

    //! @brief create a larger matrix from this matrix
    //! @param ROW_PADDING the number of additional rows to pad the matrix
    //! @param COL_PADDING the number of additional cols to pad the matrix
    //! @param FILL_VALUE assign every pad element to that value
    //! @return Matrix that has additional col on the right, and rows at the bottom, preserving this matrix as a
    //!         submatrix in the upper left corner
    template< typename t_DataType>
    Matrix< t_DataType> Matrix< t_DataType>::CreatePaddedMatrix
    (
      const size_t ROW_PADDING,
      const size_t COL_PADDING,
      const t_DataType &FILL_VALUE
    ) const
    {
      const size_t pad_num_cols( m_NumberCols + COL_PADDING);

      // create matrix
      Matrix< t_DataType> padded_matrix( m_NumberRows + ROW_PADDING, pad_num_cols, FILL_VALUE);

      // copy row by row
      for
      (
        t_DataType *data( m_Data), *data_end( m_Data + m_NumberRows * m_NumberCols), *target( padded_matrix.Begin());
        data != data_end;
        data += m_NumberCols, target += pad_num_cols
      )
      {
        std::copy
        (
          data,
          data + m_NumberCols,
          target
        );
      }

      // end
      return padded_matrix;
    }

    //! @brief create submatrix from this matrix
    //! @param NUMBER_ROWS number of rows in submatrix
    //! @param NUMBER_COLS number of cols in submatrix
    //! @param POS_ROW starting row position in matrix - default 0
    //! @param POS_COL starting col position in matrix - default 0
    //! @return Matrix that is a submatrix of this
    template< typename t_DataType>
    Matrix< t_DataType> Matrix< t_DataType>::CreateSubMatrix
    (
      const size_t NUMBER_ROWS,
      const size_t NUMBER_COLS,
      const size_t POS_ROW,
      const size_t POS_COL
    ) const
    {
      // check that there are enough cols and rows
      BCL_Assert
      (
        NUMBER_ROWS + POS_ROW <= m_NumberRows && NUMBER_COLS + POS_COL <= m_NumberCols,
        "this matrix is too small for submatrix"
      );

      // if row and cols are equal (pos have to be 0 as asserted)
      if( NUMBER_ROWS == m_NumberRows && NUMBER_COLS == m_NumberCols)
      {
        return *this;
      }

      // create matrix
      Matrix< t_DataType> sub_matrix( NUMBER_ROWS, NUMBER_COLS);

      // copy row by row
      for
      (
        t_DataType *data( m_Data + POS_ROW * m_NumberCols), *data_end( m_Data + ( POS_ROW + NUMBER_ROWS) * m_NumberCols), *target( sub_matrix.Begin());
        data != data_end;
        data += m_NumberCols, target += NUMBER_COLS
      )
      {
        std::copy
        (
          data + POS_COL,
          data + POS_COL + NUMBER_COLS,
          target
        );
      }

      // end
      return sub_matrix;
    }

    //! @brief set a submatrix to values of given matrix
    //! @param MATRIX_INTERFACE the source of the data
    //! @param POS_ROW starting row position of matrix - default 0
    //! @param POS_COL starting col position of matrix - default 0
    //! @return ref to this
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::SetSubMatrix
    (
      const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE,
      const size_t POS_ROW,
      const size_t POS_COL
    )
    {
      // check that there are enough cols and rows
      BCL_Assert
      (
        MATRIX_INTERFACE.GetNumberRows() + POS_ROW <= m_NumberRows && MATRIX_INTERFACE.GetNumberCols() + POS_COL <= m_NumberCols,
        "this matrix is too small for given submatrix"
      );

      // if cols are equal
      if( MATRIX_INTERFACE.GetNumberCols() == m_NumberCols)
      {
        std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), m_Data + POS_ROW * m_NumberCols);

        // return ref to this
        return *this;
      }

      t_DataType *target( m_Data + POS_ROW * m_NumberCols);
      // copy row by row
      for
      (
        const t_DataType *data( MATRIX_INTERFACE.Begin()), *data_end( MATRIX_INTERFACE.End());
        data != data_end;
        data += MATRIX_INTERFACE.GetNumberCols(), target += m_NumberCols
      )
      {
        std::copy
        (
          data,
          data + MATRIX_INTERFACE.GetNumberCols(),
          target + POS_COL
        );
      }

      // return ref to this
      return *this;
    }

    //! @brief discard unnecessary rows from the matrix (without reallocating it)
    //! @param NUMBER_ROWS new number of rows in matrix
    template< typename t_DataType>
    void Matrix< t_DataType>::ShrinkRows( const size_t NUMBER_ROWS)
    {
      BCL_Assert( NUMBER_ROWS <= m_NumberRows, "Can't shrink to more rows than matrix has");
      m_NumberRows = NUMBER_ROWS;
    }

    //! @brief remove the given list of rows (without reallocating memory)
    //! @param ROWS list of row ids to remove
    template< typename t_DataType>
    void Matrix< t_DataType>::RemoveRows( const storage::Vector< size_t> &ROWS)
    {
      // nothing to do if the rows are already sorted
      if( ROWS.IsEmpty())
      {
        return;
      }

      // create a sorted version of the rows vector
      storage::Vector< size_t> sorted_rows( ROWS);
      sorted_rows.Sort( std::less< size_t>());

      BCL_Assert( sorted_rows.LastElement() < m_NumberRows, "Cannot remove rows beyond end of matrix");

      // add a placeholder to the end of the matrix.  This avoids some repetitive need to check for the end of the
      // vector below
      sorted_rows.PushBack( m_NumberRows);

      // compute # of lines removed; decrement if repeated lines are given
      size_t number_lines_removed( ROWS.GetSize());

      // in-place removal of rows
      storage::Vector< size_t>::const_iterator itr_remove( sorted_rows.Begin()), itr_remove_end( sorted_rows.End());
      storage::Vector< size_t>::const_iterator itr_remove_next( itr_remove + 1);
      t_DataType *itr_insert( operator[]( *itr_remove));
      for( ; itr_remove_next != itr_remove_end; ++itr_remove, ++itr_remove_next)
      {
        // determine the next apparent line to keep
        const size_t next_to_keep( *itr_remove + 1);

        // handle consecutive line removal
        if( next_to_keep == *itr_remove_next)
        {
          continue;
        }
        else if( *itr_remove == *itr_remove_next)
        {
          --number_lines_removed;
          continue;
        }

        // copy lines between the two lines to be removed
        itr_insert =
          std::copy( m_Data + next_to_keep * m_NumberCols, m_Data + *itr_remove_next * m_NumberCols, itr_insert);
      }

      // update number of rows
      m_NumberRows -= number_lines_removed;
    }

    //! @brief remove the given list of columns (without reallocating memory)
    //! @param COLS list of column ids to remove
    template< typename t_DataType>
    void Matrix< t_DataType>::RemoveCols( const storage::Vector< size_t> &COLS)
    {
      // handle trivial case
      if( COLS.IsEmpty())
      {
        return;
      }

      // create hash vector indicating which columns to keep
      storage::Vector< int> should_keep( m_NumberCols, int( 1));
      for( auto itr( COLS.Begin()), itr_end( COLS.End()); itr != itr_end; ++itr)
      {
        should_keep( *itr) = 0;
      }

      // iterate through the matrix and copy as needed
      auto itr_place( m_Data);
      auto itr_should_keep_end( should_keep.End());
      auto itr_should_keep( should_keep.Begin());
      for( auto itr( m_Data), itr_end( End()); itr != itr_end; ++itr, ++itr_should_keep)
      {
        if( itr_should_keep == itr_should_keep_end)
        {
          itr_should_keep = should_keep.Begin();
        }
        if( *itr_should_keep)
        {
          *itr_place = *itr;
          ++itr_place;
        }
      }
      m_NumberCols = std::accumulate( should_keep.Begin(), should_keep.End(), size_t( 0));
    }

    //! @brief keep only the given list of rows (without reallocating memory)
    //! @param ROWS list of row ids to keep
    template< typename t_DataType>
    void Matrix< t_DataType>::KeepRows( const storage::Vector< size_t> &ROWS)
    {
      if( GetNumberCols() == size_t( 0))
      {
        // no columns -> empty matrix
        return;
      }

      // delete the matrix if no rows are to be keps
      if( ROWS.IsEmpty())
      {
        delete [] m_Data;
        m_Data = NULL;
        m_NumberRows = 0;
        return;
      }

      // create a sorted version of the rows vector
      storage::Vector< size_t> sorted_rows( ROWS);
      sorted_rows.Sort( std::less< size_t>());

      BCL_Assert( sorted_rows.LastElement() < m_NumberRows, "Cannot keep rows beyond end of matrix");

      // compute # of lines kept; decrement if repeated lines are given
      size_t number_lines_kept( ROWS.GetSize());

      // in-place removal of rows
      storage::Vector< size_t>::const_iterator itr_keep( sorted_rows.Begin()), itr_keep_end( sorted_rows.End());
      storage::Vector< size_t>::const_iterator itr_keep_next( itr_keep + 1);

      // do not copy any rows until the first row that will be skipped is found
      size_t insert_row( *itr_keep);
      if( insert_row == size_t( 0))
      {
        for( ; itr_keep_next != itr_keep_end; ++itr_keep, ++itr_keep_next)
        {
          // determine the next apparent line to remove
          insert_row = *itr_keep + 1;

          // ignore repeated line numbers
          if( *itr_keep == *itr_keep_next)
          {
            --number_lines_kept;
          }
          else if( insert_row != *itr_keep_next)
          {
            ++itr_keep;
            break;
          }
        }
        if( itr_keep_next == itr_keep_end)
        {
          ShrinkRows( *itr_keep + 1);
          return;
        }
      }
      else
      {
        insert_row = 0;
      }

      size_t last_line( *itr_keep);
      MatrixInterface< t_DataType>::ReplaceRow( insert_row++, MatrixInterface< t_DataType>::GetRow( *itr_keep));
      for( ++itr_keep; itr_keep != itr_keep_end; ++itr_keep, ++insert_row)
      {
        // handle consecutive line removal
        if( last_line == *itr_keep)
        {
          continue;
        }
        last_line = *itr_keep;

        // copy lines between the two lines to be removed
        MatrixInterface< t_DataType>::ReplaceRow( insert_row, MatrixInterface< t_DataType>::GetRow( *itr_keep));
      }

      // update number of rows
      m_NumberRows = insert_row;
    }

    //! @brief reorder rows in the matrix
    //! @param ROWS new order of the rows
    template< typename t_DataType>
    void Matrix< t_DataType>::ReorderRows( const storage::Vector< size_t> &ROWS)
    {
      BCL_Assert( ROWS.GetSize() == m_NumberRows, "Wrong # rows");
      if( !m_NumberRows)
      {
        // no rows, nothing to reorder
        return;
      }

      // create a vector to hold a given row of the matrix
      Vector< t_DataType> tmp_row( m_NumberCols);

      storage::Vector< size_t> rows( ROWS);

      // create the inverse mapping
      storage::Vector< size_t> inv_rows( ROWS.GetSize(), util::GetUndefined< size_t>());
      for( size_t i( 0); i < m_NumberRows; ++i)
      {
        BCL_Assert
        (
          rows( i) < m_NumberRows,
          "ReorderRows requires that all rows be maintained; problem row was: "
          + util::Format()( i) + " = " + util::Format()( rows( i))
        );
        BCL_Assert
        (
          !util::IsDefined( inv_rows( rows( i))),
          "ReorderRows requires that no rows are duplicated, but row " + util::Format()( rows( i))
          + " was"
        );
        inv_rows( rows( i)) = i;
      }
      // now rows(inv_rows(x)) = x

      // perform in-place reordering of the rows

      // Loop invariants:
      // data'[a] = data[rows[a]] for 0 <= a < leng
      // rows[inv_rows[a]] = a
      // Swap the data from position i and the position given in the map.
      // update rows and inv_rows
      // intuitively, rows is the vector that contains the row that we wish to swap here
      // inv_rows tells where to find that row

      // Example:
      // data     = {0,1,2,3,4,5}
      // rows     = {1,3,5,4,2,0}
      // inv_rows = {5,0,4,1,3,2}
      // After the first pass of the loop
      // data     = {1,0,2,3,4,5}
      // rows     = {0,3,5,4,2,1}
      // inv_rows = {0,5,4,1,3,2}

      // Step 1: swap(data[i],data[rows[i]])
      // Step 2: swap(rows[inv_rows[i]],rows[i])
      // Step 3: swap(inv_rows[i],inv_rows[rows[i]])
      // but we know that inv_rows[i'] = rows[i'] = i' for i' < i, so we can avoid
      // the second two swaps
      for( size_t i( 0); i < m_NumberRows; ++i)
      {
        if( rows( i) == i)
        {
          // row already is correct, continue
          continue;
        }
        // swap row i with rows(i)
        tmp_row = this->GetRow( i);
        this->ReplaceRow( i, this->GetRow( rows( i)));
        this->ReplaceRow( rows( i), tmp_row);

        // update the mapping vectors
        rows( inv_rows( i)) = rows( i);
        inv_rows( rows( i)) = inv_rows( i);
        inv_rows( i) = rows( i) = i;
      }
    }

    //! @brief transpose the matrix, in place if possible
    template< typename t_DataType>
    void Matrix< t_DataType>::Transpose()
    {
      if( GetNumberRows() == GetNumberCols())
      {
        // transpose in-place, if this is a square matrix
        MatrixInterface< t_DataType>::Transpose();
      }
      else
      {
        // set the matrix equal to the transpose
        *this = MatrixConstInterface< t_DataType>::Transposed();
      }
    }

    //! @brief sets all values in matrix to zero
    //! @return the matrix filled with zeroes
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::SetZero()
    {
      std::fill( m_Data, m_Data + ( m_NumberRows * m_NumberCols), t_DataType( 0));

      return *this;
    }

    //! @brief set the matrix to be the identity matrix; all offdiagonal elements == 0, all diagonal elements == 1
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::SetIdentity()
    {
      t_DataType *itr( m_Data);
      for( size_t i( 0); i < m_NumberRows; ++i)
      {
        for( size_t j( 0); j < m_NumberCols; ++j, ++itr)
        {
          *itr = t_DataType( i == j ? 1 : 0);
        }
      }
      return *this;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment from Matrix
    //! @param MATRIX the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::operator =( const Matrix< t_DataType> &MATRIX)
    {
      // copy all elements
      if( m_Data != MATRIX.m_Data)
      {
        // check that sizes match
        if( GetNumberOfElements() != MATRIX.GetNumberOfElements())
        {
          // delete m_Data
          delete[] m_Data;

          // reallocate
          m_NumberRows = MATRIX.m_NumberRows;
          m_NumberCols = MATRIX.m_NumberCols;
          m_Data = new t_DataType[ m_NumberRows * m_NumberCols];

          // check if allocation was successful
          BCL_Assert
          (
            m_Data || m_NumberRows == 0 || m_NumberCols == 0,
            "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
            GetStaticClassName< t_DataType>()
          );
        }
        else
        {
          m_NumberRows = MATRIX.m_NumberRows;
          m_NumberCols = MATRIX.m_NumberCols;
        }

        // copy elements
        std::copy( MATRIX.m_Data, MATRIX.m_Data + m_NumberRows * m_NumberCols, m_Data);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief move from Matrix
    //! @param MATRIX the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::operator =( Matrix< t_DataType> && MATRIX)
    {
      // steal data members
      if( m_Data != MATRIX.m_Data)
      {
        // delete m_Data
        delete[] m_Data;
        m_NumberRows = MATRIX.m_NumberRows;
        m_NumberCols = MATRIX.m_NumberCols;
        MATRIX.m_NumberRows = MATRIX.m_NumberCols = 0;
        m_Data = MATRIX.m_Data;
        MATRIX.m_Data = nullptr;
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from MatrixInterface
    //! @param MATRIX_INTERFACE the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
    {
      // copy all elements
      if( m_Data != MATRIX_INTERFACE.Begin())
      {
        // check that sizes match
        if( m_NumberRows != MATRIX_INTERFACE.GetNumberRows() || m_NumberCols != MATRIX_INTERFACE.GetNumberCols())
        {
          // delete m_Data
          delete[] m_Data;

          // reallocate
          m_NumberRows = MATRIX_INTERFACE.GetNumberRows();
          m_NumberCols = MATRIX_INTERFACE.GetNumberCols();
          m_Data = new t_DataType[ m_NumberRows * m_NumberCols];

          // check if allocation was successful
          BCL_Assert
          (
            m_Data || m_NumberRows == 0 || m_NumberCols == 0,
            "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
            GetStaticClassName< t_DataType>()
          );
        }

        std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), m_Data);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from value
    //! @param VALUE all elements are set to that value
    //! @return reference to this assigned Matrix
    template< typename t_DataType>
    Matrix< t_DataType> &Matrix< t_DataType>::operator =( const t_DataType &VALUE)
    {
      // set all element to given VALUE
      std::fill( m_Data, m_Data + m_NumberRows * m_NumberCols, VALUE);

      // return reference to this Vector
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Matrix< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_NumberRows, ISTREAM);
      io::Serialize::Read( m_NumberCols, ISTREAM);

      delete[] m_Data;
      m_Data = new t_DataType[ GetNumberOfElements()];

      // check if allocation was successful
      BCL_Assert
      (
        m_Data || m_NumberRows == 0 || m_NumberCols == 0,
        "unable to allocate memory for " + util::Format()( m_NumberRows * m_NumberCols) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::string templine;
      std::stringstream full_mat;
      ISTREAM >> std::ws;

      // read line-by-line into a string stream, which is then read into a matrix.
      // This is > 3x faster if the stream is coming from a compressed file than reading directly due to the peeks and
      // functions employed in reading numbers directly from a data stream, operations which are not as fast on
      // compressed streams. It is 10-20% faster on uncompressed streams as well. Only if the data is another string
      // stream would this be marginally slower
      t_DataType *ptr( Begin());
      for( size_t row( 0); row < m_NumberRows; ++row)
      {
        // read the line into a stream
        std::getline( ISTREAM, templine);
        full_mat.str( templine);
        for( t_DataType *ptr_end( ptr + m_NumberCols); ptr != ptr_end; ++ptr)
        {
          BCL_Assert( io::Serialize::Read( *ptr, full_mat), "unable to read element");
        }
        // clear the string stream -> allows for further reading after calling str()
        full_mat.clear();
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Matrix< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write NumberRows and -Cols
      io::Serialize::Write( m_NumberRows, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NumberCols, OSTREAM, INDENT) << '\n';

      //write date
      for( size_t i = 0; i < m_NumberRows;)
      {
        io::Serialize::InsertIndent( OSTREAM, INDENT);
        for( size_t j = 0; j < m_NumberCols; ++j)
        {
          io::Serialize::Write( operator()( i, j), OSTREAM) << '\t';
        }
        ++i;
        // print line brake except after last line
        if( i != m_NumberRows)
        {
          OSTREAM << '\n';
        }
      }

      //return
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_H_
