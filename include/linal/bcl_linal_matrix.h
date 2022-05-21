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

#ifndef BCL_LINAL_MATRIX_H_
#define BCL_LINAL_MATRIX_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_interface.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Matrix
    //! @brief is a general implementation of a Matrix variable rows and cols.
    //! @details The matrix is stored as an array of t_DataType, where each row is stored consecutively.
    //!
    //! @see @link example_linal_matrix.cpp @endlink
    //! @author woetzen
    //! @date Aug 18, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Matrix :
      public MatrixInterface< t_DataType>
    {

    protected:

    //////////
    // data //
    //////////

      size_t m_NumberRows; //!< number of rows in matrix
      size_t m_NumberCols; //!< number of cols in matrix
      t_DataType *m_Data;  //!< array containing all data

    public:

    //////////////
    // typedefs //
    //////////////

      //! Typedefs imported from matrix interface classes
      typedef typename MatrixConstInterface< t_DataType>::const_iterator  const_iterator;
      typedef typename MatrixConstInterface< t_DataType>::const_reference const_reference;
      typedef typename MatrixInterface< t_DataType>::iterator             iterator;
      typedef typename MatrixInterface< t_DataType>::reference            reference;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Matrix();

      //! @brief construct from dimension and possible filler
      //! @param NUMBER_ROWS number of rows in matrix
      //! @param NUMBER_COLS number of cols in matrix
      //! @param FILL_VALUE assign every element to that value
      Matrix
      (
        const size_t NUMBER_ROWS,
        const size_t NUMBER_COLS,
        const t_DataType &FILL_VALUE = t_DataType( 0)
      );

      //! @brief construct from dimension and pointer to data
      //! @param NUMBER_ROWS number of rows in matrix
      //! @param NUMBER_COLS number of cols in matrix
      //! @param DATA pointer to field of data
      Matrix
      (
        const size_t NUMBER_ROWS,
        const size_t NUMBER_COLS,
        const t_DataType *DATA
      );

      //! construct Matrix from number of rows, number of columns, and storage::Vector< t_DataType>
      Matrix( const size_t &ROWS, const size_t &COLS, const storage::Vector< t_DataType> &STORAGEVECTOR);

      //! @brief copy constructor from Matrix
      //! @param MATRIX matrix to be copied from
      Matrix( const Matrix< t_DataType> &MATRIX);

      //! @brief move constructor from Matrix
      //! @param MATRIX matrix to be copied from
      Matrix( Matrix< t_DataType> && MATRIX);

      //! @brief copy constructor from MatrixInterface
      //! @param MATRIX_INTERFACE matrix interface to be copied from
      Matrix( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE);

      template< typename t_OtherDataType>
      explicit Matrix( const MatrixConstInterface< t_OtherDataType> &MATRIX) :
        m_NumberRows( MATRIX.GetNumberRows()),
        m_NumberCols( MATRIX.GetNumberCols()),
        m_Data( new t_DataType[ MATRIX.GetNumberOfElements()])
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data || m_NumberRows == 0 || m_NumberCols == 0,
          "unable to allocate memory for " + util::Format()( MATRIX.GetNumberOfElements()) + " elements of type: " +
          GetStaticClassName< t_DataType>()
        );

        // copy data
        std::copy( MATRIX.Begin(), MATRIX.End(), m_Data);
      }

      //! @brief construct from matrix interface and add padding
      //! @param MATRIX_INTERFACE matrix interface to be copied from
      //! @param ROW_PADDING the number of additional hidden rows to pad the matrix
      //! @param COL_PADDING the number of additional hidden cols to pad the matrix
      //! @param FILL_VALUE assign every pad element to that value
      Matrix
      (
        const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE,
        const size_t ROW_PADDING,
        const size_t COL_PADDING,
        const t_DataType &FILL_VALUE = t_DataType( 0)
      );

      //! @brief Clone function
      //! @return pointer to new Matrix< t_DataType>
      Matrix< t_DataType> *Clone() const;

      //! @brief destructor
      ~Matrix()
      {
        delete[] m_Data;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get number of rows
      //! @return number of rows
      size_t GetNumberRows() const
      {
        return m_NumberRows;
      }

      //! @brief get number of columns
      //! @return number of columns
      size_t GetNumberCols() const
      {
        return m_NumberCols;
      }

      //! @brief number of elements
      //! @return total number of elements in matrix
      size_t GetNumberOfElements() const
      {
        return m_NumberRows * m_NumberCols;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Matrix
      const_iterator Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Matrix
      iterator Begin()
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      const_iterator End() const
      {
        return m_Data + m_NumberRows * m_NumberCols;
      }

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Matrix
      iterator End()
      {
        return m_Data + m_NumberRows * m_NumberCols;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a specified number of rows to the matrix
      //! @param N_ROWS the number of rows to add
      //! @param FILL_VALUE value to initialize the additional rows with
      void AddNRows( const size_t &N_ROWS, const t_DataType &FILL_VALUE = t_DataType( 0));

      //! @brief append a matrix to the bottom of the current matrix
      //! @param MATRIX matrix to add after the current matrix
      void Append( const MatrixConstInterface< t_DataType> &MATRIX);

      //! @brief append a matrix to the bottom of the current matrix
      //! @param MATRICES container with matrices that are added onto the last row of the current matrix
      //! @return matrix with addendum
      void Append( const storage::Vector< Matrix< t_DataType> > &MATRICES);

      //! @brief create a larger matrix from this matrix
      //! @param ROW_PADDING the number of additional rows to pad the matrix
      //! @param COL_PADDING the number of additional cols to pad the matrix
      //! @param FILL_VALUE assign every pad element to that value
      //! @return Matrix that has additional col on the right, and rows at the bottom, preserving this matrix as a
      //!         submatrix in the upper left corner
      Matrix< t_DataType> CreatePaddedMatrix
      (
        const size_t ROW_PADDING,
        const size_t COL_PADDING,
        const t_DataType &FILL_VALUE = t_DataType( 0)
      ) const;

      //! @brief create submatrix from this matrix
      //! @param NUMBER_ROWS number of rows in submatrix
      //! @param NUMBER_COLS number of cols in submatrix
      //! @param POS_ROW starting row position in matrix - default 0
      //! @param POS_COL starting col position in matrix - default 0
      //! @return Matrix that is a submatrix of this
      Matrix< t_DataType> CreateSubMatrix
      (
        const size_t NUMBER_ROWS,
        const size_t NUMBER_COLS,
        const size_t POS_ROW = 0,
        const size_t POS_COL = 0
      ) const;

      //! @brief set a submatrix to values of given matrix
      //! @param MATRIX_INTERFACE the source of the data
      //! @param POS_ROW starting row position of matrix - default 0
      //! @param POS_COL starting col position of matrix - default 0
      //! @return ref to this
      Matrix< t_DataType> &SetSubMatrix
      (
        const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE,
        const size_t POS_ROW = 0,
        const size_t POS_COL = 0
      );

      //! @brief discard unnecessary rows from the matrix (without reallocating it)
      //! @param NUMBER_ROWS new number of rows in matrix
      void ShrinkRows( const size_t NUMBER_ROWS);

      //! @brief remove the given list of rows (without reallocating memory)
      //! @param ROWS list of row ids to remove
      void RemoveRows( const storage::Vector< size_t> &ROWS);

      //! @brief remove the given list of columns (without reallocating memory)
      //! @param COLS list of column ids to remove
      void RemoveCols( const storage::Vector< size_t> &COLS);

      //! @brief keep only the given list of rows (without reallocating memory)
      //! @param ROWS list of row ids to keep
      void KeepRows( const storage::Vector< size_t> &ROWS);

      //! @brief reorder rows in the matrix
      //! @param ROWS new order of the rows
      void ReorderRows( const storage::Vector< size_t> &ROWS);

      //! @brief transpose the matrix, in place if possible
      void Transpose();

      //! @brief sets all values in matrix to zero
      //! @return the matrix filled with zeroes
      Matrix< t_DataType> &SetZero();

      //! @brief set the matrix to be the identity matrix; all offdiagonal elements == 0, all diagonal elements == 1
      Matrix< t_DataType> &SetIdentity();

    ///////////////
    // operators //
    ///////////////

      //! @brief return reference to changeable element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return changeable reference to the element defined bey ROW and COL number
      reference operator()( const size_t ROW, const size_t COL)
      {
        BCL_Assert
        (
          ROW < m_NumberRows && COL < m_NumberCols,
          "requested (" + util::Format()( ROW) + "," + util::Format()( COL)
          + ") but matrix has only " + util::Format()( m_NumberRows) + " rows and " + util::Format()( m_NumberCols) + " cols"
        );
        return m_Data[ ROW * m_NumberCols + COL];
      }

      //! @brief return reference to const  element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return const element defined bey ROW and COL number
      const_reference operator()( const size_t ROW, const size_t COL) const
      {
        BCL_Assert
        (
          ROW < m_NumberRows && COL < m_NumberCols,
          "requested (" + util::Format()( ROW) + "," + util::Format()( COL)
          + ") but matrix has only " + util::Format()( m_NumberRows) + " rows and " + util::Format()( m_NumberCols) + " cols"
        );
        return m_Data[ ROW * m_NumberCols + COL];
      }

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a pointer to the first member of that row
      iterator operator[]( const size_t ROW)
      {
        BCL_Assert
        (
          ROW < m_NumberRows,
          "requested [" + util::Format()( ROW) + "] but matrix has only " + util::Format()( m_NumberRows) + " rows"
        );
        return m_Data + ROW * m_NumberCols;
      }

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a constant pointer to the first member of that row
      const_iterator operator[]( const size_t ROW) const
      {
        BCL_Assert
        (
          ROW < m_NumberRows,
          "requested [" + util::Format()( ROW) + "] but matrix has only " + util::Format()( m_NumberRows) + " rows"
        );
        return m_Data + ROW * m_NumberCols;
      }

      //! @brief assignment from Matrix
      //! @param MATRIX the matrix used as source
      //! @return reference to this Matrix
      Matrix< t_DataType> &operator =( const Matrix< t_DataType> &MATRIX);

      //! @brief move from Matrix
      //! @param MATRIX the matrix used as source
      //! @return reference to this Matrix
      Matrix< t_DataType> &operator =( Matrix< t_DataType> && MATRIX);

      //! @brief assignment from MatrixInterface
      //! @param MATRIX_INTERFACE the matrix used as source
      //! @return reference to this Matrix
      Matrix< t_DataType> &operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE);

      //! @brief assignment from value
      //! @param VALUE all elements are set to that value
      //! @return reference to this assigned Matrix
      Matrix< t_DataType> &operator =( const t_DataType &VALUE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // template class Matrix

    //! @brief create an identity matrix
    //! @param SIZE size of the square identity matrix
    template< typename t_DataType>
    Matrix< t_DataType> IdentityMatrix( const size_t &SIZE)
    {
      Matrix< t_DataType> matrix( SIZE, SIZE);
      for( size_t i( 0); i < SIZE; ++i)
      {
        matrix( i, i) = t_DataType( 1);
      }
      return matrix;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_H_
