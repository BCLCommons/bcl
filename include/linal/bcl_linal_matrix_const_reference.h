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

#ifndef BCL_LINAL_MATRIX_CONST_REFERENCE_H_
#define BCL_LINAL_MATRIX_CONST_REFERENCE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_const_interface.h"
#include "bcl_linal_matrix_interface.h"
#include "bcl_linal_vector_const_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MatrixConstReference
    //! @brief This is a general implementation of a const non data-owning Matrix
    //! @details The matrix has a pointer to a continuous array of t_DataType that it doesn't own in row-major format
    //!
    //! @tparam t_DataType double, float, int, complex, etc
    //!
    //! @author woetzen, loweew
    //! @see @link example_linal_matrix_const_reference.cpp @endlink
    //! @date Aug 18, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MatrixConstReference :
      public MatrixConstInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      size_t m_NumberRows; //!< number of rows in matrix
      size_t m_NumberCols; //!< number of cols in matrix
      const t_DataType *m_Data; //!< array containing all data

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MatrixConstReference() :
        m_NumberRows( 0),
        m_NumberCols( 0),
        m_Data( NULL)
      {
      }

      //! @brief construct from dimension and pointer to data
      //! @param NUMBER_ROWS number of rows in matrix
      //! @param NUMBER_COLS number of cols in matrix
      //! @param DATA pointer to field of data
      MatrixConstReference
      (
        const size_t NUMBER_ROWS,
        const size_t NUMBER_COLS,
        const t_DataType *DATA
      ) :
        m_NumberRows( NUMBER_ROWS),
        m_NumberCols( NUMBER_COLS),
        m_Data( DATA)
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data != NULL || NUMBER_ROWS == 0 || NUMBER_COLS == 0,
          "cannot reference NULL data elements of type: " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief copy constructor from MatrixInterface
      //! @param MATRIX_INTERFACE matrix interface to be copied from
      explicit MatrixConstReference( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE) :
        m_NumberRows( MATRIX_INTERFACE.GetNumberRows()),
        m_NumberCols( MATRIX_INTERFACE.GetNumberCols()),
        m_Data( MATRIX_INTERFACE.Begin())
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data || m_NumberRows == 0 || m_NumberCols == 0,
          "cannot reference NULL data elements of type: " + util::Format()( m_NumberRows * m_NumberCols) +
          " elements of type: " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief Clone function
      //! @return pointer to new MatrixInterface< t_DataType>
      MatrixConstReference< t_DataType> *Clone() const
      {
        return new MatrixConstReference< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

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
      const t_DataType *Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      const t_DataType *End() const
      {
        return m_Data + m_NumberRows * m_NumberCols;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief return reference to const  element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return const element defined bey ROW and COL number
      const t_DataType &operator()( const size_t ROW, const size_t COL) const
      {
        return m_Data[ ROW * m_NumberCols + COL];
      }

      const t_DataType *operator[]( const size_t ROW) const
      {
        return Begin() + ROW * m_NumberCols;
      }

      //! @brief make this reference point to the given matrix
      //! @param MATRIX matrix to reference
      MatrixConstReference< t_DataType> &Reference( const MatrixConstInterface< t_DataType> &MATRIX)
      {
        m_NumberCols = MATRIX.GetNumberCols();
        m_NumberRows = MATRIX.GetNumberRows();
        m_Data = MATRIX.Begin();
        return *this;
      }

      //! @brief assignment from MatrixInterface
      //! @param MATRIX_INTERFACE the matrix used as source
      //! @return reference to this Matrix
//      MatrixConstReference< t_DataType> &operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
//      {
//        // set members
//        m_NumberRows = MATRIX_INTERFACE.GetNumberRows();
//        m_NumberCols = MATRIX_INTERFACE.GetNumberCols();
//        m_Data       = MATRIX_INTERFACE.Begin();
//
//        // check if allocation was successful
//        BCL_Assert
//        (
//          m_Data || m_NumberRows == 0 || m_NumberCols == 0,
//          "cannot reference NULL data elements of type: " + util::Format()( m_NumberRows * m_NumberCols) +
//          " elements of type: " + GetStaticClassName< t_DataType>()
//        );
//
//        return *this;
//      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "You can't read in a matrix reference because you don't have ownership!", -1);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
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

    }; // template class MatrixConstReference

    //! single instance of the Marix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> MatrixConstReference< t_DataType>::s_Instance;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX_CONST_REFERENCE_H_
