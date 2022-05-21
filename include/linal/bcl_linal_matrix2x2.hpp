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

#ifndef BCL_LINAL_MATRIX2X2_HPP_
#define BCL_LINAL_MATRIX2X2_HPP_

// include the header of this class
#include "bcl_linal_matrix2x2.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Matrix2x2< t_DataType>::Matrix2x2()
    {
      std::fill_n( &m_00a, s_NumberElements, t_DataType( 0));
    }

    //! @brief construct from filler
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Matrix2x2< t_DataType>::Matrix2x2( const t_DataType &FILL_VALUE)
    {
      // set all values to FILL_VALUE
      std::fill_n( &m_00a, s_NumberElements, FILL_VALUE);
    }

    //! @brief construct from pointer to data
    //! @param DATA pointer to field of data
    template< typename t_DataType>
    Matrix2x2< t_DataType>::Matrix2x2( const t_DataType *DATA)
    {
      // copy data
      std::copy( DATA, DATA + s_NumberElements, &m_00a);
    }

    //! @brief copy constructor
    template< typename t_DataType>
    Matrix2x2< t_DataType>::Matrix2x2( const Matrix2x2< t_DataType> &MATRIX)
    {
      std::copy( &MATRIX.m_00a, &MATRIX.m_00a + s_NumberElements, &m_00a);
    }

    //! @brief copy constructor from MatrixInterface
    //! @param MATRIX_INTERFACE matrix interface to be copied from
    template< typename t_DataType>
    Matrix2x2< t_DataType>::Matrix2x2( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
    {
      // check if dimensions match
      BCL_Assert
      (
        MATRIX_INTERFACE.GetNumberRows() == s_Dimension || MATRIX_INTERFACE.GetNumberCols() == s_Dimension,
        "can only copy from 2x2 matrix!"
      );

      std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), &m_00a);
    }

    //! @brief Clone function
    //! @return pointer to new Matrix< t_DataType>
    template< typename t_DataType>
    Matrix2x2< t_DataType> *Matrix2x2< t_DataType>::Clone() const
    {
      return new Matrix2x2< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Matrix2x2< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get number of rows
    //! @return number of rows
    template< typename t_DataType>
    size_t Matrix2x2< t_DataType>::GetNumberRows() const
    {
      return s_Dimension;
    }

    //! @brief get number of columns
    //! @return number of columns
    template< typename t_DataType>
    size_t Matrix2x2< t_DataType>::GetNumberCols() const
    {
      return s_Dimension;
    }

    //! @brief number of elements
    //! @return total number of elements in matrix
    template< typename t_DataType>
    size_t Matrix2x2< t_DataType>::GetNumberOfElements() const
    {
      return s_NumberElements;
    }

    //! @brief pointer to First Element
    //! @return const pointer to first element in range containing all elements of Matrix
    template< typename t_DataType>
    const t_DataType *Matrix2x2< t_DataType>::Begin() const
    {
      return &m_00a;
    }

    //! @brief pointer to First Element
    //! @return pointer to first element in range containing all elements of Matrix
    template< typename t_DataType>
    t_DataType *Matrix2x2< t_DataType>::Begin()
    {
      return &m_00a;
    }

    //! @brief pointer to end of range
    //! @return const pointer to address one after last element in Matrix
    template< typename t_DataType>
    const t_DataType *Matrix2x2< t_DataType>::End() const
    {
      return &m_00a + s_NumberElements;
    }

    //! @brief pointer to end of range
    //! @return pointer to address one after last element in Matrix
    template< typename t_DataType>
    t_DataType *Matrix2x2< t_DataType>::End()
    {
      return &m_00a + s_NumberElements;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is matrix a square matrix
    //! @return true if number of cols and rows are identical
    template< typename t_DataType>
    bool Matrix2x2< t_DataType>::IsSquare() const
    {
      return true;
    }

    //! @brief is matrix a diagonal matrix
    //! @return true if all but the elements in the diagonal are 0
    template< typename t_DataType>
    bool Matrix2x2< t_DataType>::IsDiagonal() const
    {
      // return true if all elements but the diagonal is 0
      const t_DataType zero( 0);
      return math::EqualWithinMachineTolerance( m_01b, zero) && math::EqualWithinMachineTolerance( m_10c, zero);
    }

    //! @brief checks if matrix is defined
    //! @return bool - true if matrix is defined
    template< typename t_DataType>
    bool Matrix2x2< t_DataType>::IsDefined() const
    {
      return math::Statistics::IsDefined( &m_00a, &m_00a + s_NumberElements);
    }

    //! check whether matrix is symmetric
    template< typename t_DataType>
    bool Matrix2x2< t_DataType>::IsSymmetric() const
    {
      return math::EqualWithinMachineTolerance( m_01b, m_10c);
    }

    //! @brief is matrix a tridiagonal matrix
    //! @return if diagonal and adjacent diagonals are filled and the rest is 0
    template< typename t_DataType>
    bool Matrix2x2< t_DataType>::IsTriDiagonal() const
    {
      // return true if all elements but the inner three diagonals are 0
      return true;
    }

    //! @brief swap elements between two rows
    //! @param ROW_A the first row
    //! @param ROW_B the second row
    template< typename t_DataType>
    void Matrix2x2< t_DataType>::SwapRows( const size_t ROW_A, const size_t ROW_B)
    {
      if( ROW_A != ROW_B && ROW_A < 2 && ROW_B < 2)
      {
        std::swap( m_00a, m_10c);
        std::swap( m_01b, m_11d);
      }
    }

    //! @brief replace row
    //! @param ROW the row to replace
    //! @param VECTOR the vector to replace the current row with
    template< typename t_DataType>
    void Matrix2x2< t_DataType>::ReplaceRow( const size_t ROW, const VectorConstInterface< t_DataType> &VECTOR)
    {
      if( ROW >= s_Dimension || VECTOR.GetSize() != s_Dimension)
      {
        return;
      }

      std::copy( VECTOR.Begin(), VECTOR.End(), operator[]( ROW));
    }

    //! @brief sort rows and given vectors (less than)
    //! @param VECTOR vector with 3 elements
    template< typename t_DataType>
    void Matrix2x2< t_DataType>::SortRowsAndVector( VectorInterface< t_DataType> &VECTOR)
    {
      // sort
      if( VECTOR( 0) < VECTOR( 1))
      {
        std::swap( VECTOR( 0), VECTOR( 1));
        SwapRows( 0, 1);
      }
    }

    //! @brief determinant of this matrix
    //! @return the determinant of the matrix
    template< typename t_DataType>
    t_DataType Matrix2x2< t_DataType>::Determinant() const
    {
      return m_00a * m_11d - m_01b * m_10c;
    }

    //! @brief trace of matrix
    //! @return the sum of the diagonal elements
    template< typename t_DataType>
    t_DataType Matrix2x2< t_DataType>::Trace() const
    {
      return m_00a + m_11d;
    }

    //! @brief transpose the matrix
    template< typename t_DataType>
    void Matrix2x2< t_DataType>::Transpose()
    {
      std::swap( m_01b, m_10c);
    }

    //! @brief eigenvalues
    //! @return the eigenvalues
    template< typename t_DataType>
    Vector< t_DataType> Matrix2x2< t_DataType>::EigenValues() const
    {
      Vector< t_DataType> eigenvalues( 2, t_DataType( 0));
      // quadratic formula equation
      const t_DataType trace( Trace());
      const t_DataType matrix_det( Determinant());

      // quadratic equation determinant
      const t_DataType det( math::Sqr( trace) - t_DataType( 4) * matrix_det);
      // detect complex roots
      BCL_Assert( det >= 0, "Complex eigenvalues detected!");
      const t_DataType rhs( math::Sqrt( det));
      eigenvalues( 0) = 0.5 * ( trace + rhs);
      eigenvalues( 1) = 0.5 * ( trace - rhs);

      // end
      return eigenvalues;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return reference to changeable element ( ROW, COL)
    //! @param ROW the row number, starting with 0
    //! @param COL the col number, starting with 0
    //! @return changeable reference to the element defined bey ROW and COL number
    template< typename t_DataType>
    t_DataType &Matrix2x2< t_DataType>::operator()( const size_t ROW, const size_t COL)
    {
      return ( ROW ? ( COL ? m_11d : m_10c) : ( COL ? m_01b : m_00a));
    }

    //! @brief return reference to const  element ( ROW, COL)
    //! @param ROW the row number, starting with 0
    //! @param COL the col number, starting with 0
    //! @return const element defined bey ROW and COL number
    template< typename t_DataType>
    const t_DataType &Matrix2x2< t_DataType>::operator()( const size_t ROW, const size_t COL) const
    {
      return ( ROW ? ( COL ? m_11d : m_10c) : ( COL ? m_01b : m_00a));
    }

    //! @brief access to a particular row
    //! @param ROW the row to access
    //! @return a pointer to the first member of that row
    template< typename t_DataType>
    t_DataType *Matrix2x2< t_DataType>::operator[]( const size_t ROW)
    {
      return ( ROW ? &m_10c : &m_00a);
    }

    //! @brief access to a particular row
    //! @param ROW the row to access
    //! @return a constant pointer to the first member of that row
    template< typename t_DataType>
    const t_DataType *Matrix2x2< t_DataType>::operator[]( const size_t ROW) const
    {
      return ( ROW ? &m_10c : &m_00a);
    }

    //! @brief assignment from Matrix2x2
    //! @param MATRIX the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix2x2< t_DataType> &Matrix2x2< t_DataType>::operator =( const Matrix2x2< t_DataType> &MATRIX)
    {
      // copy all elements
      if( &m_00a != &MATRIX.m_00a)
      {
        // copy elements
        std::copy( &MATRIX.m_00a, &MATRIX.m_00a + s_NumberElements, &m_00a);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from MatrixInterface
    //! @param MATRIX_INTERFACE the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix2x2< t_DataType> &Matrix2x2< t_DataType>::operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
    {
      // copy all elements
      if( &m_00a != MATRIX_INTERFACE.Begin())
      {
        // check if dimensions match
        BCL_Assert
        (
          MATRIX_INTERFACE.GetNumberRows() == s_Dimension || MATRIX_INTERFACE.GetNumberCols() == s_Dimension,
          "can only assign from 2x2 matrix!"
        );

        std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), &m_00a);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from value
    //! @param VALUE all elements are set to that value
    //! @return reference to this assigned Matrix
    template< typename t_DataType>
    Matrix2x2< t_DataType> &Matrix2x2< t_DataType>::operator =( const t_DataType &VALUE)
    {
      // set all element to given VALUE
      std::fill( &m_00a, &m_00a + s_NumberElements, VALUE);

      // return reference to this Vector
      return *this;
    }

    //! @brief multiply with a second matrix3x3
    //! @param MATRIX other Matrix2x2
    //! @return reference to this matrix after multiplication
    template< typename t_DataType>
    Matrix2x2< t_DataType> &Matrix2x2< t_DataType>::operator *=( const Matrix2x2< t_DataType> &MATRIX)
    {
      Matrix2x2< t_DataType> copy( *this);

      m_00a = copy.m_00a * MATRIX.m_00a + copy.m_01b * MATRIX.m_10c;
      m_01b = copy.m_00a * MATRIX.m_01b + copy.m_01b * MATRIX.m_11d;
      m_10c = copy.m_10c * MATRIX.m_00a + copy.m_11d * MATRIX.m_10c;
      m_11d = copy.m_10c * MATRIX.m_01b + copy.m_11d * MATRIX.m_11d;

      return *this;
    }

    //! @brief multiply with a second matrix3x3
    //! @param MATRIX other Matrix2x2
    //! @return the result
    template< typename t_DataType>
    Matrix2x2< t_DataType> Matrix2x2< t_DataType>::operator *( const Matrix2x2< t_DataType> &MATRIX) const
    {
      Matrix2x2< t_DataType> copy( *this);

      copy.m_00a = m_00a * MATRIX.m_00a + m_01b * MATRIX.m_10c;
      copy.m_01b = m_00a * MATRIX.m_01b + m_01b * MATRIX.m_11d;
      copy.m_10c = m_10c * MATRIX.m_00a + m_11d * MATRIX.m_10c;
      copy.m_11d = m_10c * MATRIX.m_01b + m_11d * MATRIX.m_11d;

      return copy;
    }

    //! @brief subtract second matrix
    //! @param MATRIX the matrix to subtract
    //! @return reference to this matrix
    template< typename t_DataType>
    Matrix2x2< t_DataType> &Matrix2x2< t_DataType>::operator -=( const Matrix2x2< t_DataType> &MATRIX)
    {
      std::transform( &m_00a, &m_00a + s_NumberElements, &MATRIX.m_00a, &m_00a, std::minus< t_DataType>());
      return *this;
    }

    //! @brief subtract second matrix
    //! @param MATRIX the matrix to subtract
    //! @return result matrix
    template< typename t_DataType>
    Matrix2x2< t_DataType> Matrix2x2< t_DataType>::operator -( const Matrix2x2< t_DataType> &MATRIX) const
    {
      Matrix2x2< t_DataType> copy( *this);
      std::transform( &copy.m_00a, &copy.m_00a + s_NumberElements, &MATRIX.m_00a, &copy.m_00a, std::minus< t_DataType>());
      return copy;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Matrix2x2< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read each element from ISTREAM
      for( t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
      {
        BCL_Assert( io::Serialize::Read( *ptr, ISTREAM), "unable to read element");
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Matrix2x2< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write date
      for( size_t i = 0; i < s_Dimension;)
      {
        io::Serialize::InsertIndent( OSTREAM, INDENT);
        for( size_t j = 0; j < s_Dimension; ++j)
        {
          io::Serialize::Write( operator()( i, j), OSTREAM) << '\t';
        }
        ++i;
        // print line break except after last line
        if( i != s_Dimension)
        {
          OSTREAM << '\n';
        }
      }

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! single instance of the Matrix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Matrix2x2< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Matrix2x2< t_DataType>())
    );

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX2X2_HPP_
