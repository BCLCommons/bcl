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

#ifndef BCL_LINAL_MATRIX2X2_H_
#define BCL_LINAL_MATRIX2X2_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_interface.h"
#include "bcl_linal_vector.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Matrix2x2
    //! @brief is an implementation of a matrix with 3x3 elements
    //! @details The matrix is stored as 4 elements of t_DataType of the form:
    //!          | m00 m01 |   | a b |\n
    //!          | m10 m11 | = | c d |\n
    //!          Which makes explicit implementations of formulas easier.
    //!
    //! @see @link example_linal_matrix2x2.cpp @endlink
    //! @author mendenjl
    //! @date Jun 12, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Matrix2x2 :
      public MatrixInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      static const size_t s_Dimension      = 2;                         //!< dimension of this matrix
      static const size_t s_NumberElements = s_Dimension * s_Dimension; //!< number of elements
      t_DataType m_00a;
      t_DataType m_01b;
      t_DataType m_10c;
      t_DataType m_11d;

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
      Matrix2x2();

      //! @brief construct from filler
      //! @param FILL_VALUE assign every element to that value
      explicit Matrix2x2( const t_DataType &FILL_VALUE);

      //! @brief construct from pointer to data
      //! @param DATA pointer to field of data
      explicit Matrix2x2( const t_DataType *DATA);

      //! @brief copy constructor
      Matrix2x2( const Matrix2x2< t_DataType> &MATRIX);

      //! @brief copy constructor from MatrixInterface
      //! @param MATRIX_INTERFACE matrix interface to be copied from
      Matrix2x2( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE);

      //! @brief Clone function
      //! @return pointer to new Matrix< t_DataType>
      Matrix2x2< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get number of rows
      //! @return number of rows
      size_t GetNumberRows() const;

      //! @brief get number of columns
      //! @return number of columns
      size_t GetNumberCols() const;

      //! @brief number of elements
      //! @return total number of elements in matrix
      size_t GetNumberOfElements() const;

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Matrix
      const t_DataType *Begin() const;

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Matrix
      t_DataType *Begin();

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      const t_DataType *End() const;

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Matrix
      t_DataType *End();

    ////////////////
    // operations //
    ////////////////

      //! @brief is matrix a square matrix
      //! @return true if number of cols and rows are identical
      bool IsSquare() const;

      //! @brief is matrix a diagonal matrix
      //! @return true if all but the elements in the diagonal are 0
      bool IsDiagonal() const;

      //! @brief checks if matrix is defined
      //! @return bool - true if matrix is defined
      bool IsDefined() const;

      //! check whether matrix is symmetric
      bool IsSymmetric() const;

      //! @brief is matrix a tridiagonal matrix
      //! @return if diagonal and adjacent diagonals are filled and the rest is 0
      bool IsTriDiagonal() const;

      //! @brief swap elements between two rows
      //! @param ROW_A the first row
      //! @param ROW_B the second row
      void SwapRows( const size_t ROW_A, const size_t ROW_B);

      //! @brief replace row
      //! @param ROW the row to replace
      //! @param VECTOR the vector to replace the current row with
      void ReplaceRow( const size_t ROW, const VectorConstInterface< t_DataType> &VECTOR);

      //! @brief sort rows and given vectors (less than)
      //! @param VECTOR vector with 3 elements
      void SortRowsAndVector( VectorInterface< t_DataType> &VECTOR);

      //! @brief determinant of this matrix
      //! @return the determinant of the matrix
      t_DataType Determinant() const;

      //! @brief trace of matrix
      //! @return the sum of the diagonal elements
      t_DataType Trace() const;

      //! @brief transpose the matrix
      void Transpose();

      //! @brief eigenvalues
      //! @return the eigenvalues
      Vector< t_DataType> EigenValues() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return reference to changeable element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return changeable reference to the element defined bey ROW and COL number
      t_DataType &operator()( const size_t ROW, const size_t COL);

      //! @brief return reference to const  element ( ROW, COL)
      //! @param ROW the row number, starting with 0
      //! @param COL the col number, starting with 0
      //! @return const element defined bey ROW and COL number
      const t_DataType &operator()( const size_t ROW, const size_t COL) const;

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a pointer to the first member of that row
      t_DataType *operator[]( const size_t ROW);

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a constant pointer to the first member of that row
      const t_DataType *operator[]( const size_t ROW) const;

      //! @brief assignment from Matrix2x2
      //! @param MATRIX the matrix used as source
      //! @return reference to this Matrix
      Matrix2x2< t_DataType> &operator =( const Matrix2x2< t_DataType> &MATRIX);

      //! @brief assignment from MatrixInterface
      //! @param MATRIX_INTERFACE the matrix used as source
      //! @return reference to this Matrix
      Matrix2x2< t_DataType> &operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE);

      //! @brief assignment from value
      //! @param VALUE all elements are set to that value
      //! @return reference to this assigned Matrix
      Matrix2x2< t_DataType> &operator =( const t_DataType &VALUE);

      //! @brief multiply with a second matrix3x3
      //! @param MATRIX other Matrix2x2
      //! @return reference to this matrix after multiplication
      Matrix2x2< t_DataType> &operator *=( const Matrix2x2< t_DataType> &MATRIX);

      //! @brief multiply with a second matrix3x3
      //! @param MATRIX other Matrix2x2
      //! @return the result
      Matrix2x2< t_DataType> operator *( const Matrix2x2< t_DataType> &MATRIX) const;

      //! @brief subtract second matrix
      //! @param MATRIX the matrix to subtract
      //! @return reference to this matrix
      Matrix2x2< t_DataType> &operator -=( const Matrix2x2< t_DataType> &MATRIX);

      //! @brief subtract second matrix
      //! @param MATRIX the matrix to subtract
      //! @return result matrix
      Matrix2x2< t_DataType> operator -( const Matrix2x2< t_DataType> &MATRIX) const;

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

    public:

    //////////////////////
    // helper functions //
    //////////////////////

    }; // template class Matrix

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix2x2< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Matrix2x2< double>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX2X2_H_
