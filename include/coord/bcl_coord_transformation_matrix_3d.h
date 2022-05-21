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

#ifndef BCL_COORD_TRANSFORMATION_MATRIX_3D_H_
#define BCL_COORD_TRANSFORMATION_MATRIX_3D_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_const_interface.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransformationMatrix3D
    //! @brief This class is for 3-dimensional TransformationMatrix3D
    //!
    //! @see @link example_coord_transformation_matrix_3d.cpp @endlink
    //! @author meilerj
    //! @date 21.08.2004
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class TransformationMatrix3D :
      public linal::MatrixConstInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      static const size_t s_NumberRows = 4;
      static const size_t s_NumberCols = 4;
      static const size_t s_NumberElements = s_NumberRows * s_NumberCols;

      t_DataType m_Data[ s_NumberElements];

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct TransformationMatrix3D and set to unitmatrix (empty constructor)
      TransformationMatrix3D()
      {
        SetUnit();
      }

      //! copy constructor
      TransformationMatrix3D( const TransformationMatrix3D &MATRIX)
      {
        std::copy( MATRIX.m_Data, MATRIX.m_Data + s_NumberElements, m_Data);
      }

      //! @brief constructor from a definition state
      //! @param DEFINITION_STATE defined or undefined
      TransformationMatrix3D( const util::UndefinedObject DEFINITION_STATE)
      {
        std::fill( m_Data, m_Data + s_NumberElements, util::GetUndefined< t_DataType>());
      }

      //! copy constructor
      TransformationMatrix3D *Clone() const
      {
        return new TransformationMatrix3D( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName< TransformationMatrix3D>();
      }

      //! @brief get number of rows
      //! @return number of rows
      size_t GetNumberRows() const
      {
        return s_NumberRows;
      }

      //! @brief get number of columns
      //! @return number of columns
      size_t GetNumberCols() const
      {
        return s_NumberCols;
      }

      //! @brief number of elements
      //! @return total number of elements in matrix
      size_t GetNumberOfElements() const
      {
        return s_NumberElements;
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
        return m_Data + s_NumberElements;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief is matrix a square matrix
      //! @return true if number of cols and rows are identical
      bool IsSquare() const
      {
        return true;
      }

      //! @brief is matrix a diagonal matrix
      //! @return true if all but the elements in the diagonal are 0
      bool IsDiagonal() const
      {
        return false;
      }

      //! @brief is matrix a tridiagonal matrix
      //! @return if diagonal and adjacent diagonals are filled and the rest is 0
      bool IsTriDiagonal() const
      {
        return false;
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
        return m_Data[ ROW * s_NumberCols + COL];
      }

      //! @brief access to a particular row
      //! @param ROW the row to access
      //! @return a constant pointer to the first member of that row
      const t_DataType *operator[]( const size_t ROW) const
      {
        return m_Data + ROW * s_NumberCols;
      }

      //! operator = TransformationsMatrix3D
      TransformationMatrix3D &operator =( const TransformationMatrix3D &MATRIX)
      {
        std::copy( MATRIX.m_Data, MATRIX.m_Data + s_NumberElements, m_Data);

        // end
        return *this;
      }

    ////////////////
    // operations //
    ////////////////

      //! Set Unit Matrix
      TransformationMatrix3D &SetUnit()
      {
        // all elements to 0
        std::fill( m_Data, m_Data + s_NumberElements, t_DataType( 0));

        // diagonal to 1
        for( t_DataType *ptr( m_Data), *ptr_end( m_Data + s_NumberElements); ptr < ptr_end; ptr += s_NumberCols + 1)
        {
          *ptr = t_DataType( 1);
        }

        // end
        return *this;
      }

      //! @brief returns whether this transformation is defined or not by looking at the scaling factor
      //! @return whether this transformation is defined or not
      bool IsDefined() const
      {
        return math::Statistics::IsDefined( m_Data, m_Data + s_NumberElements);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write TransformationMatrix3D to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        for( const t_DataType *row( m_Data), *row_end( m_Data + s_NumberElements); row != row_end; row += s_NumberCols)
        {
          io::Serialize::InsertIndent( OSTREAM, INDENT);
          std::copy( row, row + s_NumberCols, std::ostream_iterator< t_DataType>( OSTREAM, "\t"));
          OSTREAM << '\n';
        }

        // end
        return OSTREAM;
      }

      //! read TransformationMatrix3D from io::IFStream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read
        for( t_DataType *ptr( m_Data), *ptr_end( m_Data + s_NumberElements); ptr != ptr_end; ++ptr)
        {
          io::Serialize::Read( *ptr, ISTREAM);
        }

        // end
        return ISTREAM;
      }

    }; // class TransformationMatrix3D

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_TRANSFORMATION_MATRIX_3D_H_
