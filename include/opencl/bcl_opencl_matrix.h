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

#ifndef BCL_OPENCL_MATRIX_H_
#define BCL_OPENCL_MATRIX_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Matrix
    //! @brief opencl matrix class to ease creation of buffer objects from linal matrices while providing access to
    //!        dimensionality and padding information
    //!
    //! @tparam t_DataType can be float, double, complex, int, etc...
    //!
    //! @see @link example_opencl_matrix.cpp @endlink
    //! @author loweew
    //! @date Mar 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Matrix :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! number of rows
      size_t       m_NumberRows;

      //! number of cols
      size_t       m_NumberCols;

      //! row padding
      size_t       m_RowPadding;

      //! col padding
      size_t       m_ColPadding;

      //! command queue
      CommandQueue m_Queue;

      //! opencl buffer
      mutable Buffer       m_Data; // mutable otherwise SubMatrix does not work, since the createSubBuffer does not make a copy of the data

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Matrix() :
        m_NumberRows(),
        m_NumberCols(),
        m_RowPadding(),
        m_ColPadding(),
        m_Queue(),
        m_Data()
      {
      }

      //! @brief constructor - set all elements to 0
      //! @param ROWS, COLS, dimensions of matrix
      //! @param QUEUE command queue
      //! @param ROW_PADDING, COL_PADDING padding of matrix
      Matrix
      (
        const size_t &ROWS, const size_t &COLS, const CommandQueue &QUEUE,
        const size_t &ROW_PADDING = 0, const size_t &COL_PADDING = 0

      ) :
        m_NumberRows( ROWS + ROW_PADDING),
        m_NumberCols( COLS + COL_PADDING),
        m_RowPadding( ROW_PADDING),
        m_ColPadding( COL_PADDING),
        m_Queue( QUEUE),
        m_Data( cl::Buffer( m_Queue.GetContext(), CL_MEM_READ_WRITE, m_NumberRows * m_NumberCols * sizeof( t_DataType), NULL, NULL))
      {
        Fill( t_DataType( 0), 0, m_NumberRows, 0, m_NumberCols);
      }

      //! @brief constructor
      //! @param ROWS, COLS, dimensions of matrix
      //! @param QUEUE command queue
      //! @param ROW_PADDING, COL_PADDING padding of matrix
      //! @param ELEMENT element to initialize buffer with
      Matrix
      (
        const size_t &ROWS, const size_t &COLS, const CommandQueue &QUEUE,
        const size_t &ROW_PADDING, const size_t &COL_PADDING,
        const t_DataType &ELEMENT
      ) :
        m_NumberRows( ROWS + ROW_PADDING),
        m_NumberCols( COLS + COL_PADDING),
        m_RowPadding( ROW_PADDING),
        m_ColPadding( COL_PADDING),
        m_Queue( QUEUE),
        m_Data( cl::Buffer( m_Queue.GetContext(), CL_MEM_READ_WRITE, m_NumberRows * m_NumberCols * sizeof( t_DataType), NULL, NULL))
      {
        Fill( t_DataType( 0), 0, m_NumberRows, 0, m_NumberCols);
        Fill( ELEMENT, 0, m_NumberRows - m_RowPadding, 0, m_NumberCols - m_ColPadding);
      }

      //! @brief constructor from linal matrix
      //! @param MATRIX matrix to create buffer from
      //! @param QUEUE command queue
      //! @param ROW_PADDING, COL_PADDING the amount of padding you want to add to the matrix
      Matrix
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const CommandQueue &QUEUE,
        const size_t &ROW_PADDING = 0,
        const size_t &COL_PADDING = 0
      ) :
        m_NumberRows( MATRIX.GetNumberRows() + ROW_PADDING),
        m_NumberCols( MATRIX.GetNumberCols() + COL_PADDING),
        m_RowPadding( ROW_PADDING),
        m_ColPadding( COL_PADDING),
        m_Queue( QUEUE),
        m_Data
        (
          Buffer::CreateBufferFromMatrix< t_DataType>
          (
            linal::Matrix< t_DataType>( MATRIX, ROW_PADDING, COL_PADDING), QUEUE
          )
        )
      {
      }

      //! @brief constructor from pointer
      //! @param ROWS, COLS dimensions
      //! @param DATA pointer to data
      //! @param QUEUE command queue
      Matrix( const size_t &ROWS, const size_t &COLS, const t_DataType *DATA, const CommandQueue &QUEUE) :
        m_NumberRows( ROWS),
        m_NumberCols( COLS),
        m_RowPadding( 0),
        m_ColPadding( 0),
        m_Queue( QUEUE)
      {
        m_Queue.enqueueWriteBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * ROWS * COLS, DATA);
      }

      //! @brief copy constructor
      //! @param DATA matrix to copy
      Matrix( const Matrix< t_DataType> &DATA) :
        m_NumberRows( DATA.m_NumberRows),
        m_NumberCols( DATA.m_NumberCols),
        m_RowPadding( DATA.m_RowPadding),
        m_ColPadding( DATA.m_ColPadding),
        m_Queue( DATA.m_Queue),
        m_Data( DATA.m_Data)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Matrix
      Matrix< t_DataType> *Clone() const
      {
        return new Matrix< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get number rows
      //! @return number rows
      size_t GetNumberRows() const
      {
        return m_NumberRows;
      }

      //! @brief get number cols
      //! @return number cols
      size_t GetNumberCols() const
      {
        return m_NumberCols;
      }

      //! @brief get row padding
      //! @return row padding
      size_t GetRowPadding() const
      {
        return m_RowPadding;
      }

      //! @brief get col padding
      //! @return col padding
      size_t GetColPadding() const
      {
        return m_ColPadding;
      }

      //! @brief get number of elements
      //! @return number of elements
      size_t GetNumberOfElements() const
      {
        return m_NumberCols * m_NumberRows;
      }

      //! @brief gets buffer object
      //! @return buffer
      const Buffer &GetData() const
      {
        return m_Data;
      }

      //! @brief gets command queue associated with this buffer
      //! @return queue
      const CommandQueue &GetQueue() const
      {
        return m_Queue;
      }

      //! @brief returns linal::Matrix with user-defined padding removed
      //! @param ROW_PAD, COL_PAD the padding to remove from each dimension
      //! @return the linal::Matrix
      linal::Matrix< t_DataType> GetHostMatrix( const size_t &ROW_PAD = 0, const size_t &COL_PAD = 0) const
      {
        linal::Matrix< t_DataType> tmp( m_NumberRows, m_NumberCols);
        cl_int error_number = m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * m_NumberRows * m_NumberCols, tmp.Begin());
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));
        return tmp.CreateSubMatrix( m_NumberRows - ROW_PAD, m_NumberCols - COL_PAD);
      }

      //! @brief returns linal::Matrix with user-defined padding removed
      //! @param ROW_PAD, COL_PAD the padding to remove from each dimension
      //! @return the linal::Matrix
      linal::Matrix< t_DataType> GetHostMatrixUnpadded() const
      {
        return GetHostMatrix( m_RowPadding, m_ColPadding);
      }

      //! @brief prints the matrix to the screen
      void Print( std::string &MESSAGE)
      {
        BCL_MessageStd( MESSAGE + util::Format()( this->GetHostMatrix( m_RowPadding, m_ColPadding)));
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief sets position to value
      //! @param ROW, COL the element position to be set
      //! @param VALUE the new value
      void SetValue( const size_t ROW, const size_t COL, const t_DataType &VALUE)
      {
        m_Queue.enqueueWriteBuffer( m_Data, CL_TRUE, ( ROW * m_NumberCols + COL) * sizeof( t_DataType), sizeof( t_DataType), &VALUE);
      }

      //! @brief Fill matrix with values
      //! @param VALUE the value to set
      //! @param FIRST_ROW first row to set
      //! @param NR_ROWS number of elements per col
      //! @param FIRST_COL first column to set
      //! @param NR_COLS number of elements per row
      void Fill
      (
        const t_DataType &VALUE,
        const size_t FIRST_ROW, const size_t NR_ROWS,
        const size_t FIRST_COL, const size_t NR_COLS
      )
      {
        cl_int error_number( CL_SUCCESS);              // Error code var

        if( GetNumberOfElements() == 0)
        {
          return;
        }
        BCL_Assert
        (
          FIRST_ROW + NR_ROWS <= m_NumberRows && FIRST_COL + NR_COLS <= m_NumberCols,
          "matrix too small!"
        );
        const size_t block_size( 16);

        // Create the kernel
        cl::Kernel kernel( GetTools().GetBufferProgram( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue), "FillSubmatrix", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const cl::NDRange block_dimensions( block_size, block_size);
        const cl::NDRange offset;
        const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, NR_ROWS), Tools::RoundUp( block_size, NR_COLS));

        error_number  = kernel.setArg( 0, m_Data);
        error_number |= kernel.setArg( 1, cl_uint( m_NumberCols));
        error_number |= kernel.setArg( 2, cl_uint( FIRST_ROW));
        error_number |= kernel.setArg( 3, cl_uint( NR_ROWS));
        error_number |= kernel.setArg( 4, cl_uint( FIRST_COL));
        error_number |= kernel.setArg( 5, cl_uint( NR_COLS));
        error_number |= kernel.setArg( 6, VALUE);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      //! @brief set the diagonal elements (excluding the padding)
      //! @param VALUE the value to set to
      void SetDiagonal( const t_DataType &VALUE)
      {
        cl_int error_number( CL_SUCCESS);              // Error code var

        if( GetNumberOfElements() == 0)
        {
          return;
        }
        const size_t number_diagonal_elements( std::min( m_NumberRows - m_RowPadding, m_NumberCols - m_ColPadding));
        const size_t block_size( 256);
        const size_t stride( m_NumberCols + 1);

        // Create the kernel
        cl::Kernel kernel( GetTools().GetBufferProgram( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue), "FillStrided", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const cl::NDRange block_dimensions( block_size);
        const cl::NDRange offset;
        const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, number_diagonal_elements));

        error_number  = kernel.setArg( 0, m_Data);
        error_number |= kernel.setArg( 1, cl_uint( stride));
        error_number |= kernel.setArg( 2, cl_uint( number_diagonal_elements * m_NumberCols));
        error_number |= kernel.setArg( 3, VALUE);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment from Matrix
      //! @param MATRIX the matrix used as source
      //! @return reference to this Matrix
      Matrix< t_DataType> &operator =( const Matrix< t_DataType> &MATRIX)
      {
        m_Data = MATRIX.m_Data;
        m_NumberRows = MATRIX.m_NumberRows;
        m_NumberCols = MATRIX.m_NumberCols;
        m_RowPadding = MATRIX.m_RowPadding;
        m_ColPadding = MATRIX.m_ColPadding;
        m_Queue = MATRIX.m_Queue;

        // return reference to this Matrix
        return *this;
      }

      //! @brief returns an element in the vector
      //! @param ROW, COL position of the element
      //! @return the element
      const t_DataType operator()( const size_t ROW, const size_t COL) const
      {
        t_DataType ptr;
        m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, ( ROW * m_NumberCols + COL) * sizeof( t_DataType), sizeof( t_DataType), &ptr);
        return ptr;
      }

      //! @brief return a sub matrix from given ROW over NUMBER_ROWS
      //! @brief ROW the row to start from
      //! @param NUMBER_ROWS number of rows in the submatrix
      //! @return Matrix with subbuffer of this class
      Matrix< t_DataType> SubMatrix( const size_t ROW, const size_t NUMBER_ROWS) const
      {
        Matrix< t_DataType> new_matrix;
        cl_int error_number( CL_SUCCESS);

        const cl_buffer_region region = { ROW * sizeof( t_DataType) * m_NumberCols, NUMBER_ROWS * sizeof( t_DataType) * m_NumberCols};
        new_matrix.m_Data = Buffer( m_Data.createSubBuffer( CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, static_cast< const void *>( &region), &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "unable to create subbuffer: " + Tools::ErrorString( error_number));

        new_matrix.m_Queue      = m_Queue;
        new_matrix.m_NumberCols = m_NumberCols;
        new_matrix.m_ColPadding = m_ColPadding;
        new_matrix.m_NumberRows = NUMBER_ROWS;
        new_matrix.m_RowPadding = 0;

        // end
        return new_matrix;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_NumberRows, ISTREAM);
        io::Serialize::Read( m_NumberCols, ISTREAM);
        io::Serialize::Read( m_RowPadding, ISTREAM);
        io::Serialize::Read( m_ColPadding, ISTREAM);

        // data
        linal::Matrix< t_DataType> host_matrix;
        io::Serialize::Read( host_matrix, ISTREAM);
        m_Data = Buffer::CreateBufferFromMatrix< t_DataType>
                 (
                   linal::Matrix< t_DataType>( host_matrix, m_RowPadding, m_ColPadding), m_Queue
                 );

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_NumberRows, OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_NumberCols, OSTREAM,      0) << '\n';
        io::Serialize::Write( m_RowPadding, OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_ColPadding, OSTREAM,      0) << '\n';
        io::Serialize::Write( GetHostMatrix( m_RowPadding, m_ColPadding), OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Matrix

    //! single instance of the Marix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Matrix< t_DataType>::s_Instance;

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_MATRIX_H_ 
