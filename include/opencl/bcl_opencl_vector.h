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

#ifndef BCL_OPENCL_VECTOR_H_
#define BCL_OPENCL_VECTOR_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Vector
    //! @brief opencl vector class to ease creation of buffer objects from linal vectors while providing access to
    //!        size and padding information
    //!
    //! @tparam t_DataType can be float, double, complex, int, etc...
    //!
    //! @see @link example_opencl_vector.cpp @endlink
    //! @author loweew
    //! @date Mar 24, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Vector :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! size
      size_t       m_Size;

      //! padding
      size_t       m_Padding;

      //! command queue
      CommandQueue m_Queue;

      //! opencl buffer
      Buffer       m_Data;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Vector() :
        m_Size(),
        m_Padding(),
        m_Queue(),
        m_Data()
      {
      }

      //! @brief constructor - initializes elements with 0
      //! @param SIZE size of vector
      //! @param QUEUE command queue
      //! @param PADDING the padding
      Vector
      (
        const size_t SIZE, const CommandQueue &QUEUE,
        const size_t PADDING = 0
      ) :
        m_Size( SIZE + PADDING),
        m_Padding( PADDING),
        m_Queue( QUEUE),
        m_Data( cl::Buffer( m_Queue.GetContext(), CL_MEM_READ_WRITE, m_Size * sizeof( t_DataType), NULL, NULL))
      {
        Fill( t_DataType( 0), 0, m_Size);
      }

      //! @brief constructor - initialize to ELEMENT and padding to 0
      //! @param SIZE size of vector
      //! @param QUEUE command queue
      //! @param PADDING the padding
      //! @param ELEMENT the element to initialize buffer with
      Vector
      (
        const size_t SIZE, const CommandQueue &QUEUE,
        const size_t PADDING, const t_DataType &ELEMENT
      ) :
        m_Size( SIZE + PADDING),
        m_Padding( PADDING),
        m_Queue( QUEUE),
        m_Data( cl::Buffer( m_Queue.GetContext(), CL_MEM_READ_WRITE, m_Size * sizeof( t_DataType), NULL, NULL))
      {
        Fill( ELEMENT, 0, m_Size - PADDING);
        Fill( t_DataType( 0), m_Size - PADDING, PADDING);
      }

      //! @brief constructor from linal vector
      //! @param VECTOR vector to create buffer from
      //! @param QUEUE command queue
      //! @param PADDING the amount of padding you want to add to the vector
      Vector
      (
        const linal::VectorConstInterface< t_DataType> &VECTOR,
        const CommandQueue &QUEUE,
        const size_t PADDING = 0
      ) :
        m_Size( VECTOR.GetSize() + PADDING),
        m_Padding( PADDING),
        m_Queue( QUEUE),
        m_Data( Buffer::CreateBufferFromVector< t_DataType>( linal::Vector< t_DataType>( VECTOR, PADDING), m_Queue))
      {
      }

      //! @brief constructor from pointer
      //! @param SIZE number of elements
      //! @param DATA pointer to data
      //! @param QUEUE command queue
      Vector( const size_t SIZE, const t_DataType *DATA, const CommandQueue &QUEUE) :
        m_Size( SIZE),
        m_Padding( 0),
        m_Queue( QUEUE),
        m_Data( cl::Buffer( m_Queue.GetContext(), CL_MEM_READ_WRITE, SIZE * sizeof( t_DataType), NULL))
      {
        m_Queue.enqueueWriteBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * SIZE, DATA);
      }

      //! @brief copy constructor
      //! @param VECTOR the vector to be copied
      Vector( const Vector< t_DataType> &VECTOR) :
        m_Size( VECTOR.m_Size),
        m_Padding( VECTOR.m_Padding),
        m_Queue( VECTOR.m_Queue),
        m_Data( VECTOR.m_Data)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Vector
      Vector< t_DataType> *Clone() const
      {
        return new Vector< t_DataType>( *this);
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

      //! @brief get number of elements
      //! @return size
      size_t GetSize() const
      {
        return m_Size;
      }

      //! @brief get padding
      //! @return padding
      size_t GetPadding() const
      {
        return m_Padding;
      }

      //! @brief gives the buffer object
      //! @return buffer
      const Buffer &GetData() const
      {
        return m_Data;
      }

      //! @brief gives the queue associated with this vector
      //! @return the command queue
      const CommandQueue &GetQueue() const
      {
        return m_Queue;
      }

      //! @brief returns linal::Vector and you can choose how much padding to take off
      //! @param PADDING number of elements to remove or to unpad
      //! @return the linal::Vector
      linal::Vector< t_DataType> GetHostVector( const size_t &PADDING = 0) const
      {
        linal::Vector< t_DataType> tmp( m_Size);
        cl_int error_number = m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, 0, sizeof( t_DataType) * m_Size, tmp.Begin());
        BCL_Assert( error_number == CL_SUCCESS, "error in GetHostVector(): " + Tools::ErrorString( error_number));
        return tmp.CreateSubVector( m_Size - PADDING);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief overwrites contents of buffer object in vector without reallocating the memory
      //! @param DATA pointer to the data
      void OverWriteContents( t_DataType *DATA)
      {
        cl_int error = CL_SUCCESS;
        error = m_Queue.enqueueWriteBuffer( m_Data, CL_TRUE, 0, m_Size, DATA);
        BCL_Assert( error == CL_SUCCESS, "error in OverWriteContents: " + Tools::ErrorString( error));
      }

      //! @brief sets position to value
      //! @param POS position
      //! @param VALUE the new value
      void SetValue( const size_t POS, const t_DataType &VALUE)
      {
        Fill( VALUE, POS, 1);
      }

      //! @brief fill the vector with given value up to number of elements
      //! @param VALUE the value to set every element to
      //! @param OFFSET offset from start
      //! @param NUMBER_ELEMENTS number of elements to set
      void Fill( const t_DataType &VALUE, const size_t OFFSET, const size_t NUMBER_ELEMENTS)
      {
        cl_int error_number( CL_SUCCESS);              // Error code var

        if( m_Size == 0)
        {
          return;
        }
        BCL_Assert
        (
          OFFSET + NUMBER_ELEMENTS <= m_Size,
          "vector too short! " + util::Format()( m_Size) + '<' + util::Format()( NUMBER_ELEMENTS + OFFSET)
        );
        const size_t block_size( 256);

        // Create the kernel
        cl::Kernel kernel( GetTools().GetBufferProgram( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue), "FillBuffer", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const cl::NDRange block_dimensions( block_size);
        const cl::NDRange offset;
        const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, NUMBER_ELEMENTS));

        error_number  = kernel.setArg( 0, m_Data);
        error_number |= kernel.setArg( 1, cl_uint( OFFSET));
        error_number |= kernel.setArg( 2, cl_uint( NUMBER_ELEMENTS));
        error_number |= kernel.setArg( 3, VALUE);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment from vector
      //! @param VECTOR the vector used as source
      //! @return reference to this vector
      Vector< t_DataType> &operator =( const Vector< t_DataType> &VECTOR)
      {
        m_Data = VECTOR.m_Data;
        m_Size = VECTOR.m_Size;
        m_Padding = VECTOR.m_Padding;
        m_Queue = VECTOR.m_Queue;

        // return reference to this Matrix
        return *this;
      }

      //! @brief returns an element in the vector
      //! @param POS position of the element
      //! @return the element
      const t_DataType operator()( const size_t POS) const
      {
        t_DataType ptr;
        m_Queue.enqueueReadBuffer( m_Data, CL_TRUE, POS * sizeof( t_DataType), sizeof( t_DataType), &ptr);
        return ptr;
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
        io::Serialize::Read( m_Size   , ISTREAM);
        io::Serialize::Read( m_Padding, ISTREAM);

        // data
        linal::Vector< t_DataType> host_vector;
        io::Serialize::Read( host_vector, ISTREAM);
        m_Data = Buffer::CreateBufferFromVector( linal::Vector< t_DataType>( host_vector, m_Padding), m_Queue);

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
        io::Serialize::Write( m_Size   , OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_Padding, OSTREAM,      0) << '\n';
        io::Serialize::Write( GetHostVector( m_Padding), OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Vector

    //! single instance of the Marix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Vector< t_DataType>::s_Instance;

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_VECTOR_H_ 
