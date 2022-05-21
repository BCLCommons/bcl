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

#ifndef BCL_OPENCL_BUFFER_H_
#define BCL_OPENCL_BUFFER_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_tools.h"
#include "linal/bcl_linal_matrix_const_interface.h"
#include "linal/bcl_linal_matrix_interface.h"
#include "linal/bcl_linal_vector_interface.h"

// external includes - sorted alphabetically
#include <CL/cl.hpp>

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Buffer
    //! @brief Buffer is wrapped around cl::Buffer and provides a memory buffer object to store data
    //!
    //! @see @link example_opencl_buffer.cpp @endlink
    //! @author woetzen, loweew
    //! @date Dec 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Buffer :
      public cl::Buffer
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Buffer();

      //! @brief construct from cl::Buffer
      //! @param BUFFER the cl::Buffer
      Buffer( const cl::Buffer &BUFFER);

      //! @brief Clone function
      //! @return pointer to new Buffer
      Buffer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief size of memory in bytes
      //! @param ERROR_PTR error will be written to this location
      //! @return size_t size of buffer
      size_t GetSize( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the context for this buffer
      //! @param ERROR_PTR error will be written to this location
      //! @return Context the context
      Context GetContext( cl_int *ERROR_PTR = NULL) const;

      //! @brief the associated host pointer
      //! @param ERROR_PTR error will be written to this location
      //! @return void* the host pointer if any is associated, NULL otherwise
      void *GetHostPtr( cl_int *ERROR_PTR = NULL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief creates cl buffer on device from bcl matrix
      //! @param MATRIX matrix to create buffer from
      //! @param QUEUE the queue required for copying the data onto the buffer
      //! @param READ_ONLY bool - true if read only, false if not
      //! @return the allocated and copied to buffer
      template< typename t_DataType>
      static Buffer CreateBufferFromMatrix
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const CommandQueue &QUEUE,
        const bool READ_ONLY = false
      )
      {
        // number of elements to be copied
        const cl_uint number_elements( MATRIX.GetNumberOfElements());

        // is memory read_write or write_only
        const cl_mem_flags writable( READ_ONLY ? CL_MEM_READ_ONLY : CL_MEM_READ_WRITE);

        // error catching
        cl_int error_number = CL_SUCCESS;

        // allocate buffer and copy elements over
        Buffer buffer( cl::Buffer( QUEUE.GetContext(), writable, sizeof( t_DataType) * number_elements, NULL, &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));
        error_number = QUEUE.enqueueWriteBuffer( buffer, CL_TRUE, 0, sizeof( t_DataType) * number_elements, MATRIX.Begin());
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

        return buffer;
      }

      //! @brief creates cl buffer on device from bcl matrix
      //! @param SIZE size of buffer to allocate
      //! @param QUEUE the queue required for copying the data onto the buffer
      //! @return the allocated and copied to buffer
      template< typename t_DataType>
      static Buffer AllocateBufferOfSize( const size_t SIZE, const CommandQueue &QUEUE)
      {
        // error catching
        cl_int error_number = CL_SUCCESS;

        // allocate buffer and copy elements over
        Buffer buffer( cl::Buffer( QUEUE.GetContext(), CL_MEM_READ_WRITE, sizeof( t_DataType) * cl_uint( SIZE), NULL, &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation failed with: " + opencl::Tools::ErrorString( error_number));

        return buffer;
      }

      //! @brief creates cl buffer on device from bcl matrix
      //! @param MATRIX matrix to create buffer from
      //! @param QUEUE the queue required for copying the data onto the buffer
      //! @param READ_ONLY bool - true if read only, false if not
      //! @return the allocated and copied to buffer
      template< typename t_DataType>
      static Buffer CreateBufferFromVector
      (
        const linal::VectorConstInterface< t_DataType> &VECTOR,
        const CommandQueue &QUEUE,
        bool READ_ONLY = false
      )
      {
        // number of elements to be copied
        const cl_uint number_elements( VECTOR.GetSize());

        // is memory read_write or write_only
        const cl_mem_flags writable( READ_ONLY ? CL_MEM_READ_ONLY : CL_MEM_READ_WRITE);

        // error catching
        cl_int error_number = CL_SUCCESS;

        // allocate buffer and copy elements over
        Buffer buffer( cl::Buffer( QUEUE.GetContext(), writable, sizeof( t_DataType) * number_elements, NULL, &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));
        error_number = QUEUE.enqueueWriteBuffer( buffer, CL_TRUE, 0, sizeof( t_DataType) * number_elements, VECTOR.Begin());
        BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

        return buffer;
      }

    }; // class Buffer

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_BUFFER_H_
