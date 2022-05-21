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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opencl/bcl_opencl_buffer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Buffer::Buffer()
    {
    }

    //! @brief construct from cl::Buffer
    //! @param BUFFER the cl::Buffer
    Buffer::Buffer( const cl::Buffer &BUFFER) :
      cl::Buffer( BUFFER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Buffer
    Buffer *Buffer::Clone() const
    {
      return new Buffer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Buffer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief size of memory in bytes
    //! @param ERROR_PTR error will be written to this location
    //! @return size_t size of buffer
    size_t Buffer::GetSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief get the context for this buffer
    //! @param ERROR_PTR error will be written to this location
    //! @return Context the context
    Context Buffer::GetContext( cl_int *ERROR_PTR) const
    {
      return Context( getInfo< CL_MEM_CONTEXT>( ERROR_PTR));
    }

    //! @brief the associated host pointer
    //! @param ERROR_PTR error will be written to this location
    //! @return void* the host pointer if any is associated, NULL otherwise
    void *Buffer::GetHostPtr( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_MEM_HOST_PTR>( ERROR_PTR);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Buffer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Buffer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl
