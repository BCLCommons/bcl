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
#include "opencl/bcl_opencl_command_queue.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_context.h"
#include "opencl/bcl_opencl_device.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CommandQueue::CommandQueue()
    {
    }

    //! @brief construct from context and device
    //! @param CONTEXT the context
    //! @param DEVICE the device for this command queue
    //! @param ERROR_PTR error will be written to this location
    CommandQueue::CommandQueue( const Context &CONTEXT, const Device &DEVICE, cl_int *ERROR_PTR) :
      cl::CommandQueue( CONTEXT, DEVICE, 0, ERROR_PTR)
    {
    }

    //! @brief construct from cl::CommandQueue
    //! @param COMMAND_QUEUE the cl::CommandQueue
    CommandQueue::CommandQueue( const cl::CommandQueue &COMMAND_QUEUE) :
      cl::CommandQueue( COMMAND_QUEUE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CommandQueue
    CommandQueue *CommandQueue::Clone() const
    {
      return new CommandQueue( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CommandQueue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the context for this command queue
    //! @param ERROR_PTR error will be written to this location
    //! @return Context the context associated with this queue
    Context CommandQueue::GetContext( cl_int *ERROR_PTR) const
    {
      return Context( cl::CommandQueue::getInfo< CL_QUEUE_CONTEXT>( ERROR_PTR));
    }

    //! @brief get the devices for this command queue
    //! @param ERROR_PTR error will be written to this location
    //! @return Device the device associated with this queue
    Device CommandQueue::GetDevice( cl_int *ERROR_PTR) const
    {
      return Device( cl::CommandQueue::getInfo< CL_QUEUE_DEVICE>( ERROR_PTR));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief comparison less than
    //! @param RHS right hand side queue
    //! @return true is this queue id is smaller than the rhs queue
    bool CommandQueue::operator <( const CommandQueue &RHS) const
    {
      return cl::CommandQueue::operator ()() < RHS();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CommandQueue::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CommandQueue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
