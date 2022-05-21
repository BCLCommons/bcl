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

#ifndef BCL_OPENCL_COMMAND_QUEUE_H_
#define BCL_OPENCL_COMMAND_QUEUE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <CL/cl.hpp>

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CommandQueue
    //! @brief this is the commandqueue that wraps around the cl::CommandQueue and communicates commands to the device and
    //!        its context
    //! @details the commandqueue is the handler to copy data from the host to the device, using a buffer; enques the
    //!          the executation of kernels and tasks, handles events and copies data back to the host.
    //!
    //! @see @link example_opencl_command_queue.cpp @endlink
    //! @author woetzen
    //! @date Nov 30, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CommandQueue :
      public cl::CommandQueue,
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CommandQueue();

      //! @brief construct from context and device
      //! @param CONTEXT the context
      //! @param DEVICE the device for this command queue
      //! @param ERROR_PTR error will be written to this location
      CommandQueue( const Context &CONTEXT, const Device &DEVICE, cl_int *ERROR_PTR = NULL);

      //! @brief construct from cl::CommandQueue
      //! @param COMMAND_QUEUE the cl::CommandQueue
      CommandQueue( const cl::CommandQueue &COMMAND_QUEUE);

      //! @brief Clone function
      //! @return pointer to new CommandQueue
      CommandQueue *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the context for this command queue
      //! @param ERROR_PTR error will be written to this location
      //! @return Context the context associated with this queue
      Context GetContext( cl_int *ERROR_PTR = NULL) const;

      //! @brief get the devices for this command queue
      //! @param ERROR_PTR error will be written to this location
      //! @return Device the device associated with this queue
      Device GetDevice( cl_int *ERROR_PTR = NULL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief comparison less than
      //! @param RHS right hand side queue
      //! @return true is this queue id is smaller than the rhs queue
      bool operator <( const CommandQueue &RHS) const;

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

    }; // class CommandQueue

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_COMMAND_QUEUE_H_ 
