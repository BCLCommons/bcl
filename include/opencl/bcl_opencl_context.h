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

#ifndef BCL_OPENCL_CONTEXT_H_
#define BCL_OPENCL_CONTEXT_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

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
    //! @class Context
    //! @brief opencl context
    //! @details this class wraps around cl::Context and provides functionality to query information on the context.
    //!          a Context is constructed from an opencl device and associates commandqueues and memory objects. Multiple
    //!          contexts can be associated with a device, but commands and memory objects cannot be interchanged between
    //!          different contexts on the same device or on diferent devices.
    //!
    //! @see @link example_opencl_context.cpp @endlink
    //! @author woetzen
    //! @date Nov 19, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Context :
      public cl::Context
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Context();

      //! @brief construct from cl::Context
      //! @param CONTEXT the cl::Context to construct from
      Context( const cl::Context &CONTEXT);

      //! @brief Context from devices
      //! @param DEVICES devices to add to context
      //! @param PLATFORM platform the devices reside on
      //! @param ERROR_PTR pointer to error storage
      Context( const storage::Vector< Device> &DEVICES, cl_int *ERROR_PTR = NULL);

      //! @brief Clone function
      //! @return pointer to new Context
      Context *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

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

    private:

    }; // class Context

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_CONTEXT_H_
