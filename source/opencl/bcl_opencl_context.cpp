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
#include "opencl/bcl_opencl_context.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Context::Context()
    {
    }

    //! @brief construct from cl::Context
    //! @param CONTEXT the cl::Context to construct from
    Context::Context( const cl::Context &CONTEXT) :
      cl::Context( CONTEXT)
    {
    }

    //! @brief Context from devices
    //! @param ERROR_PTR pointer to error storage
    Context::Context( const storage::Vector< Device> &DEVICES, cl_int *ERROR_PTR)
    {
      // convert to cl::Devices
      std::vector< cl::Device> devices( DEVICES.Begin(), DEVICES.End());
      BCL_MessageVrb( "Creating context");
      cl_int error( CL_SUCCESS);
      cl::Context new_context( devices, NULL, NULL, NULL, &error);
      if( error == CL_SUCCESS)
      {
        BCL_MessageVrb( "Created new_context");
        cl::Context::operator =( new_context);
      }
      else
      {
        BCL_MessageCrt( "Error creating context: " + Tools::ErrorString( error));
        Tools::AssignError( ERROR_PTR, error);
      }
    }

    //! @brief Clone function
    //! @return pointer to new Context
    Context *Context::Clone() const
    {
      return new Context( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Context::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Context::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Context::Write( std::ostream &OSTREAM, const size_t INDENT) const
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
