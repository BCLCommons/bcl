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
#include "util/bcl_util_ptr_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief Writes out a warning message that a pointer cast failed between the two types
    //! @param PTR_TYPE the pointer type that was received
    //! @param CAST_TYPE the type the pointer was cast to
    //! @note this function provides a convenient hook for debuggers, which need only set a single breakpoint in the cpp
    void NotifyUserBadPointerCast( const std::string &PTR_TYPE, const std::string &CAST_TYPE)
    {
      BCL_MessageCrt
      (
        "was not able to cast pointer from " + PTR_TYPE + " to " + CAST_TYPE
        + " Callstack: " + CallStack().String()
      );
    }

    //! @brief get the string for an empty pointer, e.g. "NULL"
    //! @return the string for an empty pointer, e.g. "NULL"
    const std::string &GetNullDescriptor()
    {
      static const std::string s_null_descriptor( "NULL");
      return s_null_descriptor;
    }

  } // namespace util
} // namespace bcl
