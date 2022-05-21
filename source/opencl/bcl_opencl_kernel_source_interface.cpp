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
#include "opencl/bcl_opencl_kernel_source_interface.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_extensions.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief is the given Precision compatible with the given extension
    //! @param PRECISION the precision in question
    //! @param EXTENSIONS extensions of the devices
    //! @return true if the precision is supported
    bool KernelSourceInterface::PrecisionCompatibleWithExtensions( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS)
    {
      switch( PRECISION)
      {
        case util::CPPDataTypes::e_Float:
          return true;
        case util::CPPDataTypes::e_Double:
          return EXTENSIONS.Find( GetExtensions().e_khr_fp64) != EXTENSIONS.End() || EXTENSIONS.Find( GetExtensions().e_amd_fp64) != EXTENSIONS.End();
        default: break;
      }

      return false;
    }

    //! @brief get compiler flags necessary for precision
    //! @param PRECISION the precision of the kernel to be compiled
    //! @return string containing the compiler options that should be used
    std::string KernelSourceInterface::GetPrecisionCompilerOptions( const util::CPPDataTypes::Types PRECISION)
    {
      switch( PRECISION)
      {
        case util::CPPDataTypes::e_Float:
          // prevent warning
          // double-precision constant is represented as
          // single-precision constant because double is not enabled
          // PRECISION c_sub = 0.0;
          //                   ^
          return std::string( " -cl-single-precision-constant");
        default:
          return std::string();
      }

      return std::string();
    }

    //! @brief get additional compiler flags
    //! @return string containing the additional compiler options that are desired by the user
    std::string KernelSourceInterface::GetAdditionalCompilerOptions()
    {
      // kernels have to be compiled with -g to add debug information, so that it can be debugged
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        return std::string( " -g");
      }

      return std::string();
    }

  } // namespace opencl
} // namespace bcl
