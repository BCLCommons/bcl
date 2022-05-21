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

#ifndef BCL_OPENCL_KERNEL_SOURCE_INTERFACE_H_
#define BCL_OPENCL_KERNEL_SOURCE_INTERFACE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_cpp_data_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KernelSourceInterface
    //! @brief this class defines, how source for kernels are handled.
    //! @details kernel sources are defined in a way that they can be compiled for different floating point precisions.
    //!          given all device extensions, in checks if the device extensions can support the requested precision and
    //!          will return an empty source if this is not the case, otherwise will will prepend a
    //!          "#define PRECISION float" or according to the desired precision in the source, so that the program can
    //!          be compiled accordingly
    //!
    //! @remarks example unnecessary
    //!
    //! @author woetzen
    //! @date Aug 20, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KernelSourceInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new KernelSourceInterface
      virtual KernelSourceInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief identifier for kernel source, so that the progam build can be cached based on it
      //! @return identifier like the filename or the function name
      virtual const std::string &GetIdentifier() const = 0;

      //! @brief get the source for compilation
      //! @param PRECISION precision of the kernel
      //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
      //! @return source with precision set correctly
      virtual std::string GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const = 0;

      //! @brief is the given Precision compatible with the given extension
      //! @param PRECISION the precision in question
      //! @param EXTENSIONS extensions of the devices
      //! @return true if the precision is supported
      static bool PrecisionCompatibleWithExtensions( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS);

      //! @brief get compiler flags necessary for precision
      //! @param PRECISION the precision of the kernel to be compiled
      //! @return string containing the compiler options that should be used
      static std::string GetPrecisionCompilerOptions( const util::CPPDataTypes::Types PRECISION);

      //! @brief get additional compiler flags
      //! @return string containing the additional compiler options that are desired by the user
      static std::string GetAdditionalCompilerOptions();

    }; // class KernelSourceInterface

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_KERNEL_SOURCE_INTERFACE_H_ 
