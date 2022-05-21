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

#ifndef BCL_OPENCL_KERNEL_SOURCE_ALTERNATIVE_H_
#define BCL_OPENCL_KERNEL_SOURCE_ALTERNATIVE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_kernel_source_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KernelSourceAlternative
    //! @brief this class collects alternatives for Kernel sources
    //! @details It will iterate through all alternatives for GetSource, till it finds one, that does not return an
    //!          empty string. This can be used if one prefers to use a kernel source from a file, that can be modified,
    //!          but wants to have a fallback, in case the file is not available
    //!
    //! @see @link example_opencl_kernel_source_alternative.cpp @endlink
    //! @author woetzen
    //! @date Aug 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KernelSourceAlternative :
      public KernelSourceInterface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtrVector< KernelSourceInterface> m_Alternatives; //!< collection of kernel source alternatives

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from two alternatives
      //! @param ALTERNATIVE_SOURCE_1
      //! @param ALTERNATIVE_SOURCE_2
      KernelSourceAlternative
      (
        const KernelSourceInterface &ALTERNATIVE_SOURCE_1,
        const KernelSourceInterface &ALTERNATIVE_SOURCE_2
      );

      //! @brief Clone function
      //! @return pointer to new KernelSourceAlternative
      KernelSourceAlternative *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief identifier for kernel source, so that the progam build can be cached based on it
      //! @return identifier like the filename or the function name
      const std::string &GetIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the source for compilation
      //! @param PRECISION precision of the kernel
      //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
      //! @return source with precision set correctly
      std::string GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const;

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

    }; // class KernelSourceAlternative

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_KERNEL_SOURCE_ALTERNATIVE_H_ 
