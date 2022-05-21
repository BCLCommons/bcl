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

#ifndef BCL_OPENCL_KERNEL_SOURCE_STRING_H_
#define BCL_OPENCL_KERNEL_SOURCE_STRING_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_kernel_source_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KernelSourceString
    //! @brief this class handles kernel sources, that are given as a string, like hardcoded kernels in the library
    //! @details it implements the KernelSourceInterface by storing the kernel source in a string member variable
    //!
    //! @see @link example_opencl_kernel_source_interface.cpp @endlink
    //! @author woetzen
    //! @date Aug 20, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KernelSourceString :
      public KernelSourceInterface
    {

    //////////
    // data //
    //////////

      std::string m_Source; //!< actul source string for that kernel

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      KernelSourceString();

      //! @brief construct from string
      //! @param SOURCE_STRING kernel source string
      KernelSourceString( const std::string &SOURCE_STRING);

      //! @brief Clone function
      //! @return pointer to new KernelSourceString
      KernelSourceString *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief identifier for kernel source, so that the progam build can be cached based on it
      //! @return identifier like the filename or the function name
      const std::string &GetIdentifier() const;

      //! @brief get the source for compilation
      //! @param PRECISION precision of the kernel
      //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
      //! @return source with precision set correctly
      std::string GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const;

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class KernelSourceString

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_KERNEL_SOURCE_STRING_H_ 
