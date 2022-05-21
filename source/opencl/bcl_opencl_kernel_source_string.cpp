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
#include "opencl/bcl_opencl_kernel_source_string.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_extensions.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelSourceString::KernelSourceString() :
      m_Source()
    {
    }

    //! @brief construct from string
    //! @param SOURCE_STRING kernel source string
    KernelSourceString::KernelSourceString( const std::string &SOURCE_STRING) :
      m_Source( SOURCE_STRING)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceString
    KernelSourceString *KernelSourceString::Clone() const
    {
      return new KernelSourceString( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceString::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceString::GetIdentifier() const
    {
      static const std::string s_identifier;
      return s_identifier;
    }

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceString::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      std::string source;
      source += "#define PRECISION " + util::CPPDataTypes::GetCPPString( PRECISION) + '\n';
      if( PRECISION == util::CPPDataTypes::e_Double)
      {
        if( EXTENSIONS.Find( GetExtensions().e_amd_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_amd_fp64.GetName() + " : enable\n";
        }
        else if( EXTENSIONS.Find( GetExtensions().e_khr_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_khr_fp64.GetName() + " : enable\n";
        }
        else
        {
          BCL_MessageCrt( "double precision is not supported by given extensions");
          return std::string();
        }
      }
      source += m_Source + "\n";
      return source;
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
    std::istream &KernelSourceString::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Source, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &KernelSourceString::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Source, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
