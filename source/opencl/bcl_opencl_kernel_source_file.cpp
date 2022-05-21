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
#include "opencl/bcl_opencl_kernel_source_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "opencl/bcl_opencl_extensions.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from filename
    //! @param FILE_NAME filename of the kernel
    KernelSourceFile::KernelSourceFile( const std::string &FILE_NAME) :
      m_FileName( FILE_NAME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceFile
    KernelSourceFile *KernelSourceFile::Clone() const
    {
      return new KernelSourceFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceFile::GetIdentifier() const
    {
      return m_FileName;
    }

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceFile::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      // source
      std::string source;

      // open file and acquire content
      io::IFStream read;
      if( !io::File::TryOpenIFStream( read, KernelSources::AddKernelPath( m_FileName)))
      {
        BCL_MessageDbg( "cannot find kernel source file: " + KernelSources::AddKernelPath( m_FileName));
        return source;
      }

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

      // filecontent
      std::stringstream file_content;
      file_content << read.rdbuf();
      io::File::CloseClearFStream( read);

      // add file content to source
      source += file_content.str() + "\n";

      // end
      return source;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelSourceFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FileName, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KernelSourceFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FileName, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
