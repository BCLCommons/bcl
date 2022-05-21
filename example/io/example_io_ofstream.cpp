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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "io/bcl_io_ofstream.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_ofstream.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoOFStream :
    public ExampleInterface
  {
  public:

    ExampleIoOFStream *Clone() const
    {
      return new ExampleIoOFStream( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // open normal file to copy content to zipped stream
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      std::stringstream pdb;
      pdb << read.rdbuf();
      io::File::CloseClearFStream( read);

      const std::string pdb_string( pdb.str());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct ofstream
      io::OFStream write;

      const std::string base_file_name
      (
        AddExampleOutputPathToFilename( io::GetNamespaceIdentifier(), "1ubi.pdb.ofstream")
      );

      // construct and open output gz stream
      // iterate over all available compressions
      for
      (
        io::StreamBufferClasses::const_iterator
          itr( io::GetStreamBufferClasses().Begin()), itr_end( io::GetStreamBufferClasses().End());
        itr != itr_end;
        ++itr
      )
      {
        // skip uncompressed
        if( *itr == io::GetStreamBufferClasses().e_Uncompressed)
        {
          continue;
        }

        const std::string filename( base_file_name + "." + ( **itr)->GetDefaultFileExtension());

        // open file with the current compression
        if( !ExampleClass::RequestExampleFileAccess( GetClassIdentifier(), filename, std::ios::out, util::GetLogger()))
        {
          return 1;
        }
        BCL_ExampleCheck( io::File::TryOpenOFStream( write, filename), true);
        write << pdb_string;
        // close and clear
        io::File::CloseClearFStream( write);
      }

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoOFStream

  const ExampleClass::EnumType ExampleIoOFStream::s_Instance
  (
    GetExamples().AddEnum( ExampleIoOFStream())
  );

} // namespace bcl
