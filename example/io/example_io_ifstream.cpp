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
#include "io/bcl_io_ifstream.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_ifstream.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoIFStream :
    public ExampleInterface
  {
  public:

    ExampleIoIFStream *Clone() const
    {
      return new ExampleIoIFStream( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // bcl ifstream
      io::IFStream read;

      const std::string reference_pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      std::stringstream pdb;
      BCL_ExampleMustOpenInputFile( read, reference_pdb_filename);
      pdb << read.rdbuf();
      io::File::CloseClearFStream( read);

      const std::string pdb_string( pdb.str());

    /////////////////
    // data access //
    /////////////////

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

        const io::DirectoryEntry filename( reference_pdb_filename + '.' + ( **itr)->GetDefaultFileExtension());
        // check if a compressed file with the file extension exists
        if( !filename.DoesExist())
        {
          BCL_MessageCrt( "unable to find compressed file: " + filename.GetFullName());
          continue;
        }

        BCL_MessageStd( "reading in file " + filename.GetFullName());

        // open file with the current compression
        BCL_ExampleMustOpenInputFile( read, filename.GetFullName());

        // read into stringstream
        std::stringstream current_file;
        current_file << read.rdbuf();
        // close and clear
        io::File::CloseClearFStream( read);

        const std::string current_string( current_file.str());

        // compare it to the original uncompressed source
        BCL_ExampleCheck( pdb_string.size(), current_string.size());
        BCL_ExampleCheck( pdb_string, current_string);
      }

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

  }; //end ExampleIoIFStream

  const ExampleClass::EnumType ExampleIoIFStream::s_Instance
  (
    GetExamples().AddEnum( ExampleIoIFStream())
  );

} // namespace bcl
