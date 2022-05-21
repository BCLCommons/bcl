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
#include "io/bcl_io_stream_buffer_classes.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_stream_buffer_classes.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoStreamBufferClasses :
    public ExampleInterface
  {
  public:

    ExampleIoStreamBufferClasses *Clone() const
    {
      return new ExampleIoStreamBufferClasses( *this);
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
      // file to be used
      const std::string example_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // read an uncompressed file
      io::IFStream read;
      BCL_MessageStd( "standard buffer: " + read.GetCompression().GetName());

      BCL_ExampleMustOpenInputFile( read, example_file_name);

      const std::string base_file_name( AddExampleOutputPathToFilename( io::GetNamespaceIdentifier(), "1ubi.pdb"));

      // read the entire into a string stream
      std::stringstream contents_stream;
      while( !read.eof())
      {
        std::string line;

        // retrieve one line from istream
        std::getline( read, line);

        // append it to the contents stream, along with a new line
        contents_stream << line << '\n';
      }

      // store the contents of the stream in a string
      const std::string contents( contents_stream.str());

      // try to write and read compressed files
      for
      (
        io::StreamBufferClasses::const_iterator
          itr( io::GetStreamBufferClasses().Begin()),
          itr_end( io::GetStreamBufferClasses().End());
        itr != itr_end;
        ++itr
      )
      {
        // check if the stream buffer class is part of the commandline
        BCL_Example_Check
        (
          dynamic_cast< const command::Parameter &>( *io::GetStreamBufferClasses().GetFlagFileCompression()->GetFirstParameter()).IsAllowedParameter( itr->GetName(), util::GetLogger()),
          "the type of compression is not in the commandline flag: " + itr->GetName()
        );

        // skip uncompressed
        if( ( **itr)->GetDefaultFileExtension().empty())
        {
          continue;
        }

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        // construct write stream
        io::OFStream write;

      /////////////////
      // data access //
      /////////////////

        // filename of example file to be written and read uncompressed and compressed
        const std::string current_file_name( base_file_name + "." + ( **itr)->GetDefaultFileExtension());

        // get default extension for that stream
        BCL_MessageStd( "current buffer: " + ( **itr)->GetClassIdentifier());
        BCL_MessageStd( "default file extension: " + ( **itr)->GetDefaultFileExtension());

        // open write stream, the correct buffer will be chose - depending on file extension
        BCL_ExampleMustOpenOutputFile( write, current_file_name);

        // check if it is open
        BCL_ExampleIndirectCheck( write.is_open(), true, "open ofstream");

        // check that correct compression is associated
        BCL_ExampleCheck( write.GetCompression(), *itr);

      ////////////////
      // operations //
      ////////////////

        BCL_MessageStd( "write");
        // write content to file
        write << contents;
        BCL_MessageStd( "close clear");

        // clear and close
        io::File::CloseClearFStream( write);

        // construct read stream to read the written file
        io::IFStream read_written;

        // open stream for reading
        BCL_ExampleMustOpenInputFile( read_written, current_file_name);

        // clear the flags (e.g. eof) from the read stream and move it back to the start of the file
        read.clear();
        read.seekg( 0, std::ios_base::beg);

        // ensure that the files have the same strings in the same order
        BCL_ExampleIndirectCheck
        (
          io::File::StreamsMatch( read_written, read),
          true,
          "content of compressed file compared to original uncompressed file for extension: " +
          ( **itr)->GetDefaultFileExtension()
        );
      }

      // clear and close
      io::File::CloseClearFStream( read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoStreamBufferClasses

  const ExampleClass::EnumType ExampleIoStreamBufferClasses::s_Instance
  (
    GetExamples().AddEnum( ExampleIoStreamBufferClasses())
  );

} // namespace bcl

