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
#include "io/bcl_io_file.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_file.cpp
  //!
  //! @author alexanns, mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoFile :
    public ExampleInterface
  {
  public:

    ExampleIoFile *Clone() const
    { return new ExampleIoFile( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Try to register a filename for output
    //! @param FILENAME the file that the example is attempting to access
    //! @return true if the file is openable by this example
    bool RegisterOutputFilename( const std::string &FILENAME) const
    {
      std::ostringstream message_stream;
      if( !ExampleClass::RequestExampleFileAccess( GetClassIdentifier(), FILENAME, std::ios::out, message_stream))
      {
        ExampleClass::GetResults().InsertTest( GetClassIdentifier(), __LINE__, false, message_stream.str());
        BCL_MessageCrt( "ERROR: " + message_stream.str());
        return false;
      }
      return true;
    }

    //! @brief Try to register a filename for input
    //! @param FILENAME the file that the example is attempting to access
    //! @return true if the file is openable by this example
    bool RegisterInputFilename( const std::string &FILENAME) const
    {
      std::ostringstream message_stream;
      if( !ExampleClass::RequestExampleFileAccess( GetClassIdentifier(), FILENAME, std::ios::in, message_stream))
      {
        ExampleClass::GetResults().InsertTest( GetClassIdentifier(), __LINE__, false, message_stream.str());
        BCL_MessageCrt( "ERROR: " + message_stream.str());
        return false;
      }
      return true;
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a few example files
      const std::string example_output_file
      (
        AddExampleOutputPathToFilename( io::GetNamespaceIdentifier(), "example_tryopenofstream.txt")
      );
      // register the filename for this example; normally this would be done by BCL_ExampleMustOpenInput/OutputFile
      // macro, but it must be done explicitly here since we are testing functionality required by that function
      if( !RegisterOutputFilename( example_output_file))
      {
        return 1;
      }

      const std::string example_input_file( AddExampleInputPathToFilename( e_Biology, "table_example.txt"));
      if( !RegisterInputFilename( example_input_file))
      {
        return 1;
      }

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // check TryOpenOfstream function
      {
        io::OFStream ofstream;

        BCL_MessageStd( "check that TryOpenOfstream can correctly return true");
        // check that TryOpenOfstream can correctly return true
        BCL_ExampleCheck( io::File::TryOpenOFStream( ofstream, example_output_file), true);
        BCL_ExampleIndirectCheck( ofstream.is_open(), true, "TryOpenOFStream on valid file");
      }

      // check MustOpenOfstream function
      {
        // create ofstream "ofstream"
        io::OFStream ofstream;

        // give message
        BCL_MessageStd( "check that MustOpenOfstream can correctly NOT assert");
        // MustOpenOfstream should not assert
        io::File::MustOpenOFStream( ofstream, example_output_file);
        BCL_ExampleIndirectCheck( ofstream.is_open(), true, "MustOpenOFStream");
      }

      // check MustOpenIfstream function
      {
        // create ifstream "ifstream"
        io::IFStream ifstream;

        // give message
        BCL_MessageStd( "check that MustOpenIfstream can correctly NOT assert");
        // MustOpenIfstream should not assert
        io::File::MustOpenIFStream( ifstream, AddExampleInputPathToFilename( e_Biology, "table_example.txt"));

        BCL_ExampleIndirectCheck( ifstream.is_open(), true, "MustOpenIFStream");
      }

      // check TryOpenIfstream function
      {
        io::IFStream ifstream;

        BCL_MessageStd( "check that TryOpenIfstream can correctly return true");
        // check that TryOpenIfstream can correctly return true
        BCL_ExampleCheck( io::File::TryOpenIFStream( ifstream, example_input_file), true);
        BCL_ExampleIndirectCheck( ifstream.is_open(), true, "TryOpenIFStream");

        BCL_MessageStd( "check that TryOpenIfstream can correctly return false");
        // check that TryOpenIfstream can correctly return false
        BCL_ExampleCheck
        (
          io::File::TryOpenOFStream
          (
            ifstream,
            AddExampleOutputPathToFilename( util::GetNamespaceIdentifier(), "example_tryopenifstream.txt"),
            std::ios_base::in
          ),
          false
        );
      }

      // check Size function
      {
        // give message
        BCL_MessageStd( "make sure that the Size function works");

        const std::string file( AddExampleInputPathToFilename( e_Biology, "1ubiA.fasta"));
        if( !RegisterInputFilename( file))
        {
          return 1;
        }

        // make sure the file was found
        BCL_ExampleIndirectCheck( io::File::Size( file).First(), true, "file existence detection");
        BCL_ExampleIndirectCheck( io::File::Size( file).Second(), 106, "file size detection");
      }

      // check GetFullExtension function
      CheckGetFullExtensionFunction( "/home/bcl/.test.file.extension.*", "file.extension.");

      // check the GetLastExtension function to make sure it works properly
      CheckGetLastExtensionFunction( "/home/bcl/.test.file.extension.*", "");

      // check the GetLastExtension function to make sure it works properly
      CheckGetLastExtensionFunction( "/home/bcl/.test.file.extension.last~", "last");

      // check the CleanFilename function
      CheckCleanFilenameFunction( "/home/bcl/.test.file.extension.*", ".test.file.extension.");

      // check the RemoveLastExtension function
      CheckRemoveLastExtensionFunction( "/home/bcl/.test.file.extension.last~", "/home/bcl/.test.file.extension");

      // check the RemoveLastExtension function
      CheckRemoveLastExtensionFunction( "/home/bcl/.test.file.extension.*", "/home/bcl/.test.file.extension");

      // check the RemoveFullExtension function
      CheckRemoveFullExtensionFunction( "/home/bcl/.test.file.extension.last~", "/home/bcl/.test");

      // check the RemoveFullExtension function
      CheckRemoveFullExtensionFunction( "/home/bcl/.test.file.extension.*", "/home/bcl/.test");

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    // @brief CheckGetFullExtensionFunction checks the GetFullExtension function to make sure it works properly
    // @param FILENAME the file whose full extension will be gotten
    // @param CORRECT_FULL_EXTENSION the extension that should be gotten by the GetFullExtension function
    void CheckGetFullExtensionFunction( const std::string &FILENAME, const std::string &CORRECT_FULL_EXTENSION) const;

    // @brief CheckGetLastExtensionFunction checks the GetLastExtension function to make sure it works properly
    // @param FILENAME the file whose last extension will be gotten
    // @param CORRECT_LAST_EXTENSION the extension that should be gotten by the GetLastExtension function
    void CheckGetLastExtensionFunction( const std::string &FILENAME, const std::string &CORRECT_LAST_EXTENSION) const;

    // @brief CheckCleanFilenameFunction checks the CleanFilename function to make sure it works properly
    // @param FILENAME the file whose name will be cleaned
    // @param CORRECT_CLEAN_FILENAME clean filename that should be returned by the CleanFilename function
    void CheckCleanFilenameFunction( const std::string &FILENAME, const std::string &CORRECT_CLEAN_FILENAME) const;

    // @brief CheckRemoveLastExtensionFunction checks the RemoveLastExtension function to make sure it works properly
    // @param FILENAME the file whose name will have its last extension removed
    // @param CORRECT_NEW_FILENAME new filename that should be returned by the RemoveLastExtension function
    void CheckRemoveLastExtensionFunction
    (
      const std::string &FILENAME, const std::string &CORRECT_NEW_FILENAME
    ) const;

    // @brief CheckRemoveFullExtensionFunction checks the RemoveFullExtension function to make sure it works properly
    // @param FILENAME the file whose name will have its full extension removed
    // @param CORRECT_NEW_FILENAME new filename that should be returned by the RemoveFullExtension function
    void CheckRemoveFullExtensionFunction
    (
      const std::string &FILENAME, const std::string &CORRECT_NEW_FILENAME
    ) const;

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoFile

  const ExampleClass::EnumType ExampleIoFile::s_Instance
  (
    GetExamples().AddEnum( ExampleIoFile())
  );

  // @brief CheckGetFullExtensionFunction checks the GetFullExtension function to make sure it works properly
  // @param FILENAME the file whose full extension will be gotten
  // @param CORRECT_FULL_EXTENSION the extension that should be gotten by the GetFullExtension function
  void ExampleIoFile::CheckGetFullExtensionFunction
  (
    const std::string &FILENAME, const std::string &CORRECT_FULL_EXTENSION
  ) const
  {
    // create string "full_extension" and inialize with the extensions of "filename"
    const std::string full_extension( io::File::GetFullExtension( FILENAME));

    // write message
    BCL_MessageStd( "full extension of \"" + FILENAME + "\" is \"" + full_extension + "\"");

    // make sure the correct extension is gotten
    BCL_ExampleCheck( io::File::GetFullExtension( FILENAME), CORRECT_FULL_EXTENSION);
  }

  // @brief CheckGetLastExtensionFunction checks the GetLastExtension function to make sure it works properly
  // @param FILENAME the file whose last extension will be gotten
  // @param CORRECT_LAST_EXTENSION the extension that should be gotten by the GetLastExtension function
  void ExampleIoFile::CheckGetLastExtensionFunction
  (
    const std::string &FILENAME, const std::string &CORRECT_LAST_EXTENSION
  ) const
  {
    // create string "full_extension" and inialize with the extensions of "filename"
    const std::string full_extension( io::File::GetLastExtension( FILENAME));

    // write message
    BCL_MessageStd( "last extension of \"" + FILENAME + "\" is \"" + full_extension + "\"");

    // make sure the correct extension is gotten
    BCL_ExampleCheck( io::File::GetLastExtension( FILENAME), CORRECT_LAST_EXTENSION);
  }

  // @brief CheckCleanFilenameFunction checks the CleanFilename function to make sure it works properly
  // @param FILENAME the file whose name will be cleaned
  // @param CORRECT_CLEAN_FILENAME clean filename that should be returned by the CleanFilename function
  void ExampleIoFile::CheckCleanFilenameFunction
  (
    const std::string &FILENAME, const std::string &CORRECT_CLEAN_FILENAME
  ) const
  {
    // create string "cleaned_filename" and initialize with the cleaned filename of "FILENAME"
    const std::string cleaned_filename( io::File::CleanFilename( FILENAME, "~*%"));

    // write message
    BCL_MessageStd( "cleaned name of \"" + FILENAME + "\" is \"" + cleaned_filename + "\"");

    // make sure the filename was cleaned correctly
    BCL_ExampleCheck( io::File::CleanFilename( FILENAME, "~*%"), CORRECT_CLEAN_FILENAME);
  }

  // @brief CheckRemoveLastExtensionFunction checks the RemoveLastExtension function to make sure it works properly
  // @param FILENAME the file whose name will have its last extension removed
  // @param CORRECT_NEW_FILENAME new filename that should be returned by the RemoveLastExtension function
  void ExampleIoFile::CheckRemoveLastExtensionFunction
  (
    const std::string &FILENAME, const std::string &CORRECT_NEW_FILENAME
  ) const
  {
    // create string "new_filename" and initialize with the new filename of "FILENAME" given by RemoveLastExtension
    const std::string new_filename( io::File::RemoveLastExtension( FILENAME));

    // write message
    BCL_MessageStd
    (
      "name of \"" + FILENAME + "\" without last extension is \"" + new_filename + "\""
    );

    // make sure the filename correctly had its last extension removed correctly
    BCL_ExampleCheck( io::File::RemoveLastExtension( FILENAME), CORRECT_NEW_FILENAME);
  }

  // @brief CheckRemoveFullExtensionFunction checks the RemoveFullExtension function to make sure it works properly
  // @param FILENAME the file whose name will have its full extension removed
  // @param CORRECT_NEW_FILENAME new filename that should be returned by the RemoveFullExtension function
  void ExampleIoFile::CheckRemoveFullExtensionFunction
  (
    const std::string &FILENAME, const std::string &CORRECT_NEW_FILENAME
  ) const
  {
    // create string "new_filename" and initialize with the new filename of "FILENAME" given by RemoveLastExtension
    const std::string new_filename( io::File::RemoveFullExtension( FILENAME));

    // write message
    BCL_MessageStd
    (
      "name of \"" + FILENAME + "\" without full extension is \"" + new_filename + "\""
    );

    // make sure the filename correctly had its last extension removed correctly
    BCL_ExampleCheck( io::File::RemoveFullExtension( FILENAME), CORRECT_NEW_FILENAME);
  }

} // namespace bcl
