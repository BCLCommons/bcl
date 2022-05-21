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
#include "command/bcl_command_parameter_check_extensions_file_existence.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_extensions_file_existence.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckExtensionsFileExistence :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckExtensionsFileExistence *Clone() const
    {
      return new ExampleCommandParameterCheckExtensionsFileExistence( *this);
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
      // a stream for use in our test cases
      std::stringstream error_stream;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // filename that exists
      const std::string filename_one_ext( AddExampleInputPathToFilename( e_Biology, "1C1D"));
      const std::string filename_all_ext( AddExampleInputPathToFilename( e_Biology, "1ubi"));
      // list of allowed extensions
      const std::string all_extensions[] = { ".pdb", "A.fasta", "A.ascii"};

      // default and clone constructor is tested together with write and read functions

      // test ParameterCheckExtensionsFileExistence( const std::string) constructor
      util::ShPtr< command::ParameterCheckExtensionsFileExistence> parameter_check_extensions_file_existence
      (
        new command::ParameterCheckExtensionsFileExistence( ".pdb")
      );
      // test ParameterCheckExtensionsFileExistence( const storage::Vector< std::string>) constructor
      util::ShPtr< command::ParameterCheckInterface> parameter_check_extensions_file_existence_all
      (
        new command::ParameterCheckExtensionsFileExistence( storage::Vector< std::string>( 3, all_extensions))
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck
      (
        GetStaticClassName< command::ParameterCheckExtensionsFileExistence>(),
        parameter_check_extensions_file_existence->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // checks that the file exists in all of the specified forms, i.e. a .pdb, .fasta and .ascii file
      // for the specified filename
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence->IsAllowedParameter( filename_all_ext, "", error_stream),
        true
      );
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_all->IsAllowedParameter( filename_all_ext, "", error_stream),
        true
      );

      // only 1C1D.pdb exists, so first test should succeed, while the second should not
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence->IsAllowedParameter( filename_one_ext, "", error_stream),
        true
      );
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_all->IsAllowedParameter( filename_one_ext, "", error_stream),
        false
      );

      // make sure that a file that does not exist in all of the forms is not acceptable
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence->IsAllowedParameter( "abcd", "", error_stream),
        false
      );
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_all->IsAllowedParameter( "abcd", "", error_stream),
        false
      );

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_check_extensions_file_existence_all->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);
      // clear string behind the stringstream buffer
      output_stream.str( std::string());
      // test WriteListOfExtensions() and check if length > 0
      util::ShPtr< command::ParameterCheckExtensionsFileExistence>
      (
        parameter_check_extensions_file_existence_all
      )->WriteListOfExtensions( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( *parameter_check_extensions_file_existence_all);

      // test default constructor
      command::ParameterCheckExtensionsFileExistence parameter_check_extensions_file_existence_all_read;
      // test reading
      ReadBCLObject( parameter_check_extensions_file_existence_all_read);

      // test clone
      util::ShPtr< command::ParameterCheckInterface> parameter_check_extensions_file_existence_clone
      (
        parameter_check_extensions_file_existence_all_read.Clone()
      );

      // compare the two parameters
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( filename_all_ext, "", error_stream),
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( filename_all_ext, "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( filename_one_ext, "", error_stream),
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( filename_one_ext, "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( "abcd", "", error_stream),
        parameter_check_extensions_file_existence_clone->IsAllowedParameter( "abcd", "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckExtensionsFileExistence

  const ExampleClass::EnumType ExampleCommandParameterCheckExtensionsFileExistence::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckExtensionsFileExistence())
  );

} // namespace bcl

