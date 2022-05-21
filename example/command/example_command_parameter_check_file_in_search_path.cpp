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
#include "command/bcl_command_parameter_check_file_in_search_path.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_file_in_search_path.cpp
  //!
  //! @author mendenjl
  //! @date Jan 06, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckFileInSearchPath :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckFileInSearchPath *Clone() const
    {
      return new ExampleCommandParameterCheckFileInSearchPath( *this);
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
      const std::string filename_file_existent( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      // filename that does NOT exist
      const std::string filename_file_nonexistent( AddExampleInputPathToFilename( e_Biology, "bob.bob"));

      // test default constructor - that it works is tested below
      // parameter check file existence - checks that a given file exists
      util::ShPtr< command::ParameterCheckInterface> param_check_exists
      (
        new command::ParameterCheckFileInSearchPath
        (
          storage::Vector< std::string>::Create( filename_file_nonexistent, filename_file_existent)
        )
      );

      // test Clone() - that it works is tested below
      util::ShPtr< command::ParameterCheckInterface> param_check_exists_clone( param_check_exists->Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck
      (
        GetStaticClassName< command::ParameterCheckFileInSearchPath>(),
        param_check_exists_clone->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check that the file exists
      BCL_ExampleCheck( param_check_exists_clone->IsAllowedParameter( filename_file_existent, "", error_stream), true);

      // check that some nonexistent file does not exist
      BCL_ExampleCheck( param_check_exists_clone->IsAllowedParameter( filename_file_nonexistent, "", error_stream), false);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      param_check_exists_clone->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( *param_check_exists_clone);

      // read file back into different parameter
      command::ParameterCheckFileInSearchPath parameter_check_exists_read;
      ReadBCLObject( parameter_check_exists_read);

      // perform test from top with the read object
      // test with existing file
      BCL_ExampleCheck
      (
        parameter_check_exists_read.IsAllowedParameter( filename_file_existent, "", error_stream),
        param_check_exists_clone->IsAllowedParameter( filename_file_existent, "", error_stream)
      );
      // test with nonexisting file
      BCL_ExampleCheck
      (
        parameter_check_exists_read.IsAllowedParameter( filename_file_nonexistent, "", error_stream),
        param_check_exists_clone->IsAllowedParameter( filename_file_nonexistent, "", error_stream)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckFileInSearchPath

  const ExampleClass::EnumType ExampleCommandParameterCheckFileInSearchPath::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckFileInSearchPath())
  );

} // namespace bcl
