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
#include "command/bcl_command_parameter_check_extension.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_extension.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckExtension :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckExtension *Clone() const
    {
      return new ExampleCommandParameterCheckExtension( *this);
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

      // this is the extension to be checked for
      const std::string filename_extension( ".pdb");
      // filename that has pdb extension, i.e. is allowed
      const std::string filename_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      // filename that has different extension, i.e. is not allowed
      const std::string filename_fasta( AddExampleInputPathToFilename( e_Biology, "1ubi.fasta"));

      // test default constructor
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ext_undef( new command::ParameterCheckExtension());

      // construct a check extension for pdb files behind an ParameterCheck interface
      command::ParameterCheckExtension parameter_check_ext_pdb( filename_extension);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck(
        GetStaticClassName< command::ParameterCheckExtension>(),
        parameter_check_ext_pdb.GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check if pdb filename has valid extension
      BCL_ExampleCheck( parameter_check_ext_undef->IsAllowedParameter( filename_pdb, "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ext_pdb.IsAllowedParameter( filename_pdb, "", error_stream), true);

      // check if fasta file has valid extensions
      BCL_ExampleCheck( parameter_check_ext_undef->IsAllowedParameter( filename_fasta, "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ext_pdb.IsAllowedParameter( filename_fasta, "", error_stream), false);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_check_ext_pdb.WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( parameter_check_ext_pdb);

      // read file back into different parameter
      command::ParameterCheckExtension parameter_check_ext_pdb_read;
      ReadBCLObject( parameter_check_ext_pdb_read);

      // perform test from top with the read object
      BCL_ExampleCheck
      (
        parameter_check_ext_pdb_read.IsAllowedParameter( filename_pdb, "", error_stream),
        parameter_check_ext_pdb.IsAllowedParameter( filename_pdb, "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_ext_pdb_read.IsAllowedParameter( filename_fasta, "", error_stream),
        parameter_check_ext_pdb.IsAllowedParameter( filename_fasta, "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckExtension

  const ExampleClass::EnumType ExampleCommandParameterCheckExtension::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckExtension())
  );

} // namespace bcl
