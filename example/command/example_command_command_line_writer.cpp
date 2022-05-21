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
#include "command/bcl_command_command_line_writer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_command_line_writer.cpp
  //!
  //! @author mendenjl
  //! @date Jul 15, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandCommandLineWriter :
    public ExampleInterface
  {
  public:

    ExampleCommandCommandLineWriter *Clone() const
    {
      return new ExampleCommandCommandLineWriter( *this);
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

      // create sample commands to check that the created command lines is right
      storage::Map< std::string, std::string> command_to_expected_output;

      // make sure that basic args work properly
      command_to_expected_output[ "+55.34e-59abc_"] = "+55.34e-59abc_ \n";
      // make sure that arguments that need quotes get quoted
      command_to_expected_output[ "hello world"] = "'hello world' \n";
      command_to_expected_output[ "func(arg)"] = "'func(arg)' \n";
      command_to_expected_output[ "\"I am quoted\""] = "'\"I am quoted\"' \n";

      // now commands that appear different on different shells
      command_to_expected_output[ "!!"] =
        "                     TCSH / CSH / BASH / SH : \\!\\! \n"
        "Windows Powershell / Windows Command Prompt : !! \n";
      command_to_expected_output[ "#"] =
        "BASH / SH / Windows Powershell / Windows Command Prompt : '#' \n"
        "                                             TCSH / CSH : # \n";
      command_to_expected_output[ "~"] =
        "                     TCSH / CSH / BASH / SH : '~' \n"
        "Windows Powershell / Windows Command Prompt : ~ \n";
      command_to_expected_output[ "\\~*"] =
        "                                               BASH / SH : '\\\\~*' \n"
        "TCSH / CSH / Windows Powershell / Windows Command Prompt : '\\~*' \n";
      command_to_expected_output[ "`"] =
        "TCSH / CSH / BASH / SH / Windows Command Prompt : '`' \n"
        "                             Windows Powershell : `` \n";
      command_to_expected_output[ "^"] =
        "TCSH / CSH / BASH / SH / Windows Powershell : ^ \n"
        "                     Windows Command Prompt : ^^ \n";
      command_to_expected_output[ "'"] =
        "TCSH / CSH / BASH / SH : \\' \n"
        "Windows Command Prompt : ^' \n"
        "    Windows Powershell : `' \n";
      command_to_expected_output[ "~!@#$%^&*()_+`0-=p[]{}\\|l;':\"<>,.?/"] =
        "             BASH / SH : '~\\!@#$%^&*()_+`0-=p[]{}\\\\|l;\\':\"<>,.?/' \n"
        "            TCSH / CSH : \\~\\!@#$%^\\&\\*\\(\\)_+\\`0-=p[]{}\\\\\\|l\\;\\':\\\"\\<\\>,.\\?/ \n"
        "Windows Command Prompt : '~!@#$%^^&*()_+`0-=p[]{}\\|l;^':\"<>,.?/' \n"
        "    Windows Powershell : '~!@#$%^&*()_+``0-=p[]{}\\|l;`':\"<>,.?/' \n";

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      for
      (
        storage::Map< std::string, std::string>::const_iterator
          itr( command_to_expected_output.Begin()), itr_end( command_to_expected_output.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string &test_command( itr->first);
        const storage::Vector< std::string> test_command_vec( 1, test_command);
        const std::string &expected_result( itr->second);
        BCL_ExampleIndirectCheck
        (
          command::CommandLineWriter::CreateCommandLine( test_command_vec),
          expected_result,
          "CreateCommandLine(" + test_command + ")"
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandCommandLineWriter

  const ExampleClass::EnumType ExampleCommandCommandLineWriter::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandCommandLineWriter())
  );

} // namespace bcl

