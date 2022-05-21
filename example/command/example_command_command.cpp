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
#include "command/bcl_command_command.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_command.cpp
  //!
  //! @author heinzes1
  //! @date 11/06/2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandCommand :
    public ExampleInterface
  {
  public:

    ExampleCommandCommand *Clone() const
    {
      return new ExampleCommandCommand( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief called by the function Run with various command line arguments
    //! @return returns 0 upon completion
    int MainProgram( const int number_arguments, const std::string arguments[]) const
    {
      // the command object
      command::Command command_line;

      // parameterlist after .exe
      // input file parameter
      util::ShPtr< command::ParameterInterface> input_filename
      (
        new command::Parameter( "input_filename", "filename for input", command::ParameterCheckExtension( ".pdb"))
      );
      command_line.AddParameter( input_filename);

      // output file paramter
      util::ShPtr< command::ParameterInterface> output_filename
      (
        new command::Parameter( "output_filename", "filename for output", command::ParameterCheckExtension( ".bcl"))
      );
      command_line.AddParameter( output_filename);

      // simple flags used as switches
      // switch flag on
      util::ShPtr< command::FlagInterface> switch_on
      (
        new command::FlagStatic( "switch_on", "use this to switch something on")
      );

      // switch flag off
      util::ShPtr< command::FlagInterface> switch_off
      (
        new command::FlagStatic( "switch_off", "use this to switch something off")
      );

      // flag that tells this example to write itself to a stream and retrive itself from the strem
      util::ShPtr< command::FlagInterface> flag_read_write
      (
        new command::FlagStatic( "read_write_command_object", "use this to write and read command object")
      );

      // push flags to commandline
      command_line.AddFlag( switch_off);
      command_line.AddFlag( switch_on);
      command_line.AddFlag( flag_read_write);

      // flags with parameters
      // flag with static parameter list
      util::ShPtr< command::FlagStatic> blosum_matrix
      (
        new command::FlagStatic( "blosum", "this is a choice for the blosum matrix")
      );

      // create parameter with a list of allowed values
      const std::string blosum_matrix_strings[ 3] = { "BLOSUM45", "BLOSUM90", "BLOSUM120"};
      storage::Vector< std::string> allowed_blosum_matrices( 3, blosum_matrix_strings);

      util::ShPtr< command::ParameterInterface> param_blosum
      (
        new command::Parameter
        (
          "blosum_matrix",
          "a choice of possible matrices",
          command::ParameterCheckAllowed( allowed_blosum_matrices)
        )
      );
      blosum_matrix->PushBack( param_blosum);

      command_line.AddFlag( blosum_matrix);

      // flag with static parameter list
      util::ShPtr< command::FlagStatic> blosum_weight
      (
        new command::FlagStatic( "blosum_weight", "this is weight for the blosum matrix score")
      );

      // create parameter with an allowed range
      util::ShPtr< command::ParameterInterface> param_weightblosum
      (
        new command::Parameter
        (
          "blosum_score_weight",
          "weight for blosum score",
          command::ParameterCheckRanged< double>( double( 0.0), double( 1.0))
        )
      );
      blosum_weight->PushBack( param_weightblosum);

      command_line.AddFlag( blosum_weight);

      // flag with dynamic parameterlist
      // it consists of the flag name, a description for the flag, a template parameter, that shows how one of the
      // parameters should look like, including an optional parametercheck, a min and a max size
      util::ShPtr< command::FlagInterface> file_list
      (
        new command::FlagDynamic
        (
          "filelist",
          "this is list of input .fasta file",
          command::Parameter
          (
            "fasta-file",
            "any fasta file",
            command::ParameterCheckExtension( ".fasta")
          ),
          2,
          2
        )
      );

      command_line.AddFlag( file_list);

      // output help of commandline

      command_line.WriteHelp( util::GetLogger());
      BCL_MessageStd( "\n");

      // read commandline from command object
      std::stringstream error_str;
      command_line.ReadArguments
      (
        storage::Vector< std::string>( size_t( number_arguments), arguments),
        error_str,
        true
      );
      BCL_MessageStd( "this would be the output if you would use an assertion" + error_str.str());

      // write the command after it was
      // synchronized with given commandline

      command_line.WriteUserCommand( util::GetLogger());

      if( flag_read_write->GetFlag())
      {
        WriteBCLObject( command_line);
        command::Command test_command;
        ReadBCLObject( test_command);
        BCL_MessageStd( "check if the command object is working");
        test_command.WriteUserCommand( util::GetLogger());
      }

      // return 0 upon completion
      return 0;
    } // Main_program

    int Run() const
    {

      // correct commandline
      {
        // this is the given command line
        const int number_arguments = 11;
        const std::string arguments[ number_arguments] =
        {
          "test.exe", "1ubi.pdb", "1ubi.bcl", "-switch_off", "-blosum", "BLOSUM90", "-blosum_weight",
          "0.5", "-filelist", "1ubi.fasta", "1pub.fasta"
        };

        MainProgram( number_arguments, arguments);
      }

      // incorrect commandline
      {
        const int number_arguments = 12;
        const std::string arguments[ number_arguments] =
        {
          "test.exe", "1ubi.bcll", "1bla.pdd", "-switch_oof", "-blosum", "BLOSUM99", "-blosum_weight", "1.5",
          "-filelist", "1ubi.fasta", "1pub.fasta", "1uba.pdb"
        };

        MainProgram( number_arguments, arguments);
      }

      // write and retrieve Command object from a STREAM
      {
        const int number_arguments = 2;
        const std::string arguments[ number_arguments] = { "test.exe", "-read_write_command_object"};

        MainProgram( number_arguments, arguments);
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

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

  }; //end ExampleCommandCommand

  const ExampleClass::EnumType ExampleCommandCommand::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandCommand())
  );

} // namespace bcl
