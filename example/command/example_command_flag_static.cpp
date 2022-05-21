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
#include "command/bcl_command_flag_static.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_flag_static.cpp
  //!
  //! @author heinzes1
  //! @date 11/06/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandFlagStatic :
    public ExampleInterface
  {
  public:

    ExampleCommandFlagStatic *Clone() const
    {
      return new ExampleCommandFlagStatic( *this);
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

      // create flag names and descriptions
      const std::string switch_flag_name( "switch");
      const std::string switch_flag_description( "to switch some functionality on or off");
      const std::string blosum_flag_name( "blosum");
      const std::string blosum_flag_description( "this is a choice for the blosum matrix");
      const std::string pam_flag_name( "pam");
      const std::string pam_flag_description( "this is a choice for the pam matrix");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      command::FlagStatic default_constructed_flag;

      // static flag with no parameters (i.e., a switch)
      command::FlagStatic switch_flag( switch_flag_name, switch_flag_description);

      // flag with static parameter list
      // a static flag consists of the flag name, a description for the flag, and an optional parameter
      util::ShPtr< command::FlagInterface> blosum_flag
      (
        new command::FlagStatic
        (
          blosum_flag_name,
          blosum_flag_description,
          command::Parameter( "blosum_matrix", "a choice of possible matrices")
        )
      );

      // create ShPtrVector of ParameterInterfaces
      util::ShPtrVector< command::ParameterInterface> param_vector( 1, command::Parameter( "p_name", "p_description"));
      // create static flag with list of parameters
      util::ShPtr< command::FlagStatic> pam_flag
      (
        new command::FlagStatic( pam_flag_name, pam_flag_description, param_vector)
      );

      // test copy constructor and Clone()
      util::ShPtr< command::FlagInterface> cloned_flag( blosum_flag->Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_ExampleCheck( default_constructed_flag.GetClassIdentifier(), GetStaticClassName< command::FlagStatic>());

      // test GetName()
      BCL_ExampleCheck( switch_flag.GetName(), switch_flag_name);

      // test SetName()
      BCL_ExampleCheck( default_constructed_flag.GetName().empty(), true);
      default_constructed_flag.SetName( switch_flag_name);
      BCL_ExampleCheck( default_constructed_flag.GetName(), switch_flag_name);

      // test GetDescription()
      BCL_ExampleCheck( switch_flag.GetDescription(), switch_flag_description);

      // the switch flag should not yet be set
      BCL_ExampleCheck( switch_flag.GetFlag(), false);

      // now set the flag and ensure that it's now set
      switch_flag.SetFlag();
      BCL_ExampleCheck( switch_flag.GetFlag(), true);

      // now unset the flag and ensure that it's now not set
      switch_flag.UnsetFlag();
      BCL_ExampleCheck( switch_flag.GetFlag(), false);

      // test GetParameterList()
      BCL_ExampleCheck( blosum_flag->GetParameterList().GetSize(), 1);
      BCL_ExampleCheck( pam_flag->GetParameterList().GetSize(), 1);

      // test GetFirstParameter()
      BCL_ExampleCheck
      (
        blosum_flag->GetFirstParameter()->GetClassIdentifier(),
        GetStaticClassName< command::Parameter>()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // create a sample parameter for use in testing the ReadFromList function
      const std::string blosum_param_string[ 1] = { "BLOSUM60"};
      storage::Vector< std:: string> blosum_param_vector( 1, blosum_param_string);

      // make sure that it can read in a set of parameters from a vector with no trouble
      BCL_ExampleCheck( blosum_flag->ReadFromList( blosum_param_vector, error_stream), true);

      // test IsValidList()
      BCL_ExampleCheck( blosum_flag->IsValidList(), true);

      // test PushBack()
      pam_flag->PushBack
      (
        util::ShPtr< command::ParameterInterface>( new command::Parameter( "pam_matrix", "a choice of possible matrices"))
      );
      BCL_ExampleCheck( pam_flag->GetParameterList().GetSize(), 2);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      blosum_flag->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // clear string behind the stringstream buffer
      output_stream.str( std::string());
      // test WriteUserCommand() and check if length > 0
      blosum_flag->WriteUserCommand( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( switch_flag);

      // read file back into different parameter
      command::FlagStatic switch_flag_read;
      ReadBCLObject( switch_flag_read);

      // check that written and read flag are the same
      BCL_ExampleCheck
      (
        switch_flag.GetName() == switch_flag_read.GetName()
          && switch_flag.GetDescription() == switch_flag_read.GetDescription()
          && switch_flag.GetFlag() == switch_flag_read.GetFlag()
          && switch_flag.GetSize() == switch_flag_read.GetSize(),
        true
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandFlagStatic

  const ExampleClass::EnumType ExampleCommandFlagStatic::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandFlagStatic())
  );

} // namespace bcl

