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
#include "command/bcl_command_parameter.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter.cpp
  //!
  //! @author heinzes1
  //! @date 11/06/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameter :
    public ExampleInterface
  {
  public:

    //! @brief the clone function for ExampleCommandParameter
    //! @return a clone of this ExampleCommandParameter
    ExampleCommandParameter *Clone() const
    {
      return new ExampleCommandParameter( *this);
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

    //! @brief the function that is called for this example; tests the functionality of commandline parameters
    //! @return returns 0 upon completion
    int Run() const
    {
      // a stream for use in our test cases
      std::stringstream error_stream;

      // name of parameter_ptr_b
      const std::string parameter_ptr_b_name( "second_parameter");
      // description of parameter_ptr_b
      const std::string parameter_ptr_b_description( "the second test parameter");
      // default alpha numerical parameter as string
      const std::string default_parameter_string( "default parameter");
      // default numerical parameter as string and as double
      const std::string default_parameter_double_string( "3.1415");
      const double default_parameter_double( 3.1415);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct Parameter with no parameters (default constructor)
      util::ShPtr< command::ParameterInterface> parameter_ptr_a( new command::Parameter());

      // construct Parameter with name and description
      util::ShPtr< command::ParameterInterface> parameter_ptr_b
      (
        new command::Parameter( parameter_ptr_b_name, parameter_ptr_b_description)
      );

      // construct Parameter with name, description and parameter check
      util::ShPtr< command::ParameterInterface> parameter_ptr_c
      (
        new command::Parameter( "third_parameter", "the third test parameter", command::ParameterCheckDefault())
      );

      // construct Parameter with name, description, parameter check and default parameter string
      util::ShPtr< command::ParameterInterface> parameter_ptr_d
      (
        new command::Parameter
        (
          "fourth_parameter",
          "the fourth test parameter",
          command::ParameterCheckDefault(),
          default_parameter_string
        )
      );
      // same as before, but construct of type Paramter instead of ParameterInterface
      util::ShPtr< command::Parameter> parameter_ptr_d2
      (
        new command::Parameter
        (
          "fourth_parameter",
          "the fourth test parameter",
          command::ParameterCheckDefault(),
          default_parameter_string
        )
      );

      // construct Parameter with name, description and default parameter string
      util::ShPtr< command::Parameter> parameter_ptr_e
      (
        new command::Parameter( "fifth_parameter", "the fifth test parameter", "default parameter")
      );

      // test Clone()
      util::ShPtr< command::ParameterInterface> parameter_ptr_f( parameter_ptr_a->Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck( GetStaticClassName< command::Parameter>(), parameter_ptr_a->GetClassIdentifier());

      // test GetWasSetInCommandLine() for different constructors
      BCL_ExampleCheck( parameter_ptr_a->GetWasSetInCommandLine(), false);
      BCL_ExampleCheck( parameter_ptr_b->GetWasSetInCommandLine(), false);
      BCL_ExampleCheck( parameter_ptr_c->GetWasSetInCommandLine(), false);
      BCL_ExampleCheck( parameter_ptr_d->GetWasSetInCommandLine(), false);
      BCL_ExampleCheck( parameter_ptr_e->GetWasSetInCommandLine(), false);
      BCL_ExampleCheck( parameter_ptr_f->GetWasSetInCommandLine(), false);

      // test GetWasDefaultGiven() for different constructors
      BCL_ExampleCheck( parameter_ptr_a->GetWasDefaultGiven(), false);
      BCL_ExampleCheck( parameter_ptr_b->GetWasDefaultGiven(), false);
      BCL_ExampleCheck( parameter_ptr_c->GetWasDefaultGiven(), false);
      BCL_ExampleCheck( parameter_ptr_d->GetWasDefaultGiven(), true);
      BCL_ExampleCheck( parameter_ptr_e->GetWasDefaultGiven(), true);
      BCL_ExampleCheck( parameter_ptr_f->GetWasDefaultGiven(), false);

      // add the parameter name, test SetParameter()
      // test GetWasSetInCommandLine() again since it should be true now
      BCL_ExampleCheck( parameter_ptr_a->SetParameter( "first_parameter", error_stream), true);
      BCL_ExampleCheck( parameter_ptr_a->GetWasSetInCommandLine(), true);

      // test GetName()
      BCL_ExampleCheck( parameter_ptr_b->GetName(), parameter_ptr_b_name);

      // test GetDescription()
      BCL_ExampleCheck( parameter_ptr_b->GetDescription(), parameter_ptr_b_description);

      // test GetValue()
      BCL_ExampleCheck( parameter_ptr_d->GetValue(), default_parameter_string);

      // test GetDefaultValue()
      BCL_ExampleCheck( parameter_ptr_d->GetDefaultValue(), default_parameter_string);

      // test SetDefaultParameter()
      parameter_ptr_b->SetDefaultParameter( default_parameter_string);
      BCL_ExampleCheck( parameter_ptr_b->GetWasDefaultGiven(), true);
      BCL_ExampleCheck( parameter_ptr_b->GetValue(), default_parameter_string);

      // test GetNumericalValue()
      BCL_ExampleCheck( parameter_ptr_c->SetParameter( default_parameter_double_string, error_stream), true);
      BCL_ExampleCheck( parameter_ptr_c->GetValue(), default_parameter_double_string);
      BCL_ExampleCheck( parameter_ptr_c->GetNumericalValue< double>(), default_parameter_double);

      // test GetParameterCheck()
      std::string check_class( parameter_ptr_d2->GetParameterCheck().GetClassIdentifier());
      std::string check_class_expected( "bcl::command::ParameterCheckDefault");
      BCL_ExampleCheck( check_class, check_class_expected);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test IsAllowedParameter()
      BCL_ExampleCheck( parameter_ptr_d->IsAllowedParameter( "bob", error_stream), true);

      // test Reset
      parameter_ptr_a->Reset();
      BCL_ExampleCheck( parameter_ptr_a->GetValue().empty(), true);
      BCL_ExampleCheck( parameter_ptr_a->GetWasSetInCommandLine(), false);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_ptr_d->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // clear string behind the stringstream buffer
      output_stream.str( std::string());
      // test WriteUserCommand() and check if length > 0
      parameter_ptr_d->WriteUserCommand( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      WriteBCLObject( *parameter_ptr_d);
      // read file back into different parameter
      command::Parameter parameter_read;
      ReadBCLObject( parameter_read);

      // compare the two parameters
      BCL_ExampleCheck( parameter_ptr_d->GetName(), parameter_read.GetName());
      BCL_ExampleCheck( parameter_ptr_d->GetDescription(), parameter_read.GetDescription());
      BCL_ExampleCheck( parameter_ptr_d->GetValue(), parameter_read.GetValue());
      BCL_ExampleCheck( parameter_ptr_d->GetWasSetInCommandLine(), parameter_read.GetWasSetInCommandLine());
      BCL_ExampleCheck( parameter_ptr_d->GetWasDefaultGiven(), parameter_read.GetWasDefaultGiven());

    //////////////////////
    // helper functions //
    //////////////////////

      // return 0 upon completion
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCommandParameter

  const ExampleClass::EnumType ExampleCommandParameter::s_Instance( GetExamples().AddEnum( ExampleCommandParameter()));

} // namespace bcl
