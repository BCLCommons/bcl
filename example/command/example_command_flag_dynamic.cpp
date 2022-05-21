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
#include "command/bcl_command_flag_dynamic.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_flag_dynamic.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandFlagDynamic :
    public ExampleInterface
  {
  public:

    ExampleCommandFlagDynamic *Clone() const
    {
      return new ExampleCommandFlagDynamic( *this);
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

      // a stream for use in our test cases
      std::stringstream error_stream;

      // flag with dynamic parameter list
      // it consists of the flag name, a description for the flag, a template parameter, that shows how one of the
      // parameters should look like, including an optional parametercheck, a min and a max size
      // as we know, Bob Saget's jokes are not that funny... however he insists on telling at least two but no
      // more than four jokes as shown below

      // example 1: fewer than 2 jokes = bad (well really, good for us since we have to hear fewer of his jokes, but
      // bad for BCL since a minimum of 2 jokes is required!)
      command::FlagDynamic bs_jokes_1
      (
        "bob_sagets_jokes_1", // flag name
        "the titles of bob saget's jokes", // flag description
        command::Parameter // the template parameter
        (
          "bs-joke", // parameter name
          "a joke that bob saget has made", // parameter description
          command::ParameterCheckDefault() // optional Parameter Check
        ),
        2, // min # parameters for this flag
        4  // max # parameters for this flag
      );

      // example 2: just 2 jokes = A-OK!
      util::ShPtr< command::FlagInterface> bs_jokes_2
      (
        new command::FlagDynamic
        (
          "bob_sagets_jokes_2", // flag name
          "the titles of bob saget's jokes", // flag description
          command::Parameter // the template parameter
          (
            "bs-joke",
            "a joke that bob saget has made"
          ),
          2, // min # parameters for this flag
          4  // max # parameters for this flag
        )
      );

      // example 3: just 4 jokes = A-OK!
      util::ShPtr< command::FlagInterface> bs_jokes_3
      (
        new command::FlagDynamic
        (
          "bob_sagets_jokes_3", // flag name
          "the titles of bob saget's jokes", // flag description
          command::Parameter // the template parameter
          (
            "bs-joke",
            "a joke that bob saget has made"
          ),
          2, // min # parameters for this flag
          4  // max # parameters for this flag
        )
      );

      // example 4: 5 bob saget jokes = (really!) bad
      util::ShPtr< command::FlagInterface> bs_jokes_4
      (
        new command::FlagDynamic
        (
          "bob_sagets_jokes_4", // flag name
          "the titles of bob saget's jokes", // flag description
          command::Parameter // the template parameter
          (
            "bs-joke",
            "a joke that bob saget has made"
          ),
          2, // min # parameters for this flag
          4  // max # parameters for this flag
        )
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // create a sample parameter for use in testing the ReadFromList function
      const std::string one_joke_string[ 1] = { "knock-knock"};
      storage::Vector< std:: string> too_short_joke_list( 1, one_joke_string);
      // ReadFromList should return false (i.e. failure) in this case since the vector we gave it only
      // contains 1 joke whereas this flag requires a minimum of 2 jokes
      BCL_ExampleCheck( bs_jokes_1.ReadFromList( too_short_joke_list, error_stream), false);

      // create a sample parameter for use in testing the ReadFromList function
      const std::string two_joke_string[ 2] = { "knock_knock", "blonde"};
      storage::Vector< std:: string> just_right_joke_list_1( 2, two_joke_string);
      // ReadFromList should return true (i.e. success) in this case because we gave it 2 parameters
      // which is within the acceptable range (2 - 4 parameters)
      BCL_ExampleCheck( bs_jokes_2->ReadFromList( just_right_joke_list_1, error_stream), true);

      // create a sample parameter for use in testing the ReadFromList function
      const std::string four_joke_string[ 4] = { "knock_knock", "blonde", "horse_in_bar", "whos_on_first"};
      storage::Vector< std:: string> just_right_joke_list_2( 4, four_joke_string);
      // ReadFromList should return true (i.e. success) in this case because we gave it 4 parameters
      // which is within the acceptable range (2 - 4 parameters)
      BCL_ExampleCheck( bs_jokes_3->ReadFromList( just_right_joke_list_2, error_stream), true);

      // create a sample parameter for use in testing the ReadFromList function
      const std::string five_joke_string[ 5] = { "knock_knock", "blonde", "horse_in_bar", "whos_on_first", "corny"};
      storage::Vector< std:: string> too_long_joke_list( 5, five_joke_string);
      // ReadFromList should return false (i.e. failure) in this case since the vector we gave it
      // contains more than the max allowed # of jokes
      BCL_ExampleCheck( bs_jokes_4->ReadFromList( too_long_joke_list, error_stream), false);

      // print out all the error messages
      BCL_MessageStd( "those are the produced error messages\n" + error_stream.str());

    //////////////////////
    // input and output //
    //////////////////////

      // check if it can be written and read from file stream
      WriteBCLObject( bs_jokes_1);

      // read file back into different parameter
      command::FlagDynamic bs_jokes_read;
      ReadBCLObject( bs_jokes_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandFlagDynamic

  const ExampleClass::EnumType ExampleCommandFlagDynamic::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandFlagDynamic())
  );

} // namespace bcl

