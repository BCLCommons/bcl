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
#include "command/bcl_command_flag_separator.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_flag_separator.cpp
  //!
  //! @author karakam
  //! @date Nov 23, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandFlagSeparator :
    public ExampleInterface
  {
  public:

    ExampleCommandFlagSeparator *Clone() const
    {
      return new ExampleCommandFlagSeparator( *this);
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

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // initialize strings
      const std::string string_a( "TEST");
      const std::string empty_string;

      // default constructor
      command::FlagSeparator flag_default;
      BCL_ExampleCheck( flag_default.GetSeparatorText(), empty_string);

      // construct from string
      command::FlagSeparator flag( string_a);
      const command::FlagSeparator const_flag( string_a);
      BCL_ExampleCheck( flag.GetSeparatorText(), string_a);

      // clone constructor
      util::ShPtr< command::FlagSeparator> sp_flag( flag.Clone());
      BCL_ExampleCheck( sp_flag->GetSeparatorText(), flag.GetSeparatorText());

    /////////////////
    // data access //
    /////////////////

      // check static class name
      BCL_ExampleCheck( GetStaticClassName< command::FlagSeparator>(), flag.GetClassIdentifier());

      // check GetSeparatorText()
      BCL_ExampleCheck( flag.GetSeparatorText(), string_a);

      // check GetName()
      BCL_ExampleCheck( flag.GetName(), empty_string);

      // check GetDescription()
      BCL_ExampleCheck( flag.GetDescription(), empty_string);

      // check GetFlag()
      BCL_ExampleCheck( flag.GetFlag(), true);

      // check GetParameterList()
      BCL_ExampleCheck( flag.GetParameterList().GetSize(), 0);

      // check const GetParameterList()
      BCL_ExampleCheck( const_flag.GetParameterList().GetSize(), 0);

    ////////////////
    // operations //
    ////////////////

      // check ReadFromList()
      std::stringstream sstream;
      const storage::Vector< std::string> parameter_list;
      BCL_ExampleCheck( flag.ReadFromList( parameter_list, sstream), true);
      BCL_ExampleCheck( sstream.str(), empty_string);
      sstream.clear();

      // check IsValidList()
      BCL_ExampleCheck( flag.IsValidList( sstream), true);
      BCL_ExampleCheck( sstream.str(), empty_string);
      sstream.clear();

    //////////////////////
    // input and output //
    //////////////////////

      // initialize expected output string
      const std::string expected_out( '\n' + string_a + '\n');

      // check WriteHelp()
      flag.WriteHelp( sstream);
      BCL_ExampleCheck( sstream.str(), expected_out);
      sstream.clear();

      // check GetSeparatorText()
      flag.WriteUserCommand( sstream);
      BCL_ExampleCheck( sstream.str(), expected_out);

      // write object
      WriteBCLObject( flag);

      // read object
      command::FlagSeparator flag_read;
      ReadBCLObject( flag_read);
      BCL_ExampleCheck( flag_read.GetSeparatorText(), flag.GetSeparatorText());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandFlagSeparator

  const ExampleClass::EnumType ExampleCommandFlagSeparator::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandFlagSeparator())
  );
  
} // namespace bcl
