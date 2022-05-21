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
#include "command/bcl_command_parameter_check_allowed.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_allowed.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckAllowed :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckAllowed *Clone() const
    {
      return new ExampleCommandParameterCheckAllowed( *this);
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

      // test default constructor
      util::ShPtr< command::ParameterCheckAllowed> parameter_allowed_none( new command::ParameterCheckAllowed());

      //test ParameterAllowed( ALLOWED_PARAMETERS) constructor
      const std::string allowed_list[ 3] = { "dog", "cat", "horse"};
      command::ParameterCheckAllowed parameter_allowed
      (
        command::ParameterCheckAllowed( storage::Vector< std::string>( 3, allowed_list))
      );

      // test Clone()
      util::ShPtr< command::ParameterCheckAllowed> parameter_allowed_clone( parameter_allowed.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck( GetStaticClassName< command::ParameterCheckAllowed>(), parameter_allowed.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //bool IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const;
      // make sure that dogs are allowed...
      BCL_ExampleCheck( parameter_allowed_none->IsAllowedParameter( "dog", "", error_stream), false);
      BCL_ExampleCheck( parameter_allowed.IsAllowedParameter( "dog", "", error_stream), true);
      BCL_ExampleCheck( parameter_allowed_clone->IsAllowedParameter( "dog", "", error_stream), true);

      // ...whereas rabbits are not
      BCL_ExampleCheck( parameter_allowed_none->IsAllowedParameter( "rabbit", "", error_stream), false);
      BCL_ExampleCheck( parameter_allowed.IsAllowedParameter( "rabbit", "", error_stream), false);
      BCL_ExampleCheck( parameter_allowed_clone->IsAllowedParameter( "rabbit", "", error_stream), false);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_allowed.WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( parameter_allowed);

      // read file back into different parameter
      command::ParameterCheckAllowed parameter_allowed_read;
      ReadBCLObject( parameter_allowed_read);

      // compare the two parameters
      BCL_ExampleCheck
      (
        parameter_allowed.IsAllowedParameter( "dog", "", error_stream),
        parameter_allowed_read.IsAllowedParameter( "dog", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_allowed.IsAllowedParameter( "rabbit", "", error_stream),
        parameter_allowed_read.IsAllowedParameter( "rabbit", "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckAllowed

  const ExampleClass::EnumType ExampleCommandParameterCheckAllowed::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckAllowed())
  );

} // namespace bcl
