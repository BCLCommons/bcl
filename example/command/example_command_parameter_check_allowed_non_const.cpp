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
#include "command/bcl_command_parameter_check_allowed_non_const.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_allowed_non_const.cpp
  //!
  //! @author mendenjl
  //! @date Jan 18, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckAllowedNonConst :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckAllowedNonConst *Clone() const
    {
      return new ExampleCommandParameterCheckAllowedNonConst( *this);
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

      //test ParameterAllowed( ALLOWED_PARAMETERS) constructor
      storage::Vector< std::string> allowed_list( storage::Vector< std::string>::Create( "dog", "cat", "horse"));
      command::ParameterCheckAllowedNonConst parameter_allowed( allowed_list);

      // test Clone()
      util::ShPtr< command::ParameterCheckAllowedNonConst> parameter_allowed_clone( parameter_allowed.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck( GetStaticClassName< command::ParameterCheckAllowedNonConst>(), parameter_allowed.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //bool IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const;
      // make sure that dogs are allowed...
      BCL_ExampleCheck( parameter_allowed.IsAllowedParameter( "dog", "", error_stream), true);
      BCL_ExampleCheck( parameter_allowed_clone->IsAllowedParameter( "dog", "", error_stream), true);

      // ...whereas rabbits are not
      BCL_ExampleCheck( parameter_allowed.IsAllowedParameter( "rabbit", "", error_stream), false);
      BCL_ExampleCheck( parameter_allowed_clone->IsAllowedParameter( "rabbit", "", error_stream), false);

      // now add rabbit to the vector
      allowed_list.PushBack( "rabbit");

      // and check that it is now allowed
      BCL_ExampleCheck( parameter_allowed.IsAllowedParameter( "rabbit", "", error_stream), true);
      BCL_ExampleCheck( parameter_allowed_clone->IsAllowedParameter( "rabbit", "", error_stream), true);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_allowed.WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckAllowedNonConst

  const ExampleClass::EnumType ExampleCommandParameterCheckAllowedNonConst::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckAllowedNonConst())
  );

} // namespace bcl
