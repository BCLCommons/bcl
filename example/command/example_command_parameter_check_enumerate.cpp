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
#include "command/bcl_command_parameter_check_enumerate.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_enumerate.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckEnumerate :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckEnumerate *Clone() const
    {
      return new ExampleCommandParameterCheckEnumerate( *this);
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
      // a stream for use in our test cases
      std::stringstream error_stream;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      command::ParameterCheckEnumerate< ExampleClass> parameter_check_a;

      // test Clone()
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_b( parameter_check_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck
      (
        GetStaticClassName< command::ParameterCheckEnumerate< ExampleClass> >(),
        parameter_check_a.GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test IsAllowedParameter()
      const std::string allowed_parameter( ExampleClass::RemoveUnnecessaryString( this->GetClassIdentifier()));
      const std::string disallowed_parameter( "spaces are not allowed");
      BCL_ExampleCheck( parameter_check_a.IsAllowedParameter( allowed_parameter, "", error_stream), true);
      BCL_ExampleCheck( parameter_check_a.IsAllowedParameter( disallowed_parameter, "", error_stream), false);
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( allowed_parameter, "", error_stream), true);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_check_ptr_b->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( parameter_check_a);

      // read file back into different parameter
      command::ParameterCheckEnumerate< ExampleClass> parameter_check_a_read;
      ReadBCLObject( parameter_check_a_read);

      // perform test from top with the read object
      BCL_ExampleCheck( parameter_check_a_read.IsAllowedParameter( allowed_parameter, "", error_stream), true);
      BCL_ExampleCheck( parameter_check_a_read.IsAllowedParameter( disallowed_parameter, "", error_stream), false);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleParameterCheckEnumerate

  const ExampleClass::EnumType ExampleCommandParameterCheckEnumerate::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckEnumerate())
  );

} // namespace bcl
