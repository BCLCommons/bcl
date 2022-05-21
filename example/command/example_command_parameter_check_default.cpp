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
#include "command/bcl_command_parameter_check_default.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_default.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckDefault :
    public ExampleInterface
  {
  public:

    //! @brief the clone function for ExampleCommandParameterCheck
    //! @return a clone of this ExampleCommandParameterCheck
    ExampleCommandParameterCheckDefault *Clone() const
    {
      return new ExampleCommandParameterCheckDefault( *this);
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

    //! @brief the function called for this example class; tests the functionality of command line parameter checks
    //! @return returns 0 upon completion
    int Run() const
    {
      // a stream for use in our test cases
      std::stringstream error_stream;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test ParameterCheck() default constructor
      command::ParameterCheckDefault parameter_check_a;

      // test Clone()
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_b( parameter_check_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck( GetStaticClassName< command::ParameterCheckDefault>(), parameter_check_a.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test IsAllowedParameter() for parameter_check_ptr_a: everything allowed
      BCL_ExampleCheck( parameter_check_a.IsAllowedParameter( "Nate&Nils", "", error_stream), true);

      // test IsAllowedParameter() for parameter_check_ptr_b: everything allowed for cloned
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "Nate&Nils", "", error_stream), true);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;

      // check if it can be written and read from file stream
      WriteBCLObject( parameter_check_a);

      // read file back into different parameter
      command::ParameterCheckDefault parameter_check_a_read;
      ReadBCLObject( parameter_check_a_read);

      // perform test from top with the read object
      BCL_ExampleCheck
      (
        parameter_check_a.IsAllowedParameter( "Nate&Nils", "", error_stream),
        parameter_check_a_read.IsAllowedParameter( "Nate&Nils", "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // return 0 upon completion
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCommandParameterCheck

  const ExampleClass::EnumType ExampleCommandParameterCheckDefault::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckDefault())
  );

} // namespace bcl
