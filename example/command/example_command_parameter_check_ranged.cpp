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
#include "command/bcl_command_parameter_check_ranged.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_ranged.cpp
  //!
  //! @author heinzes1
  //! @date 11/06/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace command
  {
    // explicit template instantiation
    template class ParameterCheckRanged< size_t>;
  } // namespace command

  class ExampleCommandParameterCheckRanged :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckRanged *Clone() const
    {
      return new ExampleCommandParameterCheckRanged( *this);
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

      // test ParameterCheckRanged() default constructor
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_a( new command::ParameterCheckRanged< double>());

      // test ParameterCheckRanged< double>( const double MIN, const double MAX) constructor
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_b
      (
        new command::ParameterCheckRanged< double>( double( 0.0), double( 1.0))
      );

      // test ParameterCheckRanged< size_t>( const size_t MIN, const size_t MAX) constructor
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_c( new command::ParameterCheckRanged< size_t>( 1, 10));

      // test Clone()
      util::ShPtr< command::ParameterCheckInterface> parameter_check_ptr_d( parameter_check_ptr_c->Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier()
      BCL_ExampleCheck
      (
        GetStaticClassName< command::ParameterCheckRanged< double> >(),
        parameter_check_ptr_a->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test IsAllowedParameter() for parameter_check_ptr_a: ranged parameters are umerical values
      BCL_ExampleCheck( parameter_check_ptr_a->IsAllowedParameter( "Nate&Nils", "", error_stream), false);
      // test IsAllowedParameter() for parameter_check_ptr_a: if no range is given all values are allowed
      BCL_ExampleCheck( parameter_check_ptr_a->IsAllowedParameter( "-0.5", "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_a->IsAllowedParameter( "0.5",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_a->IsAllowedParameter( "1.5",  "", error_stream), true);

      // test IsAllowedParameter() for parameter_check_ptr_b: range is [0.0, 1.0]
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "-0.1", "", error_stream), false);
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "0.0",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "0.5",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "1.0",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_b->IsAllowedParameter( "1.1",  "", error_stream), false);

      // test that only parameters within the range [1, 10] are valid
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "0",  "", error_stream), false);
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "1",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "2",  "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "10", "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "11", "", error_stream), false);
      BCL_ExampleCheck( parameter_check_ptr_c->IsAllowedParameter( "EA", "", error_stream), false);

      // test IsAllowedParameter() for parameter_check_ptr_d: Clone() worked correctly
      BCL_ExampleCheck( parameter_check_ptr_d->IsAllowedParameter( "0", "", error_stream), false);
      BCL_ExampleCheck( parameter_check_ptr_d->IsAllowedParameter( "1", "", error_stream), true);
      BCL_ExampleCheck( parameter_check_ptr_d->IsAllowedParameter( "EA", "", error_stream), false);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_check_ptr_a->WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( *parameter_check_ptr_c);

      // read file back into different parameter
      command::ParameterCheckRanged< size_t> parameter_check_c_read;
      ReadBCLObject( parameter_check_c_read);

      // perform test from top with the read object
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "0", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "0", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "1", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "1", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "2", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "2", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "10", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "10", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "11", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "11", "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_check_c_read.IsAllowedParameter( "EA", "", error_stream),
        parameter_check_ptr_c->IsAllowedParameter( "EA", "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckRanged

  const ExampleClass::EnumType ExampleCommandParameterCheckRanged::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckRanged())
  );

} // namespace bcl
