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
#include "command/bcl_command_parameter_check_serializable.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_parameter_check_serializable.cpp
  //!
  //! @author mendenjl
  //! @date Jan 19, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandParameterCheckSerializable :
    public ExampleInterface
  {
  public:

    ExampleCommandParameterCheckSerializable *Clone() const
    {
      return new ExampleCommandParameterCheckSerializable( *this);
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

      // test ParameterCheckSerializable with message verbosity
      util::ShPtr< command::ParameterCheckSerializable> parameter_serializable_verbosity
      (
        new command::ParameterCheckSerializable( util::Message::MessageVerbosityEnum())
      );

      //test ParameterCheckSerializable() with message level
      command::ParameterCheckSerializable parameter_serializable_msg =
        command::ParameterCheckSerializable( util::Message::MessageLevelEnum());

      // test Clone()
      util::ShPtr< command::ParameterCheckSerializable>
        parameter_serializable_msg_clone( parameter_serializable_msg.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //bool IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const;
      // make sure that the standard messsage level is valid for the
      const std::string message_level_std( util::Message::GetLevelString( util::Message::e_Standard));
      BCL_ExampleCheck( parameter_serializable_verbosity->IsAllowedParameter( message_level_std, "", error_stream), false);
      BCL_ExampleCheck( parameter_serializable_msg.IsAllowedParameter( message_level_std, "", error_stream), true);
      BCL_ExampleCheck( parameter_serializable_msg_clone->IsAllowedParameter( message_level_std, "", error_stream), true);

      // ...whereas message verbosity (detail) is not for message levels, but is for message verbosity
      const std::string detail( util::Message::GetVerbosityString( util::Message::e_Detail));
      BCL_ExampleCheck( parameter_serializable_verbosity->IsAllowedParameter( detail, "", error_stream), true);
      BCL_ExampleCheck( parameter_serializable_msg.IsAllowedParameter( detail, "", error_stream), false);
      BCL_ExampleCheck( parameter_serializable_msg_clone->IsAllowedParameter( detail, "", error_stream), false);

    //////////////////////
    // input and output //
    //////////////////////

      // stream to write output into to test
      std::stringstream output_stream;
      // test WriteHelp() and check if length > 0
      parameter_serializable_msg.WriteHelp( output_stream);
      BCL_ExampleCheck( output_stream.str().length() > 0, true);

      // check if it can be written and read from file stream
      WriteBCLObject( parameter_serializable_msg);

      // read file back into different parameter
      command::ParameterCheckSerializable parameter_serializable_read;
      ReadBCLObject( parameter_serializable_read);

      // compare the two parameters
      BCL_ExampleCheck
      (
        parameter_serializable_msg.IsAllowedParameter( message_level_std, "", error_stream),
        parameter_serializable_read.IsAllowedParameter( message_level_std, "", error_stream)
      );
      BCL_ExampleCheck
      (
        parameter_serializable_msg.IsAllowedParameter( message_level_std, "", error_stream),
        parameter_serializable_read.IsAllowedParameter( message_level_std, "", error_stream)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandParameterCheckSerializable

  const ExampleClass::EnumType ExampleCommandParameterCheckSerializable::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandParameterCheckSerializable())
  );

} // namespace bcl
