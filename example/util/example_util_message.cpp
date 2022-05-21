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
#include "util/bcl_util_message.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_message.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilMessage :
    public ExampleInterface
  {
  public:

    ExampleUtilMessage *Clone() const
    { return new ExampleUtilMessage( *this);}

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
      // set MessageLevel to the one given in command line
      util::Message::UpdateMessengerFromCommandLine();

      BCL_MessageStd
      (
        "this is the current util::MessageLevel " +
        util::Message::GetLevelString( util::GetMessenger().GetCurrentMessageLevel())
      );

      BCL_MessageStd
      (
        "this is the current util::MessageVerbosity " +
        util::Message::GetVerbosityString( util::GetMessenger().GetMessageVerbosity())
      );

      //just util::Message evaluates current util::MessageLevel to given util::MessageLevel and prints util::Message
      util::GetMessenger().SetMessageLevel( util::Message::e_Debug);
      BCL_MessageDbg( "MessageLevel was set to Debug");
      BCL_MessageVrb( "MessageLevel was set to Debug");
      BCL_MessageStd( "MessageLevel was set to Debug");
      BCL_MessageCrt( "MessageLevel was set to Debug");
      BCL_MessageTop( "MessageLevel was set to Debug");

      BCL_MessageStd( "\n");

      util::GetMessenger().SetMessageLevel( util::Message::e_Verbose);
      BCL_MessageDbg( "MessageLevel was set to Verbose");
      BCL_MessageVrb( "MessageLevel was set to Verbose");
      BCL_MessageStd( "MessageLevel was set to Verbose");
      BCL_MessageCrt( "MessageLevel was set to Verbose");
      BCL_MessageTop( "MessageLevel was set to Verbose");

      BCL_MessageStd( "\n");

      util::GetMessenger().SetMessageLevel( util::Message::e_Standard);
      BCL_MessageDbg( "MessageLevel was set to Standard");
      BCL_MessageVrb( "MessageLevel was set to Standard");
      BCL_MessageStd( "MessageLevel was set to Standard");
      BCL_MessageCrt( "MessageLevel was set to Standard");
      BCL_MessageTop( "MessageLevel was set to Standard");

      BCL_MessageStd( "\n");

      util::GetMessenger().SetMessageLevel( util::Message::e_Critical);
      BCL_MessageDbg( "MessageLevel was set to Critical");
      BCL_MessageVrb( "MessageLevel was set to Critical");
      BCL_MessageStd( "MessageLevel was set to Critical");
      BCL_MessageCrt( "MessageLevel was set to Critical");
      BCL_MessageTop( "MessageLevel was set to Critical");

      BCL_MessageStd( "\n");

      util::GetMessenger().SetMessageLevel( util::Message::e_Silent);
      BCL_MessageDbg( "MessageLevel was set to Silent");
      BCL_MessageVrb( "MessageLevel was set to Silent");
      BCL_MessageStd( "MessageLevel was set to Silent");
      BCL_MessageCrt( "MessageLevel was set to Silent");
      BCL_MessageTop( "MessageLevel was set to Silent");

      BCL_MessageStd( "\n");

      util::Message::UpdateMessengerFromCommandLine();

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilMessage

  const ExampleClass::EnumType ExampleUtilMessage::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilMessage())
  );

} // namespace bcl
