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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "util/bcl_util_message.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief access to the only instance of the message class
    //! @return reference to MessageObject
    Message &GetMessenger()
    {
      static Message *s_message( new Message());
      return *s_message;
    }

  ///////////
  // enums //
  ///////////

    //! @brief MessageLevel as string
    //! @param MESSAGE_LEVEL the message level
    //! @return the MessageLevel as string
    const std::string &Message::GetLevelString( const Message::MessageLevel &MESSAGE_LEVEL)
    {
      static const std::string s_message_level_strings[] =
      {
        "Error", "Silent", "Critical", "Standard", "Verbose", "Debug",
        GetStaticClassName< MessageLevel>()
      };

      return s_message_level_strings[ MESSAGE_LEVEL];
    }

    //! @brief MessageVerbosity as string
    //! @param MESSAGE_VERBOSITY the message verbosity
    //! @return the MessageVerbosity as string
    const std::string &Message::GetVerbosityString( const Message::MessageVerbosity &MESSAGE_VERBOSITY)
    {
      static const std::string s_message_verbosity_strings[] =
      {
        "Summary", "Detail",
        GetStaticClassName< MessageVerbosity>()
      };

      return s_message_verbosity_strings[ MESSAGE_VERBOSITY];
    }

  //////////
  // data //
  //////////

    //! @brief command line flag to be used to set MessageLevel over the command line
    //! @return ShPtr to a FlagInterface which is used to set MessageLevel
    const ShPtr< command::FlagInterface> &Message::GetMessageLevelFlag()
    {
      static ShPtr< command::FlagInterface> s_flag_message_level;

      // first time function is called
      if( !s_flag_message_level.IsDefined())
      {
        s_flag_message_level = ShPtr< command::FlagInterface>
        (
          new command::FlagStatic
          (
            "message_level",
            "adjust the MessageLevel",
            ShPtrVector< command::ParameterInterface>::Create
            (
              GetMessageLevelParam(),
              GetMessageVerbosityParam()
            ),
            &Message::UpdateMessengerFromCommandLine
          )
        );
      }

      // end
      return s_flag_message_level;
    }

    //! @brief command line parameter to be used to set MessageLevel over the command line
    //! @return ShPtr to a parameter which is used to set MessageLevel
    const ShPtr< command::ParameterInterface> &Message::GetMessageLevelParam()
    {
      static const ShPtr< command::ParameterInterface> s_message_level_param
      (
        new command::Parameter
        (
          "level",
          "minimum level of messages that will be printed",
          command::ParameterCheckSerializable( MessageLevelEnum()),
          GetLevelString( e_Standard)
        )
      );

      // return
      return s_message_level_param;
    }

    //! @brief command line parameter to be used to set MessageVerbosity over the command line
    //! @return ShPtr to a parameter which is used to set MessageVerbosity
    const ShPtr< command::ParameterInterface> &Message::GetMessageVerbosityParam()
    {
      static const ShPtr< command::ParameterInterface> s_message_verbosity_param
      (
        new command::Parameter
        (
          "verbosity",
          "set to Detail to print the source file and line of origination for each message",
          command::ParameterCheckSerializable( MessageVerbosityEnum()),
          GetVerbosityString( e_Summary)
        )
      );

      // return
      return s_message_verbosity_param;
    }

  //////////
  // data //
  //////////

    //! @brief 3-char strings for the message level to be used when outputting messages to indicate level
    //! @return array with string for each MessageLevel enum
    const std::string *Message::GetMessageLevelTags()
    {
      //! @brief message level codes for "tagging" output messages
      static std::string *s_message_level_tags( NULL);
      if( s_message_level_tags == NULL)
      {
        s_message_level_tags = new std::string[ 6];
        s_message_level_tags[ 0] = "err";
        s_message_level_tags[ 1] = "slt";
        s_message_level_tags[ 2] = "crt";
        s_message_level_tags[ 3] = "std";
        s_message_level_tags[ 4] = "vrb";
        s_message_level_tags[ 5] = "dbg";
      }

      return s_message_level_tags;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief empty constructor - setting defaults
    Message::Message() :
      m_CurrentMessageLevel( e_Standard),
      m_CurrentMessageVerbosity( e_Summary)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief store the current MessageLevel
    //! initialize the default message level to e_Standard
    Message::MessageLevel Message::GetCurrentMessageLevel() const
    {
      return m_CurrentMessageLevel;
    }

    //! @brief stores the MessageVerbosity flag
    //! initialize the default message verbosity to e_Summary
    Message::MessageVerbosity Message::GetMessageVerbosity() const
    {
      return m_CurrentMessageVerbosity;
    }

    //! @brief set MessageLevel to a new MESSAGE_LEVEL
    //! @param MESSAGE_LEVEL level that should be the current one
    void Message::SetMessageLevel( const Message::MessageLevel MESSAGE_LEVEL)
    {
      m_CurrentMessageLevel = MESSAGE_LEVEL;
    }

    //! @brief set MessageVerbosity to a new MESSAGE_VERBOSITY
    //! @param MESSAGE_VERBOSITY level of message output: summary or detail
    void Message::SetMessageVerbosity( const Message::MessageVerbosity MESSAGE_VERBOSITY)
    {
      m_CurrentMessageVerbosity = MESSAGE_VERBOSITY;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief sets the MessageLevel to what was given in the command line
    void Message::SetMessageLevelFromCommandLineFlag()
    {
      // level from command line
      const MessageLevelEnum new_level( GetMessageLevelParam()->GetValue());

      // check if it was changed
      if( m_CurrentMessageLevel != new_level)
      {
        SetMessageLevel( new_level);
      }
    }

    //! @brief sets the MessageLevel to what was given in the command line
    void Message::SetMessageVerbosityFromCommandLineFlag()
    {
      // verbosity from command line
      const MessageVerbosityEnum new_verbosity( GetMessageVerbosityParam()->GetValue());

      // check if it was changed
      if( m_CurrentMessageVerbosity != new_verbosity)
      {
        SetMessageVerbosity( new_verbosity);
      }
    }

    //! @brief is argument MessageLevel smaller or equal current MessageLevel
    //! @param MESSAGE_LEVEL the messagelevel that should be compared
    //! @return true if the MESSAGE_LEVEL is smaller or equal than the current MessageLevel
    bool Message::IsSmallerEqualCurrentMessageLevel( const MessageLevel MESSAGE_LEVEL) const
    {
      return MESSAGE_LEVEL <= m_CurrentMessageLevel;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief updates the message level and verbosity from the command line flags for the global messenger
    void Message::UpdateMessengerFromCommandLine()
    {
      GetMessenger().SetMessageLevelFromCommandLineFlag();
      GetMessenger().SetMessageVerbosityFromCommandLineFlag();
    }

    //! @brief send a message from a namespace to the user
    //! @param NAMESPACE the name of the current namespace
    //! @param MESSAGE the user message
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser( const std::string &NAMESPACE, const std::string &MESSAGE)
    {
      // write namespace, "> " and the message
      GetLogger().LogMessage( NAMESPACE + "> " + MESSAGE);
      return GetLogger();
    }

    //! @brief a message is outputted
    //! @param MESSAGE_LEVEL is the level of severity of the message
    //! @param NAMESPACE the namespace
    //! @param MESSAGE the user Message
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      // write arrow and message
      GetLogger().LogMessage( AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE));

      return GetLogger();
    }

    //! @brief a message is output with the current filename and functionname and line
    //! @param MESSAGE_LEVEL is the level of severity of the message
    //! @param NAMESPACE the namespace
    //! @param MESSAGE the user Message
    //! @param FILE_NAME the name of the file, where the Message came from
    //! @param LINE_NUMBER the line number the message came from
    //! @param FUNCTION_NAME the name of the function where the message came from
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE,
      const char *FILE_NAME,
      const int LINE_NUMBER,
      const char *FUNCTION_NAME
    )
    {
      if( GetMessenger().GetMessageVerbosity() == Message::e_Detail)
      {
        // remove path from filename
        std::string file_name( FILE_NAME);
        file_name = file_name.substr( file_name.find_last_of( PATH_SEPARATOR) + 1);
        // output message
        GetLogger().LogMessage
        (
          AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE) +
          ", function: " + FUNCTION_NAME + " | file: " + file_name + " | line: " + Format()( LINE_NUMBER)
        );
      }
      else
      {
        GetLogger().LogMessage( AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE));
      }

      return GetLogger();
    }

    //! @brief assemble a message from Namespace and Message
    //! @param NAMESPACE the namespace the message comes from
    //! @param MESSAGE the actual message
    //! @return string of the for {namespace}=>{message}
    std::string Message::AssembleMessage
    (
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      return std::string( NAMESPACE + "=> " + MESSAGE);
    }

    //! @brief assemble a message from MessageLevel and Message
    //! @param MESSAGE_LEVEL the level for this message
    //! @param NAMESPACE the namespace the message comes from
    //! @param MESSAGE the actual message
    //! @return string of the for {level}={namespace}=>{message}
    std::string Message::AssembleMessage
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      return std::string( GetMessageLevelTags()[ MESSAGE_LEVEL] + "=" + AssembleMessage( NAMESPACE, MESSAGE));
    }

  } // namespace util
} // namespace bcl
