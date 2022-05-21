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

#ifndef BCL_UTIL_MESSAGE_H_
#define BCL_UTIL_MESSAGE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief access to the only instance of the message class
    //! @return reference to MessageObject
    BCL_API Message &GetMessenger();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Message
    //! @brief message managing class for the library
    //! @details This class is the message managing class. It contains various static functions to output different messages.\n
    //! It contains functions and members to set the MessageLevel dynamically during runtime and keeps a history\n
    //! of all previous MessageLevel changes.\n
    //! After calling the SetMessageLevel function you can call the ResetToPreviousMessagelevel function to get\n
    //! the previous in the history.\n
    //! It also provides a function so that you can ask what the current MessageLevel is and whether the MessageLevel\n
    //! is high enough for outputting a Message of a desired MessageLevel, if you are not directly using the Macros\n
    //! defined at the end of this file.
    //! You can change the stream the message gets written to, and using the e_Debug MessageLevel, you are able to\n
    //! receive a detailed Message that contains: the filename, the function name where the Message came from and the\n
    //! line number.
    //!
    //! @see @link example_util_message.cpp @endlink
    //! @author woetzen
    //! @date Nov 10, 2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Message
    {

    /////////////
    // friends //
    /////////////

      friend BCL_API Message &GetMessenger();

    public:

    ///////////
    // enums //
    ///////////

      //! custom type message level
      enum MessageLevel
      {
        e_Error,    //!< level that only sends messages that are errors
        e_Silent,   //!< level that only sends messages that can not be missed
        e_Critical, //!< level that sends silent messages and also messages that are of critical importance
        e_Standard, //!< level that sends critical messages and messages of normal importance
        e_Verbose,  //!< level that sends standard messages and additional - but not necessary - information
        e_Debug,    //!< level that send verbose and all other messages, helping to debug code (with detailed output)
        s_NumberMessageLevelTypes
      };

      //! verbosity (output detail) of messages
      enum MessageVerbosity
      {
        e_Summary, //!< output summary messages
        e_Detail,  //!< output detailed messages
        s_NumberMessageVerbosityTypes
      };

      //! @brief MessageLevel as string
      //! @param MESSAGE_LEVEL the message level
      //! @return the MessageLevel as string
      static const std::string &GetLevelString( const MessageLevel &MESSAGE_LEVEL);

      //! @brief MessageVerbosity as string
      //! @param MESSAGE_VERBOSITY the message verbosity
      //! @return the MessageVerbosity as string
      static const std::string &GetVerbosityString( const MessageVerbosity &MESSAGE_VERBOSITY);

      // wrappers for the above enums; to be used for all classes that use these enums as member data
      typedef WrapperEnum< MessageLevel, &GetLevelString, s_NumberMessageLevelTypes> MessageLevelEnum;
      typedef WrapperEnum< MessageVerbosity, &GetVerbosityString, s_NumberMessageVerbosityTypes> MessageVerbosityEnum;

    //////////
    // data //
    //////////

      //! @brief command line flag to be used to set MessageLevel over the command line
      //! @return ShPtr to a FlagInterface which is used to set MessageLevel
      static const ShPtr< command::FlagInterface> &GetMessageLevelFlag();

      //! @brief command line parameter to be used to set MessageLevel over the command line
      //! @return ShPtr to a parameter which is used to set MessageLevel
      static const ShPtr< command::ParameterInterface> &GetMessageLevelParam();

      //! @brief command line parameter to be used to set MessageVerbosity over the command line
      //! @return ShPtr to a parameter which is used to set MessageVerbosity
      static const ShPtr< command::ParameterInterface> &GetMessageVerbosityParam();

    private:

      MessageLevel               m_CurrentMessageLevel;     //!< current message level
      MessageVerbosity           m_CurrentMessageVerbosity; //!< current message verbosity

      //! @brief 3-char strings for the message level to be used when outputting messages to indicate level
      //! @return array with string for each MessageLevel enum
      static const std::string *GetMessageLevelTags();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief empty constructor - setting defaults
      Message();

      //! @brief undefined copy constructor (no message objects will exist)
      Message( const Message &);

      //! @brief undefined operator=
      Message &operator=( const Message &);

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief store the current MessageLevel
      MessageLevel GetCurrentMessageLevel() const;

      //! @brief message level verbosity controls the amount of output in messages: summary or detail
      MessageVerbosity GetMessageVerbosity() const;

      //! @brief set MessageLevel to a new MESSAGE_LEVEL
      //! @param MESSAGE_LEVEL level that should be the current one
      void SetMessageLevel( const MessageLevel MESSAGE_LEVEL);

      //! @brief set MessageVerbosity to a new MESSAGE_VERBOSITY
      //! @param MESSAGE_VERBOSITY summary or detail flag
      void SetMessageVerbosity( const MessageVerbosity MESSAGE_VERBOSITY);

      //! @brief get the string describing a given message level
      static const std::string &GetEnumDescriptor( const MessageLevel &LEVEL);

      //! @brief get the string describing a given message verbosity
      static const std::string &GetEnumDescriptor( const MessageVerbosity &LEVEL);

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the MessageLevel to what was given in the command line
      void SetMessageLevelFromCommandLineFlag();

      //! @brief sets the MessageVerbosity to what was given in the command line
      void SetMessageVerbosityFromCommandLineFlag();

      //! @brief is argument MessageLevel smaller or equal current MessageLevel
      //! @param MESSAGE_LEVEL the message level that should be compared
      //! @return true if the MESSAGE_LEVEL is smaller or equal than the current MessageLevel
      bool IsSmallerEqualCurrentMessageLevel( const MessageLevel MESSAGE_LEVEL) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief updates the message level and verbosity from the command line flags for the global messenger
      static void UpdateMessengerFromCommandLine();

      //! @brief send a message from a namespace to the user
      static std::ostream &SendToUser
      (
        const std::string &NAMESPACE,
        const std::string &MESSAGE
      );

      //! @brief a message is output
      static std::ostream &SendToUser
      (
        const MessageLevel &MESSAGE_LEVEL,
        const std::string &NAMESPACE,
        const std::string &MESSAGE
      );

      //! @brief a message is output with the current filename and function name and line
      static std::ostream &SendToUser
      (
        const MessageLevel &MESSAGE_LEVEL,
        const std::string &NAMESPACE,
        const std::string &MESSAGE,
        const char *FILE_NAME,
        const int LINE_NUMBER,
        const char *FUNCTION_NAME
      );

    private:

      //! @brief assemble a message from Namespace and Message
      //! @param NAMESPACE the namespace the message comes from
      //! @param MESSAGE the actual message
      //! @return string of the for {namespace}=>{message}
      static std::string AssembleMessage
      (
        const std::string &NAMESPACE,
        const std::string &MESSAGE
      );

      //! @brief assemble a message from MessageLevel, Namespace and Message
      //! @param MESSAGE_LEVEL the level for this message
      //! @param NAMESPACE the namespace the message comes from
      //! @param MESSAGE the actual message
      //! @return string of the for {level}={namespace}=>{message}
      static std::string AssembleMessage
      (
        const MessageLevel &MESSAGE_LEVEL,
        const std::string &NAMESPACE,
        const std::string &MESSAGE
      );

    }; // class Message

  } // namespace util

  #define \
  BCL_Message( BCL_MESSAGE_LEVEL, BCL_MESSAGE_DESCRIPTION) \
    if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( BCL_MESSAGE_LEVEL))\
    {\
      util::GetMessenger().SendToUser( BCL_MESSAGE_LEVEL, GetNamespaceIdentifier(), std::string() + ( BCL_MESSAGE_DESCRIPTION), __FILE__, __LINE__, __FUNCTION__); \
    }\

  // shortcuts for BCL_Messages of various levels
  #define BCL_MessageTop( BCL_MESSAGE_DESCRIPTION) BCL_Message( util::Message::e_Silent, BCL_MESSAGE_DESCRIPTION)
  #define BCL_MessageCrt( BCL_MESSAGE_DESCRIPTION) BCL_Message( util::Message::e_Critical, BCL_MESSAGE_DESCRIPTION)
  #define BCL_MessageStd( BCL_MESSAGE_DESCRIPTION) BCL_Message( util::Message::e_Standard, BCL_MESSAGE_DESCRIPTION)
  #define BCL_MessageVrb( BCL_MESSAGE_DESCRIPTION) BCL_Message( util::Message::e_Verbose, BCL_MESSAGE_DESCRIPTION)
  #define BCL_MessageDbg( BCL_MESSAGE_DESCRIPTION) BCL_Message( util::Message::e_Debug, BCL_MESSAGE_DESCRIPTION)

} // namespace bcl

#endif // BCL_UTIL_MESSAGE_H_
