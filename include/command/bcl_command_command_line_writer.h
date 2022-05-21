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

#ifndef BCL_COMMAND_COMMAND_LINE_WRITER_H_
#define BCL_COMMAND_COMMAND_LINE_WRITER_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CommandLineWriter
    //! @brief object that writes the command line under different environments
    //!
    //! @see @link example_command_command_line_writer.cpp @endlink
    //! @author mendenjl
    //! @date Jul 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CommandLineWriter
    {

    private:

    //////////
    // enum //
    //////////

      enum Environment
      {
        e_UnixTcsh,              //!< Unix tcsh/csh type shell environment
        e_UnixBash,              //!< Bash/sh type shell environment
        e_WindowsPowershell,     //!< Windows powershell
        e_WindowsCommandPrompt,  //!< Windows command prompt
        s_NumberEnvironments
      };

    //////////
    // data //
    //////////

      std::string m_Names[ s_NumberEnvironments];                  //!< name of each environment
      char        m_EscapeCharacters[ s_NumberEnvironments];       //!< escape character for each environment
      std::string m_SpecialCharacters[ s_NumberEnvironments];      //!< special characters for each environment

      //! characters that must be escaped even when in quotes, for each environment
      std::string m_ExtraSpecialChars[ s_NumberEnvironments];

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, private because this is a singleton class
      CommandLineWriter();

      //! @brief undefined copy constructor because this is a singleton class
      CommandLineWriter( const CommandLineWriter &);

      //! @brief undefined assignment operator because this is a singleton class
      CommandLineWriter &operator =( const CommandLineWriter &);

      //! @brief get the single instance of this class
      //! @return the single instance of this class
      static const CommandLineWriter &GetCommandLineWriter();

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief write a command line with properly-escaped arguments for known shell environments
      //! @param ARGUMENTS the arguments that were passed into the command line
      //! @return the command line with properly-escaped arguments
      static std::string CreateCommandLine( const storage::Vector< std::string> &ARGUMENTS);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the name of an environment
      //! @param ENVIRONMENT the environment for which a name is desired
      //! @return the name for that environment
      static const std::string &GetName( const Environment &ENVIRONMENT);

      //! @brief formats an argument given to main such that it can be reused on the command line
      //! @param ARG the argument to be escaped
      //! @param ENVIRONMENT the environment to escape the command with
      //! @return the argument reformatted with quotes and escaped characters if necessary for the target environment
      std::string EscapeArgument( const std::string &ARG, const Environment &ENVIRONMENT) const;

      //! @brief write the command line to the screen with properly-escapped arguments
      //! @param ARGUMENTS the arguments that were passed into the command line
      //! @param ENVIRONMENT the environment to write the command line in
      //! @return string containing the command line escaped properly for the target environment
      std::string WriteCommandLine
      (
        const storage::Vector< std::string> &ARGUMENTS,
        const Environment &ENVIRONMENT
      ) const;

    }; //class CommandLineWriter

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_COMMAND_LINE_WRITER_H_
