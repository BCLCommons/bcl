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
#include "command/bcl_command_command_line_writer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <iomanip>
#include <map>

namespace bcl
{
  namespace command
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, private because this is a singleton class
    CommandLineWriter::CommandLineWriter()
    {
      // store the names of each shell-environment
      m_Names[ size_t( e_UnixTcsh)]             = "TCSH / CSH";
      m_Names[ size_t( e_UnixBash)]             = "BASH / SH";
      m_Names[ size_t( e_WindowsPowershell)]    = "Windows Powershell";
      m_Names[ size_t( e_WindowsCommandPrompt)] = "Windows Command Prompt";

      // store the escape character for each shell/environment
      m_EscapeCharacters[ size_t( e_UnixTcsh)] = '\\';
      m_EscapeCharacters[ size_t( e_UnixBash)] = '\\';
      m_EscapeCharacters[ size_t( e_WindowsPowershell)] = '`';
      m_EscapeCharacters[ size_t( e_WindowsCommandPrompt)] = '^';

      // characters that require escape even when inside quotes for each shell
      m_ExtraSpecialChars[ size_t( e_UnixTcsh)] = "!"; // ! = last command
      m_ExtraSpecialChars[ size_t( e_UnixBash)] = "!\\'";
      m_ExtraSpecialChars[ size_t( e_WindowsPowershell)] = "`'";
      m_ExtraSpecialChars[ size_t( e_WindowsCommandPrompt)] = "^'";

      // string containing characters that must be escaped for all shells
      const std::string common_special_characters( "&()*;|<>\" \t");

      // initialize all special characters to common special characers + the extra special characters for that shell
      for( size_t environment_id( 0); environment_id < size_t( s_NumberEnvironments); ++environment_id)
      {
        m_SpecialCharacters[ environment_id] = common_special_characters + m_ExtraSpecialChars[ environment_id];
      }

      // add additional characters that have special meaning for tcsh/csh and must be escaped if outside quotes
      // ` is used to create internal commands
      // ~ expands to $HOME
      // ? means the previous command executable
      // \\ does not need escaping inside quotes, but so it is not an extra special char for tcsh, but it does need escaping outside quotes
      // ' cannot be escaped inside quotes for tcsh, but it can be outside of quotes, so it belongs here
      m_SpecialCharacters[ size_t( e_UnixTcsh)] += "`~?'\\";

      // add additional characters that have special meaning for bash/sh and must be escaped if outside quotes
      // ` is used to create internal commands
      // ~ expands to $HOME
      // # is a comment
      m_SpecialCharacters[ size_t( e_UnixBash)] += "`~#";

      // : is a special character on windows used to denote drive or volume, however, it is not interpreted directly by
      //   the shell, so no need to escape it
      // ` is also special in powershell, because it serves as the escape character, which is handled explicitly,
      //   so it does not need to go in the special characters section
      m_SpecialCharacters[ size_t( e_WindowsPowershell)] += "#";

      // : and ` are special command prompt for the same reasons as powershell, but on the command prompt the escape
      // character is neither, so both need to be added to the special characters
      m_SpecialCharacters[ size_t( e_WindowsCommandPrompt)] += "#`";

    }

    //! @brief get the single instance of this class
    //! @return the single instance of this class
    const CommandLineWriter &CommandLineWriter::GetCommandLineWriter()
    {
      // single instance of this class
      static const CommandLineWriter s_Instance;
      return s_Instance;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief write a command line with properly-escaped arguments for known shell environments
    //! @param ARGUMENTS the arguments that were passed into the command line
    //! @return the command line with properly-escaped arguments
    std::string CommandLineWriter::CreateCommandLine( const storage::Vector< std::string> &ARGUMENTS)
    {
      // create a map from escaped command line to environments that that command line is valid in
      std::map< std::string, std::string> command_line_to_env_string;

      // create the command strings from each environment
      for( size_t environment_id( 0); environment_id < size_t( s_NumberEnvironments); ++environment_id)
      {
        // get the command line for this environment
        const std::string env_cmd_line
        (
          GetCommandLineWriter().WriteCommandLine( ARGUMENTS, Environment( environment_id))
        );

        // see if the map already contains this command line
        std::map< std::string, std::string>::iterator itr( command_line_to_env_string.find( env_cmd_line));

        // if the map already contained the command line string, just append the environment name to the string
        if( itr != command_line_to_env_string.end())
        {
          itr->second += " / " + GetName( Environment( environment_id));
        }
        else // otherwise, insert the new command string, keyed to the name of the environment
        {
          command_line_to_env_string[ env_cmd_line] = GetName( Environment( environment_id));
        }
      }

      std::ostringstream output;

      // if there was only one command line, just write it out
      if( command_line_to_env_string.size() == size_t( 1))
      {
        output << command_line_to_env_string.begin()->first << '\n';
      }
      else
      {
        // we have to write the command line and the environments that it is valid for
        // so that the environments are always shown in the same order, remap the arguments by environment
        std::map< std::string, std::string> env_to_cmd_line;

        // also keep track of the longest environment name
        size_t max_env_name_length( 0);
        for
        (
          std::map< std::string, std::string>::const_iterator
            itr( command_line_to_env_string.begin()), itr_end( command_line_to_env_string.end());
          itr != itr_end;
          ++itr
        )
        {
          env_to_cmd_line[ itr->second] = itr->first;
          max_env_name_length = std::max( max_env_name_length, itr->second.size());
        }

        // now write out the environment and associated command line
        for
        (
          std::map< std::string, std::string>::const_iterator
            itr( env_to_cmd_line.begin()), itr_end( env_to_cmd_line.end());
          itr != itr_end;
          ++itr
        )
        {
          output << std::setw( max_env_name_length) << itr->first << " : " << itr->second << '\n';
        }
      }
      return output.str();
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the name of an environment
    //! @param ENVIRONMENT the environment for which a name is desired
    //! @return the name for that environment
    const std::string &CommandLineWriter::GetName( const Environment &ENVIRONMENT)
    {
      BCL_Assert( ENVIRONMENT != s_NumberEnvironments, "Invalid environment");
      return GetCommandLineWriter().m_Names[ size_t( ENVIRONMENT)];
    }

    //! @brief formats an argument given to main such that it can be reused on the command line
    //! @param ARG the argument to be escaped
    //! @param ENVIRONMENT the environment to escape the command with
    //! @return the argument reformatted with quotes and escaped characters if necessary for the target environment
    std::string CommandLineWriter::EscapeArgument( const std::string &ARG, const Environment &ENVIRONMENT) const
    {
      // determine whether any characters in the string are special; if so, the string will be quoted and/or special
      // characters will be escaped
      const bool has_special_chars
      (
        ARG.find_first_of( m_SpecialCharacters[ size_t( ENVIRONMENT)]) != std::string::npos
      );
      const bool has_non_extra_special_chars
      (
        ARG.find_first_not_of( m_ExtraSpecialChars[ size_t( ENVIRONMENT)]) != std::string::npos
      );

      // check whether the string should be quoted
      // this requires that the string contains no ' characters or that
      // the environment considers ' an extra special character (which can be escaped inside a string)
      bool should_be_quoted
      (
        has_special_chars && has_non_extra_special_chars &&
        (
          ARG.find_first_of( '\'') == std::string::npos
          || m_ExtraSpecialChars[ size_t( ENVIRONMENT)].find( '\'') != std::string::npos
        )
      );

      // if the argument is empty, it also needs to be quoted
      should_be_quoted = should_be_quoted || ARG.empty();

      // determine whether any characters will need to be escaped
      // generally this is true for any extra special char, as well as all special chars if the string cannot be quoted
      bool needs_escaped_characters
      (
        ARG.find_first_of( m_ExtraSpecialChars[ size_t( ENVIRONMENT)]) != std::string::npos
        ||
        ( has_special_chars && !should_be_quoted)
      );

      // if the string has won't be quoted and has no characters requiring escapes, just return the string
      if( !should_be_quoted && !needs_escaped_characters)
      {
        return ARG;
      }

      // shell_arg will contain the reformatted string
      std::string shell_arg;
      // reserve enough characters for the new string
      shell_arg.reserve( ( needs_escaped_characters ? 2 : 1) * ARG.size() + ( should_be_quoted ? 2 : 0));

      // if the string has characters that are special for the shell, quote the string
      // this is much faster than doing '\'' + shell_arg + '\''
      if( should_be_quoted)
      {
        shell_arg += '\'';
      }

      // it is rare to have escape characters or single quotes inside arguments, so check whether it is necessary first
      if( needs_escaped_characters)
      {
        // get the escape character for this environment
        const char escape_char( m_EscapeCharacters[ size_t( ENVIRONMENT)]);

        // make a reference to the characters that will be escaped
        const std::string &chars_to_escape
        (
          should_be_quoted
          ? m_ExtraSpecialChars[ size_t( ENVIRONMENT)]
          : m_SpecialCharacters[ size_t( ENVIRONMENT)]
        );

        // escape the characters that require escaping while building up the shell_arg
        for( int i( 0), n( ARG.size()); i < n; ++i)
        {
          const char c( ARG[ i]);
          if( chars_to_escape.find( c) != std::string::npos)
          {
            shell_arg += escape_char;
          }
          shell_arg += c;
        }
      }
      else
      {
        // add ARG to the shell argument, no characters inside the string have special meaning inside single quotes
        shell_arg += ARG;
      }

      // if the string has characters that are special for the shell, quote it
      if( should_be_quoted)
      {
        shell_arg += '\'';
      }
      return shell_arg;
    }

    //! @brief write the command line to the screen with properly-escapped arguments
    //! @param ARGUMENTS the arguments that were passed into the command line
    //! @param ENVIRONMENT the environment to write the command line in
    //! @return string containing the command line escaped properly for the target environment
    std::string CommandLineWriter::WriteCommandLine
    (
      const storage::Vector< std::string> &ARGUMENTS,
      const Environment &ENVIRONMENT
    ) const
    {
      // create a string to write the contents of the command line to
      std::string output;

      // write out each argument properly escaped
      for
      (
        storage::Vector< std::string>::const_iterator itr_arg( ARGUMENTS.Begin()), itr_arg_end( ARGUMENTS.End());
        itr_arg != itr_arg_end;
        ++itr_arg
      )
      {
        output += EscapeArgument( *itr_arg, ENVIRONMENT);
        output += ' ';
      }

      return output;
    }

  } // namespace command
} // namespace bcl
