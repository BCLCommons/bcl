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
#include "command/bcl_command_command_state.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_runtime_environments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CommandState::s_Instance( GetObjectInstances().AddInstance( new CommandState()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param DRY_RUN true if no flags should have their signal functions called
    CommandState::CommandState( const bool &DRY_RUN) :
      m_AvailableFlags(),
      m_DryRun( DRY_RUN),
      m_CurrentFlagType( e_BclCore)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to CommandState object
    CommandState *CommandState::Clone() const
    {
      return new CommandState( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CommandState::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get all flags that remain in the state
    const storage::Map< std::string, storage::Vector< std::string> > &CommandState::GetState() const
    {
      return m_AvailableFlags;
    }

    //! @brief get access to the global command state
    //! @return the command state object that will be used in int main() const
    CommandState &CommandState::GetGlobalCommandState()
    {
      static CommandState s_global_command_state;
      return s_global_command_state;
    }

    //! @brief get access to a bool that indicates whether this is the static initialization phase
    //! @return bool that indicates hether this is the static initialization phase
    bool &CommandState::IsInStaticInitialization()
    {
      static bool s_static_init( true);
      return s_static_init;
    }

    //! @brief get access to a bool that indicates whether help was requested anywhere on the command line
    //! @return bool that indicates whether help was requested anywhere on the command line
    bool &CommandState::GetWasHelpRequested()
    {
      static bool s_help_requested( false);
      return s_help_requested;
    }

    //! @brief get access to a bool that indicates whether help was given anywhere on the command line
    //! @return bool that indicates whether help was given anywhere on the command line
    bool &CommandState::GetWasHelpGiven()
    {
      static bool s_help_given( false);
      return s_help_given;
    }

    //! @brief get whether main() is currently parsing the command line
    //! @return bool that indicates whether main() is currently parsing the command line
    bool &CommandState::GetInMainCommandLineParsing()
    {
      static bool s_in_main_command_line_parsing( false);
      return s_in_main_command_line_parsing;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Use *this to handle a particular flag, if it was set; and then remove it from the state
    //! @param FLAG the flag to handle
    //! @param ERR_STREAM the stream to use for writing errors to
    //! @return true on success
    bool CommandState::Update( FlagInterface &FLAG, std::ostream &ERR_STREAM)
    {
      bool success( true);
      storage::Map< std::string, storage::Vector< std::string> >::iterator itr( m_AvailableFlags.Find( FLAG.GetName()));
      if( itr != m_AvailableFlags.End())
      {
        HandleDryRunState();
        success = FLAG.ReadFromList( itr->second, ERR_STREAM);
        if( FLAG.IsUnused())
        {
          ERR_STREAM << "Flag " << FLAG.GetName() << " is no longer used and may be removed from a future release. ";
        }
        UnsetDryRunState();
        m_AvailableFlags.RemoveElement( itr);
      }
      else if( FLAG.IsUnused())
      {
        // do nothing; no need to check validity even
      }
      else if( !m_DryRun && FLAG.GetSignal())
      {
        // call the signal function
        ( *FLAG.GetSignal())();
        success = FLAG.IsValidList( ERR_STREAM);
      }
      else
      {
        success = FLAG.IsValidList( ERR_STREAM);
      }
      FlagTypeEnum new_flag_type( GetFlagTypeFromName( FLAG.GetName()));
      if( new_flag_type > m_CurrentFlagType)
      {
        m_CurrentFlagType = new_flag_type;
      }
      return success;
    }

    //! @brief Use *this to handle a particular parameter, if it was set; and then remove it from the state
    //! @param PARAMETER the parameter to handle
    //! @param ERR_STREAM the stream to use for writing errors to
    //! @return true on success
    bool CommandState::Update( ParameterInterface &PARAMETER, std::ostream &ERR_STREAM)
    {
      bool success( true);
      // look for the parameter flag
      storage::Map< std::string, storage::Vector< std::string> >::iterator itr( m_AvailableFlags.Find( ""));

      if( itr == m_AvailableFlags.End() || itr->second.IsEmpty())
      {
        // if neither the parameter was set in command line, nor the parameter has a default
        if( !PARAMETER.GetWasDefaultGiven())
        {
          ERR_STREAM << "parameter " << PARAMETER.GetName() << ":";
          ERR_STREAM << " was not given!" << '\n';
          success = false;
        }
      }
      else
      {
        // set the parameter with the 1st element in the vector
        success = PARAMETER.SetParameter( itr->second.FirstElement(), ERR_STREAM);

        // erase the parameter from the vector
        itr->second.RemoveElements( 0, 1);

        // if the vector is now empty, remove it from the map
        if( itr->second.IsEmpty())
        {
          m_AvailableFlags.RemoveElement( itr);
        }
      }
      return success;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes CommandState to ostream OSTREAM
    //! @param OSTREAM - the outstream to write to
    //! @param INDENT indentation
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &CommandState::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AvailableFlags, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DryRun, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    //! @brief read CommandState for std::istream ISTREAM
    //! @param ISTREAM - the instream to read from
    //! @return std::istream &ISTREAM - return the stream after reading
    std::istream &CommandState::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AvailableFlags, ISTREAM);
      io::Serialize::Read( m_DryRun, ISTREAM);

      // end
      return ISTREAM;
    } // Read

    //! @brief sets the internal argument map from the given command line strings
    //! @param NUMBER_ARGUMENTS argc from the command line
    //! @param ARGUMENTS given arguments
    //! @param ERR_STREAM stream to write errors out to
    //! @return true on success
    bool CommandState::ParseArguments( const int NUMBER_ARGUMENTS, const char **ARGUMENTS, std::ostream &ERR_STREAM)
    {
      return this->ParseArguments( util::StringListFromCharacterArray( NUMBER_ARGUMENTS, ARGUMENTS), ERR_STREAM);
    }

    //! @brief sets the internal argument map from the given command line strings
    //! @param CMD_LINE command line strings
    //! @param ERR_STREAM stream to write errors out to
    //! @return true on success
    bool CommandState::ParseArguments( const storage::Vector< std::string> &CMD_LINE, std::ostream &ERR_STREAM)
    {
      HandleDryRunState();
      const std::pair< storage::Vector< std::string>, bool> argument_list_success
      (
        ExpandResponseFiles( CMD_LINE, ERR_STREAM)
      );
      UnsetDryRunState();

      const storage::Vector< std::string> &argument_list( argument_list_success.first);
      bool success_setting_flags_with_parameters( argument_list_success.second);

      std::vector< std::string>::const_iterator arg_itr( argument_list.Begin());
      const std::vector< std::string>::const_iterator arg_itr_end( argument_list.End());

      // set the name of the first param to be the application name
      storage::Vector< std::string> parameters( GetParametersTillNextFlag( arg_itr, arg_itr_end));
      if( !parameters.IsEmpty())
      {
        m_AvailableFlags[ ""].Append( parameters);
      }

      // iterate over given arguments
      while( arg_itr != arg_itr_end)
      {
        std::string flagname( *arg_itr);

        // remove leading -
        while( !flagname.empty() && flagname[ 0] == '-')
        {
          flagname = flagname.erase( 0, 1);
        }

        if( m_AvailableFlags.Has( flagname))
        {
          success_setting_flags_with_parameters = false;
          ERR_STREAM << "Duplicate flag entry for flag with name: " << flagname << '\n';
        }

        m_AvailableFlags[ flagname] = GetParametersTillNextFlag( ++arg_itr, arg_itr_end);
      }

      m_ParseArgumentsSignal.Emit( *this);
      return success_setting_flags_with_parameters;
    }

    //! @brief get the number of parameters remaining
    //! @return the number parameters that remain to be parsed
    size_t CommandState::GetNumberRemainingParameters() const
    {
      storage::Map< std::string, storage::Vector< std::string> >::const_iterator itr( m_AvailableFlags.Find( ""));
      return itr == m_AvailableFlags.End() ? size_t( 0) : itr->second.GetSize();
    }

    //! @brief get the arguments for a particular flag (empty vector if the flag was not set)
    //! @param FLAG_NAME the flag of interest, should not have the - prefix
    //! @return arguments for the given flag, if any were given
    const storage::Vector< std::string> &CommandState::GetArguments( const std::string &FLAG_NAME) const
    {
      static storage::Vector< std::string> s_empty_args;
      storage::Map< std::string, storage::Vector< std::string> >::const_iterator itr( m_AvailableFlags.Find( FLAG_NAME));
      return itr == m_AvailableFlags.End() ? s_empty_args : itr->second;
    }

    //! @brief test whether a particular flag was given
    //! @param FLAG_NAME the flag of interest, should not have the - prefix
    //! @return true if the flag was given
    bool CommandState::WasFlagGiven( const std::string &FLAG_NAME) const
    {
      return m_AvailableFlags.Has( FLAG_NAME);
    }

    //! @brief get the next flag type to consider
    //! @return the next flag type to handle
    const FlagType &CommandState::GetFlagType() const
    {
      return m_CurrentFlagType;
    }

    //! @brief set next flag type to consider
    //! @param TYPE the next type of flag to consider
    void CommandState::SetFlagType( const FlagType &TYPE)
    {
      m_CurrentFlagType = TYPE;
    }

    //! @brief handle setting the signal for all flags
    void CommandState::HandleDryRunState()
    {
      // turn signals off if dry run was set
      if( m_DryRun)
      {
        FlagInterface::SetNoSignal();
      }
      else // turn signals on if dry run was set to false
      {
        FlagInterface::SetSignal();
      }
    }

    //! @brief handle unsetting the signal for all flags
    void CommandState::UnsetDryRunState()
    {
      // turn signals back on if the dry run was set and has not been overridden by another command state
      if( m_DryRun && !FlagInterface::ShouldSignal())
      {
        FlagInterface::SetSignal();
      }
    }

    //! @brief reads in parameters until the next flag is reached
    //! @param ARG_ITR - an iterator over the arguments
    //! @param ARG_ITR_END - the end of the iterator
    //! @return storage::Vector< std::string> - the parameters
    storage::Vector< std::string> CommandState::GetParametersTillNextFlag
    (
      std::vector< std::string>::const_iterator &ARG_ITR,
      const std::vector< std::string>::const_iterator &ARG_ITR_END
    )
    {
      // initialize list of arguments
      storage::Vector< std::string> parameters;

      // pushback arguments if it is not the end and if they have no leading '-'
      // or are numerical if they have a leading '-'
      while( ARG_ITR != ARG_ITR_END && ( ARG_ITR->operator []( 0) != '-' || util::IsNumerical( *ARG_ITR)))
      {
        // store this argument
        parameters.PushBack( *ARG_ITR);
        // step to the next argument
        ++ARG_ITR;
      }

      // return collected parameters
      return parameters;
    } // GetParametersTillNextFlag

    //! @brief expand response files recursively
    //! @param ARGUMENTS list of arguments
    //! @param ERR_STREAM the stream to use for writing errors to
    //! @return arguments where arguments preceding '@' are treated as files, opened and inserted as a list of arguments
    //!         and additionally a bool indicating whether response files were read successfully
    std::pair< storage::Vector< std::string>, bool> CommandState::ExpandResponseFiles
    (
      const storage::Vector< std::string> &ARGUMENTS,
      std::ostream &ERR_STREAM
    )
    {
      bool found_response_files( false);
      // check whether any of the arguments contain response files
      for
      (
        storage::Vector< std::string>::const_iterator itr( ARGUMENTS.Begin()), itr_end( ARGUMENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->size() && ( *itr)[ 0] == s_ResponseFileChar)
        {
          found_response_files = true;
          break;
        }
      }

      // no response files found, just return
      if( !found_response_files)
      {
        return std::make_pair( ARGUMENTS, true);
      }

      storage::Vector< std::string> expanded;

      // first, look for the runtime environment flag, since it would affect interpretation of file names
      static const std::string s_runtime_env_flag_name
      (
        "-" + util::RuntimeEnvironments::GetFlagRuntimeEnvironment()->GetName()
      );

      util::SiPtr< const storage::Vector< std::string> > arguments_ptr( ARGUMENTS);

      // vector to hold arguments except the runtime environment
      storage::Vector< std::string> arguments_no_runtime_env;

      // success
      bool success( true);

      const size_t runtime_env_position( ARGUMENTS.Find( s_runtime_env_flag_name));
      storage::Vector< std::string>::const_iterator itr_runtime( ARGUMENTS.Begin() + runtime_env_position);
      storage::Vector< std::string>::const_iterator itr_runtime_end( ARGUMENTS.End());
      if( runtime_env_position < ARGUMENTS.GetSize())
      {
        // retrieve the parameters associated with that flag
        storage::Vector< std::string>::const_iterator itr_runtime_copy( itr_runtime);
        storage::Vector< std::string> runtime_env_arguments
        (
          GetParametersTillNextFlag( itr_runtime_copy, itr_runtime_end)
        );

        // handle the runtime environment flag
        success =
          util::RuntimeEnvironments::GetFlagRuntimeEnvironment()->ReadFromList( runtime_env_arguments, ERR_STREAM);

        itr_runtime_end = itr_runtime;
      }

      // iterate over all arguments
      for
      (
        storage::Vector< std::string>::const_iterator itr( ARGUMENTS.Begin()), itr_end( ARGUMENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // skip already-parsed runtime environment
        if( itr == itr_runtime)
        {
          itr = itr_runtime_end;
          --itr;
          continue;
        }
        if( itr->operator []( 0) == s_ResponseFileChar)
        {
          const std::string response_filename( itr->c_str() + 1);
          // read the arguments from the file
          io::IFStream read;
          // test whether the file exists and was opened
          if( !io::File::TryOpenIFStream( read, response_filename))
          {
            ERR_STREAM << "Response file: " << response_filename << " does not exist!\n";
            success = false;
          }
          else
          {
            // add the new response files
            std::pair< storage::Vector< std::string>, bool> arguments_success
            (
              ExpandResponseFiles( util::StringListFromIStream( read), ERR_STREAM)
            );
            expanded.Append( arguments_success.first);
            success = success && arguments_success.second;
            io::File::CloseClearFStream( read);
          }
        }
        else
        {
          // non-response file, just add it
          expanded.PushBack( *itr);
        }
      }

      // end
      return std::make_pair( expanded, success);
    }

  } // namespace command
} // namespace bcl
