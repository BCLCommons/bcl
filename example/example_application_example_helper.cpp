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
#include "example_application_example_helper.h"

// includes from bcl - sorted alphabetically
#include "example.h"
#include "app/bcl_app_apps.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_command_line_writer.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"

// external includes - sorted alphabetically

namespace bcl
{
  // note; the paths in this class do not use PATH_SEPARATOR because file handlers will take care
  // of changing the paths at runtime, and problems have been found using PATH_SEPARATOR explicitly here
  // when cross-compiling for windows and then running on windows
  std::string ApplicationExampleHelper::s_ApplicationExampleFilePath
  (
    "/dors/meilerlab/apps/bcl/app_example_files/"
  );

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

  //! @brief constructor from an application
  ApplicationExampleHelper::ApplicationExampleHelper( const app::ApplicationType &APP) :
    m_Application( APP),
    m_ApplicationInstance( m_Application->HardCopy()),
    m_Command( m_ApplicationInstance->InitializeCommand()),
    m_RunTimer(),
    m_AppDefaultFlags( new command::Command())
  {
    // add default bcl parameters
    command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *m_AppDefaultFlags);
    // save the original app default flags
    InitializeAppDefaultFlagMap();
    m_FlagParameterMap = m_AppDefaultFlagValueMap;

    // write out the default command line options as a normal command line
    ResetFlagsAndParameters();

    // reset the application parameter so that later calls can call it
    app::GetApps().GetApplicationsParameter()->Reset();

    m_RunTimer.Reset();
  }

  //! @brief destructor, resets the app default flags
  ApplicationExampleHelper::~ApplicationExampleHelper()
  {
    ResetDefaultFlags();

    // reset the application parameter so that later calls can call it
    app::GetApps().GetApplicationsParameter()->Reset();

    // reset the application parameter so that later calls can call it
    app::GetApps().GetApplicationsParameter()->SetParameter( "Examples", util::GetLogger());
  }

/////////////////
// data access //
/////////////////

  //! @brief get the application
  //! @return const access to the application that will be run internally
  const util::ShPtr< app::Interface> &ApplicationExampleHelper::GetApplication() const
  {
    return m_ApplicationInstance;
  }

  //! @brief get the command object
  //! @return const access to the command
  const command::Command &ApplicationExampleHelper::GetCommand() const
  {
    return *m_Command;
  }

  //! @brief get the stopwatch for this application example helper
  //! @return const access to the timer, which tells how long the last command took, and the total runtime of the app
  const util::Stopwatch &ApplicationExampleHelper::GetRunTimer() const
  {
    return m_RunTimer;
  }

  //! @brief get the command line
  //! @return a string representation of the current command line
  std::string ApplicationExampleHelper::GetCurrentCommandLine() const
  {
    // get the arguments and have the command line writer write them
    return command::CommandLineWriter::CreateCommandLine( GetArguments());
  }

  //! @brief get the arguments that will be passed to the application the next time it is run
  //! @return the arguments that will be passed to the application the next time it is run
  storage::Vector< std::string> ApplicationExampleHelper::GetArguments() const
  {
    return this->GetArguments( true);
  }

  //! @brief unset a particular flag
  void ApplicationExampleHelper::UnsetFlag( const std::string &FLAG)
  {
    // look for the flag in the currently set flags
    storage::Map< std::string, storage::Vector< std::string> >::iterator itr( m_FlagParameterMap.Find( FLAG));

    // if it was not set, there is nothing to do so return
    if( itr == m_FlagParameterMap.End())
    {
      return;
    }

    // look for the flag in the app default flags that were overwritten
    storage::Set< std::string>::iterator itr_app_default_flag
    (
      m_AppDefaultFlagsSetByExample.Find( FLAG)
    );

    // if the flag was found, reset it to the default value, or the value given on the command line, if applicable
    if( itr_app_default_flag != m_AppDefaultFlagsSetByExample.End())
    {
      m_AppDefaultFlagsSetByExample.RemoveElement( itr_app_default_flag);
      if( m_AppDefaultFlagsFromCommandline.Contains( FLAG))
      {
        itr->second = m_AppDefaultFlagValueMap[ FLAG];
      }
      else
      {
        m_FlagParameterMap.RemoveElement( itr);
      }
    }
    else // flag was not an app default flag, just remove it from the map
    {
      m_FlagParameterMap.RemoveElement( itr);
    }
  }

  //! @brief reset the parameters
  void ApplicationExampleHelper::ResetParameters()
  {
    UnsetFlag( "");
  }

  //! @brief reset the flags and parameters
  void ApplicationExampleHelper::ResetFlagsAndParameters()
  {
    m_FlagParameterMap.Reset();
    m_AppDefaultFlagsSetByExample.Reset();

    // add flags that were changed over the command line
    for
    (
      storage::Set< std::string>::const_iterator
        itr( m_AppDefaultFlagsFromCommandline.Begin()), itr_end( m_AppDefaultFlagsFromCommandline.End());
      itr != itr_end;
      ++itr
    )
    {
      m_FlagParameterMap[ *itr] = m_AppDefaultFlagValueMap[ *itr];
    }
  }

  //! @brief add a parameter to the list of arguments
  //! @return reference to this (so that multiple add parameters can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::AddParameter( const std::string &PARAMETER)
  {
    // check that the given parameter is not actually a flag input file
    if( !PARAMETER.empty() && PARAMETER[ 0] == command::CommandState::s_ResponseFileChar)
    {
      command::CommandState cmd_state( false);
      BCL_Assert
      (
        cmd_state.ParseArguments( storage::Vector< std::string>( size_t( 1), PARAMETER), util::GetLogger()),
        ""
      );

      const storage::Map< std::string, storage::Vector< std::string> > &flag_name_to_params( cmd_state.GetState());
      for
      (
        storage::Map< std::string, storage::Vector< std::string> >::const_iterator
          itr( flag_name_to_params.Begin()), itr_end( flag_name_to_params.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first.empty())
        {
          // parameter
          AddParameters( itr->second);
        }
        else
        {
          // flag
          SetFlag( itr->first, itr->second);
        }
      }
      return *this;
    }
    // add the parameter value to the un-flagged values
    m_FlagParameterMap[ std::string()].PushBack( PARAMETER);
    return *this;
  }

  //! @brief add several parameters, which may be response files, to the list of arguments
  //! @param PARAMETERS the parameters to add
  //! @return reference to this (so that multiple add parameters can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::AddParameters( const storage::Vector< std::string> &PARAMETERS)
  {
    for
    (
      storage::Vector< std::string>::const_iterator itr( PARAMETERS.Begin()), itr_end( PARAMETERS.End());
      itr != itr_end;
      ++itr
    )
    {
      AddParameter( *itr);
    }
    return *this;
  }

  //! @brief add a flag to the list of arguments
  //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
  //! @return reference to this (so that multiple adds can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::SetFlag( const std::string &FLAG_NAME)
  {
    // look for the flag
    storage::Map< std::string, storage::Vector< std::string> >::iterator itr( m_FlagParameterMap.Find( FLAG_NAME));

    // if it wasn't found, add the flag
    if( itr == m_FlagParameterMap.End())
    {
      m_FlagParameterMap[ FLAG_NAME] = storage::Vector< std::string>();

      // note whether this overwrites a default flag that was given a non-default value over the command line
      if( m_AppDefaultFlagValueMap.Has( FLAG_NAME))
      {
        m_AppDefaultFlagsSetByExample.Insert( FLAG_NAME);
      }
    }

    return *this;
  }

  //! @brief add a flag with one value to the list of arguments
  //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
  //! @param FLAG_VALUE the value of the flag to add
  //! @return reference to this (so that multiple adds can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::SetFlag
  (
    const std::string &FLAG_NAME,
    const std::string &FLAG_VALUE
  )
  {
    if( FLAG_NAME.empty())
    {
      ResetParameters();
      return AddParameter( FLAG_VALUE);
    }
    m_FlagParameterMap[ FLAG_NAME] = storage::Vector< std::string>( 1, FLAG_VALUE);
    // note whether this overwrites a default flag that was given a non-default value over the command line
    if( m_AppDefaultFlagValueMap.Has( FLAG_NAME))
    {
      m_AppDefaultFlagsSetByExample.Insert( FLAG_NAME);
    }
    return *this;
  }

  //! @brief add a flag with one value to the list of arguments
  //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
  //! @param FLAG_VALUES the values of the flag to add
  //! @return reference to this (so that multiple adds can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::SetFlag
  (
    const std::string &FLAG_NAME,
    const storage::Vector< std::string> &FLAG_VALUES
  )
  {
    if( FLAG_NAME.empty())
    {
      ResetParameters();
      return AddParameters( FLAG_VALUES);
    }
    m_FlagParameterMap[ FLAG_NAME] = FLAG_VALUES;
    // note whether this overwrites a default flag that was given a non-default value over the command line
    if( m_AppDefaultFlagValueMap.Has( FLAG_NAME))
    {
      m_AppDefaultFlagsSetByExample.Insert( FLAG_NAME);
    }
    return *this;
  }

  //! @brief add a parameter to an existing flag
  //! @param FLAG_NAME the flag name to add the parameter to
  //! @param FLAG_VALUE the parameter to add to the flag
  //! @return reference to this (so that multiple adds can be strung together)
  ApplicationExampleHelper &ApplicationExampleHelper::AddParameterToFlag
  (
    const std::string &FLAG_NAME,
    const std::string &FLAG_VALUE
  )
  {
    m_FlagParameterMap[ FLAG_NAME].PushBack( FLAG_VALUE);
    return *this;
  }

  //! @brief flag to change application example path
  //! @return flag to change application example path
  util::ShPtr< command::FlagInterface> &ApplicationExampleHelper::GetApplicationExamplePathFlag()
  {
    // static example path flag
    static util::ShPtr< command::FlagInterface> s_example_path_flag
    (
      new command::FlagStatic
      (
        "application_example_path",
        "change application_example_path for reading and writing application example files",
        command::Parameter
        (
          "path",
          "relative or absolute path",
          s_ApplicationExampleFilePath
        )
      )
    );

    // end
    return s_example_path_flag;
  }

  //! @brief set the application example file path
  //! @param PATH the new path to application example files
  //! The path that is set must have the same structure as the example files directory
  void ApplicationExampleHelper::SetApplicationExampleFilePath( const std::string &PATH)
  {
    s_ApplicationExampleFilePath = PATH;
  }

////////////////
// operations //
////////////////

  //! @brief test if a command line is valid
  //! @param PRINT_ERRORS whether to print errors to the screen
  //! @return true if the command line is valid for the application
  bool ApplicationExampleHelper::CheckCommandString( const bool &PRINT_ERRORS)
  {
    ResetOverriddenDefaultFlags();

    // recreate the last application
    m_ApplicationInstance = m_Application->HardCopy();

    // print out the command line being tested
    BCL_MessageStd( "Checking command line:\n" + GetCurrentCommandLine());

    // create a std::ostringstream to write errors out to
    std::ostringstream err_out;

    // test the flags with the application
    const bool result( m_Command->ReadArguments( GetArguments( false), err_out, true));

    // let the user know whether the command line succeeded
    if( result)
    {
      BCL_MessageStd( "Command line was valid");
    }
    else
    {
      BCL_MessageStd( "Command line was invalid");
    }

    // if PRINT_ERRORS was set and there were errors, write out the errors
    if( PRINT_ERRORS && !result)
    {
      BCL_MessageCrt( "Command line error: " + err_out.str());
    }

    // reset the flags
    m_Command->ResetFlagsAndParameters();
    ResetDefaultFlags();

    return result;
  }

  //! @brief run a command with an application
  //! @return return value from main of the application
  int ApplicationExampleHelper::RunCommand()
  {
    ResetOverriddenDefaultFlags();

    // recreate the application
    m_ApplicationInstance = m_Application->HardCopy();

    // create a std::ostringstream to write errors out to
    std::ostringstream err_out;

    // test the flags with the application
    const bool command_was_valid( m_Command->ReadArguments( GetArguments( false), err_out, false));

    // we should always test the command separately to ensure that it will run before
    // calling RunCommand, so here just assert that the command was okay
    BCL_Assert( command_was_valid, "Invalid command line: " + GetCurrentCommandLine() + "\nErrors:\n" + err_out.str());

    BCL_MessageStd( "Running command line: " + GetCurrentCommandLine());

    // start the stopwatch
    m_RunTimer.Start();

    // run the application
    const int return_status( m_ApplicationInstance->Main());

    BCL_MessageStd
    (
      "Run completed of " + GetCurrentCommandLine() + ", return status was: " + util::Format()( return_status)
      + "\nCommand line took " + util::Format()( m_RunTimer.GetProcessDuration().GetSecondsFractional()) + " sec to run"
    );

    // stop the stopwatch
    m_RunTimer.Stop();

    // reset the command line
    m_Command->ResetFlagsAndParameters();
    ResetDefaultFlags();

    return return_status;
  }

//////////////////////
// helper functions //
//////////////////////

  //! @brief get the application example input path
  //! @return the application example input path
  const std::string &ApplicationExampleHelper::GetApplicationExampleInputPath()
  {
    return s_ApplicationExampleFilePath;
  }

  //! @brief get the input path for this application example
  //! @return the input path for this application example
  std::string ApplicationExampleHelper::GetThisApplicationExampleInputPath() const
  {
    return s_ApplicationExampleFilePath + m_Application.GetName() + '/';
  }

  //! @brief get the arguments that will be passed to the application the next time it is run
  //! @param INCLUDE_APP_DEFAULTS whether to include flags that in the AppDefaultFlagList
  //! @return the arguments that will be passed to the application the next time it is run
  storage::Vector< std::string> ApplicationExampleHelper::GetArguments( const bool &INCLUDE_DEFAULTS) const
  {
    storage::Vector< std::string> arguments( 1, m_Application.GetName());

    for
    (
      storage::Map< std::string, storage::Vector< std::string> >::const_iterator
        itr_flag_params( m_FlagParameterMap.Begin()), itr_flag_params_end( m_FlagParameterMap.End());
      itr_flag_params != itr_flag_params_end;
      ++itr_flag_params
    )
    {
      // only include AppDefaultFlags if they were set by the example or INCLUDE_DEFAULTS was set
      if
      (
        INCLUDE_DEFAULTS
        || !m_AppDefaultFlagValueMap.Has( itr_flag_params->first)
        || m_AppDefaultFlagsSetByExample.Contains( itr_flag_params->first)
      )
      {
        // add the flag if it exists
        if( !itr_flag_params->first.empty())
        {
          arguments.PushBack( "-" + itr_flag_params->first);
        }

        // add the parameters
        arguments.Append( itr_flag_params->second);
      }
    }
    return arguments;
  }

  //! @brief create the app default flag to value map
  void ApplicationExampleHelper::InitializeAppDefaultFlagMap()
  {
    // construct a map containing the current values of all app default flags
    for
    (
      util::ShPtrVector< command::FlagInterface>::const_iterator
        itr( command::GetAppDefaultFlags().GetAllFlags().Begin()),
        itr_end( command::GetAppDefaultFlags().GetAllFlags().End());
      itr != itr_end;
      ++itr
    )
    {
      const command::FlagInterface &flag( **itr);

      // if the flag was set over the command line, note that by inserting the flag into m_AppDefaultFlagsFromCommandline
      if( flag.GetFlag())
      {
        m_AppDefaultFlagsFromCommandline.Insert( flag.GetName());
        m_AppDefaultFlagValueMap[ flag.GetName()] = storage::Vector< std::string>();
      }

      // if the flag is just either set or not set, don't add it to the map, since that would set it
      if( flag.GetParameterList().IsEmpty())
      {
        continue;
      }

      // get a reference to the arguments for that flag
      m_AppDefaultFlagValueMap[ flag.GetName()] = flag.GetStringList();
    }
  }

  //! @brief reset the default app flags to the values they had at the beginning of the example
  void ApplicationExampleHelper::ResetDefaultFlags()
  {
    for
    (
      util::ShPtrVector< command::FlagInterface>::const_iterator
        itr( m_Command->GetBclFlagsWithParams().Begin()),
        itr_end( m_Command->GetBclFlagsWithParams().End());
      itr != itr_end;
      ++itr
    )
    {
      util::ShPtr< command::FlagInterface> flag( *itr);

      // reset all given flags
      if( flag->GetFlag())
      {
        flag->ResetFlag();
      }
    }

    // reset flags in m_AppDefaultFlagsFromCommandline, then give them their old value from the m_AppDefaultFlagValueMap
    for
    (
      storage::Set< std::string>::const_iterator
        itr( m_AppDefaultFlagsFromCommandline.Begin()), itr_end( m_AppDefaultFlagsFromCommandline.End());
      itr != itr_end;
      ++itr
    )
    {
      util::SiPtr< command::FlagInterface> flag( m_AppDefaultFlags->GetFlagWithParams( *itr));

      // reread the arguments that were given over the command line
      flag->ReadFromList( m_AppDefaultFlagValueMap[ *itr], util::GetLogger());
    }
  }

  //! @brief set AppDefaultFlags that were changed by this example
  //! @param ERR_STREAM stream to write out errors to
  void ApplicationExampleHelper::ResetOverriddenDefaultFlags()
  {
    // reset flags in m_AppDefaultFlagsSetByExample
    for
    (
      storage::Set< std::string>::const_iterator
        itr( m_AppDefaultFlagsSetByExample.Begin()), itr_end( m_AppDefaultFlagsSetByExample.End());
      itr != itr_end;
      ++itr
    )
    {
      m_AppDefaultFlags->GetFlagWithParams( *itr)->ResetFlag();
    }

    // reset the application parameter so that later calls can call it
    app::GetApps().GetApplicationsParameter()->Reset();
  }

} // namespace bcl
