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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_command_line_writer.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_memory_usage.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

// leave this, otherwise BCL_Message cannot be used
// it assumes that GetNamespaceIdentifier is available, which it is not, for non-bcl namespace
using namespace bcl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @file apps.cpp
//! @author woetzen, meilerj
//!
//! @brief this is the main for bcl
//! the first argument determines the app that is started, the rest is passed to the app as the command line
//! @param NUMBER_ARGUMENTS the number of arguments passed from the shell
//! @param ARGUMENTS passed from the shell, the first (0) being the executable name, the second one should be the bcl
//!        app name; if there is only one app available (linked) then there is no need to pass an app name, but it is
//!        still possible
//! @remarks example unnecessary
//! @return the error code, 0 for successful execution
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main( const int NUMBER_ARGUMENTS, const char *ARGUMENTS[])
{
  // update whether this is static initialization time
  command::CommandState::IsInStaticInitialization() = false;

  // number of applications in bcl
  const size_t number_apps_in_bcl( app::GetApps().GetEnumCount());

  // check that some applications were linked
  if( !number_apps_in_bcl)
  {
    BCL_ExitWithoutCallstack( "no applications linked", -1);
    return 0;
  }

  // save the executable path
  app::Apps::GetExecutablePath() = std::string( ARGUMENTS[ 0]);

  // create the command state object; parse all arguments except the 1st, which is the bcl executable command
  command::CommandState &cmd_state( command::CommandState::GetGlobalCommandState());

  // set bool that we're in main's command line parsing module
  command::CommandState::GetInMainCommandLineParsing() = true;

  // create a string stream to catch all error output
  std::ostringstream errors;

  // test whether any parameters were given. If not, try to read the default command line arguments file
  bool argument_parsing_success( cmd_state.ParseArguments( NUMBER_ARGUMENTS - 1, ARGUMENTS + 1, errors));

  // if no parameters were given, and the default command line file exists, use it as a response file
  if
  (
    argument_parsing_success
    && cmd_state.GetNumberRemainingParameters() == size_t( 0)
    && io::DirectoryEntry( app::Apps::GetDefaultCommandLineArgumentsFile()).DoesExist()
  )
  {
    errors << "No application name given; trying to commands from "
           << app::Apps::GetDefaultCommandLineArgumentsFile() << '\n';
    argument_parsing_success =
      cmd_state.ParseArguments
      (
        storage::Vector< std::string>( size_t( 1), app::Apps::GetDefaultCommandLineArgumentsFile()),
        errors
      );
  }

  // handle required arguments
  bool required_flags_success( command::GetAppDefaultFlags().HandleRequiredFlags( cmd_state, errors));

  // return if bad argument parsing
  if( !argument_parsing_success)
  {
    util::GetLogger() << errors.str() << std::endl;
    BCL_ExitWithoutCallstack( "Bad command line", -1);
  }
  else if( !required_flags_success)
  {
    util::GetLogger() << errors.str() << std::endl;
    BCL_ExitWithoutCallstack( "", -1);
  }
  else if( command::AppDefaultFlags::GetHelpFlag()->GetFlag() || command::AppDefaultFlags::GetReadMeFlag()->GetFlag())
  {
    BCL_ExitWithoutCallstack( "", 0);
  }
  else if( !app::GetApps().GetApplicationsParameter()->GetWasSetInCommandLine())
  {
    app::Apps::WriteGenericHelp( util::GetLogger());
    BCL_ExitWithoutCallstack( "No application given", -1);
  }

  // application name given
  const std::string app_to_execute( app::GetApps().GetApplicationsParameter()->GetValue());

  // get the arguments as a vector
  const storage::Vector< std::string> arguments( util::StringListFromCharacterArray( NUMBER_ARGUMENTS, ARGUMENTS));

  // write commandline
  util::GetLogger() << command::CommandLineWriter::CreateCommandLine( arguments) << std::flush;

  // output information about version and build date and time
  util::GetLogger() << GetVersion().GetDetailedInfoString() << std::endl;

  // initialize the program start time
  util::Stopwatch program_timer( "bcl", util::Message::e_Silent, false);

  // find the application with the name
  app::ApplicationType selected_app( app_to_execute);

  // if there was no app with that name
  if( !selected_app.IsDefined())
  {
    return -1;
  }

  BCL_MessageStd( "executing application: " + selected_app.GetName());

  // Set the global command state object up

  // Initialize the command line object for this application
  util::ShPtr< command::Command> sp_cmd( ( *selected_app)->InitializeCommand());

  // set the remainder of the flags
  const bool valid_command_line( sp_cmd->SetFlags( cmd_state, util::GetLogger()));

  // set bool that we're done with Main's command line parsing module
  command::CommandState::GetInMainCommandLineParsing() = false;

  if( command::CommandState::GetWasHelpRequested() && command::CommandState::GetWasHelpGiven())
  {
    BCL_ExitWithoutCallstack( "", 0);
  }

  // check the commandline was valid
  BCL_UserAssert
  (
    valid_command_line,
    "incorrect command line\n"
    "\nuse -help to see all options or -readme for more information for " + selected_app.GetName()
  );

  if( !command::CommandState::GetWasHelpRequested())
  {
    // write all given parameters for checking
    sp_cmd->WriteUserCommand( util::GetLogger());
  }

  // test whether any warning messages were given
  if( errors.tellp())
  {
    BCL_MessageCrt( errors.str());
  }

  // run the selected application
  const int app_return_code( ( *selected_app)->Main());

  if( command::CommandState::GetWasHelpRequested())
  {
    // an error state is normal if help was requested
    BCL_ExitWithoutCallstack( "", 0);
  }
  util::GetLogger() << '\n' << command::Command::DefaultSectionSeparator();

  // create a message about how much time and memory the run used
  util::MemoryUsage memory;
  std::ostringstream oss;
  oss << "bcl has run for " << program_timer.GetTotalTime().GetTimeAsHourMinuteSecond();
  if( util::IsDefined( memory.GetPeakVirtualMemoryUsed()))
  {
    oss << ", peak virtual memory used: " << memory.GetPeakVirtualMemoryUsed() << " MB";
  }
  if( util::IsDefined( memory.GetPeakRAMUsed()))
  {
    oss << ", peak physical RAM used: " << memory.GetPeakRAMUsed() << " MB";
  }
  BCL_MessageTop( oss.str());

  // finalize the environment
  util::GetRuntimeEnvironment().Finalize( app_return_code);

  BCL_ExitWithoutCallstack( "", app_return_code);

  // end
  return 0;
} // main
