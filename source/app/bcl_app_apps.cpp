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
#include "command/bcl_command_parameter_check_enumerate.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Apps::Apps() :
      util::Enumerate< util::ShPtr< Interface>, Apps>( false),
      m_ReleaseApplicationsNames()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Apps::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return vector of application names
    //! @return StorageVector<std::string> with application name
    const storage::Vector< std::string> &Apps::GetApplicationNames() const
    {
      static storage::Vector< std::string> s_app_names( GetEnumStrings());
      static bool s_sorted( false);

      // if this function was called before static initialization occurred, it is possible that
      // the number of applications has increased
      if( s_app_names.GetSize() != GetEnumCount())
      {
        BCL_MessageStd
        (
          "app::Apps::GetApplicationNames() should not normally be called before static initialization is performed"
        );
        s_app_names = GetEnumStrings();
        s_sorted = false;
      }

      // sort once
      if( !s_sorted)
      {
        std::sort( s_app_names.Begin(), s_app_names.End());
        s_sorted = true;
      }

      // return
      return s_app_names;
    }

    //! @brief return vector of release application names
    //! @return StorageVector<std::string> with release application name
    const storage::Vector< std::string> &Apps::GetReleaseApplicationNames() const
    {
      // return
      return m_ReleaseApplicationsNames;
    }

    //! @brief return a parameter with all application strings
    //! @return command::Parameter with check that contains all possible application name strings
    util::ShPtr< command::ParameterInterface> &Apps::GetApplicationsParameter()
    {
      // create Parameter with all application names
      static util::ShPtr< command::ParameterInterface> s_apps_parameter
      (
        new command::Parameter
        (
          "application_group>:<application",
          "name of the bcl application to be executed",
          command::ParameterCheckEnumerate< Apps>()
        )
      );

      // return
      return s_apps_parameter;
    }

    //! @brief return the executable path
    //! @return the actual executable path
    std::string &Apps::GetExecutablePath()
    {
      static std::string s_executable_path;
      return s_executable_path;
    }

    //! @brief add new Interface derived Application to Apps
    ApplicationType &Apps::AddEnum
    (
      const Interface &APPLICATION
    )
    {
      return this->AddEnum( APPLICATION.GetNameForGroup(), util::CloneToShPtr( APPLICATION));
    }

    //! @brief add the app directly
    //! @param NAME           name of the current application
    //! @param SP_APPLICATION object to be enumerated
    ApplicationType &Apps::AddEnum( const std::string &NAME, const util::ShPtr< Interface> &SP_APPLICATION)
    {
      // add release applications to the release app vector
      if( SP_APPLICATION->IsReleaseApplication())
      {
        if( m_ReleaseApplicationsNames.Find( NAME) >= m_ReleaseApplicationsNames.GetSize())
        {
          m_ReleaseApplicationsNames.PushBack( NAME);
          m_ReleaseApplicationsNames.Sort( std::less< std::string>());
        }
      }

      return util::Enumerate< util::ShPtr< Interface>, Apps>::AddEnum( NAME, SP_APPLICATION);
    }

    //! @brief writes the list of enums
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the list was written to
    //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
    std::ostream &Apps::WriteList( std::ostream &OSTREAM) const
    {
      // defer the call to app groups
      return GetAppGroups().WriteList( OSTREAM);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Apps::WriteGenericHelp( std::ostream &OSTREAM)
    {
      if( !GetApps().GetApplicationsParameter()->GetWasSetInCommandLine())
      {
        OSTREAM << "BCL Help\n";
        OSTREAM << GetVersion().GetDetailedInfoString() << '\n';
        OSTREAM << "Usage: " << GetExecutablePath() << " <application_group>:<application> [Other parameters] [Flags] [@filenames]\n";
        GetApps().GetApplicationsParameter()->WriteHelp( OSTREAM, 1);

        command::Command cmd;
        command::GetAppDefaultFlags().AddRequiredCommandlineFlags( cmd);
        cmd.WriteHelp( OSTREAM);
      }
      else
      {
        const std::string app_name( GetApps().GetApplicationsParameter()->GetValue());
        ApplicationType app( GetApps().GetApplicationsParameter()->GetValue());
        util::ShPtr< Interface> app_copy( app->HardCopy());
        util::ShPtr< command::Command> sp_cmd( app_copy->InitializeCommand());

        OSTREAM << '\n' << command::Command::DefaultSectionSeparator() << '\n';
        OSTREAM << app.GetName() << " Help\n";
        OSTREAM << GetVersion().GetDetailedInfoString() << '\n';
        OSTREAM << "Usage: " << GetExecutablePath() << ' ' << app.GetName() << ' ';
        sp_cmd->WriteUsage( util::GetLogger());
        sp_cmd->WriteHelp( util::GetLogger());
      }
      OSTREAM
       << "additional arguments (including the application) can be loaded from files by passing e.g. @my_file.txt"
       << "\n  If no application is specified, attempts to read commands from " << GetDefaultCommandLineArgumentsFile()
       << '\n';
      return OSTREAM;
    }

    //! @brief get the default command line
    //! If no arguments are given to an application, this file will be opened and its contents will be used as a command
    const std::string &Apps::GetDefaultCommandLineArgumentsFile()
    {
      static const std::string s_default_command_line_arguments_file( "bcl_commands.txt");
      return s_default_command_line_arguments_file;
    }

    Apps &GetApps()
    {
      return Apps::GetEnums();
    }

  } // namespace app

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< app::Interface>, app::Apps>;

  } // namespace util
} // namespace bcl
