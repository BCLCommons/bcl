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
#include "command/bcl_command_app_default_flags.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "app/bcl_app_group_handler.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_enums_instances.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_runtime_environments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    //! @brief generic function to write the help
    void WriteHelp()
    {
      // test whether the help flag was set
      if( AppDefaultFlags::GetHelpFlag()->GetFlag())
      {
        app::Apps::WriteGenericHelp( util::GetLogger());
      }
    }

    //! @brief generic function to write the read me if an app was given
    void WriteReadMe()
    {
      // test whether the help flag was set
      if( AppDefaultFlags::GetReadMeFlag()->GetFlag())
      {
        if( app::GetApps().GetApplicationsParameter()->GetWasSetInCommandLine())
        {
          app::ApplicationType app( app::GetApps().GetApplicationsParameter()->GetValue());
          io::FixedLineWidthWriter writer;
          writer << ( *app)->GetReadMe();
          util::GetLogger() << writer.String() << std::endl;
        }
        else
        {
          io::FixedLineWidthWriter writer;
          writer << app::Interface::DefaultReadMe();
          util::GetLogger() << writer.String() << std::endl;
        }
      }
    }

    //! @brief default flag that writes all options and exits the program
    //! @return FlagInterface containing the help flag
    util::ShPtr< FlagInterface> &AppDefaultFlags::GetHelpFlag()
    {
      static util::ShPtr< FlagInterface> s_help_flag
      (
        new FlagStatic
        (
          "help",
          "output user help for the bcl or an application",
          &WriteHelp
        )
      );

      return s_help_flag;
    }

    //! default flag that writes readme information if one is available
    util::ShPtr< FlagInterface> &AppDefaultFlags::GetReadMeFlag()
    {
      // create static flag for readme information
      static util::ShPtr< FlagInterface> s_readme_flag
      (
        new FlagStatic
        (
          "readme",
          "output readme information for the application",
          &WriteReadMe
        )
      );

      // end
      return s_readme_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! Adds all flags that have a required order to the enum
    AppDefaultFlags::AppDefaultFlags() :
      m_FlagsByType( s_NumberFlagTypes)
    {
      // The following list of data members must be initialized in this order, hence their inclusion here

      // Message level should be set 1st so that the user sees all messages of interest
      AddDefaultFlag( util::Message::GetMessageLevelFlag(), e_BclCore);

      // Runtime environment must be set before e_Logger to enable -logger File myfile.txt, since the path to the
      // actual output directory depends on the runtime
      if( util::GetRuntimeEnvironments().GetEnumCount() > size_t( 1))
      {
        AddDefaultFlag( util::RuntimeEnvironments::GetFlagRuntimeEnvironment(), e_BclCore);
      }

      // Logger must be setup before any later options since other options may call the logger
      // This implies that any problems with setting up the earlier flags must be output to a string stream that is
      // later passed to the logger, to ensure that the desired messages always go to the correct place
      AddDefaultFlag( util::Loggers::GetFlagLogger(), e_BclCore);

      // License must be read before application parameter to determine whether to enable the various apps
      //if( GetVersion().IsLicense())
      //{
      //  AddDefaultFlag( License::GetLicenseFileNameFlag(), e_BclCore);
      //}

      // Next add flags that are (possibly) dependent on the application

      // Help flag comes checked next because its output depends on whether the application parameter was given
      AddDefaultFlag( GetHelpFlag(), e_AppGeneric);

      //! Help flag should be checked next because its output depends on whether the application parameter was given
      AddDefaultFlag( GetReadMeFlag(), e_AppGeneric);

      // even though IO is the next type, its member flags (file compression and compressed alternatives)
      // have no order dependence, so they can be setup by static initialization rather than including them here

      // List of files to read to override normal enums; must be set after sh ptr behavior
      AddDefaultFlag( util::EnumsInstances::GetFlagEnumsFiles(), e_Util);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AppDefaultFlags::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the only instance of this class
    //! @return the only instance of this class
    AppDefaultFlags &AppDefaultFlags::GetInstance()
    {
      static AppDefaultFlags s_instance;
      return s_instance;
    }

    //! @brief get all the flags
    //! @return all the flags
    const util::ShPtrVector< FlagInterface> &AppDefaultFlags::GetAllFlags() const
    {
      return m_Flags;
    }

    //! @brief detect whether a given flag is a bcl app default flag
    //! @param FLAG_NAME the flag name to test
    //! @return type of the flag
    FlagType AppDefaultFlags::GetFlagType( const std::string &FLAG_NAME)
    {
      storage::Map< std::string, FlagTypeEnum>::const_iterator itr( m_FlagToType.Find( FLAG_NAME));
      return itr != m_FlagToType.End() ? itr->second : FlagTypeEnum( e_AppSpecific);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a default flag, to the collection
    //! @param SH_PTR_FLAG shptr to flag, to be added to default commandline flags
    //! @return the enum, that was added
    const util::ShPtr< FlagInterface> &AppDefaultFlags::AddDefaultFlag
    (
      const util::ShPtr< FlagInterface> &SH_PTR_FLAG,
      const FlagType &TYPE
    )
    {
      BCL_Assert
      (
        m_FlagToType.Insert( std::make_pair( SH_PTR_FLAG->GetName(), FlagTypeEnum( TYPE))).second,
        "Tried to add duplicate default flag with name: " + SH_PTR_FLAG->GetName()
      );
      m_Flags.PushBack( SH_PTR_FLAG);
      m_FlagsByType( size_t( TYPE)).PushBack( SH_PTR_FLAG);
      return SH_PTR_FLAG;
    }

    //! @brief get the flags for a given type
    //! @param TYPE the type of flag of interest
    //! @return the flags of the requested type
    const util::ShPtrVector< FlagInterface> &AppDefaultFlags::GetDefaultFlagsOfType( const FlagType &TYPE) const
    {
      return m_FlagsByType( size_t( TYPE));
    }

    //! @brief add the bcl default flags to the Command object
    //! @param COMMAND reference to a commandline object
    //! @param TYPES types of command line flags to add.  Note that e_Required will always be added
    void AppDefaultFlags::AddDefaultCommandlineFlags
    (
      Command &COMMAND,
      const storage::Set< FlagTypeEnum> &TYPES
    )
    {
      storage::Set< FlagTypeEnum> types( TYPES);
      // always add the required types
      types.Insert( e_AppGeneric);
      types.Insert( e_BclCore);
      // iterate over each type desired
      for( storage::Set< FlagTypeEnum>::const_iterator itr( types.Begin()), itr_end( types.End()); itr != itr_end; ++itr)
      {
        // iterate over all flags of that type
        for
        (
          util::ShPtrVector< FlagInterface>::const_iterator
            itr_flag( m_FlagsByType( size_t( *itr)).Begin()), itr_flag_end( m_FlagsByType( size_t( *itr)).End());
          itr_flag != itr_flag_end;
          ++itr_flag
        )
        {
          COMMAND.AddFlag( *itr_flag);
        }
      }
    }

    //! @brief add the bcl default flags required to the Command object
    //! @param COMMAND reference to a commandline object
    //! Note that e_BclRequired and e_AppGeneric will always be added
    void AppDefaultFlags::AddRequiredCommandlineFlags( Command &COMMAND)
    {
      AddDefaultCommandlineFlags( COMMAND, storage::Set< FlagTypeEnum>());
    }

    //! @brief add the bcl default flags to the Command object
    //! @param COMMAND reference to a commandline object
    void AppDefaultFlags::AddDefaultCommandlineFlags( Command &COMMAND)
    {
      this->AddDefaultCommandlineFlags( COMMAND, GetDefaultFlagTypes());
    }

    //! @brief handle the bcl-required flags
    //! @param STATE code state; from command line
    //! @param ERR_STREAM stream for output of errors
    //! @return true on success
    bool AppDefaultFlags::HandleRequiredFlags( CommandState &STATE, std::ostream &ERR_STREAM)
    {
      // iterate over all bcl required flags
      util::ShPtrVector< FlagInterface> &required_flags( m_FlagsByType( size_t( e_BclCore)));
      bool success( true);
      for
      (
        util::ShPtrVector< FlagInterface>::iterator
          itr_flag( required_flags.Begin()), itr_flag_end( required_flags.End());
        itr_flag != itr_flag_end;
        ++itr_flag
      )
      {
        // try to set the flag
        success = STATE.Update( **itr_flag, ERR_STREAM) && success;
      }

      // test whether any parameters were given, if so, parse the application parameter, since that must be 1st
      ParameterInterface &applications_parameter( *app::GetApps().GetApplicationsParameter());
      if( STATE.GetNumberRemainingParameters() && !applications_parameter.GetWasSetInCommandLine())
      {
        success = STATE.Update( applications_parameter, ERR_STREAM) && success;

        // check whether the user is licensed to run this application
        if( success)
        {
          // create the application
          app::ApplicationType application( applications_parameter.GetValue());

          // we are always licensed to run the BCL now
//          if
//          (
//            GetVersion().IsLicense()
//            //&& !License::GetLicenseFromCommandline().IsLicensedApplication( ( *application)->GetLicensedName()) &&
//            !util::EndsWith( ( *application)->GetLicensedName(), "Help")
//          )
//          {
//            ERR_STREAM << "You are not licensed to run " << applications_parameter.GetValue() << std::endl;
//            return false;
//          }

          // test whether the application name has been deprecated
          // split the string into the group name and app name
          storage::Vector< std::string> group_app_name
          (
            util::SplitString
            (
              applications_parameter.GetValue(),
              std::string( size_t( 1), app::GroupHandler::s_GroupDelimiter)
            )
          );

          const std::string group_name( group_app_name.GetSize() == size_t( 2) ? group_app_name( 0) : std::string());

          if( ( *application)->GetNameForGroup( group_name) != applications_parameter.GetValue() || group_name.empty())
          {
            const storage::Vector< std::string> group_names
            (
              app::GetAppGroups().GetApplicationGroupsForApp( **application)
            );
            if( !group_names.IsEmpty())
            {
              io::FixedLineWidthWriter output;
              output << "WARNING: " << applications_parameter.GetValue()
                     << " is a deprecated name for this application.  Use ";
              if( group_names.GetSize() > size_t( 1))
              {
                output << "one of the following: ";
              }
              output << ( *application)->GetNameForGroup( group_names.FirstElement());
              for
              (
                storage::Vector< std::string>::const_iterator itr( group_names.Begin() + 1), itr_end( group_names.End());
                itr != itr_end;
                ++itr
              )
              {
                output << ", " << ( *application)->GetNameForGroup( *itr);
              }
              output << " instead!";
              ERR_STREAM << output.String();
            }
          }
        }

        // handle application required flags
        util::ShPtrVector< FlagInterface> &app_required_flags( m_FlagsByType( size_t( e_AppGeneric)));
        for
        (
          util::ShPtrVector< FlagInterface>::iterator
            itr_flag( app_required_flags.Begin()), itr_flag_end( app_required_flags.End());
          itr_flag != itr_flag_end;
          ++itr_flag
        )
        {
          // try to set the flag
          success = STATE.Update( **itr_flag, ERR_STREAM) && success;
        }
      }
      else
      {
        // just check for the help flag
        success = STATE.Update( *GetHelpFlag(), ERR_STREAM) && success;
      }
      return success;
    }

    //! @brief handle all default flags contained in the given command flags
    //! @param STATE code state; from command line
    //! @param FLAGS the flags to consider
    //! @param ERR_STREAM stream for output of errors
    //! @return true on success
    bool AppDefaultFlags::HandleFlags
    (
      CommandState &STATE,
      util::ShPtrVector< FlagInterface> &FLAGS,
      std::ostream &ERR_STREAM
    )
    {
      bool success( true);
      for
      (
        util::ShPtrVector< FlagInterface>::iterator itr_flag( FLAGS.Begin()), itr_flag_end( FLAGS.End());
        itr_flag != itr_flag_end;
        ++itr_flag
      )
      {
        // try to set the flag
        success = STATE.Update( **itr_flag, ERR_STREAM) && success;
      }
      return success;
    }

    //! @brief construct on access function for all AppDefaultFlags
    //! @return reference to only instances of AppDefaultFlags
    AppDefaultFlags &GetAppDefaultFlags()
    {
      return AppDefaultFlags::GetInstance();
    }

  } // namespace command

} // namespace bcl
