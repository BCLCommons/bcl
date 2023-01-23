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
#include "command/bcl_command.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_command.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Command::s_Instance( GetObjectInstances().AddInstance( new Command()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Command::Command() :
      m_Parameters
      (
        new FlagStatic
        (
          "",
          ""
        )
      )
    {
    }

    //! @brief copy constructor
    Command::Command( const Command &COMMAND) :
      m_Parameters( COMMAND.m_Parameters),
      m_FlagsWithParams( COMMAND.m_FlagsWithParams),
      m_DefaultFlags( COMMAND.m_DefaultFlags)
    {
      RegenerateFlagNameMap();
    }

    //! @brief virtual copy constructor
    //! @return pointer to Command object
    Command *Command::Clone() const
    {
      return new Command( *this);
    }

    //! @brief assignment operator
    Command &Command::operator =( const Command &COMMAND)
    {
      m_Parameters = COMMAND.m_Parameters;
      m_FlagsWithParams = COMMAND.m_FlagsWithParams;
      m_DefaultFlags = COMMAND.m_DefaultFlags;
      RegenerateFlagNameMap();
      return *this;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Command::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the parameters
    //! @return the parameters
    util::ShPtr< FlagInterface> const &Command::GetParameters() const
    {
      return m_Parameters;
    }

    //! @brief sets the parameters
    //! @param PARAMETERS parameters to set
    void Command::SetParameters( const util::ShPtr< FlagInterface> &PARAMETERS)
    {
      m_Parameters = PARAMETERS;
    }

    //! @brief returns the flags with parameters that are application-specific (e.g. not bcl-wide)
    //! @return FlagInterfaces
    const util::ShPtrVector< FlagInterface> &Command::GetAppFlagsWithParams() const
    {
      return m_FlagsWithParams;
    }

    //! @brief returns the flags with parameters that are application-specific (e.g. not bcl-wide)
    //! @return FlagInterfaces of all default flags that have been added to this command
    const util::ShPtrVector< FlagInterface> &Command::GetBclFlagsWithParams() const
    {
      return m_DefaultFlags;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief read from storage::Vector of strings
    //! @param ARGUMENT_LIST - the list of arguments
    //! @param ERROR_STREAM - the error stream to which errors should be written
    //! @param DRY_RUN - true is this is only a dry run, and no signals should be triggered (default false)
    //! @return bool - returns true on success and false if there was an error
    bool Command::ReadArguments
    (
      const storage::Vector< std::string> &ARGUMENT_LIST,
      std::ostream &ERROR_STREAM,
      const bool DRY_RUN //= false
    )
    {
      CommandState cmd_state( DRY_RUN);
      bool argument_parsing_success( cmd_state.ParseArguments( ARGUMENT_LIST, ERROR_STREAM));

      // handle required arguments
      bool required_flags_success( GetAppDefaultFlags().HandleRequiredFlags( cmd_state, ERROR_STREAM));

      // return if bad argument parsing
      if( !argument_parsing_success || !required_flags_success)
      {
        return false;
      }
      else if( AppDefaultFlags::GetHelpFlag()->GetFlag() || AppDefaultFlags::GetReadMeFlag()->GetFlag())
      {
        return true;
      }

      return SetFlags( cmd_state, ERROR_STREAM);
    } // ReadArguments

    //! @brief read from a command state
    //! @param STATE - the current command state
    //! @param ERROR_STREAM - the error stream to which errors should be written
    //! @return bool - returns true on success and false if there was an error
    bool Command::SetFlags( CommandState &STATE, std::ostream &ERROR_STREAM)
    {
      // get the current flag type to consider
      const FlagType min_flag_type( STATE.GetFlagType());
      bool success_setting_flags( true);

      // go through and set all global flags with later flag type
      for
      (
        util::ShPtrVector< FlagInterface>::iterator
          flag_itr( m_DefaultFlags.Begin()),
          flag_itr_end( m_DefaultFlags.End());
        flag_itr != flag_itr_end;
        ++flag_itr
      )
      {
        // only consider flags of type that should be considered
        if( GetFlagTypeFromName( ( *flag_itr)->GetName()) > min_flag_type)
        {
          success_setting_flags &= STATE.Update( **flag_itr, ERROR_STREAM);
        }
        else
        {
          // check that the global flag was fine with whatever value it already had
          success_setting_flags &= ( *flag_itr)->IsValidList( ERROR_STREAM);
        }
      }

      // if there are any parameters
      if( STATE.GetNumberRemainingParameters())
      {
        // reset the parameters and then set them up
        m_Parameters->ResetFlag();
        success_setting_flags &= STATE.Update( *m_Parameters, ERROR_STREAM);
      }
      else
      {
        success_setting_flags &= m_Parameters->IsValidList( ERROR_STREAM);
      }

      // next, handle all application-specific flags
      for
      (
        util::ShPtrVector< FlagInterface>::iterator
          flag_itr( m_FlagsWithParams.Begin()),
          flag_itr_end( m_FlagsWithParams.End());
        flag_itr != flag_itr_end;
        ++flag_itr
      )
      {
        success_setting_flags &= STATE.Update( **flag_itr, ERROR_STREAM);
      }

      // if there were any flags remaining that could not be set
      const storage::Map< std::string, storage::Vector< std::string> > &state( STATE.GetState());
      if( !state.IsEmpty())
      {
        success_setting_flags = false;
        storage::Vector< std::string> valid_flag_names( m_FlagNamesToFlags.GetKeysAsVector());
        // iterate over all remaining flags in the state
        for
        (
          storage::Map< std::string, storage::Vector< std::string> >::const_iterator
            itr( state.Begin()), itr_end( state.End());
          itr != itr_end;
          ++itr
        )
        {
          Guesser::GetDefaultGuesser().WriteGuesses( itr->first, valid_flag_names, ERROR_STREAM, "flag");
        }
      }
      return success_setting_flags;
    }

    //! @brief checks whether a flag is set
    //! @param STRING the flag to check
    //! @return bool - returns true if the flag is set, false if it's not
    bool Command::IsFlagSet( const std::string &STRING)
    {
      storage::Map< std::string, util::SiPtr< FlagInterface> >::const_iterator itr( m_FlagNamesToFlags.Find( STRING));
      if( itr != m_FlagNamesToFlags.End())
      {
        return itr->second->GetFlag();
      }
      return false;
    } // IsFlagSet

    //! @brief returns a flag with parameters with given NAME
    //! @param NAME - the name of the requested flag
    //! @return util::ShPtr< FlagInterface> - a pointer to the requested flag
    util::SiPtr< FlagInterface> Command::GetFlagWithParams( const std::string &NAME)
    {
      std::string flagname( NAME);

      // remove leading -
      if( flagname[ 0] == '-')
      {
        flagname = flagname.substr( 1, flagname.length() - 1);
      }
      storage::Map< std::string, util::SiPtr< FlagInterface> >::const_iterator itr( m_FlagNamesToFlags.Find( flagname));
      if( itr != m_FlagNamesToFlags.End())
      {
        return itr->second;
      }

      // return empty ShPtr
      return util::SiPtr< FlagInterface>();
    } // GetFlagWithParams

    //! @brief add an additional parameter to the commandline
    //! @param COMMANDLINE_PARAMETER adds given ParameterInterface
    void Command::AddParameter( const util::ShPtr< ParameterInterface> &COMMANDLINE_PARAMETER)
    {
      util::ShPtr< FlagStatic>( m_Parameters)->PushBack( COMMANDLINE_PARAMETER);
    }

    //! @brief add an additional FlagWithParams to the commandline
    //! @param COMMANDLINE_FLAG_PARAMS add given FlagInterface
    void Command::AddFlag
    (
      const util::ShPtr< FlagInterface> &COMMANDLINE_FLAG_PARAMS
    )
    {
      util::ShPtr< FlagInterface> flag_to_add( COMMANDLINE_FLAG_PARAMS);
      // default flags should be separated
      if( GetFlagTypeFromName( COMMANDLINE_FLAG_PARAMS->GetName()) != e_AppSpecific)
      {
        m_DefaultFlags.PushBack( flag_to_add);
        if( !m_DefaultFlags.LastElement()->GetFlag())
        {
          // calling reset flag is necessary on any unset flags in case additional setup takes place in the Clean
          // function of the parameter check interface
          m_DefaultFlags.LastElement()->ResetFlag();
        }
      }
      else
      {
        m_FlagsWithParams.PushBack( flag_to_add);
        m_FlagsWithParams.LastElement()->ResetFlag();
      }
      // only insert real flags into the map, e.g. not FlagSeparator
      if( !flag_to_add->GetName().empty())
      {
        // assert that either there was no flag with that name already in the map (e.g. it was successfully inserted)
        BCL_Assert
        (
          m_FlagNamesToFlags.Insert
          (
            std::make_pair( flag_to_add->GetName(), util::SiPtr< FlagInterface>( flag_to_add))
          ).second,
          "Tried to add duplicate flag with name " + flag_to_add->GetName()
        );
      }
    }

    //! @brief reset all flags and parameters
    void Command::ResetFlagsAndParameters()
    {
      // reset all parameters
      m_Parameters->ResetFlag();

      // reset all flags
      for
      (
        util::ShPtrVector< FlagInterface>::iterator
          itr_flagwithparams( m_FlagsWithParams.Begin()),
          itr_flagwithparams_end( m_FlagsWithParams.End());
          itr_flagwithparams != itr_flagwithparams_end;
        itr_flagwithparams++
      )
      {
        ( *itr_flagwithparams)->ResetFlag();
      }
    }

    //! @brief add additional FlagWithParams
    //! @param COMMANDLINE_FLAG_PARAMS adds all FlagInterfaces
    void Command::PushBack( const util::ShPtrVector< FlagInterface> &COMMANDLINE_FLAG_PARAMS)
    {
      // iterate over given flags
      for
      (
        util::ShPtrVector< FlagInterface>::const_iterator
          flag_itr( COMMANDLINE_FLAG_PARAMS.Begin()),
          flag_itr_end( COMMANDLINE_FLAG_PARAMS.End());
        flag_itr != flag_itr_end; ++flag_itr
      )
      {
        // pushback this flag
        AddFlag( *flag_itr);
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the string seperating sections
    //! @return line break terminated separator string
    const std::string &Command::DefaultSectionSeparator()
    {
      static const std::string s_separator
      (
        std::string( util::GetLogger().GetMaxLineWidth(), '=') + std::string( 1, '\n')
      );

      return s_separator;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Command::WriteHelp( std::ostream &OSTREAM) const
    {
      OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      if( !m_Parameters->GetParameterList().IsEmpty())
      {
        // write parameters
        OSTREAM << "PARAMETERS: arguments that immediately follow the application name\n\n";
        m_Parameters->WriteHelp( OSTREAM);
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      // write BCL flags first so that user doesn't have to scroll to see help messages
      if( !m_DefaultFlags.IsEmpty())
      {
        // write flags
        OSTREAM << "BCL FLAGS: affect general BCL functionality, but may not be relevant for all applications\n\n";
        for
        (
          util::ShPtrVector< FlagInterface>::const_iterator
            itr( m_DefaultFlags.Begin()), itr_end( m_DefaultFlags.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !( *itr)->IsUnused())
          {
            ( *itr)->WriteHelp( OSTREAM);
          }
        }
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      if( !m_FlagsWithParams.IsEmpty())
      {
        // write flags
        OSTREAM << "APPLICATION FLAGS: syntax: -flagname [flagparameter1] ... \n\n";
        util::ShPtrVector< FlagInterface> all_flags( m_DefaultFlags);
        all_flags.Append( m_FlagsWithParams);
        for
        (
          util::ShPtrVector< FlagInterface>::const_iterator
            itr( all_flags.Begin()), itr_end( all_flags.End());
          itr != itr_end;
          ++itr
        )
        {
          // get flag type
          FlagType flag_type( GetFlagTypeFromName( ( *itr)->GetName()));

          // write help for everything that is application specific
          if( flag_type == e_AppGeneric || flag_type == e_AppSpecific)
          {
            if( !( *itr)->IsUnused())
            {
              ( *itr)->WriteHelp( OSTREAM);
            }
          }
        }
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      return OSTREAM;
    } // WriteHelp

    //! @brief writes the user provided commandline
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Command::WriteUserCommand( std::ostream &OSTREAM) const
    {
      OSTREAM << '\n' << DefaultSectionSeparator() << '\n';

      if( !m_Parameters->GetParameterList().IsEmpty())
      {
        //write parameters
        OSTREAM << "PARAMETERS\n\n";
        m_Parameters->WriteUserCommand( OSTREAM);
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      // track default flags to be written out later
      util::ShPtrVector< FlagInterface> default_flags( m_DefaultFlags);
      if( !default_flags.IsEmpty())
      {
        //write bcl flags
        OSTREAM << "BCL FLAGS\n\n";
        for
        (
          util::ShPtrVector< FlagInterface>::const_iterator
            itr( default_flags.Begin()), itr_end( default_flags.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !( *itr)->IsUnused())
          {
            ( *itr)->WriteUserCommand( OSTREAM);
          }
        }
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      if( !m_FlagsWithParams.IsEmpty())
      {
        //write flags
        OSTREAM << "APPLICATION FLAGS\n\n";
        util::ShPtrVector< FlagInterface> all_flags( m_DefaultFlags);
        all_flags.Append( m_FlagsWithParams);
        for
        (
          util::ShPtrVector< FlagInterface>::const_iterator
            itr( all_flags.Begin()), itr_end( all_flags.End());
          itr != itr_end;
          ++itr
        )
        {
          // get flag type
          FlagType flag_type( GetFlagTypeFromName( ( *itr)->GetName()));

          // write help for everything that is application specific
          if( flag_type == e_AppGeneric || flag_type == e_AppSpecific)
          {
            if( !( *itr)->IsUnused())
            {
              ( *itr)->WriteUserCommand( OSTREAM);
            }
          }
        }
        OSTREAM << '\n' << DefaultSectionSeparator() << '\n';
      }

      return OSTREAM;
    } // WriteUserCommand

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Command::WriteUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_params( m_Parameters->GetNumberRequiredParameters());
      const size_t n_optional_params( m_Parameters->GetSize() - n_required_params);

      if( n_required_params)
      {
        m_Parameters->WriteRequiredUsage( OSTREAM);
      }

      if( n_optional_params)
      {
        OSTREAM << "[OPTIONAL PARAMETERS] ";
      }

      for
      (
        util::ShPtrVector< FlagInterface>::const_iterator
          itr( m_FlagsWithParams.Begin()),
          itr_end( m_FlagsWithParams.End());
          itr != itr_end;
        ++itr
      )
      {
        if( !( *itr)->IsUnused())
        {
          ( *itr)->WriteRequiredUsage( OSTREAM);
        }
      }

      OSTREAM << " [OPTIONAL FLAGS] [@FILENAMES]\n";
      return OSTREAM;
    }

    //! @brief writes Command to ostream OSTREAM
    //! @param OSTREAM - the outstream to write to
    //! @param INDENT indentation
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Command::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Parameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FlagsWithParams, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultFlags, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    //! @brief read Command for std::istream ISTREAM
    //! @param ISTREAM - the instream to read from
    //! @return std::istream &ISTREAM - return the stream after reading
    std::istream &Command::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Parameters, ISTREAM);
      io::Serialize::Read( m_FlagsWithParams, ISTREAM);
      io::Serialize::Read( m_DefaultFlags, ISTREAM);
      RegenerateFlagNameMap();

      // end
      return ISTREAM;
    } // Read

    //! @brief regenerate the flag to names map
    void Command::RegenerateFlagNameMap()
    {
      m_FlagNamesToFlags.Reset();

      util::ShPtrVector< FlagInterface> all_flags( m_DefaultFlags);
      all_flags.Append( m_FlagsWithParams);
      for
      (
        util::ShPtrVector< FlagInterface>::iterator
          itr( all_flags.Begin()), itr_end( all_flags.End());
        itr != itr_end;
        ++itr
      )
      {
        // only insert real flags into the map, e.g. not FlagSeparator
        if( !( *itr)->GetName().empty())
        {
          // assert that either there was no flag with that name already in the map (e.g. it was successfully inserted)
          BCL_Assert
          (
            m_FlagNamesToFlags.Insert
            (
              std::make_pair( ( *itr)->GetName(), util::SiPtr< FlagInterface>( *itr))
            ).second,
            "Tried to add duplicate flag with name " + ( *itr)->GetName()
          );
        }
      }
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_default_flag_types.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    //! @brief Type as string
    //! @param TYPE the flag type
    //! @return the string for the flag type
    const std::string &GetFlagTypeName( const FlagType &TYPE)
    {
      static const std::string s_names[ s_NumberFlagTypes + 1] =
      {
        "BclCore",
        "AppGeneric",
        "Io",
        "Util",
        "Random",
        "Model",
        "Score",
        "Db",
        "Opencl",
        "Pthread",
        "OpenBlas",
        "AppSpecific",
        GetStaticClassName< FlagType>()
      };
      return s_names[ TYPE];
    }

    //! @brief get the default types to be included in all applications
    const storage::Set< FlagTypeEnum> &GetDefaultFlagTypes()
    {
      // by default, include all the types, except app specific flags
      static const storage::Set< FlagTypeEnum> s_default_types
      (
        FlagTypeEnum::GetEnumVector()[ 0],
        FlagTypeEnum::GetEnumVector()[ e_AppSpecific]
      );

      return s_default_types;
    }

    //! @brief get the default flag types to be included in all applications
    //! @param FLAG_NAME name of the flag
    FlagType GetFlagTypeFromName( const std::string &FLAG_NAME)
    {
      return GetAppDefaultFlags().GetFlagType( FLAG_NAME);
    }

  } // namespace command

} // namespace bcl
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
#include "command/bcl_command_flag_dynamic.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_data_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FlagDynamic::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagDynamic())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagDynamic::FlagDynamic() :
      m_WasSetInCommandline( false)
    {
    }

    //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETER,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagDynamic::FlagDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &TEMPLATE_PARAMETER,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      FlagInterface::t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameters( 1, util::CloneToShPtr( TEMPLATE_PARAMETER)),
      m_Signal( SIGNAL)
    {
      BCL_Assert
      (
        MIN_NUMBER_PARAMETERS <= MAX_NUMBER_PARAMETERS,
        "min number parameters > max number parameters: "
          + util::Format()( MIN_NUMBER_PARAMETERS) + " > "
          + util::Format()( MAX_NUMBER_PARAMETERS)
      );
      if( TEMPLATE_PARAMETER.GetName().empty())
      {
        m_TemplateParameters.LastElement()->SetName( "1");
      }
    } // FlagDynamic

    //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETERS,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagDynamic::FlagDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const storage::Vector< Parameter> &TEMPLATE_PARAMETERS,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      FlagInterface::t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( DESCRIPTION),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameters(),
      m_Signal( SIGNAL)
    {
      BCL_Assert
      (
        MIN_NUMBER_PARAMETERS <= MAX_NUMBER_PARAMETERS,
        "min number parameters > max number parameters: "
          + util::Format()( MIN_NUMBER_PARAMETERS) + " > "
          + util::Format()( MAX_NUMBER_PARAMETERS)
      );
      BCL_Assert
      (
        ( MIN_NUMBER_PARAMETERS % TEMPLATE_PARAMETERS.GetSize()) == 0,
        "min number parameters is not a multiple of the number of internal parameters"
      );
      BCL_Assert
      (
        !util::IsDefined( MAX_NUMBER_PARAMETERS) || ( MAX_NUMBER_PARAMETERS % TEMPLATE_PARAMETERS.GetSize()) == 0,
        "max number parameters is not a multiple of the number of internal parameters"
      );

      size_t counter( 1);
      for
      (
        storage::Vector< Parameter>::const_iterator
          itr( TEMPLATE_PARAMETERS.Begin()), itr_end( TEMPLATE_PARAMETERS.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        m_TemplateParameters.PushBack( util::CloneToShPtr( *itr));
        if( itr->GetName().empty())
        {
          m_TemplateParameters.LastElement()->SetName( util::Format()( counter));
        }
      }
    } // FlagDynamic

    //! @brief virtual copy constructor
    FlagDynamic *FlagDynamic::Clone() const
    {
      return new FlagDynamic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return the name of the class as a std::string
    const std::string &FlagDynamic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagDynamic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagDynamic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagDynamic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief set flag to true if it was false
    void FlagDynamic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagDynamic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagDynamic::ResetFlag()
    {
      m_WasSetInCommandline = false;

      // iterate over parameters in current flag and reset them all
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end;
        ++flag_param_itr
      )
      {
        ( *flag_param_itr)->Reset();
      }

      m_ParameterList.Reset();

      // call the handling function, if it was defined
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagDynamic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagDynamic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagDynamic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagDynamic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagDynamic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagDynamic::GetNumberRequiredParameters() const
    {
      return m_TemplateParameters.GetSize() * m_MinNumberParameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST - the parameter list as a storage vector of strings
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false on failure
    //! does not allow you to add more than m_MaxNumberParameters parameters to the list
    //! @see SetFlag
    //! @see PushBack
    bool FlagDynamic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      // set flag
      SetFlag();

      bool success_setting_parameters( true);

      // iterate over parameter in given commandline
      for
      (
        std::vector< std::string>::const_iterator
          param_itr( PARAMETER_LIST.Begin()),
          param_itr_end( PARAMETER_LIST.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        // pushback new entries in the parameterlist
        success_setting_parameters &= PushBack( *param_itr, ERROR_STREAM);
      }

      // ensure that list is long enough - but not too long
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // call the handling function, if it was defined
      if( success_setting_parameters && m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }

      return success_setting_parameters;
    } // ReadFromList

    //! @brief checks if there are too few items in your parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true if the parameter list is of a valid size, false otherwise
    //! @see MeetSizeSpecification
    bool FlagDynamic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool is_valid = MeetSizeSpecification();

      // check if parameter list has sufficient size
      if( !is_valid)
      {
        if( m_ParameterList.GetSize() % m_TemplateParameters.GetSize())
        {
          ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has improper size: given " << GetSize()
            << " parameters, but need a multiple of " << m_TemplateParameters.GetSize() << '\n';
        }
        else
        {
          const size_t multiplicity( m_ParameterList.GetSize() / m_TemplateParameters.GetSize());
          ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has wrong size: " << multiplicity
            << " is not between: [" << m_MinNumberParameters << ".." << m_MaxNumberParameters << "]" << '\n';
        }
      }

      return is_valid;
    } // IsValidList

    //! @brief returns true if the number of parameters is between m_MinNumberParameters and m_MaxNumberParameters
    //! @return true if the parameter list meets the size specification, false otherwise
    bool FlagDynamic::MeetSizeSpecification() const
    {
      return
          ( m_ParameterList.GetSize() % m_TemplateParameters.GetSize()) == 0
          && ( GetSize() / m_TemplateParameters.GetSize() >= m_MinNumberParameters)
          && ( GetSize() / m_TemplateParameters.GetSize() <= m_MaxNumberParameters);
    }

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER - the parameter to be added to the parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false otherwise
    //! @see FlagDynamic::PushBack
    bool FlagDynamic::PushBack( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      // determine which template parameter this member of the flag belongs to
      const size_t parameter_id( m_ParameterList.GetSize() % m_TemplateParameters.GetSize());

      util::ShPtr< ParameterInterface> new_parameter
      (
        m_TemplateParameters( parameter_id).HardCopy()
      );
      // try to set parameter
      if( new_parameter->SetParameter( PARAMETER, ERROR_STREAM))
      {
        // if this was an allowed parameter => pushback
        m_ParameterList.PushBack( new_parameter);
        return true;
      }
      return false;
    } // PushBack

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the command line
    //! @param OSTREAM - the stream to write to
    //! @return the stream after you wrote to it
    std::ostream &FlagDynamic::WriteHelp( std::ostream &OSTREAM) const
    {
      io::FixedLineWidthWriter writer;
      // write name and description
      if( !m_Name.empty() || !m_Description.empty())
      {
        // write out the flag name directly, wrap the description / size part, indenting to beginning of the description
        writer << '-' << m_Name << " : ";

        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          // short description
          writer.SetIndent( writer.GetLinePosition());
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 2);
        }
        writer << m_Description << ", ";

        // get the description for the sizes
        std::string size_description;
        std::ostringstream description_stream;
        util::DataType::WriteSizeRequirements( description_stream, m_MinNumberParameters, m_MaxNumberParameters);
        size_description = description_stream.str();
        // erase with, if it starts with that
        if( util::StartsWith( size_description, " with"))
        {
          size_description.erase( 0, 5);
        }
        else
        {
          // any number
          size_description = " any number of ";
        }

        // add the size description to the writer
        writer << "This flag can be followed by " << size_description;
        writer.PopIndent();
      }

      // add a newline, if one was not already present
      if( writer.GetLinePosition())
      {
        writer.NewLine();
      }
      OSTREAM << writer.String();

      // write the writer to the output stream

      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_TemplateParameters.Begin()), param_itr_end( m_TemplateParameters.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      }

      // end
      return OSTREAM;
    } // WriteHelp

    //! @brief writes the user provided commandline
    //! @return the stream after you write to it
    std::ostream &FlagDynamic::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "-" << m_Name << ( m_WasSetInCommandline ? " set" : " not set") << '\n';

      // write WriteUserCommand for every parameter in the list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()), param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        OSTREAM << "   ";
        ( *param_itr)->WriteUserCommand( OSTREAM);
      } // for

      // end
      return OSTREAM;
    } // WriteUserCommand

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagDynamic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        OSTREAM << '-' << m_Name << ' ';
        for( size_t counter( 0); counter < m_MinNumberParameters; ++counter)
        {
          for
          (
            util::ShPtrVector< ParameterInterface>::const_iterator
              param_itr( m_TemplateParameters.Begin()),
              param_itr_end( m_TemplateParameters.End());
            param_itr != param_itr_end;
            ++param_itr
          )
          {
            OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
          }
        }
      }
      return OSTREAM;
    }

    //! @brief write Flag with Params to ostream
    //! @param OSTREAM - the stream to write to
    //! @param INDENT indentation
    //! @return the stream after you wrote to it
    std::ostream &FlagDynamic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TemplateParameters,  OSTREAM, INDENT);

      // end
      return OSTREAM;
    } // Write

    //! @brief read Flag with params from istream
    //! @param ISTREAM - the stream to read from
    //! @return the stream after you read from it
    std::istream &FlagDynamic::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);
      io::Serialize::Read( m_MinNumberParameters, ISTREAM);
      io::Serialize::Read( m_MaxNumberParameters, ISTREAM);
      io::Serialize::Read( m_TemplateParameters,  ISTREAM);

      // end
      return ISTREAM;
    } // Read

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_flag_separator.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FlagSeparator::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagSeparator())
    );

    //! static variable to undefined list of parameters
    util::ShPtrVector< ParameterInterface> FlagSeparator::s_Parameters;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagSeparator::FlagSeparator() :
      m_SeparatorText()
    {
    }

    //! @brief constructor from a separator text
    FlagSeparator::FlagSeparator( const std::string &SEPARATOR_TEXT) :
      m_SeparatorText( SEPARATOR_TEXT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FlagSeparator
    FlagSeparator *FlagSeparator::Clone() const
    {
      return new FlagSeparator( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FlagSeparator::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagSeparator::GetName() const
    {
      // initialize static empty string and return it
      static std::string s_name;
      return s_name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagSeparator::SetName( const std::string &NAME)
    {

    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagSeparator::GetDescription() const
    {
      // initialize empty string and return it
      static const std::string s_empty;
      return s_empty;
    }

    //! @brief set flag to true if it was false
    void FlagSeparator::SetFlag()
    {

    }

    //! @brief set flag to false if it was true
    void FlagSeparator::UnsetFlag()
    {

    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagSeparator::ResetFlag()
    {

    }

    //! @brief returns true if was set in command line
    //! @return whether this flag was set in the command line
    bool FlagSeparator::GetFlag() const
    {
      return true;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagSeparator::GetSignal() const
    {
      return NULL;
    }

    //! @brief returns the number of parameters
    //! @return number of parameters after flag
    size_t FlagSeparator::GetSize() const
    {
      return size_t( 0);
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagSeparator::GetParameterList() const
    {
      return s_Parameters;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagSeparator::GetParameterList()
    {
      return s_Parameters;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagSeparator::GetFirstParameter() const
    {
      // undefined first parameter
      static const util::ShPtr< ParameterInterface> s_undefined;

      // this function should never be called
      BCL_Assert( false, "this class should never have been asked for first parameter");

      // end
      return s_undefined;
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagSeparator::GetNumberRequiredParameters() const
    {
      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagSeparator::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      return true;
    }

    //! @brief checks that for every parameter was given a default value or it was passed through commandline
    //! @brief or if it is a dynamic list, if it meets the size specifications
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagSeparator::IsValidList( std::ostream &ERROR_STREAM) const
    {
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the help was written to
    std::ostream &FlagSeparator::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      OSTREAM << '\n' << m_SeparatorText << '\n';

      //end
      return OSTREAM;
    }

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &FlagSeparator::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // end
      return OSTREAM;
    }

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagSeparator::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FlagSeparator::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SeparatorText, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FlagSeparator::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SeparatorText, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_flag_static.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FlagStatic::s_Instance( GetObjectInstances().AddInstance( new FlagStatic()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagStatic::FlagStatic() :
      m_WasSetInCommandline( false),
      m_Signal( NULL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_Signal( SIGNAL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param PARAMETER the parameter associated with this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &PARAMETER,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_ParameterList( 1, PARAMETER),
      m_Signal( SIGNAL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param PARAMETERS the parameters associated with this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const util::ShPtrVector< ParameterInterface> &PARAMETERS,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_ParameterList( PARAMETERS),
      m_Signal( SIGNAL)
    {
    }

    //! @brief copy constructor
    //! @return pointer to FlagStatic object
    FlagStatic *FlagStatic::Clone() const
    {
      return new FlagStatic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return name of class as string
    //! @see GetStaticClassName()
    const std::string &FlagStatic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagStatic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagStatic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagStatic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief set flag to true if it was false
    void FlagStatic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagStatic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagStatic::ResetFlag()
    {
      m_WasSetInCommandline = false;

      // iterate over parameters in current flag and reset them all
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end;
        ++flag_param_itr
      )
      {
        ( *flag_param_itr)->Reset();
      }

      // if a signal should be emitted, emit it
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagStatic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagStatic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagStatic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagStatic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagStatic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagStatic::GetNumberRequiredParameters() const
    {
      size_t counter( m_ParameterList.GetSize());

      // find the last parameter without a default value
      for
      (
        util::ShPtrVector< ParameterInterface>::const_reverse_iterator
          itr( m_ParameterList.ReverseBegin()),
          itr_end( m_ParameterList.ReverseEnd());
        itr != itr_end;
        ++itr, --counter
      )
      {
        // if no default was given, return that #
        if( !( *itr)->GetWasDefaultGiven())
        {
          return counter;
        }
      }

      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagStatic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      // set flag
      SetFlag();

      // check to make sure they didn't pass too many parameters for the given flag
      if( PARAMETER_LIST.GetSize() > GetSize())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" expects " << GetSize() << " parameters but received "
          << PARAMETER_LIST.GetSize() << '\n';
        return false;
      }

      bool success_setting_parameters( true);

      // iterator on begin of given parameters after option
      std::vector< std::string>::const_iterator
        param_itr( PARAMETER_LIST.Begin()),
        param_itr_end( PARAMETER_LIST.End());

      // iterate parallel over parameters in current flag with parameters and given parameters in commandline
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end && param_itr != param_itr_end;
        ++flag_param_itr, ++param_itr
      )
      {
        // set parameters in static parameter list
        success_setting_parameters &= ( *flag_param_itr)->SetParameter( *param_itr, ERROR_STREAM);
      } // for

      // check if the list of parameters is valid
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // if a signal should be emitted, emit it
      if( success_setting_parameters && m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }

      // end
      return success_setting_parameters;
    } // ReadFromList

    //! @brief checks that for every parameter was given a default value or it was passed through commandline
    //! @brief or if it is a dynamic list, if it meets the size specifications
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagStatic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool valid_list( true);
      size_t counter( 0);
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr, ++counter
      )
      {
        // if neither the parameter was set in command line, nor the parameter has a default
        if( !( *param_itr)->GetWasSetInCommandLine() && !( *param_itr)->GetWasDefaultGiven())
        {
          valid_list = false;
          ERROR_STREAM << "flag: " << m_Name;

          // if the parameter has a name, use it
          if( !( *param_itr)->GetName().empty())
          {
            ERROR_STREAM << ", parameter " << ( *param_itr)->GetName() << ":";
          }
          else if( m_ParameterList.GetSize() > size_t( 1))
          {
            // multiple parameters; this parameter unnamed, use the counter to id the parameter
            ERROR_STREAM << ", parameter #" << counter;
          }

          ERROR_STREAM << " was not given!" << '\n';
        }
      } // for

      return valid_list;
    } // IsValidList

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER the parameter to be added
    void FlagStatic::PushBack( const util::ShPtr< ParameterInterface> &PARAMETER)
    {
      m_ParameterList.PushBack( PARAMETER);
      if( PARAMETER->GetName().empty())
      {
        m_ParameterList.LastElement()->SetName( util::Format()( m_ParameterList.GetSize()));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &FlagStatic::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "-" << m_Name << ( m_WasSetInCommandline ? " set" : " not set") << '\n';

      // write WriteUserCommand for every parameter in the list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        OSTREAM << "   ";
        ( *param_itr)->WriteUserCommand( OSTREAM);
      } // for

      // end
      return OSTREAM;
    } // WriteUserCommand

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the help was written to
    std::ostream &FlagStatic::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      if( !m_Name.empty() || !m_Description.empty())
      {
        io::FixedLineWidthWriter writer;
        // write out the flag name directly, wrap the description / size part, indenting to beginning of the description
        writer << '-' << m_Name << " : ";

        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          // short description
          writer.SetIndent( writer.GetLinePosition());
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 2);
        }

        writer << m_Description;
        writer.PopIndent();

        // add a newline, if one is not already present
        if( writer.GetLinePosition())
        {
          writer.NewLine();
        }
        OSTREAM << writer.String();
      }

      // static list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      } // for

      //end
      return OSTREAM;
    } // WriteHelp

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagStatic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        if( !m_Name.empty())
        {
          OSTREAM << '-' << m_Name << ' ';
        }
        size_t counter( 0);
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator param_itr( m_ParameterList.Begin());
          counter < n_required_parameters;
          ++param_itr, ++counter
        )
        {
          OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
        }
      }
      return OSTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FlagStatic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FlagStatic::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_flag_static_and_dynamic.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_data_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FlagStaticAndDynamic::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagStaticAndDynamic())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagStaticAndDynamic::FlagStaticAndDynamic() :
      m_Name(),
      m_Description(),
      m_WasSetInCommandline( false),
      m_MinNumberParameters(),
      m_MaxNumberParameters(),
      m_ParameterList(),
      m_TemplateParameter(),
      m_NumberStaticParameters( 0)
    {
    }

    //! @brief construct FlagStaticAndDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETER,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStaticAndDynamic::FlagStaticAndDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &TEMPLATE_PARAMETER,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameter( TEMPLATE_PARAMETER.Clone()),
      m_NumberStaticParameters( 0),
      m_Signal( SIGNAL)
    {
    }

    //! @brief copy constructor
    FlagStaticAndDynamic *FlagStaticAndDynamic::Clone() const
    {
      return new FlagStaticAndDynamic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return the name of the class as a std::string
    const std::string &FlagStaticAndDynamic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagStaticAndDynamic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagStaticAndDynamic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagStaticAndDynamic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagStaticAndDynamic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief set flag to true if it was false
    void FlagStaticAndDynamic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagStaticAndDynamic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagStaticAndDynamic::ResetFlag()
    {
      m_WasSetInCommandline = false;

      // iterate over parameters in current flag and reset them all
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end;
        ++flag_param_itr
      )
      {
        ( *flag_param_itr)->Reset();
      }

      // remove any dynamic parameters
      const size_t number_dynamic_parameters( m_ParameterList.GetSize() - m_NumberStaticParameters);
      if( number_dynamic_parameters > size_t( 0))
      {
        m_ParameterList.RemoveElements( m_NumberStaticParameters, number_dynamic_parameters);
      }
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagStaticAndDynamic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagStaticAndDynamic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagStaticAndDynamic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagStaticAndDynamic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagStaticAndDynamic::GetNumberRequiredParameters() const
    {
      // if dynamic parameters are required, all static parameters must be given
      if( m_MinNumberParameters)
      {
        return m_MinNumberParameters + m_NumberStaticParameters;
      }
      size_t counter( m_ParameterList.GetSize());

      // find the last parameter without a default value
      for
      (
        util::ShPtrVector< ParameterInterface>::const_reverse_iterator
          itr( m_ParameterList.ReverseBegin()),
          itr_end( m_ParameterList.ReverseEnd());
        itr != itr_end;
        ++itr, --counter
      )
      {
        // if no default was given, update the # required to the counter
        if( !( *itr)->GetWasDefaultGiven())
        {
          return counter;
        }
      }

      return 0;
    }

    //! @brief returns list of dynamic parameters
    //! @return the parameter list of dynamically added parameters
    util::SiPtrVector< const ParameterInterface> FlagStaticAndDynamic::GetDynamicParameterList() const
    {
      util::SiPtrVector< const ParameterInterface> dynamic_list;
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          itr( m_ParameterList.Begin() + m_NumberStaticParameters), itr_end( m_ParameterList.End());
        itr != itr_end;
        ++itr
      )
      {
        // insert current parameter
        dynamic_list.PushBack( *itr);
      }

      // end
      return dynamic_list;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST - the parameter list as a storage vector of strings
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false on failure
    //! does not allow you to add more than m_MaxNumberParameters parameters to the list
    //! @see SetFlag
    //! @see PushBack
    bool FlagStaticAndDynamic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      SetFlag();

      if( PARAMETER_LIST.IsEmpty())
      {
        return true;
      }

      bool success_setting_parameters( true);

      // iterators on commandline arguments
      storage::Vector< std::string>::const_iterator
        param_itr( PARAMETER_LIST.Begin()),
        param_itr_end( PARAMETER_LIST.End());

      // read static params
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          itr( m_ParameterList.Begin()), itr_end( m_ParameterList.Begin() + m_NumberStaticParameters);
        itr != itr_end && param_itr != param_itr_end;
        ++itr, ++param_itr
      )
      {
        // set parameters in static parameter list
        success_setting_parameters &= ( *itr)->SetParameter( *param_itr, ERROR_STREAM);
      }

      // iterate over dynamic parameters in given commandline
      for( ; param_itr != param_itr_end; ++param_itr)
      {
        // pushback new entries in the parameterlist
        success_setting_parameters &= PushBack( *param_itr, ERROR_STREAM);
      }

      // ensure that list is long enough - but not too long
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // call the signalling function if successful and one was given
      if( success_setting_parameters && m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }

      return success_setting_parameters;
    } // ReadFromList

    //! @brief checks if there are too few items in your parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true if the parameter list is of a valid size, false otherwise
    //! @see MeetSizeSpecification
    bool FlagStaticAndDynamic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool valid_list( true);
      size_t counter( 0);
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( GetParameterList().Begin()),
          param_itr_end( GetParameterList().End());
        param_itr != param_itr_end;
        ++param_itr, ++counter
      )
      {
        // if neither the parameter was set in command line, nor the parameter has a default
        if( !( *param_itr)->GetWasSetInCommandLine() && !( *param_itr)->GetWasDefaultGiven())
        {
          valid_list = false;
          ERROR_STREAM << "flag: " << m_Name;

          // if the parameter has a name, use it
          if( !( *param_itr)->GetName().empty())
          {
            ERROR_STREAM << ", parameter " << ( *param_itr)->GetName() << ":";
          }
          else if( m_ParameterList.GetSize() > size_t( 1))
          {
            // multiple parameters; this parameter unnamed, use the counter to id the parameter
            ERROR_STREAM << ", parameter #" << counter;
          }

          ERROR_STREAM << " was not given!" << '\n';
        }
      } // for

      // check if parameter list has sufficient size
      if( !MeetSizeSpecification())
      {
        valid_list = false;
        ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has wrong size: " << GetSize()
          << " is not between: [" << m_MinNumberParameters << ".." << m_MaxNumberParameters << "]" << '\n';
      }

      return valid_list;
    }

    //! @brief returns true if the number of parameters is between m_MinNumberParameters and m_MaxNumberParameters
    //! @return true if the parameter list meets the size specification, false otherwise
    bool FlagStaticAndDynamic::MeetSizeSpecification() const
    {
      const size_t number_dynamic_parameters( GetSize() - m_NumberStaticParameters);
      return ( number_dynamic_parameters >= m_MinNumberParameters) && ( number_dynamic_parameters <= m_MaxNumberParameters);
    }

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER - the parameter to be added to the parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false otherwise
    //! @see FlagInterface::PushBack
    bool FlagStaticAndDynamic::PushBack( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      util::ShPtr< ParameterInterface> new_parameter( m_TemplateParameter.HardCopy());
      // try to set parameter
      if( new_parameter->SetParameter( PARAMETER, ERROR_STREAM))
      {
        // if this was an allowed parameter => pushback
        m_ParameterList.PushBack( new_parameter);
        return true;
      }
      return false;
    } // PushBack

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the command line
    //! @param OSTREAM - the stream to write to
    //! @return the stream after you wrote to it
    std::ostream &FlagStaticAndDynamic::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      io::FixedLineWidthWriter writer;

      // write name and description
      if( !m_Name.empty() || !m_Description.empty())
      {
        // write out the flag name directly, wrap the description / size part, indenting to beginning of the description
        writer << '-' << m_Name << " : ";

        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          // short description
          writer.SetIndent( writer.GetLinePosition());
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 2);
        }
        writer << m_Description;
      }

      // add a newline, if one is not already present
      if( writer.GetLinePosition())
      {
        writer.NewLine();
      }
      OSTREAM << writer.String();

      // get the description for the sizes
      std::string size_description;
      {
        std::ostringstream description_stream;
        util::DataType::WriteSizeRequirements( description_stream, m_MinNumberParameters, m_MaxNumberParameters);
        size_description = description_stream.str();
        // erase with, if it starts with that
        if( util::StartsWith( size_description, " with"))
        {
          size_description.erase( 0, 5);
        }
        else
        {
          // any number
          size_description = " any number of ";
        }
      }

      // static list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( GetParameterList().Begin()),
          param_itr_end( GetParameterList().End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      } // for

      // write the template help
      OSTREAM << "  Followed by " << size_description << " of ";
      m_TemplateParameter->WriteHelp( OSTREAM, 1);

      //end
      return OSTREAM;
    }

    //! @brief writes the user provided commandline to a stream
    //! @param OSTREAM stream to write to
    //! @return the stream written to
    std::ostream &FlagStaticAndDynamic::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "-" << m_Name << ( m_WasSetInCommandline ? " set" : " not set") << '\n';

      // write WriteUserCommand for every parameter in the list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()), param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        OSTREAM << "   ";
        ( *param_itr)->WriteUserCommand( OSTREAM);
      } // for

      // end
      return OSTREAM;
    }

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagStaticAndDynamic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        OSTREAM << '-' << m_Name << ' ';
        size_t parameter_counter( 0);
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator
            param_itr( m_ParameterList.Begin()),
            param_itr_end( m_ParameterList.End());
          param_itr != param_itr_end && parameter_counter < n_required_parameters;
          ++param_itr, ++parameter_counter
        )
        {
          OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
        }
        while( parameter_counter < n_required_parameters)
        {
          OSTREAM << '<' << m_TemplateParameter->GetName() << "> ";
        }
      }
      return OSTREAM;
    }

    //! @brief write Flag with Params to ostream
    //! @param OSTREAM - the stream to write to
    //! @return the stream after you wrote to it
    std::ostream &FlagStaticAndDynamic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TemplateParameter,   OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read Flag with params from istream
    //! @param ISTREAM - the stream to read from
    //! @return the stream after you read from it
    std::istream &FlagStaticAndDynamic::Read( std::istream &ISTREAM)
    {
      // write members
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);
      io::Serialize::Read( m_MinNumberParameters, ISTREAM);
      io::Serialize::Read( m_MaxNumberParameters, ISTREAM);
      io::Serialize::Read( m_TemplateParameter,   ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_guesser.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_csi_substructure.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace command
  {
    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &Guesser::GetTypeName( const MismatchType &TYPE)
    {
      static const std::string s_descriptors[ s_NumberTypes + 1] =
      {
        "CaseOrSpace",
        "DefinedAlias",
        "FirstLetters",
        "Suffix",
        "SuffixCaseOrSpace",
        "Stems",
        "ReorderedWords",
        "ReorderedStems",
        "Explict",
        "StrongAbbreviation",
        "WeakAbbreviation",
        "SomeWordsMatch",
        "SomeStemsMatch",
        GetStaticClassName< MismatchType>()
      };
      return s_descriptors[ size_t( TYPE)];
    }
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief get the single instance of this class
    //! @return the single instance of this class
    const Guesser &Guesser::GetDefaultGuesser()
    {
      // single instance of this class, initialized with the desired replacements
      static Guesser s_instance
      (
        Guesser()
        // suffixes that should never be changed
        .ProtectSuffix( "ss")     // ss is already a depluralized suffix
        .ProtectSuffix( "is")     // e.g. this
        .ProtectSuffix( "us")     // e.g. thus
        .ProtectSuffix( "se")     // e.g. these, portugese, geese
        .ProtectSuffix( "ite")    // e.g. favorite
        // suffixes that should be changed (order is irrelevant here)
        .RegisterSuffix( "s", "")        // Rods -> Rod
        .RegisterSuffix( "ed", "")       // Checked -> Check
        .RegisterSuffix( "ies", "y")     // Properties -> Property
        .RegisterSuffix( "ied", "y")     // Tried -> Try
        .RegisterSuffix( "oes", "o")     // potatoes -> potato
        .RegisterSuffix( "xes", "x")     // foxes -> fox
        .RegisterSuffix( "ves", "fe")    // lives -> life
        .RegisterSuffix( "lves", "lf")   // wolves -> wolf, calves -> calf
        .RegisterSuffix( "rices", "rix") // matrices -> matrix
        .RegisterSuffix( "sses", "ss")   // caresses  ->  caress
        .RegisterSuffix( "shes", "sh")   // bashes  ->  bash
        .RegisterSuffix( "ches", "ch")   // patches  ->  patch
        .RegisterSuffix( "xes", "x")     // boxes  ->  box
        .RegisterSuffix( "zzes", "zz")   // buzzes  ->  buzzes
        // These are usually suffixes, though there are many exceptions that a more sophisticated word stemming
        // algorithm would need to handle.  As the name of the class suggests, this is all heuristics anyway, so
        // if some app needs more precision, please create a separate word stemming class instead.
        .RegisterSuffix( "er", "")
        .RegisterSuffix( "ced", "ce")
        .RegisterSuffix( "sed", "se")
        .RegisterSuffix( "hed", "he")
        .RegisterSuffix( "ared", "are")
        .RegisterSuffix( "ered", "ere")
        .RegisterSuffix( "ored", "ore")
        .RegisterSuffix( "aned", "ane")
        .RegisterSuffix( "ened", "ene")
        .RegisterSuffix( "ined", "ine")
        .RegisterSuffix( "oned", "one")
        .RegisterSuffix( "ured", "ure")
        .RegisterSuffix( "ers", "")
        .RegisterSuffix( "or", "")
        .RegisterSuffix( "ing", "")
        .RegisterSuffix( "ors", "")
        .RegisterSuffix( "ion", "")
        .RegisterSuffix( "est", "")
        .RegisterSuffix( "able", "")
        .RegisterSuffix( "ible", "")
        .RegisterSuffix( "iest", "y")
        .RegisterSuffix( "ize", "")
        .RegisterSuffix( "ment", "")
        .RegisterSuffix( "ments", "")
        .RegisterSuffix( "able", "")
        .RegisterSuffix( "ful", "")
        .RegisterSuffix( "ions", "")
        .RegisterSuffix( "ness", "")
        .RegisterSuffix( "ly", "")
        .RegisterSuffix( "mic", "m")
        .RegisterSuffix( "nessly", "")
        .RegisterSuffix( "ity", "e")
        .RegisterSuffix( "ative", "ate")
        .RegisterSuffix( "ysis", "yze")
      );

      // ultimately word stemming is probably not optimal to determine user intent
      // An ideal approach would use alignments, with a score that decreases by distance from the start of the word

      return s_instance;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Guess what the user intended
    //! @param ARGUMENT the argument that was passed into the command line
    //! @param EXPECTED known, valid arguments
    //! @return a map from MismatchType to the strings that fell under that mismatch type
    storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> > Guesser::Guess
    (
      const std::string &ARGUMENT,
      const storage::Vector< std::string> &EXPECTED
    ) const
    {
      // create an object for the return value
      storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> > match_type_and_closest_matches;
      Guesser::TypeEnum &match_type( match_type_and_closest_matches.First());
      storage::Vector< std::string> &closest_matches( match_type_and_closest_matches.Second());

      // Normalize argument
      const std::string normalized_arg( NormalizeWord( ARGUMENT));

      // a vector for holding the normalized expected version of each word, to avoid recalculating
      storage::Vector< std::string> normalized_expected;
      const size_t expected_size( EXPECTED.GetSize());
      normalized_expected.AllocateMemory( expected_size);
      {
        // Check for CaseOrSpace case
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string normalized_option( NormalizeWord( *itr));
          if( normalized_option == normalized_arg)
          {
            // just differed in cases or spaces
            closest_matches.PushBack( *itr);
          }
          else if( closest_matches.IsEmpty())
          {
            // cache the normalized string for further use
            normalized_expected.PushBack( normalized_option);
          }
        }
        if( !closest_matches.IsEmpty())
        {
          // options differing in only case or spaces were found, so return them
          match_type = e_CaseOrSpace;
          return match_type_and_closest_matches;
        }
      }

      // Check for DefinedAlias case
      storage::Map< std::string, std::string>::const_iterator itr( m_Aliases.Find( normalized_arg));
      if( itr != m_Aliases.End() && EXPECTED.Find( itr->second) < expected_size)
      {
        // a defined alias was found and is also present in the expected strings, so return it
        closest_matches.PushBack( itr->second);
        match_type = e_DefinedAlias;
        return match_type_and_closest_matches;
      }

      // check for first letters case
      {
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          // check whether the given word is a prefix of any options
          if( GetMatchingPrefixLength( *itr_norm, normalized_arg) == normalized_arg.size())
          {
            closest_matches.PushBack( *itr);
          }
        }

        // check for first letters case
        if( !closest_matches.IsEmpty())
        {
          // users first letters all matched the beginning of something, so return that
          match_type = e_FirstLetters;
          return match_type_and_closest_matches;
        }
      }

      // Suffix case
      // check for matches that just differ in suffix by using the word stemming algorithm
      {
        const std::string arg_stem( GetStem( ARGUMENT));
        const std::string normalized_arg_stem( NormalizeWord( arg_stem));
        match_type = e_SuffixCaseSpace;
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          // if any case sensitive matches have been found, only check for case-sensitive matching strings
          if( match_type == e_Suffix)
          {
            if( arg_stem == GetStem( *itr))
            {
              closest_matches.PushBack( *itr);
            }
          }
          else
          {
            // get the expected stem; possibly with case removed
            const std::string expected_stem( GetStem( *itr));

            // standardize case
            const std::string normalized_expected_stem( NormalizeWord( expected_stem));

            if( normalized_arg_stem == normalized_expected_stem)
            {
              // case insensitive match, check for case sensitive match
              if( arg_stem == expected_stem)
              {
                // remove inferior case insensitive matches
                match_type = e_Suffix;
                closest_matches.Reset();
              }
              closest_matches.PushBack( *itr);
            }
          }
        }

        // return matches that differed in suffix and possibly case/space
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      const storage::Vector< std::string> arg_split_into_words( SplitIntoWords( ARGUMENT));
      const storage::Vector< std::string> arg_split_into_stems( GetStems( arg_split_into_words));

      // cache each expected option split into words and stems for later use
      storage::Vector< storage::Vector< std::string> > expected_split_into_words;
      storage::Vector< storage::Vector< std::string> > expected_split_into_stems;
      expected_split_into_words.AllocateMemory( expected_size);
      expected_split_into_stems.AllocateMemory( expected_size);

      // stems case
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          expected_split_into_words.PushBack( SplitIntoWords( *itr));
          expected_split_into_stems.PushBack( GetStems( expected_split_into_words.LastElement()));
          if( arg_split_into_stems == expected_split_into_stems.LastElement())
          {
            closest_matches.PushBack( *itr);
          }
        }

        // check for same stems case
        if( !closest_matches.IsEmpty())
        {
          // words had the same stem set, so return the given words
          match_type = e_Stems;
          return match_type_and_closest_matches;
        }
      }

      // check for reordered words or stems
      {
        match_type = e_ReorderedStems;
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr_split_words( expected_split_into_words.Begin()), itr_split_stems( expected_split_into_stems.Begin());
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_split_words, ++itr_split_stems
        )
        {
          if( arg_split_into_words.GetSize() != itr_split_words->GetSize())
          {
            // skip items with wrong #s of words
            continue;
          }
          if( graph::CSISubstructure::IsContainedIn( arg_split_into_words, *itr_split_words))
          {
            if( match_type == e_ReorderedStems)
            {
              // remove the inferior reordered stems matches
              closest_matches.Reset();
              match_type = e_ReorderedWords;
            }
            // found this set of words as a subset of one of the options
            closest_matches.PushBack( *itr);
          }
          else if( match_type == e_ReorderedStems)
          {
            // no reordered words have been found, look for reordered stems instead
            if( graph::CSISubstructure::IsContainedIn( arg_split_into_stems, *itr_split_stems))
            {
              closest_matches.PushBack( *itr);
            }
          }
        }

        // check whether any reordered words or stems were found
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for explicit cases, where the ARGUMENT is a fully written out version of the option
      {
        size_t best_abbreviation_gaps( FullSplitIntoWords( ARGUMENT).GetSize());
        match_type = e_Explicit;
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr
        )
        {
          // determine the number of gaps, considering the words to be sequence alignments
          const size_t gaps_size( MatchingSequenceWithGaps( *itr, ARGUMENT));
          if( gaps_size <= best_abbreviation_gaps)
          {
            if( gaps_size < best_abbreviation_gaps)
            {
              // discard inferior matches
              best_abbreviation_gaps = gaps_size;
              closest_matches.Reset();
            }
            closest_matches.PushBack( *itr);
          }
        }
        // if case/space sensitive abbreviations were found, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for case/space-sensitive abbreviations; prefer case sensitive abbreviations
      {
        // threshold; if more than this # of non-consecutive gaps are found, then do not consider it a valid abbreviation
        size_t best_abbreviation_gaps( FullSplitIntoWords( ARGUMENT).GetSize());
        match_type = e_WeakAbbreviation;
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          // determine the number of gaps, considering the words to be sequence alignments
          const size_t gaps_size_strong( MatchingSequenceWithGaps( ARGUMENT, *itr));
          if( util::IsDefined( gaps_size_strong) && match_type == e_WeakAbbreviation)
          {
            // case sensitive abbreviation matched, prior matches were case insensitive, which is inferior
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size_strong;

            match_type = e_StrongAbbreviation;
            closest_matches.PushBack( *itr);
          }
          else if( gaps_size_strong <= best_abbreviation_gaps)
          {
            if( gaps_size_strong < best_abbreviation_gaps)
            {
              // fewer gaps than seen before, remove existing inferior matches
              closest_matches.Reset();

              // update best abbreviation gaps and add *itr as the first match
              best_abbreviation_gaps = gaps_size_strong;
            }

            // add another match
            closest_matches.PushBack( *itr);
          }

          // skip weak abbreviations, if strong abbreviations have already been found
          if( match_type == e_StrongAbbreviation)
          {
            continue;
          }

          // test gaps size, ignoring case
          const size_t gaps_size_weak( MatchingSequenceWithGaps( normalized_arg, *itr_norm));

          // ignore if gaps size is unreasonable
          if( gaps_size_weak < best_abbreviation_gaps)
          {
            // fewer gaps than seen before, remove existing inferior matches
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size_weak;
            closest_matches.PushBack( *itr);
          }
          else if( gaps_size_weak == best_abbreviation_gaps)
          {
            // add another match
            closest_matches.PushBack( *itr);
          }
        }
        // if case/space sensitive abbreviations were found, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // check for matching abbreviations, regardless of case
      {
        size_t best_abbreviation_gaps( util::GetUndefined< size_t>());
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( EXPECTED.Begin()), itr_norm( normalized_expected.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_norm
        )
        {
          const size_t gaps_size( MatchingSequenceWithGaps( normalized_arg, *itr_norm));
          if( gaps_size < best_abbreviation_gaps)
          {
            // fewer gaps than seen before, remove existing inferior matches
            closest_matches.Reset();

            // update best abbreviation gaps and add *itr as the first match
            best_abbreviation_gaps = gaps_size;
            closest_matches.PushBack( *itr);
          }
          else if( util::IsDefined( gaps_size) && gaps_size == best_abbreviation_gaps)
          {
            // add another match
            closest_matches.PushBack( *itr);
          }
        }
        if( !closest_matches.IsEmpty())
        {
          match_type = e_WeakAbbreviation;
          return match_type_and_closest_matches;
        }
      }

      // check for maximum overlapping words and stems
      // If there are more overlapping stems than words, return overlapping stems, otherwise,
      // return the overlapping words
      {
        // track the max overlap of words or stems; start at 1 so to avoid finding all cases where nothing
        // matches initially
        size_t max_overlap( 1);

        // track the type (words or stems) with the greatest overlap; in case of a tie, prefer words since
        // they consitute a better matching criteria
        match_type = e_SomeStemsMatch;

        // iterate over words and split words
        storage::Vector< storage::Vector< std::string> >::iterator
          itr_option_words( expected_split_into_words.Begin()), itr_option_stems( expected_split_into_stems.Begin());
        for
        (
          storage::Vector< std::string>::const_iterator itr( EXPECTED.Begin()), itr_end( EXPECTED.End());
          itr != itr_end;
          ++itr, ++itr_option_words, ++itr_option_stems
        )
        {
          // compute # of equal words
          const size_t word_overlap( graph::CSISubstructure::GetOverlap( arg_split_into_words, *itr_option_words));

          // if there were at least as many overlapping words as were seen before
          if( word_overlap >= max_overlap)
          {
            // # of overlapping words > than previously seen # of overlapping words -> Remove existing (inferior) matches
            if( word_overlap > max_overlap || match_type == e_SomeStemsMatch)
            {
              // type == e_SomeStemsMatch catches when there are as many overlapping words as the best # of overlapping stems,
              // Because overlapping words are better than overlapping stems, remove the existing inferior matches
              closest_matches.Reset();
              max_overlap = word_overlap;
              match_type = e_SomeWordsMatch;
            }

            // add the new match
            closest_matches.PushBack( *itr);
          }

          // now check for overlap of stems
          const size_t stem_overlap( graph::CSISubstructure::GetOverlap( arg_split_into_stems, *itr_option_stems));
          if( stem_overlap >= max_overlap)
          {
            if( stem_overlap > max_overlap)
            {
              // more matching stems than ever seen, remove all existing inferior matches
              closest_matches.Reset();
              max_overlap = stem_overlap;
              match_type = e_SomeStemsMatch;
            }
            if( match_type == e_SomeStemsMatch)
            {
              // add additional stem-based match, unless equivalent overlap has already been found with word-based
              // overlap
              closest_matches.PushBack( *itr);
            }
          }
        }
        // if there were any matches based on common stems or words, return them
        if( !closest_matches.IsEmpty())
        {
          return match_type_and_closest_matches;
        }
      }

      // no good matching mechanism, just return all the words
      match_type = s_NumberTypes;
      closest_matches = EXPECTED;
      return match_type_and_closest_matches;
    }

    //! @brief Guess what the user intended and write those guesses to STREAM
    //! @param ARGUMENT the argument that was passed into the command line
    //! @param EXPECTED known, valid arguments
    //! @param STREAM stream to write out the guesses to
    //! @param ARGUMENT_TYPE_NAME conceptual type of the argument
    void Guesser::WriteGuesses
    (
      const std::string &ARGUMENT,
      const storage::Vector< std::string> &EXPECTED,
      std::ostream &STREAM,
      const std::string &ARGUMENT_TYPE_NAME
    ) const
    {
      const storage::Pair< Guesser::TypeEnum, storage::Vector< std::string> >
        guess_type_guesses( Guess( ARGUMENT, EXPECTED));

      // write out the argument type name, if one was given, or parameter otherwise
      if( ARGUMENT_TYPE_NAME.empty())
      {
        STREAM << "Given parameter \"";
      }
      else
      {
        STREAM << "Given " << ARGUMENT_TYPE_NAME << " \"";
      }
      STREAM << ARGUMENT << "\" is unknown";

      // if no remotely good matches were found, print the full list
      if( guess_type_guesses.First() == Guesser::s_NumberTypes)
      {
        STREAM << ", here are the allowed values {";
        if( guess_type_guesses.Second().GetSize())
        {
          std::copy( EXPECTED.Begin(), EXPECTED.End() - 1, std::ostream_iterator< std::string>( STREAM, ", "));
          STREAM << EXPECTED.LastElement();
        }
        STREAM << "}";
      }
      else if( guess_type_guesses.Second().GetSize() == size_t( 1))
      {
        // write out the matches
        STREAM << ", did you mean \"" << guess_type_guesses.Second()( 0) << "\"";
      }
      else
      {
        STREAM << ", perhaps you meant one of these: { ";
        std::copy
        (
          guess_type_guesses.Second().Begin(),
          guess_type_guesses.Second().End() - 1,
          std::ostream_iterator< std::string>( STREAM, ", ")
        );
        STREAM << guess_type_guesses.Second().LastElement();
        STREAM << "}";
      }
      STREAM << '\n';
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Registers the affixes for the prefix and suffix list
    //! @param AFFIX to be normalized
    //! @param NORMALIZED affix to create stem word
    //! @return a reference to the Guesser modified
    Guesser &Guesser::RegisterSuffix( const std::string &AFFIX, const std::string &NORMALIZED)
    {
      if( AFFIX.size() < NORMALIZED.size())
      {
        return RegisterSuffix( NORMALIZED, AFFIX);
      }

      const std::string normalized_affix( NormalizeWord( AFFIX));
      const std::string normalized_word( NormalizeWord( NORMALIZED));
      if( m_Suffixes.GetSize() <= AFFIX.size())
      {
        m_Suffixes.Resize( AFFIX.size() + 1);
      }

      // add the normalized word in at the normal place
      m_Suffixes( normalized_affix.size())[ normalized_affix] = normalized_word;

      // make the suffix map to itself to prevent non-specific matching in non-trivial cases
      if( !normalized_word.empty() && normalized_affix != normalized_word)
      {
        m_Suffixes( normalized_word.size())[ normalized_word] = normalized_word;
      }

      return *this;
    }

    //! @brief Protect a suffix from being changed if it is the longest matching suffix
    //! @param SUFFIX to be normalized
    //! @return a reference to the Guesser modified
    Guesser &Guesser::ProtectSuffix( const std::string &SUFFIX)
    {
      return RegisterSuffix( SUFFIX, SUFFIX);
    }

    //! @brief Registers the aliases to be used
    //! @param ORIGINAL word
    //! @param ALIAS for original word
    //! @return a reference to the Guesser modified
    Guesser &Guesser::RegisterAlias( const std::string &ORIGINAL, const std::string &ALIAS)
    {
      m_Aliases[ NormalizeWord( ORIGINAL)] = ALIAS;
      return *this;
    }

    //! @brief Applies suffix rules to get down to the word stem
    //! @param WORD to be stemmed
    //! @return word with the suffix standardized
    std::string Guesser::GetStem( const std::string &WORD) const
    {
      // determine the word size; this may change so do not make it constant
      size_t word_size( WORD.size());

      // skip words with size <= 2, they cannot have a prefix or affix
      if( word_size <= 2)
      {
        return WORD;
      }

      // do not try to match affixes that would leave only a single letter on the word, e.g., if the affix is our,
      // do not match four.
      const size_t max_affix_size( word_size - 2);

      if( !m_Suffixes.IsEmpty())
      {
        // start by fixing the suffix, if possible
        // start with the largest suffixes
        for( size_t suffix_size( std::min( max_affix_size, m_Suffixes.GetSize() - 1)); suffix_size > 0; --suffix_size)
        {
          // get the suffix of the given size from the given word
          const std::string actual_suffix( WORD.substr( word_size - suffix_size));

          // check for it in the appropriate map
          storage::Map< std::string, std::string>::const_iterator itr( m_Suffixes( suffix_size).Find( actual_suffix));
          if( itr != m_Suffixes( suffix_size).End())
          {
            // replace the given suffix
            if( itr->first != itr->second)
            {
              return WORD.substr( 0, word_size - suffix_size) + itr->second;
            }
            else
            {
              return WORD;
            }
          }
        }
      }
      return WORD;
    }

    //! @brief GetStem for an entire vector
    //! @param WORDS words to be stemmed
    //! @return the word stems
    storage::Vector< std::string> Guesser::GetStems( const storage::Vector< std::string> &WORDS) const
    {
      storage::Vector< std::string> stems;
      stems.AllocateMemory( WORDS.GetSize());
      // call standardize suffix on each element in the vector
      for
      (
        storage::Vector< std::string>::const_iterator itr( WORDS.Begin()), itr_end( WORDS.End());
        itr != itr_end;
        ++itr
      )
      {
        stems.PushBack( Guesser::GetStem( *itr));
      }
      return stems;
    }

    // @brief Splits input string into words based on white space characters, capitalization, and underscores
    // @param STRING to be split into words
    // @return Vector of split words
    storage::Vector< std::string> Guesser::SplitIntoWords( const std::string &STRING) const
    {
      // test whether the trimmed string appears to be an object data label, if so, just take the first word
      size_t label_delimiter_pos( STRING.find_first_of( "(="));
      if( label_delimiter_pos && label_delimiter_pos != std::string::npos)
      {
        // if an object data label delimiter was found, just split the first word into words
        return SplitIntoWords( STRING.substr( 0, label_delimiter_pos));
      }

      // get the trimmed string
      const std::string trimmed( util::TrimString( STRING));

      storage::Vector< std::string> split;

      // check whether spaces or underscores are present
      if( trimmed.find_first_of( "_ \n\t\r") != std::string::npos)
      {
        // just split the string as normal
        split = util::SplitString( STRING, "_ \n\t\r");
      }
      else
      {
        // split based on capitalization
        // a word is defined by the regex: [A-Z0-9]+[a-z0-9]*
        // also, a word may not be composed entirely of numbers
        size_t pos( 0), size( trimmed.size());
        while( pos < size)
        {
          // skip non-alphanumeric characters
          if( !isalnum( trimmed[ pos]))
          {
            ++pos;
            continue;
          }

          // add caps and digits
          const size_t original_pos( pos);
          ++pos;
          while( pos < size && ( isupper( trimmed[ pos]) || isdigit( trimmed[ pos])))
          {
            ++pos;
          }

          size_t number_letters( pos - original_pos);
          // if there was only one initial letter, then we have a normal word
          // 2+ initial letters means an abbreviation sequence, followed by the end, or a normal word
          // If it was an abbreviation sequence followed by a word, then we need to backup one
          if( number_letters > size_t( 1) && pos < size)
          {
            --pos;
          }
          else
          {
            // typical case, the initial letter sequence is complete
            while( pos < size && ( islower( trimmed[ pos]) || isdigit( trimmed[ pos])))
            {
              ++pos;
            }

            // if the last letter was a digit, and we are not at the end, then that digit belongs with the next word
            if( pos < size && isdigit( trimmed[ pos - 1]))
            {
              --pos;
            }
          }

          // safety catch, if there are no letters, increment pos
          if( pos == original_pos)
          {
            ++pos;
          }
          split.PushBack( trimmed.substr( original_pos, pos - original_pos));
        }
      }
      std::for_each( split.Begin(), split.End(), ToLower);
      return split;
    }

    // @brief Normalizes a string to remove capitalization, spacing, and underscores
    // @param WORD to be normalized
    // @return normalized word
    std::string Guesser::NormalizeWord( const std::string &WORD) const
    {
      // form a new string without the spaces, _, and with all lower case letters
      std::string normalized_word;
      normalized_word.reserve( WORD.size());
      for( size_t i( 0), size( WORD.size()); i < size; ++i)
      {
        if( !isspace( WORD[ i]) && WORD[ i] != '_')
        {
          normalized_word += tolower( WORD[ i]);
        }
      }
      return normalized_word;
    }

    //! @brief Splits input string into words based on white space characters, capitalization, and underscores
    //! @param STRING to be split into words, letters of acronyms as separate entities
    //! @return Vector of split words, letters of acronyms as separate entities
    storage::Vector< std::string> Guesser::FullSplitIntoWords( const std::string &STRING)
    {
      // test whether the trimmed string appears to be an object data label, if so, just take the first word
      size_t label_delimiter_pos( STRING.find_first_of( "(="));
      if( label_delimiter_pos != std::string::npos)
      {
        // if an object data label delimiter was found, just split the first word into words
        return FullSplitIntoWords( STRING.substr( 0, label_delimiter_pos));
      }

      // first, get the trimmed string
      const std::string trimmed( util::TrimString( STRING));

      storage::Vector< std::string> split;

      // check whether spaces or underscores are present
      if( trimmed.find_first_of( "_ \n\t\r") != std::string::npos)
      {
        // just split the string as normal
        const storage::Vector< std::string> partially_split( util::SplitString( STRING, "_ \n\t\r"));

        for( size_t j( 0), size( partially_split.GetSize()); j < size; ++j)
        {
          if( !partially_split( j).empty())
          {
            // split this word into multiple words, if necessary
            split.Append( FullSplitIntoWords( partially_split( j)));
          }
        }
      }
      else
      {
        // split based on capitalization
        // a word is defined by the regex: [A-Z0-9]+[a-z0-9]*
        // also, a word may not be composed entirely of numbers
        size_t pos( 0), size( trimmed.size());
        while( pos < size)
        {
          // add caps and digits
          const size_t original_pos( pos);

          // the first letter is automatically a new word
          ++pos;

          // typical case, the initial letter sequence is complete
          while( pos < size && islower( trimmed[ pos]))
          {
            ++pos;
          }

          std::string new_string( trimmed.substr( original_pos, pos - original_pos));
          ToLower( new_string);

          // no further splitting is possible on single characters
          if( !IsObviousAcroynm( new_string))
          {
            split.PushBack( new_string);
          }
          else
          {
            // fully split the word
            std::string letter( size_t( 1), ' ');
            // words where all letters are either consonants or vowels are generally acronyms (aa, ss, pdb, sdf)
            for( size_t i( 0), new_size( new_string.size()); i < new_size; ++i)
            {
              letter[ 0] = new_string[ i];
              split.PushBack( letter);
            }
          }
        }
      }
      return split;
    }

    // @brief Just changes the capitalization of a word to lower
    // @param WORD to be lower-cased
    void Guesser::ToLower( std::string &WORD)
    {
      // form a new string without the spaces, _, and with all lower case letters
      for( size_t i( 0), size( WORD.size()); i < size; ++i)
      {
        WORD[ i] = tolower( WORD[ i]);
      }
    }

    //! @brief Get the number of characters that match at the start of the string
    //! @param A, B the two strings of interest
    //! @return the number of letters at the start of the strings that match
    size_t Guesser::GetMatchingPrefixLength( const std::string &A, const std::string &B)
    {
      if( A.size() <= B.size())
      {
        return std::mismatch( A.begin(), A.end(), B.begin()).first - A.begin();
      }
      return std::mismatch( B.begin(), B.end(), A.begin()).first - B.begin();
    }

    //! @brief Test whether one string is the abbreviation of another
    //! @param A, B the two strings of interest
    //! @return true if A is an abbreviation prefix of B or B is an abbreviation of A
    bool Guesser::IsAbbreviation( const std::string &A, const std::string &B)
    {
      if( A.size() <= B.size())
      {
        return std::mismatch( A.begin(), A.end(), B.begin()).first == A.end();
      }
      return std::mismatch( B.begin(), B.end(), A.begin()).first == B.end();
    }

    //! @brief Determine whether a given string is obviously an acronym
    //! @param A the string of interest
    //! @return true if A is > 1 letter, and contains all vowels, or all consonants
    //! @note if y's are present, they count as a vowel only if all the rest of the word is consonants
    //! @note if the A is just a sequences of Y's, it counts as an acronym
    bool Guesser::IsObviousAcroynm( const std::string &A)
    {
      static const std::string s_vowels( "aeiou");

      // test for 1 letter words
      if( A.size() <= size_t( 1))
      {
        return false;
      }

      size_t vowels_count( 0);
      // match the first letter independent of case
      if( s_vowels.find( tolower( A[ 0])) != std::string::npos)
      {
        ++vowels_count;
      }

      // count vowels for the rest of the string
      for( std::string::const_iterator itr( A.begin() + 1), itr_end( A.end()); itr != itr_end; ++itr)
      {
        if( s_vowels.find( *itr) != std::string::npos)
        {
          ++vowels_count;
        }
      }

      // if there appeared to be no vowels, look for a y
      // ignore the first character, since Y is always a consonant when it is the first letter
      if( !vowels_count && A.find( 'y', 1) != std::string::npos)
      {
        // make sure that the word was not all y's
        if( tolower( A[ 0]) != 'y' || A.find_first_not_of( 'y', 1) != std::string::npos)
        {
          ++vowels_count;
        }
      }

      // return true if there were no vowels or no consonants
      return !vowels_count || vowels_count == A.size();
    }

    //! @brief Test whether one string is an exact abbreviation / substring of another, with arbitrary number of gaps
    //! @param SUBSTR the posited substring
    //! @param STRING the string to test whether SUBSTR is a subsequence of
    //! @param START position to start at
    //! @return # of non-consecutive gaps, or undefined if there is no matching sub-sequence
    size_t Guesser::MatchingSequenceWithGaps
    (
      const std::string &WITH_GAPS,
      const std::string &STRING,
      const size_t &START
    )
    {
      // check that the alignment could exist
      if( STRING.size() < START + WITH_GAPS.size() || STRING.empty())
      {
        return util::GetUndefined< size_t>();
      }
      else if( WITH_GAPS.empty())
      {
        return 0;
      }

      // find the first matching letter
      const size_t first_match( STRING.find( WITH_GAPS[ 0], START));

      // no match: return
      if( first_match == std::string::npos || STRING.size() < first_match + WITH_GAPS.size())
      {
        return util::GetUndefined< size_t>();
      }

      bool last_was_gap( false);
      size_t number_non_consecutive_gaps( 0);

      // track # of available gaps left in string
      size_t remaining_gaps( STRING.size() - first_match - WITH_GAPS.size());

      // get iterators onto both strings
      std::string::const_iterator
        itr_gap( WITH_GAPS.begin()), itr_gap_end( WITH_GAPS.end()), itr_full( STRING.begin() + first_match);

      while( itr_gap != itr_gap_end)
      {
        if( *itr_gap != *itr_full)
        {
          // test whether any gaps remain
          if( !remaining_gaps)
          {
            return util::GetUndefined< size_t>();
          }
          // do not allow gaps to traverse upper letters except at the beginning and end of the string
          else if( itr_gap != WITH_GAPS.begin() && isupper( *itr_full))
          {
            // look for a match further on in the string
            return MatchingSequenceWithGaps( WITH_GAPS, STRING, first_match + 1);
          }
          --remaining_gaps;
          if( !last_was_gap)
          {
            ++number_non_consecutive_gaps;
          }
          last_was_gap = true;
        }
        else
        {
          ++itr_gap;
          last_was_gap = false;
        }
        ++itr_full;
      }

      if( !number_non_consecutive_gaps)
      {
        return 0;
      }
      return std::min( number_non_consecutive_gaps, MatchingSequenceWithGaps( WITH_GAPS, STRING, first_match + 1));
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_default.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Parameter::s_Instance
    (
      GetObjectInstances().AddInstance( new Parameter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Parameter::Parameter() :
      m_Name(),
      m_Description(),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from description string, parameter check function
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from description string, parameter check function
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterCheckInterface &PARAMETER_CHECK
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( PARAMETER_CHECK.Clone()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from name, description string and default parameter string
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param DEFAULT_PARAMETER Default string value of the parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const std::string &DEFAULT_PARAMETER
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_Parameter( DEFAULT_PARAMETER),
      m_WasSetInCommandline( false),
      m_DefaultGiven( true),
      m_DefaultParameter( DEFAULT_PARAMETER)
    {
    }

    //! @brief construct from description string, parameter check, and default parameter string
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
    //! @param DEFAULT_PARAMETER Default string value of the parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterCheckInterface &PARAMETER_CHECK,
      const std::string &DEFAULT_PARAMETER
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( PARAMETER_CHECK.Clone()),
      m_Parameter( DEFAULT_PARAMETER),
      m_WasSetInCommandline( false),
      m_DefaultGiven( true),
      m_DefaultParameter( DEFAULT_PARAMETER)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new Parameter object
    Parameter *Parameter::Clone() const
    {
      return new Parameter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return std::string - the class name as const ref std::string
    const std::string &Parameter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns Name
    //! @return std::string - name
    const std::string &Parameter::GetName() const
    {
      return m_Name;
    }

    //! @brief set the name
    //! @param NAME name for the parameter
    void Parameter::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

    //! @brief returns description
    //! @return std::string-description
    const std::string &Parameter::GetDescription() const
    {
      return m_Description;
    }

    //! @brief checks if PARAMETER is allowed and returns
    //! @param PARAMETER the parameter to set
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return if the parameter was set (parameters are not set if they are not allowed parameters)
    bool Parameter::SetParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      // get the result of the parameter check
      io::ValidationResult result( m_ParameterCheck->Check( PARAMETER, m_Name, ERROR_STREAM));

      if( !result.IsInvalid())
      {
        // clean the parameter, if it was
        m_Parameter = PARAMETER;
        m_WasSetInCommandline = true;
      }
      return result.IsAllowed();
    } // Parameter::SetParameter

    //! @brief returns m_Parameter
    //! @return value of the parameter
    const std::string &Parameter::GetValue() const
    {
      return m_Parameter;
    }

    //! @brief returns m_DefaultParameter
    //! @return value of the default parameter
    const std::string &Parameter::GetDefaultValue() const
    {
      return m_DefaultParameter;
    }

    //! @brief set the default parameter
    //! @param PARAMETER default value for the parameter
    void Parameter::SetDefaultParameter( const std::string &PARAMETER)
    {
      if( !m_WasSetInCommandline)
      {
        m_DefaultParameter = PARAMETER;
        m_Parameter        = PARAMETER;
        m_DefaultGiven     = true;
      }
    }

    //! @brief returns if parameter was set in the command line
    //! @return m_WasSetInCommandline
    bool Parameter::GetWasSetInCommandLine() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief return whether default value was given
    //! @return bool - m_DefaultGiven
    bool Parameter::GetWasDefaultGiven() const
    {
      return m_DefaultGiven;
    }

    //! @brief acces to the parameter check
    //! @return the ParameterCheck
    const ParameterCheckInterface &Parameter::GetParameterCheck() const
    {
      return *m_ParameterCheck;
    }

    //! @brief set the parameter check
    //! @param CHECK the new parameter check
    void Parameter::SetParameterCheck( const ParameterCheckInterface &CHECK)
    {
      m_ParameterCheck = util::ShPtr< ParameterCheckInterface>( CHECK.Clone());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. the object behind the ParameterCheckInterface allows it
    //! @param PARAMETER the parameter to check
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool Parameter::IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const
    {
      return m_ParameterCheck->IsAllowedParameter( PARAMETER, m_Name, ERROR_STREAM);
    }

    //! @brief reset into original state
    void Parameter::Reset()
    {
      m_WasSetInCommandline = false;

      // default was given, set the parameter to it
      if( m_DefaultGiven)
      {
        m_Parameter = m_DefaultParameter;

        // clean the default parameter.  This is necessary after calling reset since the clean function may have
        // runtime dependencies
        m_ParameterCheck->Clean( m_Parameter);
      }
      // otherwise clear
      else
      {
        m_Parameter.clear();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT number of indentations
    //! @return the given stream to which the help was written to
    std::ostream &Parameter::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::FixedLineWidthWriter writer;
      writer.SetBclIndent( INDENT);

      if( !m_Name.empty())
      {
        // write name
        writer << '<' << m_Name << "> ";
      }

      if( !m_Description.empty())
      {
        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          writer.SetIndent( writer.GetLinePosition());
          writer << m_Description << ", ";
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 4);
          writer << m_Description << ", ";
        }
      }
      else
      {
        writer.AddIndent( 2);
      }

      // if default parameter was given, also write this
      if( !m_WasSetInCommandline && m_DefaultGiven)
      {
        if( m_Parameter.empty())
        {
          writer << "optional, ";
        }
        else
        {
          writer.WriteOnOneLine( "default: \"" + m_Parameter + "\", ");
        }
      }
      writer.PopIndent();

      // write ParameterCheck
      std::ostringstream output;
      m_ParameterCheck->WriteHelp( output, 0);
      const std::string parameter_check( util::TrimString( output.str()));

      // parameter checks that have new lines look best with an extra autoindent 4
      if( parameter_check.find( '\n') != std::string::npos)
      {
        writer.SetAutoIndentExtra( 4);
      }

      // write parameter check position;
      if( parameter_check.size() < writer.GetRemainingSpaceOnLine())
      {
        // short parameter check, write it on the same line
        writer.SetIndent( writer.GetLinePosition());
      }
      else if( parameter_check.size() < writer.GetEffectiveLineWidth() - 2)
      {
        // just indent by 2 and write the parameter check on the next line
        // it is too long to fit in a fully-indented line
        writer.AddIndent( 2);
        writer.NewLine();
      }
      else if( parameter_check.size() < 2 * writer.GetRemainingSpaceOnLine())
      {
        // parameter check will look best with its own block
        writer.SetIndent( writer.GetLinePosition());
      }
      else
      {
        // long parameter check, just write with a standard indent of 2, starting on the same
        // line as the description
        writer.AddIndent( 2);
      }

      // write out the parameter check
      writer << parameter_check;
      const std::string full_help( writer.String());
      // write the whole description, skipping terminal commas, spaces, and new lines
      OSTREAM << full_help.substr( 0, full_help.find_last_not_of( ", \n") + 1) << '\n';

      // end
      return OSTREAM;
    } // Parameter::WriteHelp

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &Parameter::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "<" << m_Name << "> " << m_Parameter;

      if( !m_WasSetInCommandline && m_DefaultGiven)
      {
        OSTREAM << " (default)";
      }

      OSTREAM << '\n';

      // end
      return OSTREAM;
    } // Parameter::WriteUserCommand

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Parameter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Name               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description        , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterCheck     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Parameter          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultGiven       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultParameter   , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Parameter::Read( std::istream &ISTREAM)
    {
      //read members
      io::Serialize::Read( m_Name               , ISTREAM);
      io::Serialize::Read( m_Description        , ISTREAM);
      io::Serialize::Read( m_ParameterCheck     , ISTREAM);
      io::Serialize::Read( m_Parameter          , ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_DefaultGiven       , ISTREAM);
      io::Serialize::Read( m_DefaultParameter   , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief StringVectorFromFilenameParameter gives all the strings in a file as a Vector
    //! @param PARAMETER the name of the file
    //! @return return a Vector which has all of the strings contained within the file denoted by "FLAG"
    storage::Vector< std::string> StringVectorFromFilenameParameter( const ParameterInterface &PARAMETER)
    {
      // initialize write and read stream object
      io::IFStream read;

      // open pdb list file
      read.open( PARAMETER.GetValue().c_str());

      // make sure the file opened
      BCL_Assert( read.is_open(), "unable to open file: " + PARAMETER.GetValue());

      // create "string_vector" initialize with vector of strings returned by function util::StringListFromIStream
      const storage::Vector< std::string> string_vector( util::StringListFromIStream( read));

      // close and clear "read"
      io::File::CloseClearFStream( read);

      // return Vector filled with the strings in the pdb list
      return string_vector;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_allowed.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckAllowed::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckAllowed())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckAllowed::ParameterCheckAllowed()
    {
    }

    //! @brief construct from allowed parameters
    //! @param ALLOWED_PARAMETERS vector of allowed strings
    ParameterCheckAllowed::ParameterCheckAllowed( const storage::Vector< std::string> &ALLOWED_PARAMETERS) :
      m_AllowedParameters( ALLOWED_PARAMETERS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckAllowed object
    ParameterCheckAllowed *ParameterCheckAllowed::Clone() const
    {
      return new ParameterCheckAllowed( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckAllowed::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckAllowed::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // iterate over all parameters in m_AllowedParameters
      if( m_AllowedParameters.Find( PARAMETER) < m_AllowedParameters.GetSize())
      {
        return true;
      }

      // try fuzzy matching
      Guesser::GetDefaultGuesser().WriteGuesses( PARAMETER, m_AllowedParameters, ERROR_STREAM, PARAMETER_NAME);

      // return false because parameter was not found in list
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckAllowed::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write a short description
      OSTREAM << "allowed values: ";

      // write all entries in the list
      WriteList( OSTREAM);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckAllowed::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AllowedParameters, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckAllowed::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AllowedParameters, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param OSTREAM - the stream to which you should write
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckAllowed::WriteList( std::ostream &OSTREAM) const
    {
      //write opening parentheses
      OSTREAM << "{";

      // copy to output stream
      if( !m_AllowedParameters.IsEmpty())
      {
        std::copy
        (
          m_AllowedParameters.Begin(),
          m_AllowedParameters.End() - 1,
          std::ostream_iterator< std::string>( OSTREAM, ", ")
        );
        OSTREAM << m_AllowedParameters.LastElement();
      }

      // write closing parentheses
      OSTREAM << "}";

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_allowed_non_const.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace command
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from allowed parameters
    //! @param ALLOWED_PARAMETERS vector of allowed strings
    ParameterCheckAllowedNonConst::ParameterCheckAllowedNonConst( const storage::Vector< std::string> &ALLOWED_PARAMETERS) :
      m_AllowedParameters( &ALLOWED_PARAMETERS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckAllowedNonConst object
    ParameterCheckAllowedNonConst *ParameterCheckAllowedNonConst::Clone() const
    {
      return new ParameterCheckAllowedNonConst( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckAllowedNonConst::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckAllowedNonConst::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // iterate over all parameters in m_AllowedParameters
      if( m_AllowedParameters->Find( PARAMETER) < m_AllowedParameters->GetSize())
      {
        return true;
      }

      // try fuzzy matching
      Guesser::GetDefaultGuesser().WriteGuesses( PARAMETER, *m_AllowedParameters, ERROR_STREAM, PARAMETER_NAME);

      // return false because parameter was not found in list
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckAllowedNonConst::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write a short description
      OSTREAM << "allowed values: ";

      // write all entries in the list
      WriteList( OSTREAM);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckAllowedNonConst::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AllowedParameters, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckAllowedNonConst::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AllowedParameters, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param OSTREAM - the stream to which you should write
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckAllowedNonConst::WriteList( std::ostream &OSTREAM) const
    {
      //write opening parentheses
      OSTREAM << "{";

      // copy to output stream
      if( !m_AllowedParameters->IsEmpty())
      {
        std::copy
        (
          m_AllowedParameters->Begin(),
          m_AllowedParameters->End() - 1,
          std::ostream_iterator< std::string>( OSTREAM, ", ")
        );
        OSTREAM << m_AllowedParameters->LastElement();
      }

      // write closing parentheses
      OSTREAM << "}";

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_default.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_validation_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckDefault::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckDefault())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckDefault::ParameterCheckDefault()
    {
    }

    //! @brief virtual copy constructor
    ParameterCheckDefault *ParameterCheckDefault::Clone() const
    {
      return new ParameterCheckDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief check the parameter, returning CheckResult
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns CheckResult
    io::ValidationResult ParameterCheckDefault::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      if( PARAMETER == io::ValidationResult::GetHelpString() && CommandState::GetInMainCommandLineParsing())
      {
        io::ValidationResult val( io::e_Help);
        // in this case, no help can be given, so the application will need to handle the request for help
        CommandState::GetWasHelpGiven() = false;
      }
      return io::ValidationResult( true);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the outstream to write to
    //! @param INDENT the amount to indent each new line after the first
    //! @return return the outstream after you write to it
    std::ostream &ParameterCheckDefault::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckDefault::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_extension.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckExtension::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckExtension())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckExtension::ParameterCheckExtension()
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSION file extension
    ParameterCheckExtension::ParameterCheckExtension( const std::string &EXTENSION) :
      m_Extension( EXTENSION)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to ParameterCheckExtension object
    ParameterCheckExtension *ParameterCheckExtension::Clone() const
    {
      return new ParameterCheckExtension( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckExtension::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. if filename has allowed extension
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckExtension::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // determine the position of the Extension in the given PARAMETER
      const size_t pos( PARAMETER.rfind( m_Extension, PARAMETER.length()));

      // if extension was found and is at the end of PARAMETER return true
      if( util::IsDefined( pos) && ( pos + m_Extension.length() == PARAMETER.length()))
      {
        return true;
      }

      // otherwise it was not found or is in the wrong position, return false
      ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" does not have the extension: " << m_Extension << '\n';
      return false;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM output stream to be written to
    //! @param INDENT the amount to indent each new line after the first
    //! @return ostream which was written to
    std::ostream &ParameterCheckExtension::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write extension
      OSTREAM << "file pattern <*" << m_Extension << "> ";

      //end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckExtension::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Extension, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckExtension::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Extension, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_extensions_file_existence.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckExtensionsFileExistence::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckExtensionsFileExistence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence()
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSIONS vector of file extensions
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence( const storage::Vector< std::string> &EXTENSIONS) :
      m_Extensions( EXTENSIONS)
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSION single file extension
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence( const std::string &EXTENSION) :
      m_Extensions()
    {
      m_Extensions.PushBack( EXTENSION);
    }

    //! @brief virtual copy constructor
    //! @return pointer to ParameterCheckExtensionsFileExistence object
    ParameterCheckExtensionsFileExistence *ParameterCheckExtensionsFileExistence::Clone() const
    {
      return new ParameterCheckExtensionsFileExistence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckExtensionsFileExistence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is an allowed (i.e. it exists & has the correct extension
    //! @param PARAMETER - the parameter in question
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM - write errors here
    //! @return true if the parameter is allowed, false otherwise
    bool ParameterCheckExtensionsFileExistence::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // iterate over given extensions to check if those files exists
      for
      (
        std::vector< std::string>::const_iterator itr( m_Extensions.Begin()), itr_end( m_Extensions.End());
        itr != itr_end; ++itr
      )
      {
        const std::string name( PARAMETER + *itr);

        // if this file was not found
        if( !io::DirectoryEntry( name).DoesExist())
        {
          ERROR_STREAM << "File \"" << name << "\" does not exist!" << '\n';
          return false;
        } // if
      } // for

      return true;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM output stream to be written to
    //! @param INDENT the amount to indent each new line after the first
    //! @return ostream which was written to
    std::ostream &ParameterCheckExtensionsFileExistence::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write list of extensions
      OSTREAM << "file pattern <";
      WriteListOfExtensions( OSTREAM);
      OSTREAM << ">";

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param &OSTREAM - the stream to which you write
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckExtensionsFileExistence::WriteListOfExtensions( std::ostream &OSTREAM) const
    {
      // write opening parantheses
      OSTREAM << "\"";

      // write all entries in the list
      for
      (
        std::vector< std::string>::const_iterator itr( m_Extensions.Begin()), itr_end( m_Extensions.End());
        itr != itr_end;
        ++itr
      )
      {
        OSTREAM << "*" << *itr;
        if( itr != m_Extensions.End() - 1)
        {
          OSTREAM << ", ";
        }
      }
      // write closing parentheses
      OSTREAM << "\"";

      // end
      return OSTREAM;
    } // WriteListOfExtensions

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckExtensionsFileExistence::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Extensions, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckExtensionsFileExistence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Extensions, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_file_existence.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckFileExistence::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckFileExistence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new ParameterCheckFileExistence
    ParameterCheckFileExistence *ParameterCheckFileExistence::Clone() const
    {
      return new ParameterCheckFileExistence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckFileExistence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns true if PARAMETER is an allowed parameters. i.e. true if the file exists
    //! @param PARAMETER - the parameter in question. This is the filename.
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM - the ostream to which errors are written
    //! @return bool - returns true if the filename given by PARAMETER exists, false otherwise
    bool ParameterCheckFileExistence::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // if this file was not found
      if( !io::DirectoryEntry( PARAMETER).DoesExist())
      {
        ERROR_STREAM << "File \"" << PARAMETER << "\" does not exist!" << '\n';
        return false;
      }

      // file was found
      return true;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckFileExistence::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write list of allowed parameters. i.e. any filename
      OSTREAM << "any existent file";

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckFileExistence::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckFileExistence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_file_in_search_path.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "app/bcl_app_apps.h"
#include "io/bcl_io_directory_entry.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckFileInSearchPath::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckFileInSearchPath( storage::Vector< std::string>()))
    );

    //! @brief get a string variable that, when used with this class, will be replaced by the directory containing the bcl executable
    //! @return a string whose value will be replaced with the directory containing the bcl executable
    const std::string &ParameterCheckFileInSearchPath::GetExecutableDirectoryVariable()
    {
      static const std::string s_bcl_var( "$BCL_EXE$");
      return s_bcl_var;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new ParameterCheckFileInSearchPath
    ParameterCheckFileInSearchPath *ParameterCheckFileInSearchPath::Clone() const
    {
      return new ParameterCheckFileInSearchPath( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief default constructor
    //! @param TARGET_NAME target file or directory name to search for in the default search path
    //! use the path the bcl was started from or up to 3 directories above it, before resorting to
    //! the directory on blue, or the current directory.  This allows both installed bcl executables to
    //! be able to find the models directory without relying on the often relative and often moved install path.
    //! It also allows modification and use of files in checkouts without needing to set this flag
    //! @param LAST_RESORT last file or directory to try before giving up; this is usually a meilerlab-specific path;
    //! The last resort filename allows MeilerLab users to copy locally-built executables around on the filesystem
    //! without needing to copy all associated files, unless they have in fact been modified
    //! @param TYPE the actual type of object to search for; be it a file or directory
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const std::string &TARGET_NAME,
      const std::string &LAST_RESORT,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths
      (
        GetVersion().IsLicense()
        // release/licensed copies are installed, so everything is usually together with the executable
        // but use the putative install directory as a backup -- its rarely correct because users can move the bcl
        // folder post-installation; moreover its just whatever the user types as the install directory, which is usually
        // a relative path.
        // The absolute last resort in all cases is the current directory
        ? storage::Vector< std::string>::Create
          (
            GetExecutableDirectoryVariable(),
            LAST_RESORT,
            ""
          )
        // for locally-built (unlicensed) versions, look in the super-directories of the bcl-exectuable 1st
        // because the build-folder often has the same folders as the working directory, so e.g. we don't want to get
        // build/linux64_release/histogram, rather, we want build/linux64_release/../../histogram == histogram
        // indeed, in this case, the folder the bcl is located in should be the last one searched due to the likelyhood
        // of finding similarly named folders in the build directory
        : storage::Vector< std::string>::Create
          (
            GetExecutableDirectoryVariable() + "/../..",
            GetExecutableDirectoryVariable() + "/../../..",
            GetExecutableDirectoryVariable() + "/..",
            GetExecutableDirectoryVariable() + "/bcl", // needed for ligand-design gui
            LAST_RESORT,
            "",
            GetExecutableDirectoryVariable()
          )
      )
    {
      SetDirectories();
      for( storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End()); itr != itr_end; ++itr)
      {
        if( *itr != LAST_RESORT)
        {
          *itr += TARGET_NAME;
        }
      }
      TYPE == io::Directory::e_File ? SetFiles() : SetDirectories();
    }

    //! @brief constructor
    //! @param TARGET_NAME target file or directory name to search for in the given paths
    //! @param PATHS ordered paths to search for the target file or directory
    //! @param TYPE the actual type of object to search for; be it a file or directory
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const std::string &TARGET_NAME,
      const storage::Vector< std::string> &PATHS,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths( PATHS)
    {
      if( !TARGET_NAME.empty())
      {
        // force all existing path components to be directories
        SetDirectories();

        // add or remove / depending on the desired entry type TYPE
        std::string to_append( TARGET_NAME);
        if( TARGET_NAME[ TARGET_NAME.size() - 1] == '/')
        {
          // remove / if type was supposed to be file
          if( TYPE == io::Directory::e_File)
          {
            to_append.erase( to_append.size() - 1);
          }
        }
        else if( TYPE == io::Directory::e_Dir)
        {
          // append / if it was supposed to be a directory
          to_append += '/';
        }

        // add the target name to each path
        for
        (
          storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
          itr != itr_end;
          ++itr
        )
        {
          *itr += to_append;
        }
      }
      else
      {
        GetSearchType() == io::Directory::e_File ? SetFiles() : SetDirectories();
      }
    }

    //! @brief constructor
    //! @param PATHS ordered paths to search for a valid entry of the given type
    //! @param TYPE the actual type of object to search for; be it a file or directory
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const storage::Vector< std::string> &PATHS,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths( PATHS)
    {
      TYPE == io::Directory::e_File ? SetFiles() : SetDirectories();
    }

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckFileInSearchPath::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns true if PARAMETER is an allowed parameters. i.e. true if the file exists
    //! @param PARAMETER - the parameter in question. This is the filename.
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM - the ostream to which errors are written
    //! @return bool - returns true if the filename given by PARAMETER exists, false otherwise
    bool ParameterCheckFileInSearchPath::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {

      // this function is only called when the user has provided a path in which case the default search path is ignored
      io::DirectoryEntry entry( PARAMETER);
      if( entry.DoesExist() && entry.GetType() == GetSearchType())
      {
        return true;
      }
      ERROR_STREAM << PARAMETER << " is not a valid "
                   << ( GetSearchType() == io::Directory::e_Dir ? "directory" : "file")
                   << "for parameter " << PARAMETER_NAME;
      return false;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckFileInSearchPath::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // determine default
      const std::string default_path( FindFile( ""));
      // write list of allowed parameters. i.e. any filename
      OSTREAM << "any " << ( GetSearchType() == io::Directory::e_Dir ? "directory" : "file")
              << ", if not provided, defaults to " << default_path
              << "; search path is: {" << GetCandidatePathsString() << "}\n";

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckFileInSearchPath::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Paths, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckFileInSearchPath::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Paths, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief find the file in the search paths; if it exists, otherwise, return an empty string
    //! @param PARAMETER the first path to check; may be a default, blank, or may be something entered over the command line
    //! @return the first occurrance of the target file/dir in the search paths; empty if it is not in the search path
    std::string ParameterCheckFileInSearchPath::FindFile( const std::string &PARAMETER) const
    {
      if( PARAMETER.size())
      {
        io::DirectoryEntry entry( PARAMETER);
        if( entry.DoesExist() && entry.GetType() == GetSearchType())
        {
          return PARAMETER;
        }
      }

      // create a string replacer to replace the special string with the bcl executable path
      const util::StringReplacement replacer
      (
        util::StringReplacement::e_Any,
        GetExecutableDirectoryVariable(),
        io::DirectoryEntry( app::Apps::GetExecutablePath()).GetDirectory().GetPath()
      );

      // get the type of entry being searched for
      io::Directory::EntryType entry_type( GetSearchType());
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string candidate_path( *itr);
        replacer.ReplaceEachIn( candidate_path);
        io::DirectoryEntry candidate( candidate_path);
        if( candidate.DoesExist() && candidate.GetType() == entry_type)
        {
          return SimplifyPath( candidate.GetFullName() + ( entry_type == io::Directory::e_File ? "" : "/"));
        }
      }

      return std::string();
    }

    //! @brief clean the parameter value; that is, make it usable by code that is independent of the parameter check
    //! Some checks may include or-type functionality, such as checking for a file among a list of paths
    //! The clean function updates the parameter value with the correct path
    //! @param PARAMETER the parameter to clean
    void ParameterCheckFileInSearchPath::Clean( std::string &PARAMETER) const
    {
      std::string file( FindFile( PARAMETER));
      if( !file.empty())
      {
        PARAMETER = file;
      }
    }

    //! @brief return the type of directory entry being searched for
    io::Directory::EntryType ParameterCheckFileInSearchPath::GetSearchType() const
    {
      return
        !m_Paths.IsEmpty() && m_Paths.FirstElement()[ m_Paths.FirstElement().size() - 1] == '/'
        ? io::Directory::e_Dir
        : io::Directory::e_File;
    }

    //! @brief get a string representation of all paths that will be searched
    std::string ParameterCheckFileInSearchPath::GetCandidatePathsString() const
    {
      // create a string replacer to replace the special string with the bcl executable path
      const util::StringReplacement replacer
      (
        util::StringReplacement::e_Any,
        GetExecutableDirectoryVariable(),
        io::DirectoryEntry( app::Apps::GetExecutablePath()).GetDirectory().GetPath()
      );

      std::string candidates;
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string candidate_path( *itr);
        replacer.ReplaceEachIn( candidate_path);
        candidates += SimplifyPath( candidate_path);
        if( itr + 1 != itr_end)
        {
          // add colon for next possibility
          candidates += ":";
        }
      }
      return candidates;
    }

    //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
    void ParameterCheckFileInSearchPath::SetDirectories()
    {
      // force all paths in m_Paths to have a '/' at the end
      for
      (
        storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string &str( *itr);
        if( !str.empty() && str[ str.size() - 1] != '/')
        {
          str += '/';
        }
      }
    }

    //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
    void ParameterCheckFileInSearchPath::SetFiles()
    {
      // force all paths in m_Paths to not have a '/' at the end
      for
      (
        storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string &str( *itr);
        while( !str.empty() && str[ str.size() - 1] == '/')
        {
          str.erase( str.size() - 1);
        }
      }
    }

    //! @brief function to simplify a given path by removing the pattern X/../
    std::string ParameterCheckFileInSearchPath::SimplifyPath( const std::string &PATH)
    {
      // remove duplicate slashes, also ./ as a directory or in the middle of the path
      util::StringReplacement replace_double_slash( util::StringReplacement::e_Any, "//", "/");
      util::StringReplacement replace_slash_dot_slash( util::StringReplacement::e_Any, "/./", "/");
      std::string temp_string( PATH);
      replace_double_slash.ReplaceAllIn( temp_string);
      replace_slash_dot_slash.ReplaceAllIn( temp_string);
      while( util::StartsWith( temp_string, "./"))
      {
        temp_string.erase( 0, 2);
      }

      if( temp_string.size() < size_t( 4))
      {
        return temp_string;
      }

      // walk through the string, find anytime /../ appears
      size_t last_found_pos( 1);
      size_t prev_directory_pos( temp_string.find( "/../", last_found_pos));
      while( prev_directory_pos < temp_string.size())
      {
        size_t start_pos( temp_string.rfind( '/', prev_directory_pos - 1));
        if( start_pos == std::string::npos)
        {
          start_pos = 0;
        }
        else
        {
          ++start_pos;
        }
        // remove the unnecessary part of the string
        temp_string.erase( start_pos, prev_directory_pos + 4 - start_pos);
        last_found_pos = start_pos ? start_pos - 1 : start_pos;
        prev_directory_pos = temp_string.find( "/../", last_found_pos);
      }
      return temp_string;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_interface.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_validation_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    //! @brief process the parameter, returning ValidationResult, and cleaning the parameter if it was valid
    //! @param PARAMETER the parameter to process
    //! @param PARAMETER_NAME the name of the parameter being processed
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns ValidationResult
    io::ValidationResult ParameterCheckInterface::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // by default, check whether the parameter is equal to the global help string, if so, write the parameter's help
      // and return true
      if( PARAMETER == io::ValidationResult::GetHelpString())
      {
        ERROR_STREAM << "Requirements for <" << PARAMETER_NAME << ">:\n";
        WriteHelp( ERROR_STREAM);
        return io::e_Help;
      }

      // check whether the parameter was allowed
      return io::ValidationResult( IsAllowedParameter( PARAMETER, PARAMETER_NAME, ERROR_STREAM));
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_or.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_validation_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckOr::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckOr())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckOr::ParameterCheckOr()
    {
    }

    //! @brief construct from two checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B)
        )
      )
    {
    }

    //! @brief construct from three checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    //! @param CHECK_C the third parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B,
      const ParameterCheckInterface &CHECK_C
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B),
          util::CloneToShPtr( CHECK_C)
        )
      )
    {
    }

    //! @brief construct from four checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    //! @param CHECK_C the third parameter check that will be considered
    //! @param CHECK_D the fourth parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B,
      const ParameterCheckInterface &CHECK_C,
      const ParameterCheckInterface &CHECK_D
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B),
          util::CloneToShPtr( CHECK_C),
          util::CloneToShPtr( CHECK_D)
        )
      )
    {
    }

    //! @brief construct from vector of allowed checks
    //! @param CHECKS vector of parameter checks, any of which can be satisfied for IsAllowedParameter to return true
    ParameterCheckOr::ParameterCheckOr( const util::ShPtrVector< ParameterCheckInterface> &CHECKS) :
      m_Choices( CHECKS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckOr object
    ParameterCheckOr *ParameterCheckOr::Clone() const
    {
      return new ParameterCheckOr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckOr::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckOr::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // create a temporary error stream; whose results will be written to ERROR_STREAM only if none of
      // the parameter checks is satisfied
      std::ostringstream temp_err_stream;

      // iterate over all possible choices
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin()), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        std::ostringstream stream;

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        // if this parameter check was valid, return true; errors for other checks will be disregarded
        if( parameter_check.IsAllowedParameter( PARAMETER, PARAMETER_NAME, stream))
        {
          return true;
        }

        temp_err_stream << stream.str();
      }

      ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is none of the following:\n";
      WriteList( ERROR_STREAM, 1);
      ERROR_STREAM << "Errors for each check: " << temp_err_stream.str() << '\n';

      // return false because parameter was not found in list
      return false;
    }

    //! @brief check the parameter, returning CheckResult
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns CheckResult
    io::ValidationResult ParameterCheckOr::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // create a temporary error stream; whose results will be written to ERROR_STREAM only if none of
      // the parameter checks is satisfied
      std::ostringstream temp_err_stream;

      io::ValidationResult overall_result( io::e_Allowed);

      // iterate over all possible choices
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin()), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        std::ostringstream stream;

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        // get the check result
        io::ValidationResult result( parameter_check.Check( PARAMETER, PARAMETER_NAME, stream));

        // handle trivial result
        if( result.IsAllowed())
        {
          return result;
        }

        // handle help request
        if( result.IsHelp())
        {
          if( overall_result.IsInvalid())
          {
            // clear the current temp error stream
            temp_err_stream.str( std::string());
          }
          overall_result = result;
        }
        else if( overall_result.IsAllowed()) // first invalid parameter
        {
          overall_result = result;
        }

        if( result == overall_result)
        {
          temp_err_stream << stream.str();
        }
      }

      if( overall_result.IsInvalid())
      {
        ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is none of the following:\n";
        WriteList( ERROR_STREAM, 1);
        ERROR_STREAM << "Errors for each check: " << temp_err_stream.str() << '\n';
      }
      else // help
      {
        ERROR_STREAM << temp_err_stream.str() << '\n';
      }
      // return false because parameter was not found in list
      return overall_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckOr::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write a short description
      OSTREAM << "Any of the following:\n";

      // write all entries in the list
      WriteList( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckOr::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Choices, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckOr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Choices, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT number of indentations
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckOr::WriteList( std::ostream &OSTREAM, const size_t &INDENT) const
    {
      // no choices, nothing to print
      if( m_Choices.IsEmpty())
      {
        return OSTREAM;
      }

      std::ostringstream temp_ostream;

      // write help for the first parameter check
      m_Choices.FirstElement()->WriteHelp( temp_ostream, INDENT);

      // get the string from the output stream
      std::string temp_string( temp_ostream.str());
      if( !temp_string.empty() && temp_string[ temp_string.length() - 1] != '\n')
      {
        temp_string += '\n';
      }
      OSTREAM << temp_string;

      // iterate over remaining options
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin() + 1), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        // write Or to separate the next set of options
        io::Serialize::InsertIndent( OSTREAM, INDENT) << "Or ";

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        temp_ostream.str( std::string());
        // write help for this parameter check
        parameter_check.WriteHelp( temp_ostream, INDENT);
        temp_string = temp_ostream.str();
        if( !temp_string.empty() && temp_string[ temp_string.length() - 1] != '\n')
        {
          temp_string += '\n';
        }

        OSTREAM << temp_string;
      }

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_ranged.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    template class BCL_API ParameterCheckRanged< double>;

  } // namespace command
} // namespace bcl
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
#include "command/bcl_command_parameter_check_serializable.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckSerializable::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckSerializable())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckSerializable::ParameterCheckSerializable()
    {
    }

    //! @brief construct from serializable object
    //! @param OBJECT vector of allowed strings
    ParameterCheckSerializable::ParameterCheckSerializable( const io::SerializationInterface &OBJECT) :
      m_SerializableInstance( OBJECT.Clone())
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckSerializable object
    ParameterCheckSerializable *ParameterCheckSerializable::Clone() const
    {
      return new ParameterCheckSerializable( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckSerializable::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckSerializable::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // try to read the parameter into an object data label
      util::ObjectDataLabel internal;
      if( !internal.TryAssign( PARAMETER, ERROR_STREAM))
      {
        return false;
      }

      // if no object was given, the parameter is allowed
      if( !m_SerializableInstance.IsDefined())
      {
        return true;
      }

      // copy the pointer / clone the object, to allow for read attempt
      util::OwnPtr< io::SerializationInterface> fresh_instance( m_SerializableInstance);

      // try to read the variable; return true on success
      return fresh_instance->TryRead( internal, ERROR_STREAM);
    }

    //! @brief check the parameter, returning CheckResult
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns CheckResult
    io::ValidationResult ParameterCheckSerializable::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // by default, check whether the parameter is equal to the global help string, if so, write the parameter's help
      // and return e_Help
      if( PARAMETER == io::ValidationResult::GetHelpString())
      {
        ERROR_STREAM << "Help for <" << PARAMETER_NAME << ">:\n";
        WriteHelp( ERROR_STREAM);
        return io::e_Help;
      }
      // create an object data label out of the parameter
      util::ObjectDataLabel label;
      if( !label.TryAssign( PARAMETER, ERROR_STREAM))
      {
        return io::e_Invalid;
      }

      // copy the pointer / clone the object, to allow for read attempt
      util::OwnPtr< io::SerializationInterface> fresh_instance( m_SerializableInstance);
      return fresh_instance->ValidateRead( label, ERROR_STREAM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckSerializable::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // call write help o the internal object
      m_SerializableInstance->WriteHelp( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckSerializable::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SerializableInstance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckSerializable::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SerializableInstance, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
