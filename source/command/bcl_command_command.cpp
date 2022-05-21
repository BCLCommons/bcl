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
