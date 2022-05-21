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
#include "util/bcl_util_loggers.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_logger_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! command line flag to be used to set Logger over the command line
    ShPtr< command::FlagInterface> &Loggers::GetFlagLogger()
    {
      static ShPtr< command::FlagInterface> s_logger_flag
      (
        new command::FlagStatic
        (
          "logger",
          "change the logger this executable uses",
          ShPtrVector< command::ParameterInterface>::Create
          (
            ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "",
                "",
                command::ParameterCheckEnumerate< Loggers>(),
                GetLoggers().e_Default.GetName()
              )
            ),
            ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "identifier",
                "define a logger identifier - for file, it is the filename to be opened",
                ""
              )
            )
          ),
          &Loggers::UpdateCurrentLoggerFromCommandLineFlag
        )
      );

      return s_logger_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Default constructor
    Loggers::Loggers() :
      Enumerate< ShPtr< LoggerInterface>, Loggers>( false),
      e_Default( AddEnum( "Default", ShPtr< LoggerInterface>( new LoggerDefault()))),
      m_CurrentLogger( new LoggerDefault()),
      m_CurrentLoggerEnum( e_Default)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Loggers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Initialize the logger from the command line flag
    void Loggers::UpdateCurrentLoggerFromCommandLineFlag()
    {
      const Logger new_logger( GetFlagLogger()->GetFirstParameter()->GetValue());
      const std::string &new_identifier( GetFlagLogger()->GetParameterList()( 1)->GetValue());

      Loggers &loggers( GetEnums());
      if( new_logger == loggers.m_CurrentLoggerEnum && new_identifier == loggers.m_CurrentLogger->GetLogIdentifier())
      {
        // no need to update the current logger
        return;
      }

      // switching the current logger
      BCL_MessageStd( "change Logger to: " + new_logger.GetName() + " with identifier: " + new_identifier);
      loggers.m_CurrentLogger = ShPtr< LoggerInterface>( ( *new_logger)->Clone());
      loggers.m_CurrentLoggerEnum = new_logger;

      // set the identifier of the current logger
      BCL_Assert
      (
        loggers.m_CurrentLogger->SetLogIdentifier( new_identifier),
        "unable to initialize logger: " + new_logger.GetName() + " with identifier: " + new_identifier
      );

      // switch was successful
      BCL_MessageStd( "Logger was changed to: " + new_logger.GetName() + " with identifier: " + new_identifier);
    }

    //! @brief Initialize and return the current logger
    //! @return one and only reference to one of the loggers
    LoggerInterface &Loggers::GetCurrentLogger()
    {
      // return the logger
      return *m_CurrentLogger;
    }

    //! @brief get enumerated list of Loggers
    Loggers &GetLoggers()
    {
      return Loggers::GetEnums();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< LoggerInterface>, Loggers>;

  } // namespace util
} // namespace bcl
