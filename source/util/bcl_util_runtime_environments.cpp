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
#include "util/bcl_util_runtime_environments.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_runtime_environment_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! command line flag to be used to set RuntimeEnvironment over the command line
    ShPtr< command::FlagInterface> &RuntimeEnvironments::GetFlagRuntimeEnvironment()
    {
      static ShPtr< command::FlagInterface> s_runtime_environment_flag
      (
        new command::FlagStatic
        (
          "runtime_environment",
          "change the runtime environment this executable runs in",
          command::Parameter
          (
            "environment",
            "choice of environment",
            command::ParameterCheckEnumerate< RuntimeEnvironments>(),
            "Default"
          ),
          &RuntimeEnvironments::UpdateCurrentRuntimeEnvironmentFromCommandLineFlag
        )
      );

      return s_runtime_environment_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RuntimeEnvironments::RuntimeEnvironments() :
      Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>( false),
      e_Default( AddEnum( "Default", ShPtr< RuntimeEnvironmentInterface>( new RuntimeEnvironmentDefault()))),
      m_CurrentRuntimeEnvironment( *e_Default)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RuntimeEnvironments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Initialize the runtime environment from the command line flag
    void RuntimeEnvironments::UpdateCurrentRuntimeEnvironmentFromCommandLineFlag()
    {
      // name of runtime environment from commandline
      const std::string &re_name( GetFlagRuntimeEnvironment()->GetFirstParameter()->GetValue());

      // construct the environment enum from the name given in the commandline
      ShPtr< RuntimeEnvironmentInterface> new_runtime_environment( *RuntimeEnvironment( re_name));

      // the current environment is correct
      if( new_runtime_environment == GetEnums().m_CurrentRuntimeEnvironment)
      {
        return;
      }

      // update the current runtime environment to the commandline selected one
      GetEnums().m_CurrentRuntimeEnvironment = new_runtime_environment;

      // message to user
      BCL_MessageTop( "updating runtime environment to " + re_name);

      // initialize environment
      BCL_Assert
      (
        GetEnums().m_CurrentRuntimeEnvironment->Initialize(),
        "unable to initialize the runtime environment " + re_name
      );
    }

    //! @brief Initialize and return the current environment
    //! @return one and only reference to one of the environments
    const RuntimeEnvironmentInterface &RuntimeEnvironments::GetCurrentEnvironment() const
    {
      // return the current runtime environment
      return *m_CurrentRuntimeEnvironment;
    }

    //! @brief returns RuntimeEnvironments
    //! @return RuntimeEnvironments
    RuntimeEnvironments &GetRuntimeEnvironments()
    {
      return RuntimeEnvironments::GetEnums();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>;

  } // namespace util
} // namespace bcl
