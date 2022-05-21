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
