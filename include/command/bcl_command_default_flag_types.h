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

#ifndef BCL_COMMAND_DEFAULT_FLAG_TYPES_H_
#define BCL_COMMAND_DEFAULT_FLAG_TYPES_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_command_default_flag_types
    //! Contains the flag type enum
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Feb 21, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! all possible types of flags
    //! Flags will always be added to commands ordered by type
    enum FlagType
    {
      e_BclCore,      //!< All flags that are required for core capabilities: runtime, logger, message_level
      e_AppGeneric,   //!< All flags that are required to be declared for every application
      e_Io,           //!< File compression and alternatives
      e_Util,         //!< combine_sh_ptr_by_address, enums_files
      e_Random,       //!< Random seed
      e_Model,        //!< model_path
      e_Score,        //!< histogram path
      e_Db,           //!< mysql_databank
      e_Opencl,       //!< opencl
      e_Pthread,      //!< scheduler
      e_OpenBlas,     //!< OpenBlas
      e_AppSpecific,  //!< Application-specific flags
      s_NumberFlagTypes
    };

    //! @brief Type as string
    //! @param TYPE the flag type
    //! @return the string for the flag type
    const std::string &GetFlagTypeName( const FlagType &TYPE);

    //! @brief Type enum I/O helper
    typedef util::WrapperEnum< FlagType, &GetFlagTypeName, s_NumberFlagTypes> FlagTypeEnum;

    //! @brief get the default flag types to be included in all applications
    const storage::Set< FlagTypeEnum> &GetDefaultFlagTypes();

    //! @brief get the default flag types to be included in all applications
    //! @param FLAG_NAME name of the flag
    FlagType GetFlagTypeFromName( const std::string &FLAG_NAME);

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_DEFAULT_FLAG_TYPES_H_
