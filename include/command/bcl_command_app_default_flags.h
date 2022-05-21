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

#ifndef BCL_COMMAND_APP_DEFAULT_FLAGS_H_
#define BCL_COMMAND_APP_DEFAULT_FLAGS_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_command_default_flag_types.h"
#include "bcl_command_flag_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AppDefaultFlags
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Aug 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AppDefaultFlags
    {
    public:

    ///////////
    // enums //
    ///////////

      //! @brief default flag that writes all options and exits the program
      //! @return FlagInterface containing the help flag
      static util::ShPtr< FlagInterface> &GetHelpFlag();

      //! default flag that writes readme information if one is available
      static util::ShPtr< FlagInterface> &GetReadMeFlag();

    private:

    //////////
    // data //
    //////////

      //! All flags
      util::ShPtrVector< FlagInterface> m_Flags;

      //! Flags, indexed by type
      storage::Vector< util::ShPtrVector< FlagInterface> > m_FlagsByType;

      //! Flag name to type
      storage::Map< std::string, FlagTypeEnum> m_FlagToType;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! Adds all flags that have a required order to the enum
      AppDefaultFlags();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the only instance of this class
      //! @return the only instance of this class
      static AppDefaultFlags &GetInstance();

      //! @brief get all the flags
      //! @return all the flags
      const util::ShPtrVector< FlagInterface> &GetAllFlags() const;

      //! @brief detect whether a given flag is a bcl app default flag
      //! @param FLAG_NAME the flag name to test
      //! @return type of the flag
      FlagType GetFlagType( const std::string &FLAG_NAME);

    ////////////////
    // operations //
    ////////////////

      //! @brief add a default flag, to the collection
      //! @param SH_PTR_FLAG shptr to flag, to be added to default commandline flags
      //! @param TYPE type of flag that will be added
      //! @return the flag that was added
      const util::ShPtr< FlagInterface> &AddDefaultFlag
      (
        const util::ShPtr< FlagInterface> &SH_PTR_FLAG,
        const FlagType &TYPE
      );

      //! @brief get the flags for a given type
      //! @param TYPE the type of flag of interest
      //! @return the flags of the requested type
      const util::ShPtrVector< FlagInterface> &GetDefaultFlagsOfType( const FlagType &TYPE) const;

      //! @brief add the bcl default flags of specified type to the Command object
      //! @param COMMAND reference to a commandline object
      //! @param TYPES types of command line flags to add.
      //! Note that e_BclRequired and e_AppGeneric will always be added
      void AddDefaultCommandlineFlags
      (
        Command &COMMAND,
        const storage::Set< FlagTypeEnum> &TYPES
      );

      //! @brief add the bcl default flags required to the Command object
      //! @param COMMAND reference to a commandline object
      //! Note that e_BclRequired and e_AppGeneric will always be added
      void AddRequiredCommandlineFlags( Command &COMMAND);

      //! @brief add the bcl default flags to the Command object
      //! @param COMMAND reference to a commandline object
      void AddDefaultCommandlineFlags( Command &COMMAND);

      //! @brief handle the bcl-required flags
      //! @param STATE code state; from command line
      //! @param ERR_STREAM stream for output of errors
      //! @return true on success
      bool HandleRequiredFlags( CommandState &STATE, std::ostream &ERR_STREAM);

      //! @brief handle all default flags contained in the given command flags
      //! @param STATE code state; from command line
      //! @param FLAGS the flags to consider
      //! @param ERR_STREAM stream for output of errors
      //! @return true on success
      bool HandleFlags( CommandState &STATE, util::ShPtrVector< FlagInterface> &FLAGS, std::ostream &ERR_STREAM);

    }; // class AppDefaultFlags

    //! @brief construct on access function for all AppDefaultFlags
    //! @return reference to only instances of AppDefaultFlags
    BCL_API
    AppDefaultFlags &GetAppDefaultFlags();

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_APP_DEFAULT_FLAGS_H_
