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

#ifndef EXAMPLE_APPLICATION_EXAMPLE_HELPER_H_
#define EXAMPLE_APPLICATION_EXAMPLE_HELPER_H_

// include forward header of this class

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @class ApplicationExampleHelper
  //! @brief helper class for examples of applications
  //! All application examples should begin by initializing an ApplicationExampleHelper
  //!
  //! @author mendenjl
  //! @date May 23, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ApplicationExampleHelper
  {

  private:

  //////////
  // data //
  //////////

    const app::ApplicationType     m_Application;             //!< Application enum that we are running
    util::ShPtr< app::Interface>   m_ApplicationInstance;     //!< ShPtr to instance of app run previously
    util::ShPtr< command::Command> m_Command;                 //!< Last command that was run
    util::Stopwatch                m_RunTimer;                //!< Application run duration (both last run time and total run time)
    storage::Map< std::string, storage::Vector< std::string> > m_FlagParameterMap; //!< Arguments to run the application with
    storage::Set< std::string>     m_AppDefaultFlagsSetByExample; //!< app default flags that were set by this example

    util::ShPtr< command::Command> m_AppDefaultFlags;                //!< command containing all app default flags
    storage::Set< std::string>     m_AppDefaultFlagsFromCommandline; //!< set of flags given on the command line

    //! AppDefaultFlags values; both those with defaults  automatically added to commands that are
    //! run, though the application example may change them for specific runs.  When destructed, the helper will reset
    //! them to the original values
    storage::Map< std::string, storage::Vector< std::string> > m_AppDefaultFlagValueMap;

    //! base path for all application example files
    static std::string s_ApplicationExampleFilePath;

  public:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from an application
    explicit ApplicationExampleHelper( const app::ApplicationType &APP);

    //! @brief destructor, resets the app default flags
    ~ApplicationExampleHelper();

  /////////////////
  // data access //
  /////////////////

    //! @brief get the application
    //! @return const access to the application that will be run internally
    const util::ShPtr< app::Interface> &GetApplication() const;

    //! @brief get the command object
    //! @return const access to the command
    const command::Command &GetCommand() const;

    //! @brief get the command line
    //! @return a string representation of the current command line
    std::string GetCurrentCommandLine() const;

    //! @brief get the stopwatch for this application example helper
    //! @return const access to the timer, which tells how long the last command took, and the total runtime of the app
    const util::Stopwatch &GetRunTimer() const;

    //! @brief get the arguments that will be passed to the application the next time it is run
    //! @return the arguments that will be passed to the application the next time it is run
    storage::Vector< std::string> GetArguments() const;

    //! @brief unset a particular flag
    void UnsetFlag( const std::string &FLAG);

    //! @brief reset the parameters
    void ResetParameters();

    //! @brief reset the flags and parameters
    void ResetFlagsAndParameters();

    //! @brief add a parameter, which may be a response file, to the list of arguments
    //! @param PARAMETER the parameter to add
    //! @return reference to this (so that multiple add parameters can be strung together)
    ApplicationExampleHelper &AddParameter( const std::string &PARAMETER);

    //! @brief add several parameters, which may be response files, to the list of arguments
    //! @param PARAMETERS the parameters to add
    //! @return reference to this (so that multiple add parameters can be strung together)
    ApplicationExampleHelper &AddParameters( const storage::Vector< std::string> &PARAMETERS);

    //! @brief add a flag to the list of arguments
    //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
    //! @return reference to this (so that multiple adds can be strung together)
    ApplicationExampleHelper &SetFlag( const std::string &FLAG_NAME);

    //! @brief add a flag with one value to the list of arguments
    //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
    //! @param FLAG_VALUE the value of the flag to add
    //! @return reference to this (so that multiple adds can be strung together)
    ApplicationExampleHelper &SetFlag( const std::string &FLAG_NAME, const std::string &FLAG_VALUE);

    //! @brief add a flag with one value to the list of arguments
    //! @param FLAG_NAME the flag name to add (e.g. remove_h, do not include the -)
    //! @param FLAG_VALUES the values of the flag to add
    //! @return reference to this (so that multiple adds can be strung together)
    ApplicationExampleHelper &SetFlag( const std::string &FLAG_NAME, const storage::Vector< std::string> &FLAG_VALUES);

    //! @brief add a parameter to an existing flag
    //! @param FLAG_NAME the flag name to add the parameter to
    //! @param FLAG_VALUE the parameter to add to the flag
    //! @return reference to this (so that multiple adds can be strung together)
    ApplicationExampleHelper &AddParameterToFlag( const std::string &FLAG_NAME, const std::string &FLAG_VALUE);

    //! @brief flag to change application example path
    //! @return flag to change application example path
    static util::ShPtr< command::FlagInterface> &GetApplicationExamplePathFlag();

    //! @brief set the application example file path
    //! @param PATH the new path to application example files
    static void SetApplicationExampleFilePath( const std::string &PATH);

  ////////////////
  // operations //
  ////////////////

    //! @brief test if a command line is valid
    //! @param PRINT_ERRORS whether to print errors to the screen
    //! @return true if the command line is valid for the application
    bool CheckCommandString( const bool &PRINT_ERRORS);

    //! @brief run a command with an application
    //! @return return value from main of the application
    int RunCommand();

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the application example input path
    //! @return the application example input path
    static const std::string &GetApplicationExampleInputPath();

    //! @brief get the input path for this application example
    //! @return the input path for this application example
    std::string GetThisApplicationExampleInputPath() const;

  private:

    //! @brief get the arguments that will be passed to the application the next time it is run
    //! @param INCLUDE_APP_DEFAULTS whether to include flags that in the AppDefaultFlagList
    //! @return the arguments that will be passed to the application the next time it is run
    storage::Vector< std::string> GetArguments( const bool &INCLUDE_DEFAULTS) const;

    //! @brief initialize m_AppDefaultFlagValueMap to the values and defaults given on the command line
    void InitializeAppDefaultFlagMap();

    //! @brief reset the default app flags to the values they had at the beginning of the example
    void ResetDefaultFlags();

    //! @brief reset AppDefaultFlags that will be changed by this example
    void ResetOverriddenDefaultFlags();

  }; // ApplicationExampleHelper

} // namespace bcl

#endif // EXAMPLE_APPLICATION_EXAMPLE_HELPER_H_
