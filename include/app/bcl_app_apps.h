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

#ifndef BCL_APP_APPS_H_
#define BCL_APP_APPS_H_

// include the namespace header
#include "bcl_app.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "bcl_app_groups.h"
#include "bcl_app_interface.h"
#include "bcl_app_interface_release.h"
#include "command/bcl_command_parameter.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Apps
    //! @brief Class collecting all applications
    //!
    //! @author woetzen
    //! @see @link example_apps.cpp @endlink
    //! @date Mar 24, 2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Apps :
      public util::Enumerate< util::ShPtr< Interface>, Apps>
    {
      friend class util::Enumerate< util::ShPtr< Interface>, Apps>;
      friend class GroupHandler;

    private:

    //////////
    // data //
    //////////

      storage::Vector< std::string> m_ReleaseApplicationsNames; //!< names of all release apps

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Apps();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief add new Interface derived Application to Apps
      //! @param APPLICATION the application to be enumerated
      //! @return the newly created enum
      ApplicationType &AddEnum( const Interface &APPLICATION);

      //! @brief add the app directly
      //! @param NAME           name of the current application
      //! @param SP_APPLICATION object to be enumerated
      ApplicationType &AddEnum( const std::string &NAME, const util::ShPtr< Interface> &SP_APPLICATION);

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return vector of application names
      //! @return StorageVector<std::string> with application name
      const storage::Vector< std::string> &GetApplicationNames() const;

      //! @brief return vector of release application names
      //! @return StorageVector<std::string> with release application name
      const storage::Vector< std::string> &GetReleaseApplicationNames() const;

      //! @brief return a parameter with all application strings
      //! @return command::Parameter with check that contains all possible application name strings
      util::ShPtr< command::ParameterInterface> &GetApplicationsParameter();

      //! @brief return the executable path
      //! @return the first argument to the application == the executable path
      static std::string &GetExecutablePath();

    ////////////////
    // operations //
    ////////////////

      //! @brief writes the list of enums
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the list was written to
      //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
      virtual std::ostream &WriteList( std::ostream &OSTREAM) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      static std::ostream &WriteGenericHelp( std::ostream &OSTREAM);

      //! @brief get the default command line
      //! If no arguments are given to an application, this file will be opened and its contents will be used as a command
      static const std::string &GetDefaultCommandLineArgumentsFile();

    }; // class Apps

    BCL_API Apps &GetApps();

  } // namespace app
} // namespace bcl

#endif // BCL_APP_APPS_H_
