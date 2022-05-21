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

#ifndef BCL_APP_INTERFACE_H_
#define BCL_APP_INTERFACE_H_

// include the namespace header
#include "bcl_app.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Interface
    //! @brief for interfacing applications
    //! @remarks every main application needs to be derived from this interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Mar 23, 2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Interface :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      static const std::string &GetSupportEmailCommercial(); //!< e-mail address for commercial license support
      static const std::string &GetSupportEmailAcademic();   //!< e-mail address for academic license support

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Interface
      virtual Interface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief return the Command object
      virtual util::ShPtr< command::Command> InitializeCommand() const = 0;

      //! @brief the expiration interval from last compilation
      //! @return the time interval from the last compilation, the license is valid - by default, it is around 20 years
      virtual const util::Time &GetLicenseExpirationTime() const;

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      virtual std::string GetBCLScopedName() const;

      //! @brief return the original name of the application, as it was used in license files
      //! @return string for the bcl::commons name of that application
      //! This is necessary so that, if release application names change, licenses will continue to work
      virtual std::string GetLicensedName() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      virtual int Main() const = 0;

      //! @brief returns readme information
      //! @return string containing information about application
      virtual const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      virtual const std::string &GetWebText() const;

      //! @brief return flag whether application is release application. By default the return value is false.
      //!        Release applications have to overwrite the method to return true and indicate the release status.
      //! @return true if application is release application, false otherwise
      virtual bool IsReleaseApplication() const
      {
        // relase applications have to overwrite this method to change return value!
        return false;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the name of the application, excluding the group name (if present
      //! @param GROUP_NAME name of the group to exclude from the beginning of teh applications name
      //! @return the application name if it were part of the group in question
      virtual std::string GetNameForGroup( const std::string &GROUP_NAME = "") const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      virtual std::string GetDescription() const
      {
        return "No description provided";
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      virtual storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>();
      }

      //! @brief returns the string that can be used as default for the technical support section in the readme
      //! it adds the BCLCommons name in the appropriate places
      //! @return the line break terminated support string
      std::string DefaultTechnicalSupportString() const;

      //! @brief returns the string that can be used as default for the terms of use section in the readme
      //! @return line break terminated terms of use string
      static const std::string &DefaultTermsOfUseString();

      //! @brief returns the string that can be used as default for the installation procedure section in the readme
      //! @return line break terminate installation procedure
      static const std::string &DefaultInstallationProcedure();

      //! @brief returns the string seperating sections within a readme
      //! @return line break terminated separator string
      static const std::string &DefaultSectionSeparator();

      //! @brief returns the default read me
      //! @return the default read me
      static const std::string &DefaultReadMe();

    }; // Interface

  } // namespace app
} // namespace bcl

#endif // BCL_APP_INTERFACE_H_
