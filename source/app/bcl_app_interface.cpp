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
#include "app/bcl_app_interface.h"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_group_handler.h"
#include "command/bcl_command_command.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    //! e-mail address for commercial license support
    const std::string &Interface::GetSupportEmailCommercial()
    {
      static const std::string s_support_email_commercial( "bcl-support-commercial@meilerlab.org");
      return s_support_email_commercial;
    }

    //! e-mail address for academic license support
    const std::string &Interface::GetSupportEmailAcademic()
    {
      static const std::string s_support_email_academic(   "bcl-support-academic@meilerlab.org");
      return s_support_email_academic;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief the expiration interval from last compilation
    //! @return the time interval from the last compilation, the license is valid - by default, it is around 20 years
    const util::Time &Interface::GetLicenseExpirationTime() const
    {
      // every 14 months which turns out to be about 425 days
      static const util::Time s_default_license_expiration( 425 * 24 * 60 * 60, 0);

      // by default, it returns 425 days
      return s_default_license_expiration;
    }

    //! @brief return the bcl::commons name
    //! @return string for the bcl::commons name of that application
    std::string Interface::GetBCLScopedName() const
    {
      return "BCL" + this->GetClassIdentifier().substr( 8);
    }

    //! @brief return the original name of the application, as it was used in license files
    //! @return string for the bcl::commons name of that application
    //! This is necessary so that, if release application names change, licenses will continue to work
    std::string Interface::GetLicensedName() const
    {
      // by default, return the app name
      return GetNameForGroup();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns ReadMe information
    //! @return string containing information about application
    const std::string &Interface::GetReadMe() const
    {
      // read me message
      return DefaultReadMe();
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &Interface::GetWebText() const
    {
      return GetReadMe();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the name of the application, excluding the group name (if present
    //! @param GROUP_NAME name of the group to exclude from the beginning of teh applications name
    //! @return the application name if it were part of the group in question
    std::string Interface::GetNameForGroup( const std::string &GROUP_NAME) const
    {
      // generate the string that needs to be removed e.g. "const-bcl::app::"
      static const std::string s_obsolete_string( GetNamespaceIdentifier() + "::");

      // get the full app name
      std::string full_app_name_no_group( GetClassIdentifier().substr( s_obsolete_string.length()));

      // return it if group name was empty
      if( GROUP_NAME.empty())
      {
        return full_app_name_no_group;
      }

      // otherwise, check whether the application name and group name start off the same, ignoring case
      if( util::StartsWith( util::ToLower( full_app_name_no_group), util::ToLower( GROUP_NAME)))
      {
        // erase the group name from the full app name
        full_app_name_no_group.erase( 0, GROUP_NAME.size());
      }

      // return the group:full_app_name_no_group
      return GROUP_NAME + GroupHandler::s_GroupDelimiter + full_app_name_no_group;
    }

    //! @brief returns the string that can be used as default for the technical support section in the readme
    std::string Interface::DefaultTechnicalSupportString() const
    {
      return
        "Technical support is available throughout the license period.  For scientific or technical questions related "
        "to " + GetBCLScopedName() + ", commercial entities should contact us at " + GetSupportEmailCommercial() + " while "
        "academic institutions should contact us at " + GetSupportEmailAcademic() + ".  If there is a major bug fix or "
        "significant improvements to the application, entities with current licenses will be notified via email, and "
        "will be given access to the updated software.\n";
    }

    //! @brief returns the string that can be used as default for the terms of use section in the readme
    //! @return line break terminated terms of use string
    const std::string &Interface::DefaultTermsOfUseString()
    {
      static const std::string s_terms_of_use
      (
        "Licensees are bound by the license agreement in place at the time this software was obtained, unless legally "
        "modified by the parties. The following brief summary does not replace the original agreement. It is expected "
        "that the software will be used for its intended purpose and that appropriate citations and acknowledgments "
        "will be given in connection to this software.  The license provides access to the software, updates, and "
        "reasonable technical support for one year, after which time renewal will be required for continued use and "
        "support.  The license allows the licensee to have the software installed only on site(s) designated at the "
        "time the license was obtained, usually a single personal computer (PC). The license does not allow the "
        "software to be passed to third parties or to be modified in any way.\n"
      );

      return s_terms_of_use;
    }

    //! @brief returns the string that can be used as default for the installation procedure section in the readme
    //! @return line break terminate installation procedure
    const std::string &Interface::DefaultInstallationProcedure()
    {
      static const std::string s_installation_procedure
      (
        "Run the installer and follow the directions. On Linux, if a \"lib\" folder is present in the installation "
        "directory, it will be necessary to add that with its full path to the LD_LIBRARY_PATH environment variable, "
        "in order for the dynamic linker (ldd) is able to find all libraries to run the application. Please be also "
        "aware of the licenses for the external libraries, which are installed as well.\n"
        "Uninstallers are provided for some operating systems.\n"
        "If there are major bug fixes or revisions, an email will be sent out to entities with current licenses so "
        "that the newest application can be installed.\n"
      );

      return s_installation_procedure;
    }

    //! @brief returns the string seperating sections within a readme
    //! @return line break terminated separator string
    const std::string &Interface::DefaultSectionSeparator()
    {
      return command::Command::DefaultSectionSeparator();
    }

    //! @brief returns the default read me
    //! @return the default read me
    const std::string &Interface::DefaultReadMe()
    {
      static const std::string s_read_me( "no ReadMe information is provided for this application");
      return s_read_me;
    }

  } // namespace app
} // namespace bcl
