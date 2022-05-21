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

#ifndef BCL_APP_GENERATE_LICENSE_FILE_H_
#define BCL_APP_GENERATE_LICENSE_FILE_H_

// include the namespace header
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically

// external libraries

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateLicenseFile
    //! @brief The application generates a BCL license file which gives controlled access of applications in a
    //!        BCL executable file. This bcl license file is encrypted and used a by VUInnovations as a mechanism to
    //!        unlock specific bcl applications for a given license time frame.
    //!
    //! @author butkiem1
    //! @date Aug 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateLicenseFile :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! filename of license file
      util::ShPtr< command::FlagInterface> m_LicenseFileName;

      //! name of licensee
      util::ShPtr< command::FlagInterface> m_NameOfLicensee;

      //! name of organization
      util::ShPtr< command::FlagInterface> m_NameOfOrganization;

      //! email address of licensee
      util::ShPtr< command::FlagInterface> m_Email;

      //! vuinnovations id
      util::ShPtr< command::FlagInterface> m_VuInnovationsId;

      //! expiration duration of license
      util::ShPtr< command::FlagInterface> m_ExpirationDuration;

      //! supported operating systems
      util::ShPtr< command::FlagInterface> m_OperatingSystems;

      //! all available release applications
      util::ShPtr< command::FlagInterface> m_BclApplications;

      //! file name for license file to be decrypted
      util::ShPtr< command::FlagInterface> m_DecryptLicenseFile;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GenerateLicenseFile();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      GenerateLicenseFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! Main
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

      static const ApplicationType GenerateLicenseFile_Instance;

    }; // class GenerateLicenseFile

  } // namespace app
} // namespace bcl
#endif // BCL_APP_GENERATE_LICENSE_FILE_H_
