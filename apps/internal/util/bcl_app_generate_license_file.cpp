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
#include "bcl_app_generate_license_file.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_allowed_non_const.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "storage/bcl_storage_template_instantiations.h"

// external libraries

namespace bcl
{
  namespace app
  {

    //! @brief default constructor
    GenerateLicenseFile::GenerateLicenseFile() :
      m_LicenseFileName
      (
        new command::FlagStatic
        (
          "license_file_name",
          "file name of license file",
          command::Parameter
          (
            "license_file_name",
            "file name of license file",
            "bcl.license"
          )
        )
      ),
      m_NameOfLicensee
      (
        new command::FlagDynamic
        (
          "name",
          "name of licensee",
          command::Parameter
          (
            "name",
            "name of licensee",
            ""
          ),
          1
        )
      ),
      m_NameOfOrganization
      (
        new command::FlagDynamic
        (
          "organization",
          "name of organization",
          command::Parameter
          (
            "organization",
            "name of organization",
            ""
          ),
          1
        )
      ),
      m_Email
      (
        new command::FlagStatic
        (
          "email",
          "email address of licensee",
          command::Parameter
          (
            "email",
            "email address of licensee",
            ""
          )
        )
      ),
      m_VuInnovationsId
      (
        new command::FlagStatic
        (
          "vu_innovations_id",
          "vu innovations customer id",
          command::Parameter
          (
            "vu_innovations_id",
            "vu innovations customer id",
            ""
          )
        )
      ),
      m_ExpirationDuration
      (
        new command::FlagStatic
        (
          "duration",
          "duration until license expires in days, eg. 365 = 1 year",
          command::Parameter
          (
            "duration",
            "duration until license expires [in days], eg. 365 = 1 year",
            command::ParameterCheckRanged< size_t>( 0, 365 * 100),
            "365"
          )
        )
      ),
      m_OperatingSystems
      (
        new command::FlagDynamic
        (
          "operating_system",
          "name of operating system",
          command::Parameter
          (
            "operating_system",
            "name of operating system",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create( "Linux", "Windows", "Darwin")
            ),
            "Linux"
          ),
          1
        )
      ),
      m_BclApplications
      (
        new command::FlagDynamic
        (
          "applications",
          "unlocked applications for this license",
          command::Parameter
          (
            "applications",
            "unlocked applications for this license",
            command::ParameterCheckAllowedNonConst( GetApps().GetReleaseApplicationNames()),
            ""
          ),
          1
        )
      ),
      m_DecryptLicenseFile
      (
        new command::FlagStatic
        (
          "decrypt_license",
          "file name of license file to decrypt content",
          command::Parameter
          (
            "decrypt_license",
            "file name of license file to decrypt content",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    GenerateLicenseFile *GenerateLicenseFile::Clone() const
    {
      return new GenerateLicenseFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GenerateLicenseFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateLicenseFile::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input tag of sequence
      sp_cmd->AddFlag( m_LicenseFileName);
      sp_cmd->AddFlag( m_NameOfLicensee);
      sp_cmd->AddFlag( m_NameOfOrganization);
      sp_cmd->AddFlag( m_Email);
      sp_cmd->AddFlag( m_VuInnovationsId);
      sp_cmd->AddFlag( m_ExpirationDuration);
      sp_cmd->AddFlag( m_OperatingSystems);
      sp_cmd->AddFlag( m_BclApplications);
      sp_cmd->AddFlag( m_DecryptLicenseFile);

      // skip default bcl parameters; they are not applicable to this application

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string GenerateLicenseFile::GetDescription() const
    {
      return "Write a license file that enables use of specified apps in a bcl release";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GenerateLicenseFile::GetReadMe() const
    {
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::GenerateLicenseFile, terms of use, installation "
        "procedures, execution, and technical support.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::GenerateLicenseFile?\n"
        "BCL::GenerateLicenseFile is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::GenerateLicenseFile is a utility for generating"
        "an encrypted license file containing license information that is used by the bcl to enable licensed features\n"
        "over a specified time period\n"
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "V. RUNNING BCL::GenerateLicenseFile.\n"
        "Running BCL::GenerateLicenseFile\n"
        "\nTo create a license file named bcl.license for John Doe, with user id 567, for 90 days, for Jufo and SSPred,\n"
        "for Linux and Windows systems, the following command line would be used\n"
        "bcl.exe GenerateLicenseFile -license_file_name bcl.license -name \"John Doe\" \n"
        "        -email john.doe@domain.com -organization \"Vanderbilt University\" \n"
        "        -vu_innovations_id 567 -duration 90 -operating_system Linux -applications Jufo SSPred"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe GenerateLicenseFile -help\n"
        "\n"
        "For more general information about the product, type bcl.exe GenerateLicenseFile -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator()
      );
      return s_readme;
    }

    //! Main
    int GenerateLicenseFile::Main() const
    {
      // if flag for license file decryption is set
      if( m_DecryptLicenseFile->GetFlag())
      {
        BCL_MessageStd( "Print decrypted license file ..");

        // read in encrypted license file
        // License license_read( m_DecryptLicenseFile->GetFirstParameter()->GetValue());

        // write decrypted license file content
        // license_read.WriteUnencryptedMembers( util::GetLogger());

        return 0;
      }

      // duration in seconds, conversion from days into seconds, x days * 24h * 60 min * 60 sec
      const size_t duration_seconds
      (
        m_ExpirationDuration->GetFirstParameter()->GetNumericalValue< size_t>() * 24 * 60 * 60
      );

      // bcl version as string
      std::string version( "");
      version += util::Format()( GetVersion().GetMajor());
      version += ".";
      version += util::Format()( GetVersion().GetMinor());

      // available release set_applications
      const storage::Vector< ApplicationType> release_apps( m_BclApplications->GetObjectList< ApplicationType>());

      // get the license name of each application passed over the command line
      storage::Set< std::string> set_applications;
      for
      (
        storage::Vector< ApplicationType>::const_iterator itr( release_apps.Begin()), itr_end( release_apps.End());
        itr != itr_end;
        ++itr
      )
      {
        set_applications.Insert( ( **itr)->GetLicensedName());
      }

      // available operating systems
      const storage::Vector< std::string> operating_systems( m_OperatingSystems->GetObjectList< std::string>());
      storage::Set< std::string> set_operating_systems( operating_systems.Begin(), operating_systems.End());

      // construct license file
//      License license
//      (
//        util::Join( " ", m_NameOfLicensee->GetStringList()),
//        util::Join( " ", m_NameOfOrganization->GetStringList()),
//        m_Email->GetFirstParameter()->GetValue(),
//        m_VuInnovationsId->GetFirstParameter()->GetValue(),
//        version,
//        util::Time::GetCurrent(),
//        util::Time::GetCurrent() + util::Time( duration_seconds, 0),
//        storage::Set< std::string>::Create( "x86_64", "x86", "ppc"),
//        set_operating_systems,
//        set_applications
//      );

      // write license file
//      license.WriteEncryptedMembersToFile( m_LicenseFileName->GetFirstParameter()->GetValue());

      BCL_MessageStd
      (
        "Encrypted license file written to " + util::Format()( m_LicenseFileName->GetFirstParameter()->GetValue())
      );

      // end
      return 0;
    }

    const ApplicationType GenerateLicenseFile::GenerateLicenseFile_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateLicenseFile(), GetAppGroups().e_Bcl)
    );

  } // namespace app
} // namespace bcl
