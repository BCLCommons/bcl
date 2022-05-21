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

#ifndef BCL_VERSION_H_
#define BCL_VERSION_H_

// include the namespace header
#include "bcl.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{

  //! @brief reference to static instance of class Version
  BCL_API const Version &GetVersion();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @class Version
  //! @brief Class used for version and licensing information
  //!
  //! @see @link example_bcl_version.cpp @endlink
  //! @author woetzen
  //! @date Feb 20, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class BCL_API Version
  {
  /////////////
  // friends //
  /////////////

    friend const Version &GetVersion();

  private:

  //////////
  // data //
  //////////

    size_t      m_Major;           //! major version
    size_t      m_Minor;           //! minor version
    size_t      m_Patch;           //! patch version
    std::string m_Revision;        //! svn revision
    bool        m_IsRelease;       //! is built as release version (not by compiler flags, but intended for distribution
    bool        m_IsLicense;       //! is built with license information (expiration time)
    std::string m_InstallPrefix;   //! the folder in which everything is installed
    std::string m_Architecture;    //! architecture this application was compiled for
    std::string m_OperatingSystem; //! operating system this application was compiled for

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Version();

  public:

  /////////////////
  // data access //
  /////////////////

    //! @brief access to major version
    //! @return the major version number
    const size_t GetMajor() const
    {
      return m_Major;
    }

    //! @brief access to minor version
    //! @return the minor version
    const size_t GetMinor() const
    {
      return m_Minor;
    }

    //! @brief access to patch version
    //! @return the patch version
    const size_t GetPatch() const
    {
      return m_Patch;
    }

    //! @brief access to trunk revision number
    //! @return the svn revision number
    const std::string GetRevision() const
    {
      return m_Revision;
    }

    //! @brief is release
    //! @return true if released version
    const bool   IsRelease() const
    {
      return m_IsRelease;
    }

    //! @brief is license
    //! @return true if licensed version
    const bool   IsLicense() const
    {
      return m_IsLicense;
    }

    //! @brief return the install prefix that can be used to find installed components
    //! @brief the install prefix
    const std::string &GetInstallPrefix() const
    {
      return m_InstallPrefix;
    }

    //! @brief return the architecture this application was compiled for
    //! @return architecture string
    const std::string &GetArchitecture() const
    {
      return m_Architecture;
    }

    //! @brief return the operating system this application was compiled for
    //! @return operating system string
    const std::string &GetOperatingSystem() const
    {
      return m_OperatingSystem;
    }

    //! @brief the version string
    //! @return the version string of the form major.minor.patch
    const std::string GetString() const;

    //! @brief get the complete version string
    //! @return the version string of the form major.minor.patch_svnrevision
    const std::string GetCompleteString() const;

    //! @brief return compilation time
    //! @return date and time of compilation
    static const util::Time &GetCompilationTime();

    //! @brief get all version / build info as a string
    //! @return BCL v{major}.{minor}, r_{revision}, built {Compilation Date & Time}
    //! For releases, svn revision is omitted
    const std::string &GetDetailedInfoString() const;

  }; // class Version

} // namespace bcl

#endif // BCL_VERSION_H_
