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
#include "bcl.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl //! global namespace of the biochemistry library
{

  //! @brief identifier for the name space
  //! @return the name of the namespace
  const std::string &GetNamespaceIdentifier()
  {
    // please do leave that line, until the visual studio linker error is fixed. The ExtractNamespaceIdentifier function
    // is a BCL_API function and bcl.cpp is in the bcl.dll project and cannot import symbols, so you would get linker error
    static const std::string *s_namespace_name( new std::string( "bcl"));
//    static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
    return *s_namespace_name;
  }

  //! @brief get the copyright
  //! @return the copyright
  const std::string &GetCopyright()
  {
    static const std::string s_copyright( "BioChemistry Library (BCL), Copyright Meiler Lab 2006-2014, www.meilerlab.org");
    return s_copyright;
  }

} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_version.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_format.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////
  // ATTENTION!!
  // changes to bcl_version.cpp has no effect when using cmake!
  // only change version.cpp.in
  /////////////////////////////////////////////////////////////
  
  //! @brief reference to static instance of class Version
  const Version &GetVersion()
  {
    static const Version s_version;
    return s_version;
  }

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

  //! @brief default constructor
  Version::Version() :
    m_Major(           4),
    m_Minor(           3),
    m_Patch(           0),
    m_Revision(        "0"),
    m_IsRelease(       false),
    m_IsLicense(       false),
    m_InstallPrefix(   "./"),
    m_Architecture(    "x86_64"),
    m_OperatingSystem( "Linux")
  {
  }

/////////////////
// data access //
/////////////////

  //! @brief the version string
  //! @brief @return the version string of the form major.minor.patch
  const std::string Version::GetString() const
  {
    return util::Format()( m_Major) + '.' + util::Format()( m_Minor) + '.' + util::Format()( m_Patch);
  }
  
  //! @brief get the complete version string
  //! @return the version string of the form major.minor.patch_svnrevision
  const std::string Version::GetCompleteString() const
  {
    return GetString() + "_r" + m_Revision;
  }

  //! @brief return compilation time
  //! @return date and time of compilation
  const util::Time &Version::GetCompilationTime()
  {
    // statically initialized time of compilation
    static const util::Time s_compilation_time( util::Time::CreateTimeFromCompilerMacro( __DATE__, __TIME__));

    // end
    return s_compilation_time;
  }
  
  //! @brief get all version / build info as a string
  //! @return BCL v{major}.{minor}, r_{revision}, built {Compilation Date & Time}
  //! For releases, svn revision is omitted
  const std::string &Version::GetDetailedInfoString() const
  {
    static const std::string s_detailed_info
    (
      "BCL v" + GetString()
      + ( m_IsRelease ? std::string() : ", r" + m_Revision)
      + ", compiled on " + GetVersion().GetCompilationTime().GetTimeAsDate()
    );
    
    return s_detailed_info;
  }

} // namespace bcl
