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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "bcl_version.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_bcl_version.cpp
  //! @brief this examples presents how the version object can be used, with no further checking
  //!
  //! @author woetzen
  //! @date Apr 8, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBclVersion :
    public ExampleInterface
  {
  public:

    ExampleBclVersion *Clone() const
    {
      return new ExampleBclVersion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // no public constructor

    /////////////////
    // data access //
    /////////////////

      // version number
      BCL_MessageStd( "major: " + util::Format()( GetVersion().GetMajor()));
      BCL_MessageStd( "minor: " + util::Format()( GetVersion().GetMinor()));
      BCL_MessageStd( "patch: " + util::Format()( GetVersion().GetPatch()));

      // svn revision
      BCL_MessageStd( "svn revision: " + util::Format()( GetVersion().GetRevision()));

      // is release
      BCL_MessageStd( "is release: " + util::Format()( GetVersion().IsRelease()));

      // is license
      BCL_MessageStd( "is license: " + util::Format()( GetVersion().IsLicense()));

      // install prefix
      BCL_MessageCrt( "install prefix: " + util::Format()( GetVersion().GetInstallPrefix()));

      // architecture
      BCL_MessageCrt( "architecture: " + util::Format()( GetVersion().GetArchitecture()));

      // operating system
      BCL_MessageCrt( "operating system: " + util::Format()( GetVersion().GetOperatingSystem()));

      // check that the major revision number is at least two
      BCL_ExampleCheck( GetVersion().GetMajor() >= size_t( 2), true);

    ////////////////
    // operations //
    ////////////////

      // compilation time

      // version string

      // complete version string

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBclVersion

  const ExampleClass::EnumType ExampleBclVersion::s_Instance
  (
    GetExamples().AddEnum( ExampleBclVersion())
  );

} // namespace bcl
