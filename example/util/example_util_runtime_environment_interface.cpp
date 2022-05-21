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
#include "util/bcl_util_runtime_environment_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example  example_util_runtime_environment_interface.cpp
  //!
  //! @author woetzen
  //! @date Nov 10, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilRuntimeEnvironmentInterface :
    public ExampleInterface
  {
  public:

    ExampleUtilRuntimeEnvironmentInterface *Clone() const
    {
      return new ExampleUtilRuntimeEnvironmentInterface( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////
      BCL_MessageStd
      (
        "Entering: " + GetStaticClassName< ExampleUtilRuntimeEnvironmentInterface>()
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////
      BCL_MessageStd( "Calling the runtime environment function ResolveFileName");

      const std::string resolved_filename
      (
        util::GetRuntimeEnvironment().ResolveFileName( "example.cpp")
      );

      BCL_MessageStd
      (
        "Resolved example.cpp to " + resolved_filename +
        " was successful?: " + util::Format()( !resolved_filename.empty())
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilRuntimeEnvironmentInterface

  const ExampleClass::EnumType ExampleUtilRuntimeEnvironmentInterface::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilRuntimeEnvironmentInterface())
  );

} // namespace bcl
