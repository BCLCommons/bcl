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
#include "opencl/bcl_opencl_extension_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_extension_data.cpp
  //!
  //! @author loweew
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclExtensionData :
    public ExampleInterface
  {
  public:

    ExampleOpenclExtensionData *Clone() const
    {
      return new ExampleOpenclExtensionData( *this);
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

      const std::string sample_extension_string1( "cl_amd_fp64");
      const std::string sample_extension_string2( "cl_khr_global_int32_extended_atomics");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::ExtensionData extension_default;

      // form sample string
      opencl::ExtensionData extension1( sample_extension_string1);
      opencl::ExtensionData extension2( sample_extension_string2);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( extension1.GetVendor(), "amd");
      BCL_ExampleCheck( extension1.GetName(), "fp64");
      BCL_ExampleCheck( extension2.GetVendor(), "khr");
      BCL_ExampleCheck( extension2.GetName(), "global_int32_extended_atomics");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclExtensionData

  const ExampleClass::EnumType ExampleOpenclExtensionData::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclExtensionData())
  );

} // namespace bcl
