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
#include "restraint/bcl_restraint_sas_data_parameters.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_data_parameters.cpp
  //! @details This example tests the create of a saxs data parameters object
  //!          1) Q_value 2) Sasa value 3) c1 - excluded volume 4) c2 - hydration shell
  //!
  //! @author putnamdk
  //! @date Oct 22, 1014
  //! @remarks status complete
  //! @remarks reviewed by putnamdk on Oct 22, 2014
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSaxsDataParameters :
    public ExampleInterface
  {
  public:

      ExampleRestraintSaxsDataParameters *Clone() const
    {
      return new ExampleRestraintSaxsDataParameters( *this);
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

      // test initialization constructor

      // test setting and returning the value
      restraint::SasDataParameters initialized_point( false, 1.0, 0.83, 1.1, 3.2, 0.0);

      // Test the q-value
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetQValue(),
        1.0,
        0.001,
        "q-value"
      );

      // Test Sasa
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetSasaValue(),
        0.83,
        0.001,
        "Intensity"
      );

      // Test c1
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetExcludedVolume(),
        1.1,
        0.001,
        "Error"
      );

      // Test c2
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetHydrationShell(),
        3.2,
        0.001,
        "Error"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxsScatteringPoint

  const ExampleClass::EnumType ExampleRestraintSaxsDataParameters::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSaxsDataParameters())
  );

} // namespace bcl
