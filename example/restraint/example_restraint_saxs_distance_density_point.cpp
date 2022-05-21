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
#include "restraint/bcl_restraint_sas_distance_density_point.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_distance_density_point.cpp
  //! @details This example tests the create of a single saxs distance density point with three values:
  //!          1) disance in angstroms (r) 2) density (P) 3) Experimental Error
  //!
  //! @author putnamdk
  //! @date May 23, 2014
  //! @remarks status complete
  //! @remarks reviewed by putnamdk on May 23, 2014
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSaxsDistanceDensityPoint :
    public ExampleInterface
  {
  public:

    ExampleRestraintSaxsDistanceDensityPoint *Clone() const
    {
      return new ExampleRestraintSaxsDistanceDensityPoint( *this);
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
      restraint::SasDistanceDensityPoint initialized_point( 0.152, 15423, 0.54);

      // Test the r-value
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetRvalue(),
        0.152,
        0.001,
        "r-value"
      );

      // Test density
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetDensity(),
        15423,
        0.001,
        "Density"
      );

      // Test error
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetError(),
        0.54,
        0.001,
        "Error"
      );

      // Test the mutators

      initialized_point.SetDensity( 4.2);
      initialized_point.SetError( 5.1);

      // Test intensity
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetDensity(),
        4.2,
        0.001,
        "Density"
      );

      // Test error
      BCL_ExampleIndirectCheckWithinTolerance
      (
        initialized_point.GetError(),
        5.1,
        0.001,
        "Error"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxsDistanceDensityPoint

  const ExampleClass::EnumType ExampleRestraintSaxsDistanceDensityPoint::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSaxsDistanceDensityPoint())
  );

} // namespace bcl
