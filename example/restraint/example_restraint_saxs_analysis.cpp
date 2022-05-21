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
#include "restraint/bcl_restraint_sas_analysis.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_sas_density_data.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_analysis.cpp
  //!
  //! @author putnamdk
  //! @date May 28, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSaxsAnalysis :
    public ExampleInterface
  {
  public:

    ExampleRestraintSaxsAnalysis *Clone() const
    {
      return new ExampleRestraintSaxsAnalysis( *this);
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
      // create object to perform analysis functions on
      restraint::SasScatteringData saxs_data;
      restraint::SasDensityData saxs_data_density;

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "A10_2.mccd.gnom"));
      saxs_data.ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      io::IFStream read_density;
      BCL_ExampleMustOpenInputFile( read_density, AddExampleInputPathToFilename( e_Biology, "A10_2.mccd.gnom"));
      saxs_data_density.ReadFromDataFile( read_density);
      io::File::CloseClearFStream( read_density);

    /////////////////////////////
    // Saxs Analysis Functions //
    /////////////////////////////

      // Test the qmax value
      double qmax( restraint::SasAnalysis::ComputeQmax( saxs_data));
      BCL_ExampleIndirectCheckWithinTolerance( qmax, 0.323157, 0.001, "qmax");

      // Test the dmax value
      double dmax( restraint::SasAnalysis::ComputeDmax( saxs_data_density));
      BCL_ExampleIndirectCheckWithinTolerance( dmax, 100, 0.001, "dmax");

      // Test scaled data function
      restraint::SasScatteringData scaled_data( restraint::SasAnalysis::ScaleData( saxs_data, 10));

      BCL_ExampleIndirectCheckWithinTolerance
      (
        scaled_data.GetScatteringData().FirstElement().GetQvalue(),
        0.014898,
        0.001,
        "scaling q-value"
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        scaled_data.GetScatteringData().FirstElement().GetIntensity(),
        1.28377,
        0.001,
        "scaling intensity"
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        scaled_data.GetScatteringData().FirstElement().GetError(),
        29.2959,
        0.001,
        "scaling error"
      );

      // Test log10 data function
      restraint::SasScatteringData log_data( restraint::SasAnalysis::Log10( saxs_data));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        log_data.GetScatteringData().FirstElement().GetQvalue(),
        0.014898,
        0.001,
        "log q-value"
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        log_data.GetScatteringData().FirstElement().GetIntensity(),
        -0.891513,
        0.001,
        "log intensity"
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        log_data.GetScatteringData().FirstElement().GetError(),
        1.37695,
        0.001,
        "log error"
      );

      // Test Derivative data function
      restraint::SasScatteringData derivative_data( restraint::SasAnalysis::Derivative( saxs_data));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        derivative_data.GetScatteringData().FirstElement().GetQvalue(),
        0.014898,
        0.001,
        "Derivative q-value"
      );

      // the derivative should be 0 because the first two points are zero (so y' should be 0) and we specified a natural
      // boundary condition so y'' should also be zero at the end points
      BCL_ExampleIndirectCheckWithinTolerance
      (
        derivative_data.GetScatteringData().FirstElement().GetIntensity(),
        0.0,
        0.001,
        "Derivative intensity"
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        derivative_data.GetScatteringData().FirstElement().GetError(),
        2.92959,
        0.001,
        "Derivative error"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSaxsAnalysis

  const ExampleClass::EnumType ExampleRestraintSaxsAnalysis::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSaxsAnalysis())
  );

} // namespace bcl
