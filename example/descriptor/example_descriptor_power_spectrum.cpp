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
#include "descriptor/bcl_descriptor_fourier_analysis_window_creator.h"
#include "descriptor/bcl_descriptor_power_spectrum.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically
namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_power_spectrum.cpp
  //!
  //! @author mendenjl
  //! @date Sep 30, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPowerSpectrum :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPowerSpectrum *Clone() const
    {
      return new ExampleDescriptorPowerSpectrum( *this);
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

      // retrieve a window of 2 characters around a string.  Use the reflecting instance so that we don't have to worry
      // about undefined characters
      util::Implementation< descriptor::Base< char, float> > periodiogram_periods
      (
        "ReflectingPowerSpectrum(size=4,AlphabeticNumber,alignment=Center)"
      );

      BCL_ExampleAssert( periodiogram_periods.IsDefined(), true);

      BCL_ExampleCheck( periodiogram_periods->GetSizeOfFeatures(), 5);

      // try a string with a strong period of 2.0
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "ababababab", 2),
        util::Repeat( "1.49 0.47 0.32 0.03 0.04 ; 1.51 0.53 0.32 0.03 0.04 ; ", 5)
      );

      // try a string with a period closer to 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabc", 2),
        "2.02 0.06 0.45 0.98 0.54 ; 2.00 0.04 0.49 0.97 0.61 ; 1.98 0.02 0.45 0.98 0.54 ; "
        "2.00 0.04 0.49 0.97 0.61 ; 2.02 0.06 0.45 0.98 0.54 ; "
      );

      // try something with some mixed periods
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabeaeaeae", 2),
        "2.10 0.14 0.32 1.11 0.23 ; 2.08 0.08 0.45 1.17 0.59 ; 2.25 0.39 0.89 1.36 1.04 ; 2.45 0.69 1.54 1.44 1.14 ; "
        "2.75 1.28 1.85 1.13 0.78 ; 2.90 1.60 1.75 0.59 0.54 ; 3.01 1.93 1.50 0.18 0.13 ; 2.98 1.92 1.31 0.11 0.09 ; "
        "3.04 2.06 1.28 0.14 0.16 ; 2.96 1.94 1.28 0.14 0.16 ; 3.04 2.06 1.28 0.14 0.16 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPowerSpectrum

  const ExampleClass::EnumType ExampleDescriptorPowerSpectrum::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPowerSpectrum())
  );
} // namespace bcl
