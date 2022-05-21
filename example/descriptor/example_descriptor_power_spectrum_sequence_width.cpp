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
#include "descriptor/bcl_descriptor_power_spectrum_sequence_width.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_power_spectrum_sequence_width.cpp
  //!
  //! @author mendenjl
  //! @date Nov 07, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPowerSpectrumSequenceWidth :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPowerSpectrumSequenceWidth *Clone() const
    {
      return new ExampleDescriptorPowerSpectrumSequenceWidth( *this);
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
        "PowerSpectrumSequenceWidth(AlphabeticNumber,alignment=Center, number bins=3,frequencies per bin=3,min period=2)"
      );

      BCL_ExampleAssert( periodiogram_periods.IsDefined(), true);

      BCL_ExampleCheck( periodiogram_periods->GetSizeOfFeatures(), 2);

      // try a string with a strong period of 2.0
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "ababababab", 2),
        util::Repeat( "1.50 0.18 ; ", 10)
      );

      // try a string with a period closer to 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabc", 2),
        "2.02 0.99 ; 2.00 0.98 ; 1.98 0.99 ; 2.00 0.98 ; 2.02 0.99 ; "
      );

      // try something with some mixed periods
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabeaeaeae", 2),
        "2.49 1.15 ; 2.47 1.17 ; 2.54 1.12 ; 2.55 0.98 ; 2.65 0.78 ; 2.69 0.64 ; 2.79 0.46 ; 2.83 0.30 ; 2.87 0.14 ; "
        "2.91 0.06 ; 2.93 0.01 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPowerSpectrumSequenceWidth

  const ExampleClass::EnumType ExampleDescriptorPowerSpectrumSequenceWidth::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPowerSpectrumSequenceWidth())
  );
} // namespace bcl
