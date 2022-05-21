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
#include "descriptor/bcl_descriptor_band_pass_filter.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_band_pass_filter.cpp
  //!
  //! @author mendenjl
  //! @date Oct 16, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorBandPassFilter :
    public ExampleInterface
  {
  public:

    ExampleDescriptorBandPassFilter *Clone() const
    {
      return new ExampleDescriptorBandPassFilter( *this);
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
      util::Implementation< descriptor::Base< char, float> > highpass_filter
      (
        "ReflectingBandPassFilter("
          "size=4,"                // # of characters in the window
          "AlphabeticNumber,"      // descriptor
          "alignment=Center,"      // how the window is aligned to the central object / character
          "maxp=4"                 // maximum period of interest (makes this a high pass filter)
        ")"
      );

      BCL_ExampleAssert( highpass_filter.IsDefined(), true);

      BCL_ExampleCheck( highpass_filter->GetSizeOfFeatures(), 1);

      // try a string with a strong period of 2.0
      BCL_ExampleIndirectCheck
      (
        descriptor::StringSequence::WriteIterations( highpass_filter, "ababababab", 2),
        util::Repeat( "-0.46 ; 0.46 ; ", 5),
        "Filtered value should alternate for this simple window"
      );

      // try a string with a period = 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( highpass_filter, "cbabc", 2),
        "0.49 ; 0.00 ; -0.49 ; -0.00 ; 0.49 ; "
      );

      // try something that changes frequency over time
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( highpass_filter, "cbabeaeaeae", 2),
        "0.46 ; 0.20 ; -0.64 ; -0.25 ; 1.65 ; -2.24 ; 1.96 ; -1.84 ; 1.83 ; -1.83 ; 1.83 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorBandPassFilter

  const ExampleClass::EnumType ExampleDescriptorBandPassFilter::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorBandPassFilter())
  );
} // namespace bcl
