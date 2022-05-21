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
#include "descriptor/bcl_descriptor_periodogram.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_periodogram.cpp
  //!
  //! @author mendenjl
  //! @date Mar 12, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPeriodogram :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPeriodogram *Clone() const
    {
      return new ExampleDescriptorPeriodogram( *this);
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
        "ReflectingPeriodogram("
          "size=4,"             // # of characters in the window
          "AlphabeticNumber,"   // descriptor
          "alignment=Center,"   // how the window is aligned to the central object / character
          "bins([1],[2],[2.5],[3],[3.5],[4])," // periods of interest
          "bin resolution(0,0,0,0,0,0)"        // resolution
        ")"
      );

      BCL_ExampleAssert( periodiogram_periods.IsDefined(), true);

      BCL_ExampleCheck( periodiogram_periods->GetSizeOfFeatures(), 6);

      // try a string with a strong period of 2.0
      BCL_ExampleIndirectCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "ababababab", 2),
        util::Repeat( "1.50 0.50 0.00 0.00 0.00 0.00 ; ", 10),
        "period amplitude should not change by position in the window"
      );

      // try a string with a period = 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabc", 2),
        util::Repeat( "2.00 0.00 0.00 0.00 0.00 1.00 ; ", 5)
      );

      // try something that changes frequency over time
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( periodiogram_periods, "cbabeaeaeae", 2),
        "2.34 0.21 0.14 0.60 0.40 1.46 ; 2.15 0.14 0.57 0.74 0.00 1.36 ; 2.44 0.60 1.35 0.55 1.17 0.79 ; "
        "2.44 0.68 1.28 0.96 0.07 0.70 ; 2.62 1.03 1.39 0.00 0.09 0.84 ; 2.66 1.17 1.21 0.00 0.34 0.78 ; "
        "2.97 1.06 2.21 1.70 0.00 0.51 ; 3.07 1.88 0.35 0.36 0.00 0.21 ; "
        + util::Repeat( "3.00 2.00 0.00 0.00 0.00 0.00 ; ", 3)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPeriodogram

  const ExampleClass::EnumType ExampleDescriptorPeriodogram::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPeriodogram())
  );
} // namespace bcl
