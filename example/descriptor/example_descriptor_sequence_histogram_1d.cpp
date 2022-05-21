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
#include "descriptor/bcl_descriptor_sequence_histogram_1d.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_sequence_histogram_1d.cpp
  //!
  //! @author mendenjl
  //! @date May 03, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSequenceHistogram1D :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSequenceHistogram1D *Clone() const
    {
      return new ExampleDescriptorSequenceHistogram1D( *this);
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
      util::Implementation< descriptor::Base< char, float> > histogram
      (
        "StringHistogram1D(NElements,min=1,max=9,bin size=2,catchall=False)"
      );

      BCL_ExampleAssert( histogram.IsDefined(), true);

      BCL_ExampleCheck( histogram->GetSizeOfFeatures(), 4);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( histogram, "abce", 1),
        "0.0 1.0 0.0 0.0 ; "
      );

      util::Implementation< descriptor::Base< char, float> > smoothed_histogram
      (
        "StringHistogram1D(NElements,min=1,max=9,bin size=2,catchall=True,smoothing=1.0)"
      );

      BCL_ExampleAssert( smoothed_histogram.IsDefined(), true);

      BCL_ExampleCheck( smoothed_histogram->GetSizeOfFeatures(), 4);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( smoothed_histogram, "abcedefg", 1),
        "0.0 0.1 0.2 1.5 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorSequenceHistogram1D

  const ExampleClass::EnumType ExampleDescriptorSequenceHistogram1D::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSequenceHistogram1D())
  );
} // namespace bcl
