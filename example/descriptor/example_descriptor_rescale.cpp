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
#include "descriptor/bcl_descriptor_rescale.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_rescale.cpp
  //!
  //! @author mendenjl
  //! @date Oct 02, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorRescale :
    public ExampleInterface
  {
  public:

    ExampleDescriptorRescale *Clone() const
    {
      return new ExampleDescriptorRescale( *this);
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

      // create a descriptor that will return undefined values for some numbers
      util::Implementation< descriptor::Base< char, float> > alphabetic_number( "AlphabeticNumber");

      // create a hybrid descriptor that rescales the underlying descriptor
      descriptor::Rescale< char> rescaled_alphabetic_number( alphabetic_number);

      // test size of features; should be identical to the size of the core descriptor
      BCL_ExampleAssert( rescaled_alphabetic_number.GetSizeOfFeatures(), alphabetic_number->GetSizeOfFeatures());

      // test that the descriptor returns the expected values
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alphabetic_number, "abc", 1),
        "1.0 ; 2.0 ; 3.0 ; "
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations
        (
          util::Implementation< descriptor::Base< char, float> >( rescaled_alphabetic_number),
          "cab",
          1
        ),
        "1.2 ; -1.2 ; 0.0 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorRescale

  const ExampleClass::EnumType ExampleDescriptorRescale::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorRescale())
  );
} // namespace bcl
