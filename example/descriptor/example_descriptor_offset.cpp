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
#include "descriptor/bcl_descriptor_offset.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_offset.cpp
  //!
  //! @author mendenjl
  //! @date Dec 09, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorOffset :
    public ExampleInterface
  {
  public:

    ExampleDescriptorOffset *Clone() const
    {
      return new ExampleDescriptorOffset( *this);
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
      util::Implementation< descriptor::Base< char, float> > next_letter
      (
        "ReflectingOffset("
          "offset=1,"             // Distance away from the central AA to retrieve the value
          "AlphabeticNumber"      // descriptor
        ")"
      );

      BCL_ExampleAssert( next_letter.IsDefined(), true);

      BCL_ExampleCheck( next_letter->GetSizeOfFeatures(), 1);

      // try a string with a strong period of 2.0
      BCL_ExampleIndirectCheck
      (
        descriptor::StringSequence::WriteIterations( next_letter, "ababababab", 2),
        util::Repeat( "2.00 ; 1.00 ; ", 5),
        "Offset value should alternate for this simple window"
      );

      // try a string with a period = 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( next_letter, "cbabc", 2),
        "2.00 ; 1.00 ; 2.00 ; 3.00 ; 2.00 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorOffset

  const ExampleClass::EnumType ExampleDescriptorOffset::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorOffset())
  );
} // namespace bcl
