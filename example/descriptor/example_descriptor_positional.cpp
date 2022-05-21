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
#include "descriptor/bcl_descriptor_positional.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_positional.cpp
  //!
  //! @author mendenjl
  //! @date Mar 08, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPositional :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPositional *Clone() const
    {
      return new ExampleDescriptorPositional( *this);
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

      // retrieve the first character of a multi-character descriptor
      util::Implementation< descriptor::Base< char, char> > first_char( "1st(Character)");
      util::Implementation< descriptor::Base< char, char> > second_char( "2nd(Character)");
      util::Implementation< descriptor::Base< char, char> > both_chars( "Character");

      // a descriptor for the first char requires only a single char to be valid, so the dimension should remain 1
      BCL_ExampleCheck( first_char->GetType().GetDimension(), 1);

      first_char->SetDimension( 2);
      second_char->SetDimension( 2);
      both_chars->SetDimension( 2);

      BCL_ExampleAssert( first_char.IsDefined(), true);
      BCL_ExampleAssert( second_char.IsDefined(), true);

      BCL_ExampleCheck( first_char->GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( second_char->GetSizeOfFeatures(), 1);

      // a descriptor for the second char requires at least two chars to be valid, so the dimension should be 2
      BCL_ExampleCheck( second_char->GetType().GetDimension(), 2);

      // try calling set dimension with several different values; if this fails, it will probably assert
      first_char->SetDimension( 5);
      BCL_ExampleIndirectCheck
      (
        first_char->GetSizeOfFeatures(),
        1,
        "Changing dimension of Positional descriptor should not increase the size, since "
        "the number of object considered remains the same"
      );
      first_char->SetDimension( 2);

      // ensure that the characters descriptor is working properly
      BCL_ExampleIndirectAssert
      (
        descriptor::StringSequence::WriteIterations( both_chars, "abcd"),
        "ab; ac; ad; bc; bd; cd; ",
        "Simple characters descriptor"
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( first_char, "abcd"),
        "a; a; a; b; b; c; "
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( second_char, "abcd"),
        "b; c; d; c; d; d; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPositional

  const ExampleClass::EnumType ExampleDescriptorPositional::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPositional())
  );
} // namespace bcl
