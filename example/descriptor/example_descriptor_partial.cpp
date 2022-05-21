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
#include "descriptor/bcl_descriptor_partial.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_partial.cpp
  //!
  //! @author mendenjl
  //! @date Feb 13, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPartial :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPartial *Clone() const
    {
      return new ExampleDescriptorPartial( *this);
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
      // Then use the partial to only get the +1'st and the -2 character
      util::Implementation< descriptor::Base< char, char> > partial_window_char
      (
        "Partial(ReflectingWindow(size=2,Character),indices(1,4))"
      );

      BCL_ExampleAssert( partial_window_char.IsDefined(), true);

      BCL_ExampleCheck( partial_window_char->GetSizeOfFeatures(), 2);

      // try calling set dimension with several different values; if this fails, it will probably assert
      partial_window_char->SetDimension( 5);
      BCL_ExampleIndirectCheck( partial_window_char->GetSizeOfFeatures(), 2, "SetDimensionHook");
      partial_window_char->SetDimension( 2);
      partial_window_char->SetDimension( 1);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( partial_window_char, "abcdefg"),
        "bc; ad; be; cf; dg; ef; fe; "
      );
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( partial_window_char, "abc"), "bc; ab; ba; ");

      // now try to read in a partial where the indices are higher than the dimension of the descriptor; for elementwise
      // descriptors (like window), this should work, but just automatically adjust the type
      BCL_ExampleIndirectCheck
      (
        partial_window_char.TryRead
        (
          util::ObjectDataLabel( "Partial(ReflectingWindow(size=2,Character),indices(1,7))"),
          util::GetLogger()
        ),
        true,
        "Partial's handling of elementwise descriptors"
      );

      // make sure the dimension is correct
      BCL_ExampleCheck( partial_window_char->GetSizeOfFeatures(), 2);
      BCL_ExampleCheck( partial_window_char->GetType().GetDimension(), 2);
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( partial_window_char, "abcdefg"),
        "bb; ba; bb; bc; bd; be; aa; ab; ac; ad; ae; bb; bc; bd; be; cc; cd; ce; dd; de; ee; "
      );

      // try calling set dimension a few times on it as well
      partial_window_char->SetDimension( 3);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPartial

  const ExampleClass::EnumType ExampleDescriptorPartial::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPartial())
  );
} // namespace bcl
