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
#include "descriptor/bcl_descriptor_for_each.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_for_each.cpp
  //!
  //! @author mendenjl
  //! @date Feb 02, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorForEach :
    public ExampleInterface
  {
  public:

    ExampleDescriptorForEach *Clone() const
    {
      return new ExampleDescriptorForEach( *this);
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

      // create an object data label that should be able to initialized a combine
      util::ObjectDataLabel label
      (
        "ForEach("
        "template=ReflectingWindow(size=1,Desc),"          // size = 3
        "variable=Desc,"
        "descriptors(Character,Partial(ReflectingWindow(size=1,Character),indices(1)))"                       // size = 5
        ")"
      );

      // retrieve a window of 2 characters around a string.  Use the reflecting instance so that we don't have to worry
      // about undefined characters
      // Then use the partial to only get the +1'st and the -2 character
      util::Implementation< descriptor::Base< char, char> > combine_window_seq_size;

      // assert that it can be read
      BCL_ExampleIndirectAssert
      (
        combine_window_seq_size.TryRead( label, util::GetLogger()),
        true,
        "Reading in " + label.ToString() + " to an implementation to create a ForEach object"
      );

      BCL_ExampleCheck( combine_window_seq_size->GetSizeOfFeatures(), 6);
      BCL_ExampleCheck( combine_window_seq_size->GetType().GetDimension(), 1);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( combine_window_seq_size, "abcd"),
        "abbbaa; bacabb; cbdbac; dcccbb; "
      );

      // now try to increase the dimension; this should modify the sizes of the element-wise descriptors, since
      // they must be performed for each part of the iterator, but not for the sequence-based descriptors, or the partial
      combine_window_seq_size->SetDimension( 5);
      BCL_ExampleIndirectCheck( combine_window_seq_size->GetSizeOfFeatures(), 30, "SetDimensionHook");
      combine_window_seq_size->SetDimension( 2);

      // ensure that the output of the combine is as expected in light of the change in dimension
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( combine_window_seq_size, "abc"),
        "abbbacbaaabb; abbcbbbaabaa; baccbbabbbaa; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorForEach

  const ExampleClass::EnumType ExampleDescriptorForEach::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorForEach())
  );
} // namespace bcl
