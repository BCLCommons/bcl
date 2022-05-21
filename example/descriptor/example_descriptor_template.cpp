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
#include "descriptor/bcl_descriptor_named_template.h"
#include "descriptor/bcl_descriptor_template.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_template.cpp
  //!
  //! @author mendenjl
  //! @date Jan 30, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorTemplate :
    public ExampleInterface
  {
  public:

    ExampleDescriptorTemplate *Clone() const
    {
      return new ExampleDescriptorTemplate( *this);
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
      util::Implementation< descriptor::Base< char, float> > alphabetic_number_else_five_impl
      (
        "DefineNaN(AlphabeticNumber,replacement=Constant(0))"
      );

      // create a hybrid descriptor that replaces undefined values with a different descriptor
      descriptor::Template< char, float> alphabetic_number_else_five_define
      (
        util::ObjectDataLabel( "AlphaNumberElse(Replacement)"),
        util::ObjectDataLabel( "DefineNaN(AlphabeticNumber,replacement=Constant(Replacement))")
      );

      util::Implementation< descriptor::Base< char, float> > alphabetic_number_else_five_impl_b
      (
        util::ObjectDataLabel
        (
          "Template("
          "  signature=AlphaNumberElse(Replacement),"
          "  DefineNaN(AlphabeticNumber,replacement=Constant(Replacement))"
          ")"
        )
      );

      util::Implementation< descriptor::Base< char, float> > alphabetic_number_else_five_impl_via_define;
      BCL_ExampleAssert
      (
        alphabetic_number_else_five_impl_via_define.TryRead
        (
          util::ObjectDataLabel( "AlphaNumberElse(5.0)"),
          util::GetLogger()
        ),
        true
      );

      // test that the descriptors created are defined
      BCL_ExampleAssert( alphabetic_number_else_five_impl.IsDefined(), true);

      // test size of features; should be identical to the size of the core descriptor
      BCL_ExampleAssert
      (
        alphabetic_number_else_five_impl->GetSizeOfFeatures(),
        alphabetic_number_else_five_impl_via_define->GetSizeOfFeatures()
      );

      // test that the descriptor returns the expected values
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( *alphabetic_number_else_five_impl_via_define, "a.b", 1),
        "1.0 ; 5.0 ; 2.0 ; "
      );

      BCL_ExampleCheck( alphabetic_number_else_five_impl_via_define->GetAlias(), "AlphaNumberElse");

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations
        (
          util::Implementation< descriptor::Base< char, float> >( alphabetic_number_else_five_impl_via_define),
          "a.b  ",
          1
        ),
        "1.0 ; 5.0 ; 2.0 ; 5.0 ; 5.0 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorTemplate

  const ExampleClass::EnumType ExampleDescriptorTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorTemplate())
  );
} // namespace bcl
