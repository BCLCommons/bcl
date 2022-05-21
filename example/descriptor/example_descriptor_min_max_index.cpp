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
#include "descriptor/bcl_descriptor_min_max_index.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_min_max_index.cpp
  //!
  //! @author mendenjl
  //! @date Mar 06, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMinMaxIndex :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMinMaxIndex *Clone() const
    {
      return new ExampleDescriptorMinMaxIndex( *this);
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

      // test constructors
      descriptor::MinMaxIndex< char> max_index( true), min_index( false);

      BCL_ExampleCheck( max_index.GetAlias(), "MaxIndex");
      BCL_ExampleCheck( min_index.GetAlias(), "MinIndex");
      BCL_ExampleCheck( max_index.GetProperty().IsDefined(), false);

      // create a test descriptor which will be the alphabetic index of the character, and then a constant 4.5
      // the max index should always be 0 or 1, depending on whether a given character is above 4
      util::Implementation< descriptor::Base< char, float> > test_descriptor( "Combine(AlphabeticNumber,Constant(4.5))");
      max_index.SetProperty( test_descriptor);
      min_index.SetProperty( test_descriptor);
      BCL_ExampleIndirectCheck( max_index.GetProperty(), test_descriptor, "SetProperty");

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( max_index, "abcdefgA", 1),
        // a     b     c     d     e     f     g     A
        "1.0 ; 1.0 ; 1.0 ; 1.0 ; 0.0 ; 0.0 ; 0.0 ; 1.0 ; "
      );
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( min_index, "abcdefgA", 1),
        // a     b     c     d     e     f     g     A
        "0.0 ; 0.0 ; 0.0 ; 0.0 ; 1.0 ; 1.0 ; 1.0 ; 0.0 ; "
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMinMaxIndex

  const ExampleClass::EnumType ExampleDescriptorMinMaxIndex::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMinMaxIndex())
  );

} // namespace bcl
