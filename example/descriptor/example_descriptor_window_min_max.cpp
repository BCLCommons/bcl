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
#include "descriptor/bcl_descriptor_window_min_max.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_window_min_max.cpp
  //!
  //! @author mendenjl
  //! @date May 21, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWindowMinMax :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWindowMinMax *Clone() const
    {
      return new ExampleDescriptorWindowMinMax( *this);
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
      util::Implementation< descriptor::Base< char, float> >
        window_ave( "ReflectingWindowMinMax(size=2,AlphabeticNumber)");

      BCL_ExampleAssert( window_ave.IsDefined(), true);

      BCL_ExampleCheck( window_ave->GetSizeOfFeatures(), 2);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_ave, "abcdefg", 2),
        "1.00 3.00 ; 1.00 4.00 ; 1.00 5.00 ; 2.00 6.00 ; 3.00 7.00 ; 4.00 7.00 ; 5.00 7.00 ; "
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_ave, "abc", 2),
        "1.00 3.00 ; 1.00 3.00 ; 1.00 3.00 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWindowMinMax

  const ExampleClass::EnumType ExampleDescriptorWindowMinMax::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWindowMinMax())
  );
} // namespace bcl
