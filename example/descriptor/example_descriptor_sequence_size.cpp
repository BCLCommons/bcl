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
#include "descriptor/bcl_descriptor_sequence_size.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_sequence_size.cpp
  //!
  //! @author mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSequenceSize :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSequenceSize *Clone() const
    {
      return new ExampleDescriptorSequenceSize( *this);
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

      // default constructor
      descriptor::SequenceSize< char, float> seq_size_as_float;

      // constructor from a description
      descriptor::SequenceSize< char, char> seq_size_as_characters;

    /////////////////
    // data access //
    /////////////////

      // test the GetNormalSizeOfFeatures function
      BCL_ExampleCheck( seq_size_as_float.GetNormalSizeOfFeatures(), 1);

      // test the GetNormalSizeOfFeatures function
      BCL_ExampleCheck( seq_size_as_characters.GetNormalSizeOfFeatures(), 5);

    ///////////////
    // operators //
    ///////////////

      // initialize strings to be passed
      const descriptor::StringSequence string_a( "asdf"), string_1234as( std::string( size_t( 1234), 'a'));

      // test the operators

      // check the result of default constructor
      BCL_ExampleCheck( seq_size_as_float( string_a)( 0), float( string_a.GetSize()));
      BCL_ExampleCheck( seq_size_as_float( string_1234as)( 0), float( string_1234as.GetSize()));
      BCL_ExampleCheck( std::string( seq_size_as_characters( string_a).Begin(), size_t( 5)), "4    ");
      BCL_ExampleCheck( std::string( seq_size_as_characters( string_1234as).Begin(), size_t( 5)), "1234 ");

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorSequenceSize

  const ExampleClass::EnumType ExampleDescriptorSequenceSize::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSequenceSize())
  );

} // namespace bcl
