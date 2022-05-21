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
#include "descriptor/bcl_descriptor_within_range_smooth.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_within_range_smooth.cpp
  //!
  //! @author geanesar, mendenjl
  //! @date Feb 05, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWithinRangeSmooth :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWithinRangeSmooth *Clone() const
    {
      return new ExampleDescriptorWithinRangeSmooth( *this);
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

      // Test the default implementation
      descriptor::WithinRangeSmooth< char> within_range_smooth_def;

      // Test a normal range
      descriptor::WithinRangeSmooth< char> within_range_smooth_normal( 3, 4, 2, 2, true);

      // Test a zero-width range on the left
      descriptor::WithinRangeSmooth< char> within_range_smooth_noleft( 3, 4, 0, 2, true);

      // Test a zero-width range on the right
      descriptor::WithinRangeSmooth< char> within_range_smooth_noright( 3, 4, 2, 0, true);

      // Test a zero-width range on both sides with non-inclusive bounds
      descriptor::WithinRangeSmooth< char> within_range_smooth_noninclusive( 3, 4, 0, 0, false);

      // test Clone function
      util::ShPtr< descriptor::WithinRangeSmooth< char> > sp_within_range_smooth
      (
        within_range_smooth_normal.Clone()
      );
      BCL_ExampleCheck( sp_within_range_smooth.IsDefined(), true);

      // Destroy the new object
      sp_within_range_smooth = util::ShPtr< descriptor::WithinRangeSmooth< char> >();

    /////////////////
    // data access //
    /////////////////

      // Check the class name is correct
      BCL_ExampleCheck
      (
        GetStaticClassName< descriptor::WithinRangeSmooth< char> >(),
        within_range_smooth_normal.GetClassIdentifier()
      );
      BCL_MessageStd( "class name: " + within_range_smooth_normal.GetClassIdentifier());

      // The property used to test the range
      util::Implementation< descriptor::Base< char, float> >
              str_descriptor( "AlphabeticNumber");

      within_range_smooth_normal.SetDescriptor( str_descriptor);
      within_range_smooth_noleft.SetDescriptor( str_descriptor);
      within_range_smooth_noright.SetDescriptor( str_descriptor);
      within_range_smooth_noninclusive.SetDescriptor( str_descriptor);

      BCL_ExampleCheck( within_range_smooth_normal.GetDescriptor(), str_descriptor);

    //////////////////////
    // input and output //
    //////////////////////

      // String used to test the descriptor
      std::string test_string( "ABCDEF");

      // Test the normal range object
      std::string correct_normal_output( "0.00 ; 0.50 ; 1.00 ; 1.00 ; 0.50 ; 0.00 ; ");
      std::string normal_output( descriptor::StringSequence::WriteIterations( within_range_smooth_normal, test_string, 2));
      BCL_ExampleCheck( normal_output, correct_normal_output);

      // Test the noleft range object
      std::string correct_noleft_output( "0.00 ; 0.00 ; 1.00 ; 1.00 ; 0.50 ; 0.00 ; ");
      std::string noleft_output( descriptor::StringSequence::WriteIterations( within_range_smooth_noleft, test_string, 2));
      BCL_ExampleCheck( noleft_output, correct_noleft_output);

      // Test the noright range object
      std::string correct_noright_output( "0.00 ; 0.50 ; 1.00 ; 1.00 ; 0.00 ; 0.00 ; ");
      std::string noright_output( descriptor::StringSequence::WriteIterations( within_range_smooth_noright, test_string, 2));
      BCL_ExampleCheck( noright_output, correct_noright_output);

      // Test the noninclusive range object
      std::string correct_noninclusive_output( "0.00 ; 0.00 ; 0.00 ; 0.00 ; 0.00 ; 0.00 ; ");
      std::string noninclusive_output( descriptor::StringSequence::WriteIterations( within_range_smooth_noninclusive, test_string, 2));
      BCL_ExampleCheck( noninclusive_output, correct_noninclusive_output);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWithinRangeSmooth

  const ExampleClass::EnumType ExampleDescriptorWithinRangeSmooth::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWithinRangeSmooth())
  );

} // namespace bcl
