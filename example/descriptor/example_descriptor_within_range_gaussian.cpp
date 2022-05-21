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
#include "descriptor/bcl_descriptor_within_range_gaussian.h"

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
  //! @example example_descriptor_within_range_gaussian.cpp
  //!
  //! @author geanesar
  //! @date Jan 28, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWithinRangeGaussian :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWithinRangeGaussian *Clone() const
    {
      return new ExampleDescriptorWithinRangeGaussian( *this);
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
      descriptor::WithinRangeGaussian< char> within_range_gauss_def;

      // Test a normal range
      descriptor::WithinRangeGaussian< char> within_range_gauss_normal( 0.0, 1.0);

      // Test a zero-width range on the left
      descriptor::WithinRangeGaussian< char> within_range_gauss_broad( 0.0, 2.0);

      // Test a zero-width range on the right
      descriptor::WithinRangeGaussian< char> within_range_gauss_shift( 1.0, 1.0);

      // test Clone function
      util::ShPtr< descriptor::WithinRangeGaussian< char> > sp_within_range_gauss
      (
        within_range_gauss_normal.Clone()
      );
      BCL_ExampleCheck( sp_within_range_gauss.IsDefined(), true);

      // Destroy the new object
      sp_within_range_gauss = util::ShPtr< descriptor::WithinRangeGaussian< char> >();

    /////////////////
    // data access //
    /////////////////

      // Check the class name is correct
      BCL_ExampleCheck
      (
        GetStaticClassName< descriptor::WithinRangeGaussian< char> >(),
        within_range_gauss_normal.GetClassIdentifier()
      );
      BCL_MessageStd( "class name: " + within_range_gauss_normal.GetClassIdentifier());

      // The property used to test the range
      util::Implementation< descriptor::Base< char, float> >
              str_descriptor( "AlphabeticNumber");

      within_range_gauss_normal.SetDescriptor( str_descriptor);
      within_range_gauss_broad.SetDescriptor( str_descriptor);
      within_range_gauss_shift.SetDescriptor( str_descriptor);

      BCL_ExampleCheck( within_range_gauss_normal.GetDescriptor(), str_descriptor);

    //////////////////////
    // input and output //
    //////////////////////

      // String used to test the descriptor
      std::string test_string( "ABCDEF");

      // Test the normal range object
      std::string correct_normal_output( "0.61 ; 0.14 ; 0.01 ; 0.00 ; 0.00 ; 0.00 ; ");
      std::string normal_output( descriptor::StringSequence::WriteIterations( within_range_gauss_normal, test_string, 2));
      BCL_ExampleCheck( normal_output, correct_normal_output);

      // Test the broad range object
      std::string correct_broad_output( "0.14 ; 0.00 ; 0.00 ; 0.00 ; 0.00 ; 0.00 ; ");
      std::string broad_output( descriptor::StringSequence::WriteIterations( within_range_gauss_broad, test_string, 2));
      BCL_ExampleCheck( broad_output, correct_broad_output);

      // Test the shift range object
      std::string correct_shift_output( "1.00 ; 0.61 ; 0.14 ; 0.01 ; 0.00 ; 0.00 ; ");
      std::string shift_output( descriptor::StringSequence::WriteIterations( within_range_gauss_shift, test_string, 2));
      BCL_ExampleCheck( shift_output, correct_shift_output);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWithinRangeGaussian

  const ExampleClass::EnumType ExampleDescriptorWithinRangeGaussian::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWithinRangeGaussian())
  );

} // namespace bcl
