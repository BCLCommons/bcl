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
#include "descriptor/bcl_descriptor_sigmoid.h"

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
  //! @example example_descriptor_sigmoid.cpp
  //!
  //! @author mendenjl
  //! @date Jun 21, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSigmoid :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSigmoid *Clone() const
    {
      return new ExampleDescriptorSigmoid( *this);
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
      descriptor::Sigmoid< char> sigmoid_def;

      // Test a normal range
      descriptor::Sigmoid< char> sigmoid_normal( 1.0, 0.0, 0.0, 1.0, 0.73);

      // Test a zero-width range on the left
      descriptor::Sigmoid< char> sigmoid_mod( 5.0, 9.0, 1.0, 0.0, 2.0);

    /////////////////
    // data access //
    /////////////////

      // The property used to test the range
      util::Implementation< descriptor::Base< char, float> >
              str_descriptor( "AlphabeticNumber");

      sigmoid_def.SetDescriptor( str_descriptor);
      sigmoid_normal.SetDescriptor( str_descriptor);
      sigmoid_mod.SetDescriptor( str_descriptor);

      BCL_ExampleCheck( sigmoid_normal.GetDescriptor(), str_descriptor);

    //////////////////////
    // input and output //
    //////////////////////

      // String used to test the descriptor
      std::string test_string( "ABCDEFGHIJKLP");

      // Test the normal range object
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sigmoid_normal, "ABCDEF", 2),
        "0.73 ; 0.88 ; 0.95 ; 0.98 ; 0.99 ; 1.00 ; "
      );

      // Test the noleft range object
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sigmoid_def, "ABCDEF", 2),
        "0.73 ; 0.88 ; 0.95 ; 0.98 ; 0.99 ; 1.00 ; "
      );

      // Test the noright range object
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sigmoid_mod, "ABCDEFGHIJKLP", 2),
        "2.13 ; 2.27 ; 2.42 ; 2.58 ; 2.75 ; 2.93 ; 3.12 ; 3.31 ; 3.50 ; 3.69 ; 3.88 ; 4.07 ; 4.73 ; "
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorSigmoid

  const ExampleClass::EnumType ExampleDescriptorSigmoid::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSigmoid())
  );

} // namespace bcl
