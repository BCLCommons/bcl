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
#include "descriptor/bcl_descriptor_binary_operation.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_binary_operation.cpp
  //!
  //! @author mendenjl
  //! @date Feb 07, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorBinaryOperation :
    public ExampleInterface
  {
  public:

    ExampleDescriptorBinaryOperation *Clone() const
    {
      return new ExampleDescriptorBinaryOperation( *this);
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

      // because the class will always be used behind
      // an implementation (it lacks conventional constructors), it makes more sense to test that functionality
      util::Implementation< descriptor::Base< char, float> > alpha_plus_three( "Add(AlphabeticNumber,Constant(3))");
      util::Implementation< descriptor::Base< char, float> > alpha_squared_directly( "Multiply(AlphabeticNumber,AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > alpha_squared_via_exponentiate( "Exponentiate(lhs=AlphabeticNumber,rhs=Constant(2))");

    /////////////////
    // data access //
    /////////////////

      // ensure that all implementations could be created
      BCL_ExampleIndirectAssert
      (
        alpha_plus_three.IsDefined() && alpha_squared_directly.IsDefined()
        && alpha_squared_via_exponentiate.IsDefined(),
        true,
        "Implementation constructor"
      );

      // test feature size and dimension functions
      BCL_ExampleCheck( alpha_plus_three->GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( alpha_squared_directly->GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      // initialize strings to be passed
      const std::string asdf( "asdf");

      // test the operators

      // check the results on each string

      // add 3
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_plus_three, asdf, 1),
        "4.0 ; 22.0 ; 7.0 ; 9.0 ; "
      );

      // square via multiply
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_squared_directly, asdf, 1),
        "1.0 ; 361.0 ; 16.0 ; 36.0 ; "
      );

      // square via exponentiate should be identical
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_squared_via_exponentiate, asdf, 2),
        descriptor::StringSequence::WriteIterations( alpha_squared_directly, asdf, 2)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleBinaryOperation

  const ExampleClass::EnumType ExampleDescriptorBinaryOperation::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorBinaryOperation())
  );
} // namespace bcl
