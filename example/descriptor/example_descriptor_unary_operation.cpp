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
#include "descriptor/bcl_descriptor_unary_operation.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_unary_operation.cpp
  //!
  //! @author mendenjl
  //! @date Mar 11, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorUnaryOperation :
    public ExampleInterface
  {
  public:

    ExampleDescriptorUnaryOperation *Clone() const
    {
      return new ExampleDescriptorUnaryOperation( *this);
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
      util::Implementation< descriptor::Base< char, float> > alpha_square( "Sqr(AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > alpha_square_sqrt( "Sqrt(Sqr(AlphabeticNumber))");
      util::Implementation< descriptor::Base< char, float> > alpha_negative( "Negative(AlphabeticNumber)");

    /////////////////
    // data access //
    /////////////////

      // ensure that all implementations could be created
      BCL_ExampleIndirectAssert
      (
        alpha_square.IsDefined() && alpha_square_sqrt.IsDefined()
        && alpha_negative.IsDefined(),
        true,
        "Implementation constructor"
      );

      // test feature size and dimension functions
      BCL_ExampleCheck( alpha_square->GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( alpha_negative->GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      // initialize strings to be passed
      const std::string asdf( "asdf");

      // test the operators

      // square
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_square, asdf, 1),
        "1.0 ; 361.0 ; 16.0 ; 36.0 ; "
      );

      // Square/Sqrt
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_square_sqrt, asdf, 1),
        "1.0 ; 19.0 ; 4.0 ; 6.0 ; "
      );

      // double negative
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alpha_negative, asdf, 1),
        "-1.0 ; -19.0 ; -4.0 ; -6.0 ; "
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleUnaryOperation

  const ExampleClass::EnumType ExampleDescriptorUnaryOperation::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorUnaryOperation())
  );
} // namespace bcl
