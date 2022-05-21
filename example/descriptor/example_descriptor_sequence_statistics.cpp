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
#include "descriptor/bcl_descriptor_sequence_statistics.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_sequence_statistics.cpp
  //!
  //! @author mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSequenceStatistics :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSequenceStatistics *Clone() const
    {
      return new ExampleDescriptorSequenceStatistics( *this);
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

      // the template parameter for this class are relatively long; and because the class will always be used behind
      // an implementation, it makes more sense to test that functionality
      util::Implementation< descriptor::Base< char, float> > mean_alpha( "StringMean(AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > min_alpha(  "StringMin(AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > max_alpha(  "StringMax(AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > std_alpha(  "StringStandardDeviation(AlphabeticNumber)");
      util::Implementation< descriptor::Base< char, float> > sum_alpha(  "StringSum(AlphabeticNumber)");

    /////////////////
    // data access //
    /////////////////

      // ensure that all implementations could be created
      BCL_ExampleIndirectAssert
      (
        mean_alpha.IsDefined() && min_alpha.IsDefined()
        && max_alpha.IsDefined() && std_alpha.IsDefined()
        && sum_alpha.IsDefined(),
        true,
        "Implementation constructor"
      );

      // test feature size and dimension functions
      BCL_ExampleCheck( mean_alpha->GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( sum_alpha->GetType().GetDimension(), 0);

    ///////////////
    // operators //
    ///////////////

      // initialize strings to be passed
      const std::string string_a( "asdf"), string_1234as( size_t( 1234), 'a');

      // test the operators

      // check the results on each string

      // mean
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( mean_alpha, string_a, 2), "7.50 ; ");
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( mean_alpha, string_1234as, 2), "1.00 ; ");

      // sum
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( sum_alpha, string_a, 2), "30.00 ; ");
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( sum_alpha, string_1234as, 2), "1234.00 ; ");

      // min
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( min_alpha, string_a, 2), "1.00 ; ");
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( min_alpha, string_1234as, 2), "1.00 ; ");

      // max
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( max_alpha, string_a, 2), "19.00 ; ");
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( max_alpha, string_1234as, 2), "1.00 ; ");

      // std
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( std_alpha, string_a, 2), "6.87 ; ");
      BCL_ExampleCheck( descriptor::StringSequence::WriteIterations( std_alpha, string_1234as, 2), "0.00 ; ");

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

  }; //end ExampleDescriptorSequenceStatistics

  const ExampleClass::EnumType ExampleDescriptorSequenceStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSequenceStatistics())
  );

} // namespace bcl
