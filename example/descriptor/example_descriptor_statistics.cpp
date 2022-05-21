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
#include "descriptor/bcl_descriptor_statistics.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_statistics.cpp
  //!
  //! @author mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorStatistics :
    public ExampleInterface
  {
  public:

    ExampleDescriptorStatistics *Clone() const
    {
      return new ExampleDescriptorStatistics( *this);
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

      // create a string for the window descriptor that will be used internally
      const std::string window_desc( "(ReflectingWindow( AlphabeticNumber, size = 1))");

      // the template parameter for this class are relatively long; and because the class will always be used behind
      // an implementation, it makes more sense to test that functionality
      util::Implementation< descriptor::Base< char, float> > mean_alpha( "DescriptorMean" + window_desc);
      util::Implementation< descriptor::Base< char, float> > min_alpha( "DescriptorMin" + window_desc);
      util::Implementation< descriptor::Base< char, float> > max_alpha( "DescriptorMax" + window_desc);
      util::Implementation< descriptor::Base< char, float> > std_alpha( "DescriptorStandardDeviation" + window_desc);
      util::Implementation< descriptor::Base< char, float> > sum_alpha( "DescriptorSum"+ window_desc);

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
      BCL_ExampleCheck( sum_alpha->GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      // test the operators

      // check the results on each string

      // mean
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( mean_alpha, "asdf", 2),
        "13.00 ; 8.00 ; 9.67 ; 4.67 ; "
      );

      // sum
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sum_alpha, "asdf", 2),
        "39.00 ; 24.00 ; 29.00 ; 14.00 ; "
      );

      // min
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( min_alpha, "asdf", 2),
        "1.00 ; 1.00 ; 4.00 ; 4.00 ; "
      );

      // max
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( max_alpha, "asdf", 2),
        "19.00 ; 19.00 ; 19.00 ; 6.00 ; "
      );

      // std
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( std_alpha, "asdf", 2),
        "8.49 ; 7.87 ; 6.65 ; 0.94 ; "
      );

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

  }; //end ExampleDescriptorStatistics

  const ExampleClass::EnumType ExampleDescriptorStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorStatistics())
  );

} // namespace bcl
