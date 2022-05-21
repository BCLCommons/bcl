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
#include "descriptor/bcl_descriptor_window_weights.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_window_weights.cpp
  //!
  //! @author mendenjl
  //! @date Mar 15, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWindowWeights :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWindowWeights *Clone() const
    {
      return new ExampleDescriptorWindowWeights( *this);
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

      // create an instance of the weights class
      descriptor::WindowWeights rectangular_window( descriptor::WindowWeights::e_Rectangle);
      descriptor::WindowWeights triangular_window( descriptor::WindowWeights::e_Triangular);
      descriptor::WindowWeights welch_window( descriptor::WindowWeights::e_Welch);
      descriptor::WindowWeights gaussian_window( descriptor::WindowWeights::e_Gaussian);
      descriptor::WindowWeights hamming_window( descriptor::WindowWeights::e_Hamming);
      descriptor::WindowWeights hann_window( descriptor::WindowWeights::e_Hann);

      BCL_ExampleCheck( rectangular_window( 5).Sum(), float( 5.0));
      BCL_ExampleAssert( rectangular_window( 5).GetSize(), 5);
      BCL_ExampleCheck( rectangular_window( 5)( 4), float( 1.0));

      BCL_ExampleCheckWithinTolerance( triangular_window( 5).Sum(), float( 3.0), 0.0001);
      BCL_ExampleAssert( triangular_window( 5).GetSize(), 5);
      BCL_ExampleCheckWithinTolerance( triangular_window( 5)( 4), float( 0.2), 0.0001);

      BCL_ExampleCheckWithinTolerance( welch_window( 5).Sum(), float( 3.8), 0.0001);
      BCL_ExampleAssert( welch_window( 5).GetSize(), 5);
      BCL_ExampleCheckWithinTolerance( welch_window( 5)( 4), float( 9.0) / float( 25.0), 0.0001);

      BCL_ExampleCheckWithinTolerance( gaussian_window( 5).Sum(), float( 3.17746), 0.0001);
      BCL_ExampleAssert( gaussian_window( 5).GetSize(), 5);
      BCL_ExampleCheckWithinTolerance( gaussian_window( 5)( 4), float( 0.2), 0.0001);

      BCL_ExampleCheckWithinTolerance( hamming_window( 5).Sum(), float( 2.71739), 0.0001);
      BCL_ExampleAssert( hamming_window( 5).GetSize(), 5);
      BCL_ExampleCheckWithinTolerance( hamming_window( 5)( 4), float( 4.0) / float( 46.0), 0.0001);

      BCL_ExampleCheckWithinTolerance( hann_window( 5).Sum(), float( 2.5), 0.0001);
      BCL_ExampleAssert( hann_window( 5).GetSize(), 5);
      BCL_ExampleCheckWithinAbsTolerance( hann_window( 5)( 4), float( 0.0), 0.0001);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWindowWeights

  const ExampleClass::EnumType ExampleDescriptorWindowWeights::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWindowWeights())
  );
} // namespace bcl
