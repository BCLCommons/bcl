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
#include "contact/bcl_contact_calculate_correlations_mi.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_calculate_correlations_mi.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author teixeipl
  //! @date Aug 17, 2011
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactCalculateCorrelationsMI :
    public ExampleInterface
  {
  public:

    ExampleContactCalculateCorrelationsMI *Clone() const
    {
      return new ExampleContactCalculateCorrelationsMI( *this);
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

      // Test default constructor
      contact::CalculateCorrelationsMI test_calculator;

      // Test Clone further down

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // Create alignments to run on CalculateCorrelationsMI

      // Create correlation matrix

      // Print correlation matrix

      // Check correlation matrix values

    //////////////////////
    // input and output //
    //////////////////////

      // Test Write and Read

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactCalculateCorrelationsMI

  const ExampleClass::EnumType ExampleContactCalculateCorrelationsMI::s_Instance
  (
    GetExamples().AddEnum( ExampleContactCalculateCorrelationsMI())
  );

} // namespace bcl
