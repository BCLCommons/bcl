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
#include "math/bcl_math_comparisons.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_comparisons.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathComparisons :
    public ExampleInterface
  {
  public:

    ExampleMathComparisons *Clone() const
    {
      return new ExampleMathComparisons( *this);
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

    /////////////////
    // data access //
    /////////////////

      const util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
        greater( *math::Comparisons< double>::GetEnums().e_Greater);

      const util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
        greater_equal( *math::Comparisons< double>::GetEnums().e_GreaterEqual);

      const util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
        less( *math::Comparisons< double>::GetEnums().e_Less);

      const util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
        less_equal( *math::Comparisons< double>::GetEnums().e_LessEqual);

      BCL_Example_Check
      (
        math::Comparisons< double>::GetEnums().GetClassIdentifier() == "bcl::math::Comparisons<double>",
        "GetClassIdentifier() should give back bcl::math::Comparisons but gives back "
        + math::Comparisons< double>::GetEnums().GetClassIdentifier()
      );

      BCL_Example_Check
      (
        !math::Comparisons< double>::GetUndefinedData().IsDefined(),
        "GetUndefinedData() should give undefined data."
      );

    ///////////////
    // operators //
    ///////////////

      // test 'greater' binary function
      BCL_Example_Check
      (
        ( *greater)( 5, 3),
        "5 should be greater than 3."
      );

      BCL_Example_Check
      (
        !( *greater)( 3, 3),
        "3 should not be greater than 3."
      );

      BCL_Example_Check
      (
        !( *greater)( 3, 5),
        "3 should not be greater than 5."
      );

      // test 'greater equal' binary function
      BCL_Example_Check
      (
        ( *greater_equal)( 5, 3),
        "5 should be greater or equal than 3."
      );

      BCL_Example_Check
      (
        ( *greater_equal)( 3, 3),
        "3 should be greater or equal than 3."
      );

      BCL_Example_Check
      (
        !( *greater_equal)( 3, 5),
        "3 should not be greater or equal than 5."
      );

      // test 'less' binary function
      BCL_Example_Check
      (
        ( *less)( 3, 5),
        "3 should be less than 5."
      );

      BCL_Example_Check
      (
        !( *less)( 3, 3),
        "3 should not be less than 3."
      );

      BCL_Example_Check
      (
        !( *less)( 5, 3),
        "5 should not be less than 3."
      );

      // test 'less equal' binary function
      BCL_Example_Check
      (
        ( *less_equal)( 3, 5),
        "3 should be less or equal than 5."
      );

      BCL_Example_Check
      (
        ( *less_equal)( 3, 3),
        "3 should  be less or equal than 3."
      );

      BCL_Example_Check
      (
        !( *less_equal)( 5, 3),
        "5 should not be less or equal than 3."
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

  }; //end ExampleMathComparisons

  const ExampleClass::EnumType ExampleMathComparisons::s_Instance
  (
    GetExamples().AddEnum( ExampleMathComparisons())
  );

} // namespace bcl
