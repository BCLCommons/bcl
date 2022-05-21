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
#include "math/bcl_math_sum_function_mixin.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_quadratic_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_sum_function_mixin.cpp
  //!
  //! @author mendenjl
  //! @date Sep 19, 2016
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathSumFunctionMixin :
    public ExampleInterface
  {

  public:

    ExampleMathSumFunctionMixin *Clone() const
    {
      return new ExampleMathSumFunctionMixin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////
    // data preparation //
    //////////////////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct an example function
      math::SumFunction< double, double> function( "test");
      function.NewOperand( math::LinearFunction( 1.0, 5.0), 1.0);
      function.NewOperand( math::QuadraticFunction( 1.0, 5.0, 2.0), 2.0);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

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

  }; //end ExampleMathSumFunctionMixin

  const ExampleClass::EnumType ExampleMathSumFunctionMixin::s_Instance
  (
    GetExamples().AddEnum( ExampleMathSumFunctionMixin())
  );

} // namespace bcl

