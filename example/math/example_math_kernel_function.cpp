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
#include "math/bcl_math_kernel_function.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_kernel_function.cpp
  //!
  //! @author mueller
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathKernelFunction :
    public ExampleInterface
  {
  public:

    ExampleMathKernelFunction *Clone() const
    { return new ExampleMathKernelFunction( *this);}

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

      // default constructor
      math::KernelFunction kernel_default;

      // construct from range
      math::KernelFunction kernel( math::Range< double>( 0.0, 1.0));

      BCL_Example_Check
      (
        kernel.GetKernel().GetMin() == 0
        && kernel.GetKernel().GetMax() == 1,
        "Constructor didn't initialize the kernel correctly to [0,1]."
      );

      // clone
      util::ShPtr< math::KernelFunction> kernel_cloned( kernel.Clone());

      BCL_Example_Check
      (
        kernel_cloned->GetKernel().GetMin() == 0
        && kernel_cloned->GetKernel().GetMax() == 1,
        "Clone didn't initialize the kernel correctly to [0,1]."
      );

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + kernel_cloned->GetClassIdentifier());
      BCL_Example_Check
      (
        GetStaticClassName< math::KernelFunction>() == "bcl::math::KernelFunction"
        && kernel_cloned->GetClassIdentifier() == GetStaticClassName< math::KernelFunction>(),
        "incorrect class name"
      );

      // set range of kernel
      kernel_default.SetKernel( math::Range< double>( math::RangeBorders::e_LeftClosed, 0, 1, math::RangeBorders::e_RightOpen));

      // get range of kernel
      BCL_MessageStd( "class name: " + util::Format()( kernel_default.GetKernel()));
      BCL_Example_Check
      (
        kernel_default.GetKernel().GetMin() == 0
        && kernel_default.GetKernel().GetMax() == 1,
        "Range of the kernel should be [0, 1) but is " + util::Format()( kernel_default.GetKernel())
      );

    ///////////////
    // operators //
    ///////////////

      //! calculate f(x)
      const double y_1( kernel_default( 1));
      BCL_MessageStd
      (
        "function value of a kernel at x = 1: " + util::Format()( y_1)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_1, 0),
        "incorrect result for x = 1: " + util::Format()( y_1) + " != 0."
      );

      const double y_05( kernel_default( 0.5));
      BCL_MessageStd
      (
        "function value of a kernel at x = 0.5: "
        + util::Format()( y_05)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_05, 1),
        "incorrect result for x = 0.5: " + util::Format()( y_05) +
        " != 1."
      );

      const double y_0( kernel_default( 0));
      BCL_MessageStd
      (
        "function value of a kernel at x = 0: "
        + util::Format()( y_0)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_0, 1),
        "incorrect result for x = 0: " + util::Format()( y_0) +
        " != 1."
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathKernelFunction

  const ExampleClass::EnumType ExampleMathKernelFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathKernelFunction())
  );

} // namespace bcl
