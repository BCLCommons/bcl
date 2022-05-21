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
#include "math/bcl_math_running_sum.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_running_sum.cpp
  //!
  //! @author mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRunningSum :
    public ExampleInterface
  {
  public:

    ExampleMathRunningSum *Clone() const
    {
      return new ExampleMathRunningSum( *this);
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

      // initialize a running sum of doubles
      math::RunningSum< double> running_sum_double;

      // initialize a running sum of integers, but sum in the integers as doubles to prevent roundoff
      math::RunningSum< int> running_sum_ints;

      // element-wise min/max over a series of vector 3d's
      math::RunningSum< linal::Vector3D> running_sum_vector_3d( linal::Vector3D( 0.0, 0.0, 0.0));

      // element-wise min/max a series of vectors with dimension two
      math::RunningSum< linal::Vector< double> > running_var_ave_vector_2d( linal::Vector< double>( 2, 0.0));

    /////////////////
    // data access //
    /////////////////

      // test the default constructors
      BCL_ExampleCheck( running_sum_double.GetSum(), 0.0);
      BCL_ExampleCheck( running_sum_ints.GetSum(), 0);
      BCL_ExampleCheck( running_sum_vector_3d.GetSum(), linal::Vector3D( 0.0));
      BCL_ExampleCheck( running_var_ave_vector_2d.GetSum(), linal::Vector< double>( 2, 0.0));

    ///////////////
    // operators //
    ///////////////

      // test the operators on running sum of doubles
      // test += and -=
      BCL_ExampleCheck( ( running_sum_double += 5.1).GetSum(), 5.1);
      BCL_ExampleCheck( running_sum_double.GetSum(), 5.1);
      BCL_ExampleCheck( ( running_sum_double += 10.0).GetSum(), 15.1);
      BCL_ExampleCheck( running_sum_double.GetSum(), 15.1);
      BCL_ExampleCheck( ( running_sum_double += -5.0).GetSum(), 10.1);

      // sum 1, 2, 3, 4, 5
      running_sum_double.Reset();
      BCL_ExampleCheck( running_sum_double.GetSum(), 0.0);
      for( size_t i( 1); i <= 5; ++i)
      {
        running_sum_double += double( i);
      }
      BCL_ExampleIndirectCheck
      (
        running_sum_double.GetSum(),
        15.0,
        "running_sum_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );

      // create some vector 3d's to sum
      linal::Vector3D x1( 1.0, -1.5, 0.1), y2z3( 0.0, 2.0, 3.0), yn2zn3( 0.0, -2.0, -3.0);
      BCL_ExampleCheck
      (
        ( running_sum_vector_3d += x1).GetSum(),
        x1
      );
      BCL_ExampleCheck
      (
        ( running_sum_vector_3d += x1).GetSum(),
        x1 * 2.0
      );

      // reset the running sum and add some other vectors
      running_sum_vector_3d.Reset();
      running_sum_vector_3d += y2z3;
      running_sum_vector_3d += y2z3;
      running_sum_vector_3d -= y2z3;
      BCL_ExampleIndirectCheck
      (
        running_sum_vector_3d.GetSum(),
        y2z3,
        "Reset"
      );

      // reset the running status and use the full set of vectors
      running_sum_vector_3d.Reset();
      running_sum_vector_3d += y2z3;
      running_sum_vector_3d += x1;
      running_sum_vector_3d += yn2zn3;
      BCL_ExampleIndirectCheck
      (
        running_sum_vector_3d.GetSum(),
        yn2zn3 + x1 + y2z3,
        "Adding vectors to running sum"
      );

      // create some vectors
      linal::Vector< double> x2( linal::MakeVector< double>( 2.0, -1.0));
      linal::Vector< double> y1( linal::MakeVector< double>( 0.0, 1.0));
      linal::Vector< double> yn1( linal::MakeVector< double>( 0.0, -1.0));
      BCL_ExampleCheck
      (
        ( running_var_ave_vector_2d += x2).GetSum(),
        x2
      );
      BCL_ExampleCheck
      (
        ( running_var_ave_vector_2d += x2).GetSum(),
        x2 * 2.0
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( running_sum_double, math::RunningSum< double>()),
        true,
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleMathRunningSum

  const ExampleClass::EnumType ExampleMathRunningSum::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRunningSum())
  );
} // namespace bcl

