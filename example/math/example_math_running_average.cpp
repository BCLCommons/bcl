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
#include "math/bcl_math_running_average.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_running_average.cpp
  //!
  //! @author mendenjl
  //! @date Dec 04, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRunningAverage :
    public ExampleInterface
  {
  public:

    ExampleMathRunningAverage *Clone() const
    {
      return new ExampleMathRunningAverage( *this);
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

      // initialize a running average of doubles
      math::RunningAverage< double> running_average_double;

      // initialize a running average of integers, but average in the integers as doubles to prevent roundoff
      math::RunningAverage< double> running_average_ints;

      // average a series of vector 3d's
      math::RunningAverage< linal::Vector3D> running_average_vector_3d;

      // average a series of vectors with dimension two
      math::RunningAverage< linal::Vector< double> > running_average_vector_2d;

    /////////////////
    // data access //
    /////////////////

      // test the default constructors
      BCL_ExampleCheck( running_average_double.GetAverage(), 0.0);
      BCL_ExampleCheck( running_average_double.GetWeight(), 0.0);
      BCL_ExampleCheck( running_average_vector_3d.GetAverage(), linal::Vector3D( 0.0, 0.0, 0.0));
      BCL_ExampleCheck( running_average_vector_2d.GetAverage(), linal::Vector< double>());

    ///////////////
    // operators //
    ///////////////

      // test the operators on running average of doubles
      // test += and -=
      BCL_ExampleCheck( running_average_double += 5.0, 5.0);
      BCL_ExampleCheck( running_average_double += 10.0, 7.5);
      BCL_ExampleCheck( running_average_double -= 5.0, 10.0);

      // average 1, 2, 3, 4, 5
      running_average_double.Reset();
      for( size_t i( 1); i <= 5; ++i)
      {
        running_average_double += double( i);
      }
      BCL_ExampleIndirectCheck
      (
        running_average_double.GetAverage(),
        3.0,
        "running_average_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );
      BCL_ExampleIndirectCheck
      (
        running_average_double.GetWeight(),
        5.0,
        "running_average_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );

      // remove 1, so now the average is of 2, 3, 4, 5
      BCL_ExampleCheck( running_average_double -= 1.0, 3.5);
      BCL_ExampleIndirectCheck( running_average_double.GetWeight(), 4.0, "operator-=");

      // add 5.0 with a weight of 2.5, so now the average is of 2, 3, 4, 5 x 3.5 times
      BCL_ExampleCheck( running_average_double.AddWeightedObservation( 5.0, 2.5), ( 2.0 + 3.0 + 4.0 + 5.0 * 3.5) / 6.5);
      BCL_ExampleIndirectCheck( running_average_double.GetWeight(), 6.5, "AddWeightedObservation");

      // do operator+= with ints
      BCL_ExampleCheck( running_average_ints += 5, 5);
      BCL_ExampleCheck( running_average_ints += 10, 7.5);
      BCL_ExampleCheck( running_average_ints -= 5, 10.0);

      // create some vector 3d's to average
      linal::Vector3D x1( 1.0, 0.0, 0.0), y2z3( 0.0, 2.0, 3.0), yn2zn3( 0.0, -2.0, -3.0);
      BCL_ExampleCheck( linal::Vector3D( running_average_vector_3d += x1), x1);

      // reset the running average and average some other vectors
      running_average_vector_3d.Reset();
      running_average_vector_3d += y2z3;
      running_average_vector_3d += y2z3;
      running_average_vector_3d += y2z3;
      BCL_ExampleIndirectCheck
      (
        running_average_vector_3d.GetAverage(),
        y2z3,
        "Reset"
      );

      // reset the running average and average some other vectors
      running_average_vector_3d.Reset();
      running_average_vector_3d += y2z3;
      running_average_vector_3d += x1;
      running_average_vector_3d += yn2zn3;
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          running_average_vector_3d.GetAverage(),
          x1 / 3.0,
          1e-7
        ),
        true
      );

      // create some vectors
      linal::Vector< double> x2( linal::MakeVector< double>( 2.0, 0.0));
      linal::Vector< double> y1( linal::MakeVector< double>( 0.0, 1.0));
      linal::Vector< double> yn1( linal::MakeVector< double>( 0.0, -1.0));
      BCL_ExampleCheck( linal::Vector< double>( running_average_vector_2d += x2), x2);

      // test reset, which should clear the weight (but not necessarily the average)
      running_average_vector_2d.Reset();
      BCL_ExampleIndirectCheck( linal::Vector< double>( running_average_vector_2d += y1), y1, "Reset");

      running_average_vector_2d.Reset();
      running_average_vector_2d += x2;
      running_average_vector_2d += y1;
      running_average_vector_2d += yn1;
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          running_average_vector_2d.GetAverage(),
          linal::Vector< double>( x2 / 3.0),
          1e-7
        ),
        true
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( running_average_double, math::RunningAverage< double>()),
        true,
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleMathRunningAverage

  const ExampleClass::EnumType ExampleMathRunningAverage::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRunningAverage())
  );
} // namespace bcl

