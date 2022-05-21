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
#include "math/bcl_math_running_average_sd.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_running_average_sd.cpp
  //!
  //! @author mendenjl
  //! @date Apr 09, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRunningAverageSD :
    public ExampleInterface
  {
  public:

    ExampleMathRunningAverageSD *Clone() const
    {
      return new ExampleMathRunningAverageSD( *this);
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
      math::RunningAverageSD< double> running_var_ave_double;

      // initialize a running average of integers, but average in the integers as doubles to prevent roundoff
      math::RunningAverageSD< double> running_var_ave_ints;

      // average a series of vector 3d's
      math::RunningAverageSD< linal::Vector3D> running_var_ave_vector_3d;

      // average a series of vectors with dimension two
      math::RunningAverageSD< linal::Vector< double> > running_var_ave_vector_2d;

    /////////////////
    // data access //
    /////////////////

      // test the default constructors
      BCL_ExampleCheck( running_var_ave_double.GetAverage(), 0.0);
      BCL_ExampleCheck( running_var_ave_double.GetWeight(), 0.0);
      BCL_ExampleCheck( running_var_ave_vector_3d.GetAverage(), linal::Vector3D( 0.0, 0.0, 0.0));
      BCL_ExampleCheck( running_var_ave_vector_2d.GetAverage(), linal::Vector< double>());
      BCL_ExampleCheck( running_var_ave_vector_2d.GetVariance(), linal::Vector< double>());

    ///////////////
    // operators //
    ///////////////

      // test the operators on running average of doubles
      // test += and -=
      BCL_ExampleCheck( ( running_var_ave_double += 5.0).GetAverage(), 5.0);
      BCL_ExampleCheck( running_var_ave_double.GetVariance(), 0.0);
      BCL_ExampleCheck( ( running_var_ave_double += 10.0).GetAverage(), 7.5);
      BCL_ExampleCheck( running_var_ave_double.GetVariance(), 6.25);
      BCL_ExampleCheck( ( running_var_ave_double -= 5.0).GetAverage(), 10.0);
      BCL_ExampleCheck( running_var_ave_double.GetVariance(), 0.0);

      // average 1, 2, 3, 4, 5
      running_var_ave_double.Reset();
      for( size_t i( 1); i <= 5; ++i)
      {
        running_var_ave_double += double( i);
      }
      BCL_ExampleIndirectCheck
      (
        running_var_ave_double.GetAverage(),
        3.0,
        "running_var_ave_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );
      BCL_ExampleIndirectCheck
      (
        running_var_ave_double.GetVariance(),
        2.0,
        "running_var_ave_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );

      // remove 1, so now the average is of 2, 3, 4, 5
      BCL_ExampleCheck( ( running_var_ave_double -= 1.0).GetAverage(), 3.5);
      BCL_ExampleIndirectCheck( running_var_ave_double.GetVariance(), 1.25, "operator-=");
      BCL_ExampleIndirectCheck( running_var_ave_double.GetWeight(), 4.0, "operator-=");

      // add 5.0 with a weight of 2.5, so now the average is of 2, 3, 4, 5 x 3.5 times
      BCL_ExampleCheck
      (
        running_var_ave_double.AddWeightedObservation( 5.0, 2.5).GetAverage(),
        ( 2.0 + 3.0 + 4.0 + 5.0 * 3.5) / 6.5
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          running_var_ave_double.GetVariance(),
          1.3017751479,
          1e-7
        ),
        true,
        "AddWeightedObservation"
      );

      // create some vector 3d's to average
      linal::Vector3D x1( 1.0, 0.0, 0.0), y2z3( 0.0, 2.0, 3.0), yn2zn3( 0.0, -2.0, -3.0);
      BCL_ExampleCheck
      (
        ( running_var_ave_vector_3d += x1).GetVariance(),
        linal::Vector3D( 0.0)
      );

      // reset the running average and average some other vectors
      running_var_ave_vector_3d.Reset();
      running_var_ave_vector_3d += y2z3;
      running_var_ave_vector_3d += y2z3;
      running_var_ave_vector_3d += y2z3;
      BCL_ExampleIndirectCheck
      (
        running_var_ave_vector_3d.GetAverage(),
        y2z3,
        "Reset"
      );
      BCL_ExampleIndirectCheck
      (
        running_var_ave_vector_3d.GetVariance(),
        linal::Vector3D( 0.0),
        "Reset"
      );

      // reset the running average and average some other vectors
      running_var_ave_vector_3d.Reset();
      running_var_ave_vector_3d += y2z3;
      running_var_ave_vector_3d += x1;
      running_var_ave_vector_3d += yn2zn3;
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          running_var_ave_vector_3d.GetAverage(),
          x1 / 3.0,
          1e-7
        ),
        true
      );
      const linal::Vector3D expected_variance
      (
        linal::Vector3D( double( 2.0) / double( 9.0), double( 8.0) / double( 3.0), double( 6.0))
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          running_var_ave_vector_3d.GetVariance(),
          expected_variance,
          1e-7
        ),
        true,
        "Adding up number preserves variances; got " + util::Format()( running_var_ave_vector_3d.GetVariance())
        + " expected: " + util::Format()( expected_variance)
      );

      // create some vectors
      linal::Vector< double> x2( linal::MakeVector< double>( 2.0, 0.0));
      linal::Vector< double> y1( linal::MakeVector< double>( 0.0, 1.0));
      linal::Vector< double> yn1( linal::MakeVector< double>( 0.0, -1.0));
      BCL_ExampleCheck
      (
        linal::Vector< double>( ( running_var_ave_vector_2d += x2).GetVariance()),
        linal::MakeVector< double>( 0.0, 0.0)
      );

      // test reset, which should clear the weight (but not necessarily the average)
      running_var_ave_vector_2d.Reset();
      BCL_ExampleIndirectCheck
      (
        linal::Vector< double>( ( running_var_ave_vector_2d += y1).GetVariance()),
        linal::MakeVector< double>( 0.0, 0.0),
        "Reset"
      );

      running_var_ave_vector_2d.Reset();
      running_var_ave_vector_2d += x2;
      running_var_ave_vector_2d += y1;
      running_var_ave_vector_2d += yn1;
      const linal::Vector< double> expected_variances_2d
      (
        linal::MakeVector( double( 8.0) / double( 9.0), double( 2.0) / double( 3.0))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          running_var_ave_vector_2d.GetAverage(),
          linal::Vector< double>( x2 / 3.0),
          1e-7
        ),
        true
      );
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          running_var_ave_vector_2d.GetVariance(),
          expected_variances_2d,
          1e-7
        ),
        true,
        "Adding up number preserves variances; got " + util::Format()( running_var_ave_vector_2d.GetVariance())
        + " expected: " + util::Format()( expected_variances_2d)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( running_var_ave_double, math::RunningAverageSD< double>()),
        true,
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleMathRunningAverageSD

  const ExampleClass::EnumType ExampleMathRunningAverageSD::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRunningAverageSD())
  );
} // namespace bcl

