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
#include "math/bcl_math_running_min_max.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_running_min_max.cpp
  //!
  //! @author mendenjl
  //! @date Aug 17, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRunningMinMax :
    public ExampleInterface
  {
  public:

    ExampleMathRunningMinMax *Clone() const
    {
      return new ExampleMathRunningMinMax( *this);
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
      math::RunningMinMax< double> running_min_max_double;

      // initialize a running average of integers, but average in the integers as doubles to prevent roundoff
      math::RunningMinMax< int> running_min_max_ints;

      // element-wise min/max over a series of vector 3d's
      math::RunningMinMax< linal::Vector3D> running_min_max_vector_3d( linal::Vector3D( 0.0, 0.0, 0.0));

      // element-wise min/max a series of vectors with dimension two
      math::RunningMinMax< linal::Vector< double> > running_var_ave_vector_2d( linal::Vector< double>( 2, 0.0));

    /////////////////
    // data access //
    /////////////////

      // test the default constructors
      BCL_ExampleCheck( running_min_max_double.GetMin(), std::numeric_limits< double>::infinity());
      BCL_ExampleCheck( running_min_max_double.GetMax(), -std::numeric_limits< double>::infinity());
      BCL_ExampleCheck( running_min_max_ints.GetMin(), std::numeric_limits< int>::max());
      BCL_ExampleCheck( running_min_max_ints.GetMax(), std::numeric_limits< int>::min());
      BCL_ExampleCheck( running_min_max_vector_3d.GetMin(), linal::Vector3D( 0.0, 0.0, 0.0));
      BCL_ExampleCheck( running_min_max_vector_3d.GetMin(), running_min_max_vector_3d.GetMax());
      BCL_ExampleCheck( running_var_ave_vector_2d.GetMin(), linal::Vector< double>( 2, 0.0));

    ///////////////
    // operators //
    ///////////////

      // test the operators on running average of doubles
      // test += and -=
      BCL_ExampleCheck( ( running_min_max_double += 5.1).GetMax(), 5.1);
      BCL_ExampleCheck( running_min_max_double.GetMin(), 5.1);
      BCL_ExampleCheck( ( running_min_max_double += 10.0).GetMax(), 10.0);
      BCL_ExampleCheck( running_min_max_double.GetMin(), 5.1);
      BCL_ExampleCheck( ( running_min_max_double += -5.0).GetMax(), 10.0);
      BCL_ExampleCheck( running_min_max_double.GetMin(), -5.0);

      // average 1, 2, 3, 4, 5
      running_min_max_double.Reset();
      BCL_ExampleCheck( running_min_max_double.GetMin(), std::numeric_limits< double>::infinity());
      BCL_ExampleCheck( running_min_max_double.GetMax(), -std::numeric_limits< double>::infinity());
      for( size_t i( 1); i <= 5; ++i)
      {
        running_min_max_double += double( i);
      }
      BCL_ExampleIndirectCheck
      (
        running_min_max_double.GetMin(),
        1.0,
        "running_min_max_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );
      BCL_ExampleIndirectCheck
      (
        running_min_max_double.GetMax(),
        5.0,
        "running_min_max_double += 1.0 += 2.0 += 3.0 += 4.0 += 5.0"
      );

      // create some vector 3d's to average
      linal::Vector3D x1( 1.0, -1.5, 0.1), y2z3( 0.0, 2.0, 3.0), yn2zn3( 0.0, -2.0, -3.0);
      BCL_ExampleCheck
      (
        ( running_min_max_vector_3d += x1).GetMin(),
        linal::Vector3D( 0.0, -1.5, 0.0)
      );
      BCL_ExampleCheck
      (
        ( running_min_max_vector_3d += x1).GetMax(),
        linal::Vector3D( 1.0, 0.0, 0.1)
      );

      // reset the running average and average some other vectors
      running_min_max_vector_3d.Reset();
      running_min_max_vector_3d += y2z3;
      running_min_max_vector_3d += y2z3;
      running_min_max_vector_3d += y2z3;
      BCL_ExampleIndirectCheck
      (
        running_min_max_vector_3d.GetMin(),
        y2z3,
        "Reset"
      );
      BCL_ExampleIndirectCheck
      (
        running_min_max_vector_3d.GetMax(),
        y2z3,
        "Reset"
      );

      // reset the running status and use the full set of vectors
      running_min_max_vector_3d.Reset();
      running_min_max_vector_3d += y2z3;
      running_min_max_vector_3d += x1;
      running_min_max_vector_3d += yn2zn3;
      BCL_ExampleIndirectCheck
      (
        running_min_max_vector_3d.GetMin(),
        yn2zn3,
        "Adding vectors to running min/max"
      );
      BCL_ExampleIndirectCheck
      (
        running_min_max_vector_3d.GetMax(),
        linal::Vector3D( 1.0, 2.0, 3.0),
        "Adding vectors to running min/max"
      );

      // create some vectors
      linal::Vector< double> x2( linal::MakeVector< double>( 2.0, -1.0));
      linal::Vector< double> y1( linal::MakeVector< double>( 0.0, 1.0));
      linal::Vector< double> yn1( linal::MakeVector< double>( 0.0, -1.0));
      BCL_ExampleCheck
      (
        ( running_var_ave_vector_2d += x2).GetMin(),
        linal::MakeVector< double>( 0.0, -1.0)
      );
      BCL_ExampleCheck
      (
        ( running_var_ave_vector_2d += x2).GetMax(),
        linal::MakeVector< double>( 2.0, 0.0)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( running_min_max_double, math::RunningMinMax< double>()),
        true,
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleMathRunningMinMax

  const ExampleClass::EnumType ExampleMathRunningMinMax::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRunningMinMax())
  );
} // namespace bcl

