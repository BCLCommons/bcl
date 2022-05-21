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
#include "math/bcl_math_statistics.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_statistics.cpp
  //!
  //! @author butkiem1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathStatistics :
    public ExampleInterface
  {
  public:

    ExampleMathStatistics *Clone() const
    { return new ExampleMathStatistics( *this);}

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
      BCL_MessageStd( "Univariate Statistics");

      linal::Matrix< double> ma( 4, 5);
      ma.AsVector().SetRand( 0.0, 10.0);
      BCL_MessageStd( " This is a random 4x5 matrix ma:\n" + util::Format()( ma));
      BCL_MessageStd
      (
        " This is the minimal value of matrix ma: "
        + util::Format()( math::Statistics::MinimumValue( ma.Begin(), ma.End()))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::MinimumValue( ma.Begin(), ma.End()), double( 0.192711), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is also the minimal value of matrix ma: "
        + util::Format()( *math::Statistics::PointerToMinimumValue( ma.Begin(), ma.End()))
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          *math::Statistics::PointerToMinimumValue( ma.Begin(), ma.End()), double( 0.192711), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is the maximal value of matrix ma: "
        + util::Format()( math::Statistics::MaximumValue( ma.Begin(), ma.End()))
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::MaximumValue( ma.Begin(), ma.End()), double( 9.46668), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is also the maximal value of matrix ma: "
        + util::Format()( *math::Statistics::PointerToMaximumValue( ma.Begin(), ma.End()))
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          *math::Statistics::PointerToMaximumValue( ma.Begin(), ma.End()), double( 9.46668), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is the sum of the elements of matrix ma: "
        + util::Format()( math::Statistics::Sum( ma.Begin(), ma.End(), 0))
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::Sum( ma.Begin(), ma.End(), 0), double( 90.6887), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is the sum of the elements of matrix ma: "
        + util::Format()( ma.Sum())
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          ma.Sum(), double( 90.6887), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        " This is the square norm of the elements of matrix ma: "
        + util::Format()( math::Statistics::SquareNorm( ma.Begin(), ma.End()))
      );

      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::SquareNorm( ma.Begin(), ma.End()), double( 545.491), double( 0.001)
        ),
        true
      );

      linal::Vector< double> vector( 2, 2.0);
      BCL_ExampleCheckWithinAbsTolerance( math::Statistics::SquareNorm( vector.Begin(), vector.End()), 8.0, 0.0001);
      BCL_ExampleCheckWithinAbsTolerance( ma.AsVector().SquareNorm(), 545.491, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( math::Statistics::StandardDeviation( ma.Begin(), ma.End()), 2.59104, 0.0001);

      math::Statistics::Normalize( ma.Begin(), ma.End());
      BCL_ExampleIndirectCheckWithinTolerance( math::Statistics::Norm( ma.Begin(), ma.End()), 1, 0.001, "Normalize");

      linal::Vector< double> v1( 10);
      v1.SetRand( 0.0, 20.0);

      BCL_MessageStd( "This is a random vector: " + util::Format()( v1));

      BCL_MessageStd( "sum: " + util::Format()( math::Statistics::Sum( v1.Begin(), v1.End())));
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::Sum( v1.Begin(), v1.End()), double( 92.2477), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd( "mean: " + util::Format()( math::Statistics::Mean( v1.Begin(), v1.End())));
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::Mean( v1.Begin(), v1.End()), double( 9.22477), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd( "variance: " + util::Format()( math::Statistics::Variance( v1.Begin(), v1.End())));
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::Variance( v1.Begin(), v1.End()), double( 41.4324), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd( "standard deviation: " + util::Format()( math::Statistics::StandardDeviation( v1.Begin(), v1.End())));
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::StandardDeviation( v1.Begin(), v1.End()), double( 6.4368), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd( "\n\nBivariate Statistics\n");
      linal::Vector< double> v2( 10);
      linal::Vector< double> v3( 10);
      linal::Vector< double> v4( 10);

      v2.SetRand( 0.0, 20.0);
      v4.SetRand( -2, 2);
      v3 = v2 + v4;

      BCL_MessageStd( "data set x: " + util::Format()( v2));
      BCL_MessageStd( "data set y: " + util::Format()( v3));

      BCL_MessageStd
      (
        "r (Pearson): "
        + util::Format()( math::Statistics::CorrelationCoefficient( v2.Begin(), v2.End(), v3.Begin(), v3.End()))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::CorrelationCoefficient( v2.Begin(), v2.End(), v3.Begin(), v3.End()), double( 0.968378), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        "r (Spearman): "
        + util::Format()( math::Statistics::CorrelationSpearman( v2.Begin(), v2.End(), v3.Begin(), v3.End()))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::CorrelationSpearman( v2.Begin(), v2.End(), v3.Begin(), v3.End()), double( 0.987879), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        "r_squared: "
        + util::Format()( math::Statistics::RSquared( v2.Begin(), v2.End(), v3.Begin(), v3.End()))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::RSquared( v2.Begin(), v2.End(), v3.Begin(), v3.End()), double( 0.937755), double( 0.0001)
        ),
        true
      );

      BCL_MessageStd
      (
        "covariance: "
        + util::Format()( math::Statistics::Covariance( v2.Begin(), v2.End(), v3.Begin(), v3.End()))
      );
      BCL_ExampleCheck
      (
        math::EqualWithinAbsoluteTolerance
        (
          math::Statistics::Covariance( v2.Begin(), v2.End(), v3.Begin(), v3.End()), double( 23.0125), double( 0.0001)
        ),
        true
      );

//      storage::Pair< double, double> lsl( math::LSL(v2, v3));
//      BCL_MessageStd( "least squares line: " + util::Format()( lsl.First()) +" * x + " + util::Format()( lsl.Second()) + " = y");

      // CumulativeEuclidian
      {
        const double x_arr[] = {1.0, 0.0, 4.0, 9.0};
        const double y_arr[] = {3.0, 5.0, 2.0, 7.0};
        const linal::Vector< double> x( 4, x_arr), y( 4, y_arr);
        const double result( math::Statistics::CumulativeEuclidian( x.Begin(), x.End(), y.Begin(), y.End()));
        BCL_MessageDbg( "CumulativeEuclidian is " + util::Format()( result));
        const double expected_result( 9.327379);
        BCL_MessageDbg( "expected_result is " + util::Format()( expected_result));
        BCL_ExampleCheckWithinTolerance( result, expected_result, 0.00001);
      }

      // Kolmogorov-Smirnov
      {
        linal::Vector< double> a( size_t( 5), double( 1.0));
        double b_arr[] = { 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19};
        linal::Vector< double> b( 10, b_arr);
        double c_arr[] = { 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26 };
        linal::Vector< double> c( 10, c_arr);
        double d_arr[] = { 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35 }; 
        linal::Vector< double> d( 10, d_arr);

        std::pair< double, double> res_aa( math::Statistics::KolmogorovSmirnov( a.Begin(), a.End(), a.Begin(), a.End()));
        BCL_ExampleCheckWithinTolerance( res_aa.first, 0, 0.000001);
        BCL_ExampleCheckWithinTolerance( res_aa.second, 0.86014, 0.000001);
        std::pair< double, double> res_ab( math::Statistics::KolmogorovSmirnov( a.Begin(), a.End(), b.Begin(), b.End()));
        BCL_ExampleCheckWithinTolerance( res_ab.first, 1, 0.000001);
        BCL_ExampleCheckWithinTolerance( res_ab.second, 0.744903, 0.000001);
        std::pair< double, double> res_bc( math::Statistics::KolmogorovSmirnov( b.Begin(), b.End(), c.Begin(), c.End()));
        BCL_ExampleCheckWithinTolerance( res_bc.first, 0.7, 0.000001);
        BCL_ExampleCheckWithinTolerance( res_bc.second, 0.60821, 0.000001);
        std::pair< double, double> res_cd( math::Statistics::KolmogorovSmirnov( c.Begin(), c.End(), d.Begin(), d.End()));
        BCL_ExampleCheckWithinTolerance( res_cd.first, 0.9, 0.000001);
        BCL_ExampleCheckWithinTolerance( res_cd.second, 0.60821, 0.000001);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathStatistics

  const ExampleClass::EnumType ExampleMathStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleMathStatistics())
  );

} // namespace bcl
