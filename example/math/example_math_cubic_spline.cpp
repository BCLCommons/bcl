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
#include "math/bcl_math_cubic_spline.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_cubic_spline.cpp
  //!
  //! @author mueller
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathCubicSpline :
    public ExampleInterface
  {
  public:

    ExampleMathCubicSpline *Clone() const
    {
      return new ExampleMathCubicSpline( *this);
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
      //init, uses same results vector for every spline type
      const double values[] =
      {
        26, 3, 1, 2, 1, 3, 6, 3, 8,
         2, 7, 8, 3, 4, 2, 1, 2, 5,
        30, 0, 2, 4, 6, 3, 4, 3, 3,
        4, 11, 5, 8, 5, 2, 0, 2, 2
      };

      const linal::Vector< double> input_values( 36, values);

      //define the output formats
      util::Format fmt, fmt2;
      fmt.W( 9).FFP( 6); fmt.S();
      fmt2.W( 4).FFP( 2); fmt2.S();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::CubicSpline natural_spline;

      // Clone
      util::ShPtr< math::CubicSpline> sp_natural_spline( natural_spline.Clone());

    ////////////////
    // operations //
    ////////////////

      const double lnx[] = { 0.6931, 1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026};

      const linal::Vector< double> lnx_values( 9, lnx);

      math::CubicSpline bcl_lnx;
      bcl_lnx.Train( math::e_Natural, 2, 1, lnx_values, storage::Pair< double, double>( 0, 0));

      double check = bcl_lnx.dF( 3.3);

      BCL_Example_Check
       (
          math::EqualWithinTolerance( check, 0.306105),
         "expected natural spline function value 0.306105 does not match computed function value: "
          + util::Format()( check)
       );

      // natural means, the second order derivative at beginning and end is 0
      // this leads to a rather slow changing first order derivative meaning
      // "nearly linear" behavior of the spline
      natural_spline.Train( math::e_Natural, -180, 10, input_values, storage::Pair< double, double>( 10, 10));

      BCL_MessageStd( "Example for a natural spline (f''(x_0)=f''(x_dim-1)=0):");

      BCL_MessageStd( "Lower bound: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -180.5); x < -179.5; x += 0.1)
      {
        BCL_MessageStd( fmt2( x) + fmt( natural_spline( x)) + fmt( natural_spline.dF( x)));
      }

      BCL_Example_Check
      (
        math::EqualWithinTolerance( natural_spline( -180.30), 26.850871)
        && math::EqualWithinTolerance( natural_spline.dF( -180.30), -2.836237)
        && math::EqualWithinTolerance( natural_spline( -179.70), 25.149274)
        && math::EqualWithinTolerance( natural_spline.dF( -179.70), -2.834789)
        && math::EqualWithinTolerance( natural_spline( 179.70), 1.207596)
        && math::EqualWithinTolerance( natural_spline.dF( 179.70), -0.081691)
        && math::EqualWithinTolerance( natural_spline( 180.30),  1.158582)
        && math::EqualWithinTolerance( natural_spline.dF( 180.30), -0.081691),
        "expected natural spline function values do not match actual function values"
      );

      BCL_MessageStd( "Between bounds: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( natural_spline( x)) + fmt( natural_spline.dF( x)));
      }

      BCL_MessageStd( "Upper bound: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( 179.5); x < 180.5; x += 0.1)
      {
        BCL_MessageStd( fmt2( x) + fmt( natural_spline( x)) + fmt( natural_spline.dF( x)));
      }

      BCL_MessageDbg( "The natural spline itself:");

      // output the spline
      BCL_MessageDbg( util::Format()( natural_spline));

      // writing spline to file
      WriteBCLObject( natural_spline);
      math::CubicSpline natural_spline_read;
      ReadBCLObject( natural_spline_read);

      // output should be the same like the last output
      BCL_MessageDbg( "The same natural spline after writing and reading in:");
      BCL_MessageDbg( util::Format()( natural_spline_read));

      // firstder sets the values of the first derivative on start and end, shows the possibility to
      // "continue linearly" with a non-periodic spline over its range
      math::CubicSpline first_spline;
      first_spline.Train( math::e_FirstDer, -180, 10, input_values, storage::Pair< double, double>( -4, 17));
      BCL_MessageStd( "Example for a spline with set values for f'(x_0) and f'(x_dim-1):");
      BCL_MessageDbg( util::Format()( first_spline));

      BCL_MessageStd( "Lower bound: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -180.1); x < -179.5; x += 0.1)
      {
        BCL_MessageStd( fmt2( x) + fmt( first_spline( x)) + fmt( first_spline.dF( x)));
      }

      BCL_MessageStd( "Between bounds: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( first_spline( x)) + fmt( first_spline.dF( x)));
      }

      BCL_MessageStd( "Upper bound: ");
      for( double x( 179.8); x < 180.2; x += 0.1)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( first_spline.FdF( x).First())
        );
      }

      BCL_Example_Check
      (
        math::EqualWithinTolerance( first_spline( -180.10), 26.4)
        && math::EqualWithinTolerance( first_spline.dF( -180.10), -4.0)
        && math::EqualWithinTolerance( first_spline( -179.70), 24.818056)
        && math::EqualWithinTolerance( first_spline.dF( -179.70), -3.879911)
        && math::EqualWithinTolerance( first_spline( 179.80), 168.6)
        && math::EqualWithinTolerance( first_spline.dF( 179.80), 17)
        && math::EqualWithinTolerance( first_spline( 180.10), 173.7)
        && math::EqualWithinTolerance( first_spline.dF( 180.10), 17),
        "expected first spline function values do not match actual function values"
      );

      // periodic repeats the spline over and over with the same values in a range of equal length
      // automatically converts arguments to the correct range
      math::CubicSpline periodic_spline;
      periodic_spline.TrainWithPreprocessing
      (
        math::e_Periodic, -180, 10, input_values, storage::Pair< double, double>( 1, 1), 0.5
      );

      BCL_MessageStd
      (
        "Example for use of FdF (trained spline): FdF( -400)="
        + util::Format()( storage::Pair< double, double>( periodic_spline.FdF( -400)))
      );

      BCL_MessageStd( "Example for  a periodic spline:");

      BCL_MessageStd( " x      f( x)     f'( x)    f( x+360) f'( x+360) ");
      for( double x( -180.1); x < -179.9; x += 0.02)
      {
        BCL_MessageStd
        (
          fmt2( x)  + fmt( periodic_spline( x))          + fmt( periodic_spline.dF( x))
                    + fmt( periodic_spline( x + 360))    + fmt( periodic_spline.dF( x + 360))
        );
      }

      BCL_MessageStd( "Between bounds: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( periodic_spline( x)) + fmt( periodic_spline.dF( x)));
      }

      BCL_Example_Check
      (
        math::EqualWithinTolerance( periodic_spline( -180.10), 14.448139)
        && math::EqualWithinTolerance( periodic_spline.dF( -180.10), 0.536993)
        && math::EqualWithinTolerance( periodic_spline( -179.96), 14.519708)
        && math::EqualWithinTolerance( periodic_spline.dF( -179.96), 0.485281)
        && math::EqualWithinTolerance( periodic_spline( 179.90), 14.448139)
        && math::EqualWithinTolerance( periodic_spline.dF( 179.90), 0.536993)
        && math::EqualWithinTolerance( periodic_spline( 180.04), 14.519708)
        && math::EqualWithinTolerance( periodic_spline.dF( 180.04), 0.485281),
        "expected periodic spline function values do not match actual function values"
      );

      const std::string gnuplot_filename( AddExampleOutputPathToFilename( natural_spline, "natural_spline.gnuplot"));
      io::OFStream gnuplot_stream;

      BCL_ExampleMustOpenOutputFile( gnuplot_stream, gnuplot_filename);

      math::GnuplotHeatmap natural_spline_heatmap;
      natural_spline_heatmap.SetFromCubicSpline( natural_spline, false, true);
      natural_spline_heatmap.SetTitleAndLabel( "Natural spline", "", "", "");
      natural_spline_heatmap.WriteScript( gnuplot_stream);
      io::File::CloseClearFStream( gnuplot_stream);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( GetStaticClassName< math::CubicSpline>(), sp_natural_spline->GetClassIdentifier());

      BCL_ExampleCheck( natural_spline.GetDsecox().GetSize(), 36);

      BCL_ExampleCheck( natural_spline.GetStart(), -180);

      BCL_ExampleCheck( natural_spline.GetDelta(), 10);

      BCL_ExampleCheck( natural_spline.GetValues().GetSize(), 36);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathCubicSpline

  const ExampleClass::EnumType ExampleMathCubicSpline::s_Instance
  (
    GetExamples().AddEnum( ExampleMathCubicSpline())
  );
} // namespace bcl
