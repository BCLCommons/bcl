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
#include "math/bcl_math_cubic_spline_damped.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_variable_delta.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_cubic_spline_damped.cpp
  //!
  //! @author mendenjl
  //! @date Feb 11, 2016
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathCubicSplineDamped :
    public ExampleInterface
  {
  public:

    ExampleMathCubicSplineDamped *Clone() const
    {
      return new ExampleMathCubicSplineDamped( *this);
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

      // Test the Spline with constant Delta

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
      math::CubicSplineDamped natural_spline;

      // Clone
      util::ShPtr< math::CubicSplineDamped> sp_natural_spline( natural_spline.Clone());

    ////////////////
    // operations //
    ////////////////

      // Now Run the Same tests on same system with varying deltas between abscissa values

      const double abscissa_vd[] =
      {
        -180.00, -172.00, -157.00, -151.00, -145.00, -130.00, -125.00, -114.00, -93.00,
         -90.00,  -81.00,  -76.00,  -60.00,  -53.00,  -47.00,  -30.00,  -27.00,  -15.00,
           0.00,    7.00,   19.00,   25.00,   40.00,   50.00,   59.00,   68.00,   81.00,
          86.00,   95.00,  110.00,  119.00,  121.00,  139.00,  150.00,  160.00,  170.00
      };

      const double ordinate[] =
      {
        26.000000, 6.055637, 1.614554, 2.041015, 1.446470, 3.000000, 5.217905, 3.606333, 7.643058,
         2.000000, 6.354713, 8.598536, 3.000000, 3.582321, 3.807709, 1.000000, 1.620007, 0.276776,
        30.000000, 9.479543, 1.077421, 3.892720, 6.000000, 3.000000, 3.914016, 3.215717, 2.922658,
         2.698211, 8.290267, 5.000000, 7.752718, 8.093626, 2.269422, 0.000000, 2.000000, 2.050000
      };

      const linal::Vector< double> abscissa_values( 36, abscissa_vd);
      const linal::Vector< double> ordinate_values( 36, ordinate);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::CubicSplineDamped natural_spline_ne;
      math::CubicSplineVariableDelta natural_spline_vd;

      natural_spline_ne.Train( abscissa_values, ordinate_values);
      natural_spline_vd.Train( math::e_Natural, abscissa_values, ordinate_values);

      BCL_MessageStd( "Example for a natural spline variable delta (f''(x_0)=f''(x_dim-1)=0):");

      BCL_MessageStd( "Lower / Upper  bound: ");
      BCL_MessageStd( " x       f_NE(x)      f_NE'(x)  f_NormCubicSpline(x)      f_NormCubicSpline'(x)");

      linal::Vector< double> yvals( 351);
      size_t i( 0);
      for( double x( -180.0); x <= 170; x += 1.0, ++i)
      {
        yvals( i) = natural_spline_ne( x);
      }
      math::CubicSplineDamped csd;
      csd.Train( -180.0, 1.0, yvals);

      for( double x( -185.0); x <= 175.0; x += 0.5, ++i)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( natural_spline_ne( x)) + fmt( natural_spline_ne.dF( x))
          + fmt( csd( x)) + fmt( csd.dF( x))
          //+ fmt( natural_spline_vd( x)) + fmt( natural_spline_vd.dF( x))
        );
      }

      BCL_ExampleCheckWithinTolerance( natural_spline_ne( -180.00), 26.000000, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne.dF( -180.00), -3.34732, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne( -175.00), 10.7065, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne.dF( -175.00), -2.36124, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne( -110.00), 3.98991, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne.dF( -110.00), 0.177841, 0.001);

      // test there is no ringing / new maxima introduced at 165.0, which there is in the conventional spline
      BCL_ExampleCheckWithinTolerance( natural_spline_ne( 165.0), 2.0375, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne.dF( 165.0), 0.005, 0.001);

      BCL_Message( util::Message::e_Standard, "Between bounds: ");
      BCL_Message( util::Message::e_Standard, " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_Message( util::Message::e_Standard, fmt2( x) + fmt( natural_spline_ne( x)) + fmt( natural_spline_ne.dF( x)));
      }

      BCL_ExampleCheckWithinTolerance( natural_spline_ne( -80.1), 6.8258, 0.001);
      BCL_ExampleCheckWithinTolerance( natural_spline_ne.dF( -80.1), 0.571427, 0.001);

    ///////////////////////////////////////////////////////
    // Test Derivative function of spline with this data //
    ///////////////////////////////////////////////////////

      const double x_vd[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10};
      const double lnx[] = { 0.6931, 1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026};

      const linal::Vector< double> x_values( 9, x_vd);
      const linal::Vector< double> lnx_values( 9, lnx);

      math::CubicSpline bcl_lnx;
      bcl_lnx.Train( math::e_Natural, 2, 1, lnx_values, storage::Pair< double, double>( 0, 0));

      math::CubicSplineDamped vd_lnx;
      vd_lnx.Train( x_values, lnx_values);
      BCL_ExampleCheckWithinTolerance( vd_lnx( 2.5), 0.910575, 0.001);
      // real derivative of ln(x) is of course 1/x, so 0.4 is the true value.
      BCL_ExampleCheckWithinTolerance( vd_lnx.dF( 2.5), 0.4055, 0.001);
      BCL_ExampleCheckWithinTolerance( vd_lnx( 1.5), 0.9253, 0.001);
      BCL_ExampleCheckWithinTolerance( vd_lnx.dF( 1.5), 0.4644, 0.001);

      BCL_MessageStd( "Ln(x) with BCL Splines and Variable Delta Implementation: ");
      BCL_MessageStd( " x     f(x)_BCL  f'(x)_BCL  f(x)_VD  f'(x)_VD");

      // This is a sample of the trained spline.  If you request a point outside of the spline, an error will be given
      for( double x( 2.0); x <= 10.0; x += 0.1)
      {
        BCL_MessageStd
        (
          fmt2( x) +
          fmt( bcl_lnx( x)) +
          fmt( bcl_lnx.dF( x)) +
          fmt( vd_lnx( x)) +
          fmt( vd_lnx.dF( x))
        );
      }

      // writing spline to file
      WriteBCLObject( natural_spline_ne);
      math::CubicSplineDamped natural_spline_vd_read;
      ReadBCLObject( natural_spline_vd_read);

      // output should be the same like the last output
      BCL_Message( util::Message::e_Debug, "The same natural spline after writing and reading in:");
      BCL_Message( util::Message::e_Debug, util::Format()( natural_spline_vd_read));

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( GetStaticClassName< math::CubicSplineDamped>(), sp_natural_spline->GetClassIdentifier());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathCubicSpline

  const ExampleClass::EnumType ExampleMathCubicSplineDamped::s_Instance
  (
    GetExamples().AddEnum( ExampleMathCubicSplineDamped())
  );
} // namespace bcl
