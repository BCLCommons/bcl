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
#include "math/bcl_math_cubic_spline_variable_delta.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_cubic_spline_variable_delta.cpp
  //!
  //! @author putnamdk
  //! @date June 12, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathCubicSplineVariableDelta :
    public ExampleInterface
  {
  public:

    ExampleMathCubicSplineVariableDelta *Clone() const
    {
      return new ExampleMathCubicSplineVariableDelta( *this);
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
      math::CubicSplineVariableDelta natural_spline;
      math::CubicSpline natural_spline_constant;

      // Clone
      util::ShPtr< math::CubicSplineVariableDelta> sp_natural_spline( natural_spline.Clone());

    ////////////////
    // operations //
    ////////////////

      natural_spline.Train( math::e_Natural, -180, 10, input_values, storage::Pair< double, double>( 10, 10));

      BCL_Message( util::Message::e_Standard, "Example for a natural spline (f''(x_0)=f''(x_dim-1)=0):");

      BCL_Message( util::Message::e_Standard, "Lower / Upper  bound: ");
      BCL_Message( util::Message::e_Standard, " x       f(x)      f'(x)");

      // This is a sample of the trained spline.  If you request a point outside of the spline, an error will be given
      for( double x( -180.0); x <= -179.5; x += 0.1)
      {
        BCL_Message( util::Message::e_Standard, fmt2( x) + fmt( natural_spline( x)) + fmt( natural_spline.dF( x)));
      }

      double value( natural_spline( -180.00));

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 26.000000),
        "expected natural spline function values 26.000000 does not match computed function values: "
         + util::Format()( value)
      );

      value = natural_spline.dF( -180.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -2.836237),
        "expected natural spline derivative value -2.836237 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline( -179.70);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 25.149274),
        "expected natural spline function value 25.149274 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline.dF( -179.70);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -2.834789),
        "expected natural spline derivative value -2.834789 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline( 169.70);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 2.024485),
        "expected natural spline function value 2.024485 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline.dF( 169.70);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -0.081471),
        "expected natural spline derivative value -0.081471 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline( 170.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 2.000000),
        "expected natural spline function value 2.000000 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline.dF( 170.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -0.081691),
        "expected natural spline derivative value -0.081691 does not match computed derivative value: "
         + util::Format()( value)
      );

      BCL_MessageStd( "Between bounds: ");
      BCL_MessageStd( " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( natural_spline( x)) + fmt( natural_spline.dF( x)));
      }

      value = natural_spline( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 6.939413),
        "expected natural spline function value 6.939413 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline.dF( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.610792),
        "expected natural spline derivative value 0.610792 does not match computed derivative value: "
         + util::Format()( value)
      );

      BCL_MessageDbg( "The natural spline itself:");

      // output the spline
      BCL_MessageDbg( util::Format()( natural_spline));

      // writing spline to file
      WriteBCLObject( natural_spline);
      math::CubicSplineVariableDelta natural_spline_read;
      ReadBCLObject( natural_spline_read);

      // output should be the same like the last output
      BCL_MessageDbg( "The same natural spline after writing and reading in:");
      BCL_MessageDbg( util::Format()( natural_spline_read));

      // Now Run the Same tests on same system with varying deltas between abscissa values

      const double abscissa_vd[] =
      {
        -180.00, -172.00, -157.00, -151.00, -145.00, -130.00, -125.00, -114.00, -102.00,
         -90.00,  -81.00,  -76.00,  -60.00,  -53.00,  -47.00,  -30.00,  -27.00,  -15.00,
           0.00,    7.00,   19.00,   25.00,   40.00,   50.00,   59.00,   68.00,   81.00,
          86.00,   95.00,  110.00,  119.00,  121.00,  139.00,  150.00,  160.00,  170.00
      };

      const double ordinate[] =
      {
        26.000000, 6.055637, 1.614554, 2.041015, 1.446470, 3.000000, 5.217905, 3.606333, 7.643058,
         2.000000, 6.354713, 8.598536, 3.000000, 3.582321, 3.807709, 1.000000, 1.620007, 0.276776,
        30.000000, 9.479543, 1.077421, 3.892720, 6.000000, 3.000000, 3.914016, 3.215717, 2.922658,
         2.698211, 8.290267, 5.000000, 7.752718, 8.093626, 2.269422, 0.000000, 2.000000, 2.000000
      };

      const linal::Vector< double> abscissa_values( 36, abscissa_vd);
      const linal::Vector< double> ordinate_values( 36, ordinate);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::CubicSplineVariableDelta natural_spline_vd;

      // Clone
      util::ShPtr< math::CubicSplineVariableDelta> sp_natural_spline_vd( natural_spline.Clone());

      natural_spline_vd.Train( math::e_Natural, abscissa_values, ordinate_values);

      BCL_MessageStd( "Example for a natural spline variable delta (f''(x_0)=f''(x_dim-1)=0):");

      BCL_MessageStd( "Lower / Upper  bound: ");
      BCL_MessageStd( " x       f(x)      f'(x)");

      for( double x( -180.0); x <= 170.0; x += 20.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( natural_spline_vd( x)) + fmt( natural_spline_vd.dF( x)));
      }

      value = natural_spline_vd( -180.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 26.000000),
        "expected natural spline function values 26.000000 does not match computed function values: "
         + util::Format()( value)
      );

      value = natural_spline_vd.dF( -180.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -2.897439),
        "expected natural spline derivative value -2.897439 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline_vd( -175.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 12.302637),
        "expected natural spline function value 12.302637 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline_vd.dF( -175.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -2.423540),
        "expected natural spline derivative value -2.423540 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline_vd( -110.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 4.928556),
        "expected natural spline function value 4.928556 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline_vd.dF( -110.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.486123),
        "expected natural spline derivative value 0.486123 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = natural_spline_vd( 165.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 2.284686),
        "expected natural spline function value 2.284686 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline_vd.dF( 165.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -0.018979),
        "expected natural spline derivative value -0.018979 does not match computed derivative value: "
         + util::Format()( value)
      );

      BCL_Message( util::Message::e_Standard, "Between bounds: ");
      BCL_Message( util::Message::e_Standard, " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_Message( util::Message::e_Standard, fmt2( x) + fmt( natural_spline_vd( x)) + fmt( natural_spline_vd.dF( x)));
      }

      value = natural_spline_vd( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 6.91889),
        "expected natural spline function value 6.91889 does not match computed function value: "
         + util::Format()( value)
      );

      value = natural_spline_vd.dF( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.60097),
        "expected natural spline derivative value 0.60097 does not match computed derivative value: "
         + util::Format()( value)
      );

      // Not a Knot Boundary Condition

      // default constructor
      math::CubicSplineVariableDelta not_a_knot;

      not_a_knot.Train( math::e_NotAKnot, abscissa_values, ordinate_values);

      BCL_MessageStd( "Example for a NotAKnot spline variable delta (f''(x_0)=f''(x_dim-1)=0):");

      BCL_MessageStd( "Lower / Upper  bound: ");
      BCL_MessageStd( " x       f(x)      f'(x)");

      for( double x( -180.0); x <= 170.0; x += 20.0)
      {
        BCL_MessageStd( fmt2( x) + fmt( not_a_knot( x)) + fmt( not_a_knot.dF( x)));
      }

      value = not_a_knot( -180.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 26.000000),
        "expected natural spline function values 26.000000 does not match computed function values: "
         + util::Format()( value)
      );

      value = not_a_knot.dF( -180.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -3.77879),
        "expected natural spline derivative value -3.77879 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = not_a_knot( -175.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 11.3366),
        "expected natural spline function value 11.3366 does not match computed function value: "
         + util::Format()( value)
      );

      value = not_a_knot.dF( -175.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, -2.15743),
        "expected natural spline derivative value -2.15743 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = not_a_knot( -110.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 4.928556),
        "expected natural spline function value 4.928556 does not match computed function value: "
         + util::Format()( value)
      );

      value = not_a_knot.dF( -110.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.486123),
        "expected natural spline derivative value 0.486123 does not match computed derivative value: "
         + util::Format()( value)
      );

      value = not_a_knot( 165.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 2.72339),
        "expected natural spline function value 2.72339 does not match computed function value: "
         + util::Format()( value)
      );

      value = not_a_knot.dF( 165.00);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.0315596),
        "expected natural spline derivative value 0.0315596 does not match computed derivative value: "
         + util::Format()( value)
      );

      BCL_Message( util::Message::e_Standard, "Between bounds: ");
      BCL_Message( util::Message::e_Standard, " x       f(x)      f'(x)");
      for( double x( -80.1); x < -71.0; x += 3.0)
      {
        BCL_Message( util::Message::e_Standard, fmt2( x) + fmt( not_a_knot( x)) + fmt( not_a_knot.dF( x)));
      }

      value = not_a_knot( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 6.91889),
        "expected natural spline function value 6.91889 does not match computed function value: "
         + util::Format()( value)
      );

      value = not_a_knot.dF( -80.10);

      BCL_Example_Check
      (
         math::EqualWithinTolerance( value, 0.60097),
        "expected natural spline derivative value 0.60097 does not match computed derivative value: "
         + util::Format()( value)
      );

    ///////////////////////////////////////////////////////
    // Test Derivative function of spline with this data //
    ///////////////////////////////////////////////////////

      const double x_vd[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10};
      const double lnx[] = { 0.6931, 1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026};

      const linal::Vector< double> x_values( 9, x_vd);
      const linal::Vector< double> lnx_values( 9, lnx);

      math::CubicSpline bcl_lnx;
      bcl_lnx.Train( math::e_Natural, 2, 1, lnx_values, storage::Pair< double, double>( 0, 0));

      math::CubicSplineVariableDelta vd_lnx;
      vd_lnx.Train( math::e_Natural, x_values, lnx_values, storage::Pair< double, double>( 0, 0));

      value = vd_lnx.dF( 3.3);

      BCL_Example_Check
       (
          math::EqualWithinTolerance( value, 0.306105),
         "expected natural spline function value 0.306105 does not match computed function value: "
          + util::Format()( value)
       );

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
      WriteBCLObject( natural_spline_vd);
      math::CubicSplineVariableDelta natural_spline_vd_read;
      ReadBCLObject( natural_spline_vd_read);

      // output should be the same like the last output
      BCL_Message( util::Message::e_Debug, "The same natural spline after writing and reading in:");
      BCL_Message( util::Message::e_Debug, util::Format()( natural_spline_vd_read));

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( GetStaticClassName< math::CubicSplineVariableDelta>(), sp_natural_spline->GetClassIdentifier());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathCubicSpline

  const ExampleClass::EnumType ExampleMathCubicSplineVariableDelta::s_Instance
  (
    GetExamples().AddEnum( ExampleMathCubicSplineVariableDelta())
  );
} // namespace bcl
