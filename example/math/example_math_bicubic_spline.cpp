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
#include "math/bcl_math_bicubic_spline.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_bicubic_spline.cpp
  //!
  //! @author mueller
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathBicubicSpline :
    public ExampleInterface
  {
  public:

    ExampleMathBicubicSpline *Clone() const
    { return new ExampleMathBicubicSpline( *this);}

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
      math::BicubicSpline default_spline;

      // Clone
      util::ShPtr< math::BicubicSpline> sp_spline( default_spline.Clone());

    /////////////////
    // data access //
    /////////////////

    //////////////////////////
    // set up training data //
    //////////////////////////

      // init, uses same results vector for every spline type
      const double values[] =
      {
        26,  3,   1,  2, 1, 3, 6, 3, 8, 2, 7, 8, 3, 4, 2, 1, 2,  5,  30,  0,  2,  4, 6, 3, 4, 3, 3, 4, 11, 5, 8, 5, 2, 0, 2,  2,
        2,   0,   0,  1, 2, 0, 2, 1, 1, 3, 3, 0, 1, 0, 0, 1, 0,  0,  4,   0,  0,  0, 1, 3, 0, 5, 0, 0, 1,  1, 1, 0, 0, 3, 0,  0,
        1,   0,   0,  0, 0, 0, 1, 0, 0, 1, 2, 3, 0, 0, 0, 1, 0,  0,  76,  0,  0,  0, 0, 0, 0, 3, 2, 2, 0,  0, 0, 0, 0, 0, 0,  0,
        476, 107, 19, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 9, 3, 10, 14, 447, 14, 13, 4, 9, 5, 1, 3, 0, 0, 0,  0, 0, 0, 1, 5, 14, 104,
        1,   1,   4,  2, 1, 2, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 0,  1,  35,  0,  0,  0, 4, 0, 0, 1, 4, 3, 0,  0, 0, 2, 2, 1, 1,  1
      };

      const linal::Matrix< double> input_values( 5, 36, values);
      util::Format fmt, fmt2;
      fmt.W( 11).FFP( 8).S();
      fmt2.W( 7).FFP( 5).S();

      // natural means, the second order derivative at beginning and end is 0
      // this leads to a rather slow changing first order derivative meaning
      // "nearly linear" behavior of the spline
      // in y-direction the spline is made periodic
      math::BicubicSpline naturalspline;
      math::SplineBorderType behavior[] = { math::e_Natural, math::e_Periodic};
      const double start[ 2] = { 10, -180};
      const double delta[ 2] = { 10, 10};
      const bool lin_cont[ 2] = { true, true};
      const storage::Pair< double, double> first_be[ 2] = { storage::Pair< double, double>( 10, 10), storage::Pair< double, double>( 10, 10)};
      naturalspline.Train( behavior, start, delta, input_values, lin_cont, first_be);

    /////////////////
    // data access //
    /////////////////

      BCL_Example_Check
      (
        naturalspline.GetDsecox().GetNumberCols() == 36 && naturalspline.GetDsecox().GetNumberRows() == 5,
        "Dsecox doesn't have right dimensions."
      );

      BCL_Example_Check
      (
        naturalspline.GetDsecoy().GetNumberCols() == 36 && naturalspline.GetDsecoy().GetNumberRows() == 5,
        "Dsecoy doesn't have right dimensions."
      );

      BCL_Example_Check
      (
        naturalspline.GetDsecoxy().GetNumberCols() == 36 && naturalspline.GetDsecoxy().GetNumberRows() == 5,
        "Dsecoxy doesn't have right dimensions."
      );

    ////////////////
    // operations //
    ////////////////

      BCL_Example_Check
      (
        naturalspline.F( linal::MakeVector( 3.0, 0.0)) == 33.5,
        "F( 3, 0) should return ... but returns " + util::Format()( naturalspline.F( linal::MakeVector( 3, 0)))
      );

      BCL_Example_Check
      (
        naturalspline.dFdx( linal::MakeVector( 3, 0)) == -0.5,
        "dFdx( 3, 0) should return ... but returns " + util::Format()( naturalspline.dFdx( linal::MakeVector( 3, 0)))
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( naturalspline.dFdy( linal::MakeVector( 3, 0)), -0.391435),
        "dFdy( 3, 0) should return ... but returns " + util::Format()( naturalspline.dFdy( linal::MakeVector( 3, 0)))
      );

      BCL_Example_Check
      (
        naturalspline.FdF( linal::MakeVector( 3, 0)).First() == 33.5
        && math::EqualWithinTolerance( naturalspline.FdF( linal::MakeVector( 3, 0)).Second().Norm(), 0.634997),
        "FdF( 3, 0) returns wrong result: " + util::Format()( naturalspline.FdF( linal::MakeVector( 3, 0)))
      );

      BCL_MessageStd
      (
        "Example for a spline with natural behavior in x-direction (F_xx(x_0, y)=F_xx(x_dim-1, y)=0)"
        " and periodic behavior in y-direction (F(x, y + n * dimy * deltay)=F(x, y) :"
      );

      BCL_MessageStd( " x        F(x, -170)       F_x(x, -170)");
      for( double x( 49.99); x <= 50.01; x += 0.001)
//      for( double x( 9.99); x <= 10.01; x += 0.001)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( naturalspline.F( linal::MakeVector( x, -170.0))) + "  "
          + fmt( naturalspline.dFdx( linal::MakeVector( x, -170.0)))
        );
      }

      BCL_MessageStd( " y           F(10, y)      F(10, y+360)  F_y(10, y)");
      for( double y( -190); y <= -168; y += 0.5)
      {
        BCL_MessageStd
        (
          fmt2( y) + fmt( naturalspline.F( linal::MakeVector( 10.0, y))) + "  " +
          fmt( naturalspline.F( linal::MakeVector( 10.0, y + 360.0))) + "  " +
          fmt( naturalspline.dFdy( linal::MakeVector( 10.0, y)))
        );
      }

      //writing spline to file
      WriteBCLObject( naturalspline);
      math::BicubicSpline naturalspline_read;
      ReadBCLObject( naturalspline_read);

      //firstder means, that you give the desired value of the first order derivative in x_0/y_0 and x_dimx-1/y_dimy-1
      //the data will be preprocessed in this case to smooth it
      math::BicubicSpline firstderspline;
      math::SplineBorderType behavior2[ 2] = { math::e_FirstDer, math::e_FirstDer};
      const double start2[ 2]        = { 3.5, -180};
      const double delta2[ 2]        = { 10, 10};
      const bool lin_cont2[ 2]       = { true, true};
      const storage::Pair< double, double> first_be2[ 2] = { storage::Pair< double, double>( 2, -2), storage::Pair< double, double>( 3, 0)};
      firstderspline.TrainWithPreprocessing( behavior2, start2, delta2, input_values, lin_cont2, first_be2, 0.5);

      BCL_MessageStd( "Example for a spline with F_x=2/-2 and F_y=3/0 at start end  :");

      BCL_MessageStd( " x        F(x, -170)       F_x(x, -170)");
      for( double x( 3); x <= 4; x += 0.1)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( firstderspline.F( linal::MakeVector( x, -165.0))) + "  "
          + fmt( firstderspline.dFdx( linal::MakeVector( x, -165.0)))
        );
      }

      BCL_MessageStd( " x        F(x, -170)       F_x(x, -170)");
      for( double x( 43); x <= 44; x += 0.1)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( firstderspline.F( linal::MakeVector( x, -165.0))) + "  "
          + fmt( firstderspline.dFdx( linal::MakeVector( x, -165.0)))
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathBicubicSpline

  const ExampleClass::EnumType ExampleMathBicubicSpline::s_Instance
  (
    GetExamples().AddEnum( ExampleMathBicubicSpline())
  );

} // namespace bcl

