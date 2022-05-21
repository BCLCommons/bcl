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
#include "math/bcl_math_tricubic_spline.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_tricubic_spline.cpp
  //!
  //! @author mueller
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathTricubicSpline :
    public ExampleInterface
  {
  public:

    ExampleMathTricubicSpline *Clone() const
    { return new ExampleMathTricubicSpline( *this);}

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
      //the layers describe the x-axis, the rows of each layer are y,
      //the columns are z

      const double values[] =
      {
        //layer 0 (x0)
        //z->
        26,3, 1, 2, 1, 3, //y
        6, 3, 8, 2, 7, 8,
        3, 4, 2, 1, 2, 5,
        30,0, 2, 4, 6, 3,
        4, 3, 3, 4, 11,5,
        8, 5, 2, 0, 2, 2,
        //layer 1 (x1)
        //z->
        2, 0, 0, 1, 2, 0, //y
        2, 1, 1, 3, 3, 0,
        1, 0, 0, 1, 0, 0,
        4, 0, 0, 0, 1, 3,
        0, 5, 0, 0, 1, 1,
        1, 0, 0, 3, 0, 0,
        //layer 2 (x2)
        //z->
        1, 0, 0, 0, 0, 0, //y
        1, 0, 0, 1, 2, 3,
        0, 0, 0, 1, 0, 0,
        76,0, 0, 0, 0, 0,
        0, 3, 2, 2, 0, 0,
        0, 0, 0, 0, 0, 0,
        //layer 3 (x3)
        //z->
        476, 107, 19, 2, 2, 0, //y
        0,     0,  0, 0, 0, 0,
        0,    1, 9, 3, 10, 14,
        447, 14, 13, 4, 9,  5,
          1,  3,  0, 0, 0,  0,
          0, 0, 1, 5, 14, 104,
        //layer 4 (x4)
        //z->
        1, 1, 4, 2, 1, 2, //y
        2, 1, 1, 0, 1, 0,
        2, 1, 0, 0, 0, 1,
        35,0, 0, 0, 4, 0,
        0, 1, 4, 3, 0, 0,
        0, 2, 2, 1, 1, 1
      };

      //the tensor describes layers, rows, columns, source of input
      const math::Tensor< double> input_tensor( 5, 6, 6, values);

      math::SplineBorderType border[ 3] = { math::e_Natural, math::e_Natural, math::e_Periodic};

      //these vectors are used to input the starting points and the
      //grid width delta of every dimension (x, y, z) into the spline
      const double start[] = {10, 10, 10};
      const double delta[] = {10, 10, 10};

      const bool lin_cont[ 3] = { true, true, true};

      //this vector controls the behavior of the spline at the beginning and
      //end of every dimension, only has impact for SplineBorderType FIRSTDER

      //every pair describes the value of the first order derivative at
      //start and end
      const storage::Pair< double, double> first_be[ 3] =
      {
        storage::Pair< double, double>( 10, 10),
        storage::Pair< double, double>( 10, 10),
        storage::Pair< double, double>( 10, 10)
      };

      //setting up the formatting of the output using util::Format
      util::Format fmt, fmt2;
      fmt.W( 11).FFP( 8); fmt.S(); fmt2.W( 7).FFP( 5); fmt2.S();

      // natural means, the second order derivative at beginning and end is 0
      // this leads to a rather slow changing first order derivative meaning
      // "nearly linear" behavior of the spline
      // in z-direction the spline is made periodic
      math::TricubicSpline naturalspline;
      naturalspline.Train( border, start, delta, input_tensor, lin_cont, first_be);

      BCL_MessageStd
      (
        "Example for a spline with natural behavior in x- and y-direction (F_xx(x_0, y, z)=F_xx(x_dim-1, y, z)="
        "F_yy(x, y_0, z)=F_yy(x, y_dim-1, z)=0) "
        "and periodic behavior in z-direction (F(x, y, z + n * dimz * deltaz)=F(x, y, z)) :"
      );
      BCL_MessageStd( "To show continuous behavior at the end of a cell");
      BCL_MessageStd( " x        F(x, 10, 10)       F_x(x, 10, 10)");
      for( double x( 19.9); x < 20.1; x += 0.01)
      {
        BCL_MessageStd( fmt2( x) + fmt( naturalspline.F( x, 10, 10)) + fmt( naturalspline.dFdx( x, 10, 10)));
      }

      BCL_MessageStd( "Behavior at the end of the defined region");
      for( double x( 49); x < 50.5; x += 0.1)
      {
        BCL_MessageStd
        (
          fmt2( x) + fmt( naturalspline.F( x, 10, 10)) + fmt( naturalspline.dFdx( x, 10, 10))
        );
      }
      BCL_MessageStd( "To show periodic behavior in z-direction");
      BCL_MessageStd( " z        F(10, 10, z)       F(10, 10, z+60)");
      for( double z( 20); z < 60; z += 4)
      {
        BCL_MessageStd
        (
          fmt2( z) + fmt( naturalspline.F( 10, 10, z)) + fmt( naturalspline.F( 10, 10, z + 60))
        );
      }

      //this example describes a function f(x, y, z), where (x, y, z) resembles a helix
      //the spline is periodic in every direction, but this is not shown here
      math::TricubicSpline ts;

      BCL_MessageStd( "Training new purely periodic spline");

      border[ coord::GetAxes().e_X] = math::e_Periodic;
      border[ coord::GetAxes().e_Y] = math::e_Periodic;

      ts.Train( border, start, delta, input_tensor, lin_cont, first_be);

      //writing spline to standard output
      BCL_MessageStd( "The periodic spline itself:");
      BCL_MessageStd( util::Format()( ts));

      //writing spline to file
      WriteBCLObject( ts);
      //reading spline from file
      ReadBCLObject( ts);

      //writing spline to standard output, should be the same like above
      BCL_MessageStd( "The periodic spline after writing to file and reading back in:");
      BCL_MessageStd( util::Format()( ts));

      //last function example
      BCL_MessageStd( "Function values along a helix");
      BCL_MessageStd
      (
        "x FdF( 10+x, 10+10*std::cos( math::g_Pi / 18 * x), 10+10*std::sin( math::g_Pi / 18 * x))"
      );
      for( double x( 0); x <= 36; ++x)
      {
        BCL_MessageStd
        (
          util::Format()( x) + " "
          + util::Format()
          (
            ts.FdF( 10 + x, 10 + 10 * std::cos( math::g_Pi / 18 * x), 10 + 10 * std::sin( math::g_Pi / 18 * x))
          )
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathTricubicSpline

  const ExampleClass::EnumType ExampleMathTricubicSpline::s_Instance
  (
    GetExamples().AddEnum( ExampleMathTricubicSpline())
  );

} // namespace bcl
