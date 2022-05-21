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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_bicubic_spline.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_smooth_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BicubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new BicubicSpline())
    );

  //////////////
  // operator //
  //////////////

    //! @brief operator
    //! @param X argument x to spline
    //! @param Y argument y to spline
    //! @return the interpolate value at (x,y)
    double BicubicSpline::operator()( const double &X, const double &Y) const
    {
      return F( linal::MakeVector< double>( X, Y));
    }

  ////////////////
  // operations //
  ////////////////

    //! train BicubicSpline
    void BicubicSpline::Train( const SplineBorderType BORDER[ 2], const double START[ 2], const double DELTA[ 2], const linal::Matrix< double> &RESULTS, const bool LINCONT[ 2], const storage::Pair< double, double> FIRSTBE[ 2])
    {
      //check, if the points are given in positive direction
      BCL_Assert( DELTA[ coord::GetAxes().e_X] > 0, "deltax <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Y] > 0, "deltay <= 0 not supported");

      //determine values for all dimensions
      const size_t dimx( RESULTS.GetNumberRows());
      const size_t dimy( RESULTS.GetNumberCols());

      //assigning values
      m_Border[ coord::GetAxes().e_X] = BORDER[ coord::GetAxes().e_X];
      m_Border[ coord::GetAxes().e_Y] = BORDER[ coord::GetAxes().e_Y];
      m_Start[ coord::GetAxes().e_X]  = START[ coord::GetAxes().e_X];
      m_Start[ coord::GetAxes().e_Y]  = START[ coord::GetAxes().e_Y];
      m_Delta[ coord::GetAxes().e_X]  = DELTA[ coord::GetAxes().e_X];
      m_Delta[ coord::GetAxes().e_Y]  = DELTA[ coord::GetAxes().e_Y];
      m_DeltaNorm[ coord::GetAxes().e_X]  = Sqr( DELTA[ coord::GetAxes().e_X]) / 6.0; // cache for speed
      m_DeltaNorm[ coord::GetAxes().e_Y]  = Sqr( DELTA[ coord::GetAxes().e_Y]) / 6.0; // cache for speed

      BCL_Assert( m_Border[ 0] != e_NotAKnot && m_Border[ 1] != e_NotAKnot, "The specified boundary condition is not yet supported");

      m_Values  = RESULTS;
      m_Dsecox  = RESULTS;
      m_Dsecoy  = RESULTS;
      m_Dsecoxy = RESULTS;

      m_LinCont[ coord::GetAxes().e_X] = LINCONT[ coord::GetAxes().e_X];
      m_LinCont[ coord::GetAxes().e_Y] = LINCONT[ coord::GetAxes().e_Y];
      m_FirstBe[ coord::GetAxes().e_X] = FIRSTBE[ coord::GetAxes().e_X];
      m_FirstBe[ coord::GetAxes().e_Y] = FIRSTBE[ coord::GetAxes().e_Y];

      //train three times for fxx, fyy, fxxyy
      //reduction to Spline by training only rows/columns at the same time
      for( size_t row( 0); row < dimx; ++row)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_Y], START[ coord::GetAxes().e_Y], DELTA[ coord::GetAxes().e_Y], RESULTS.GetRow( row), FIRSTBE[ coord::GetAxes().e_Y]);
        m_Dsecoy.ReplaceRow( row, cs.GetDsecox());
      }

      for( size_t col( 0); col < dimy; ++col)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], RESULTS.GetCol( col), FIRSTBE[ coord::GetAxes().e_X]);
        m_Dsecox.ReplaceCol( col, cs.GetDsecox());
      }

      for( size_t row( 0); row < dimx; ++row)
      {
        CubicSpline cs;
        cs.Train( BORDER[ coord::GetAxes().e_Y], START[ coord::GetAxes().e_Y], DELTA[ coord::GetAxes().e_Y], m_Dsecox.GetRow( row), FIRSTBE[ coord::GetAxes().e_Y]);
        m_Dsecoxy.ReplaceRow( row, cs.GetDsecox());
      }
      return;
    }

    void BicubicSpline::TrainWithPreprocessing( const SplineBorderType BORDER[ 2], const double START[ 2], const double DELTA[ 2], const linal::Matrix< double> &RESULTS, const bool LINCONT[ 2], const storage::Pair< double, double> FIRSTBE[ 2], const double PREPROC)
    {
      // passing of smoothed data set
      Train( BORDER, START, DELTA, SmoothData::SmoothMatrix( RESULTS, PREPROC, true), LINCONT, FIRSTBE);
    }

    //! return value at certain (x, y)
    double BicubicSpline::F( const linal::Vector< double> &ARGUMENTS) const
    {
      // check that there are two argument values given
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS( coord::GetAxes().e_X));
      const double y( ARGUMENTS( coord::GetAxes().e_Y));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      // check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
      {
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (X < start : "
                             + util::Format()( x)
                             + " < "
                             + util::Format()( m_Start[ coord::GetAxes().e_X])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X]) * m_FirstBe[ coord::GetAxes().e_X].First();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (X < start : "
                            + util::Format()( x)
                            + " < "
                            + util::Format()( m_Start[ coord::GetAxes().e_X])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X]) * dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
      {
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (X > end : "
                             + util::Format()( x)
                             + " > "
                             + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * m_FirstBe[ coord::GetAxes().e_X].Second();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (X > end : "
                            + util::Format()( x)
                            + " > "
                            + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y)) + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
      {
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (y < start : "
                             + util::Format()( y)
                             + " < "
                             + util::Format()( m_Start[ coord::GetAxes().e_Y])
                             + "), using linear continuation"
                           );
                           return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y]) * m_FirstBe[ coord::GetAxes().e_Y].First();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (y < start : "
                            + util::Format()( y)
                            + " < "
                            + util::Format()( m_Start[ coord::GetAxes().e_Y])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y]) * dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
      {
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                           BCL_MessageDbg
                           (
                             "argument out of range (y > end : "
                             + util::Format()( y)
                             + " > "
                             + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                             + "), using linear continuation"
                           );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * m_FirstBe[ coord::GetAxes().e_Y].Second();

          case e_Natural: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
                          BCL_MessageDbg
                          (
                            "argument out of range (y > end : "
                            + util::Format()( y)
                            + " > "
                            + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                            + "), using linear continuation"
                          );
                          return F( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])) + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]));

          case e_Periodic: break;
          case e_NotAKnot: break;
        }
      }

      // determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X]) + 1));

      // determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y]) + 1));

      // see mkbicubw.f in /home/mj/pspline/pspline
      // slightly corrected formula, since the last two lines of the expression weren't symmetric

//C  on input:  f(1,i,j) = f(x(i),y(j))
//C  on output:  f(1,i,j) unchanged
//C              f(2,i,j) = d2f/dx2(x(i),y(j))
//C              f(3,i,j) = d2f/dy2(x(i),y(j))
//C              f(4,i,j) = d4f/dx2dy2(x(i),y(j))
//C
//C  and the interpolation formula for (x,y) in (x(i),x(i+1))^(y(j),y(j+1))
//C  is:
//C        hx = x(i+1)-x(i)   hy = y(j+1)-y(j)
//C        dxp= (x-x(i))/hx   dxm= 1-dxp     dxp,dxm in (0,1)
//C        dyp= (x-x(i))/hx   dym= 1-dyp     dyp,dym in (0,1)
//C        dx3p = dxp**3-dxp  dx3m = dxm**3-dxm     dxp3,dxm3 in (0,1)
//C
//C   finterp = dxm*(dym*f(1,i,j)+dyp*f(1,i,j+1))
//C            +dxp*(dym*f(1,i+1,j)+dyp*f(1,i+1,j+1))
//C     +1/6*hx**2*
//C            dx3m*(dym*f(2,i,j)+dyp*f(2,i,j+1))
//C           +dx3p*(dym*f(2,i+1,j)+dyp*f(2,i+1,j+1))
//C     +1/6*hy**2*
//C            dxm*(dy3m*f(3,i,j)+dy3p*f(3,i,j+1))
//C           +dxp*(dy3m*f(3,i+1,j)+dy3p*f(3,i+1,j+1))
//C     +1/36*hx**2*hy**2*
//C            dx3m*(dy3m*f(4,i,j)+dy3p*f(4,i,j+1))
//C           +dx3p*(dy3m*f(4,i+1,j)+dy3p*f(4,i+1,j+1))

      const double dxp( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X] - floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X]));
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp*dxp*dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]); // =0 at the grid points, adds cubic part of the spline
      const double dx3m( ( dxm*dxm*dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]); // =0 at the grid points, adds cubic part of the spline

      const double dyp( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y] - floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y]));
      const double dym( ( 1 - dyp));
      const double dy3p( ( dyp*dyp*dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]); // =0 at the grid points, adds cubic part of the spline
      const double dy3m( ( dym*dym*dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]); // =0 at the grid points, adds cubic part of the spline

      // generate positive values to prevent some problems with the indices
      while( i < 1) i += dimx;
      while( j < 1) j += dimy;

      return
           dxm * (  dym * m_Values( (  i - 1) % dimx, ( j - 1) % dimy) +  dyp * m_Values( (  i - 1) % dimx, j % dimy))
        +  dxp * (  dym * m_Values(    i      % dimx, ( j - 1) % dimy) +  dyp * m_Values(    i      % dimx, j % dimy))
        + dx3m * (  dym * m_Dsecox( (  i - 1) % dimx, ( j - 1) % dimy) +  dyp * m_Dsecox( (  i - 1) % dimx, j % dimy))
        + dx3p * (  dym * m_Dsecox(    i      % dimx, ( j - 1) % dimy) +  dyp * m_Dsecox(    i      % dimx, j % dimy))
        +  dxm * ( dy3m * m_Dsecoy( (  i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy( (  i - 1) % dimx, j % dimy))
        +  dxp * ( dy3m * m_Dsecoy(    i      % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy(    i      % dimx, j % dimy))
        + dx3m * ( dy3m * m_Dsecoxy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy( ( i - 1) % dimx, j % dimy))
        + dx3p * ( dy3m * m_Dsecoxy(   i      % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy(   i      % dimx, j % dimy))
      ;
    }

    //! return partial derivative at certain (x, y) for x
    double BicubicSpline::dFdx( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS(0));
      const double y( ARGUMENTS(1));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_X].First();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_X].Second();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      // determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) { i++;}

      // determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) { j++;}

      //see F(x, y) for a short explanation of the values
      const double delta_aktx( x-m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y-m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]);
      const double dy3m( ( dym * dym * dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]);

      //generate positive values to prevent some problems with the indices
      while( i < 1){ i += dimx;}
      while( j < 1){ j += dimy;}
      BCL_MessageCrt( "Final i orig: " + util::Format()( i));

      return
        -( dym * m_Values( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * m_Values( ( i - 1) % dimx, j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        +( dym * m_Values( i % dimx      , ( j - 1) % dimy) + dyp * m_Values( i % dimx        , j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        - ( 3 * dxm * dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6 * ( dym * m_Dsecox( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * m_Dsecox( ( i - 1) % dimx, j % dimy))
        + ( 3 * dxp * dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6 *( dym*m_Dsecox( i % dimx      , ( j - 1) % dimy) + dyp * m_Dsecox( i % dimx      , j%dimy))
        - ( dy3m * m_Dsecoy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoy( ( i - 1) % dimx, j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        +( dy3m * m_Dsecoy( i % dimx       , ( j - 1) % dimy) + dy3p * m_Dsecoy( i % dimx     , j % dimy)) / m_Delta[ coord::GetAxes().e_X]
        - ( 3 * dxm * dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6 * ( dy3m * m_Dsecoxy( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * m_Dsecoxy( ( i - 1) % dimx, j % dimy))
        + ( 3 * dxp * dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6 *( dy3m*m_Dsecoxy( i % dimx    , ( j - 1) % dimy) + dy3p * m_Dsecoxy( i % dimx      , j % dimy))
      ;
    }

    //! return partial derivative at certain (x, y) for y
    double BicubicSpline::dFdy( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      const double x( ARGUMENTS( 0));
      const double y( ARGUMENTS( 1));

      const int dimx( m_Values.GetNumberRows());
      const int dimy( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if( x < m_Start[ coord::GetAxes().e_X])
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)
        switch( m_Border[ coord::GetAxes().e_X])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //check if argument y is in range for non-periodic splines
      if( y < m_Start[ coord::GetAxes().e_Y])
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_Y].First();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      if( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)
        switch( m_Border[ coord::GetAxes().e_Y])
        {
          case e_FirstDer: BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return m_FirstBe[ coord::GetAxes().e_Y].Second();
          case e_Natural : BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");BCL_MessageVrb( "argument out of range, using linear continuation");return dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          case e_Periodic: break;
          case e_NotAKnot: break;
        }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) { i++;}
      if( !i){
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] > x) { i--;}
        i++;
      }

      //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) { j++;}
      if( !j){
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] > y) { j--;}
        j++;
      }

      //see F(x, y) for a short explanation of the values
      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]);
      const double dx3m( ( dxm * dxm * dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);

      //generate positive values to prevent some problems with the indices
      while( i < 1){ i += dimx;}
      while( j < 1){ j += dimy;}

      return
          dxm *( -m_Values( ( i-1)%dimx  , (j-1)%dimy)+m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        + dxp *( -m_Values( i%dimx      , (j-1)%dimy)+m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        +dx3m *( -m_Dsecox( ( i-1)%dimx  , (j-1)%dimy)+m_Dsecox( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        +dx3p *( -m_Dsecox( i%dimx      , (j-1)%dimy)+m_Dsecox(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
        + dxm *( -( 3 * dym * dym - 1) * m_Dsecoy( ( i-1)%dimx , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoy( ( i-1)%dimx , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        + dxp *( -( 3 * dym * dym - 1) * m_Dsecoy( i%dimx     , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoy( i%dimx     , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        +dx3m *( -( 3 * dym * dym - 1) * m_Dsecoxy( ( i-1)%dimx, ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoxy( ( i-1)%dimx, j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
        +dx3p *( -( 3 * dym * dym - 1) * m_Dsecoxy( i%dimx    , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * m_Dsecoxy( i%dimx    , j % dimy)) * m_Delta[ coord::GetAxes().e_Y]/ 6
      ;
    }

    //!return value and derivative at certain (x, y)
    storage::Pair< double, linal::Vector< double> > BicubicSpline::FdF( const linal::Vector< double> &ARGUMENTS) const
    {
      BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

      double x = ARGUMENTS(0);
      double y = ARGUMENTS(1);

      int dimx = m_Values.GetNumberRows();
      int dimy = m_Values.GetNumberCols();

      //auxiliary variables for the function value and the derivatives
      double fvalue( 0), dfdxvalue( 0), dfdyvalue( 0);

      //check if argument is in range for non-periodic splines
      if( ( ( m_Border[ coord::GetAxes().e_X] != e_Periodic) && ( x < m_Start[ coord::GetAxes().e_X] || m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
        || ( ( m_Border[ coord::GetAxes().e_Y] != e_Periodic) && ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y)))
      {
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (X < start : "
            + util::Format()( x)
            + " < "
            + util::Format()( m_Start[ coord::GetAxes().e_X])
            + "), using linear continuation"
          );
          fvalue    = F( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y))+(x-m_Start[ coord::GetAxes().e_X])*dFdx( linal::MakeVector(m_Start[ coord::GetAxes().e_X], y));
          dfdxvalue = dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
          dfdyvalue = dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X], y));
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (X > end : "
            + util::Format()( x)
            + " > "
            + util::Format()( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y))+(x-m_Start[ coord::GetAxes().e_X] - ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X])*dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X], y));
          dfdxvalue = dFdx( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y));
          dfdyvalue = dFdy( linal::MakeVector( m_Start[ coord::GetAxes().e_X] + ( dimx-1 ) * m_Delta[ coord::GetAxes().e_X] , y));
        }
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (Y < start : "
            + util::Format()( y)
            + " < "
            + util::Format()( m_Start[ coord::GetAxes().e_Y])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]))+(y-m_Start[ coord::GetAxes().e_Y])*dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          dfdxvalue = dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
          dfdyvalue = dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y]));
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          BCL_MessageDbg
          (
            "argument out of range (Y > end : "
            + util::Format()( y)
            + " > "
            + util::Format()( m_Start[ coord::GetAxes().e_Y] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_Y])
            + "), using linear continuation"
          );
          fvalue    = F(    linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]))+(y-m_Start[ coord::GetAxes().e_Y] - ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y])*dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          dfdxvalue = dFdx( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
          dfdyvalue = dFdy( linal::MakeVector( x, m_Start[ coord::GetAxes().e_Y] + ( dimy-1 ) * m_Delta[ coord::GetAxes().e_Y]));
        }
      }else

      {
        //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
        int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x)i++;
        if( !i) {
          while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] > x)i--;
          i++;
        }

        //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
        int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y)j++;
        if( !j) {
          while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] > y)j--;
          j++;
        }

        //see method F(x,y) for detailed formula

        double delta_aktx = x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X];
        double delta_akty = y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y];

        double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
        double dxm( 1 - dxp);
        double dx3p( ( dxp*dxp*dxp - dxp) * m_DeltaNorm[ coord::GetAxes().e_X]);
        double dx3m( ( dxm*dxm*dxm - dxm) * m_DeltaNorm[ coord::GetAxes().e_X]);

        double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
        double dym( 1 - dyp);
        double dy3p( ( dyp*dyp*dyp - dyp) * m_DeltaNorm[ coord::GetAxes().e_Y]);
        double dy3m( ( dym*dym*dym - dym) * m_DeltaNorm[ coord::GetAxes().e_Y]);

        fvalue =

            dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Values( (i-1)%dimx  , j%dimy))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy)+dyp*m_Values(i%dimx      , j%dimy))
          +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy)+dyp*m_Dsecox(i%dimx      , j%dimy))
          + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy)+dy3p*m_Dsecoy(i%dimx     , j%dimy))
          +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy)+dy3p*m_Dsecoxy(i%dimx    , j%dimy))
        ;

        dfdxvalue =
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy)+dyp*m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy)+dyp*m_Dsecox(i%dimx      , j%dimy))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy)+dy3p*m_Dsecoy(i%dimx     , j%dimy))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy)+dy3p*m_Dsecoxy(i%dimx    , j%dimy))
        ;

        dfdyvalue =
          dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy)+m_Values( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy)+m_Values(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy)+m_Dsecox( (i-1)%dimx  , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy)+m_Dsecox(i%dimx      , j%dimy))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy))* m_Delta[ coord::GetAxes().e_Y]/ 6
        ;
      };

      double dfvalues[] = { dfdxvalue, dfdyvalue};
      linal::Vector< double> dfvector( 2, dfvalues);

      return storage::Pair< double, linal::Vector< double> >( fvalue, dfvector);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write BicubicSpline into std::ostream
    std::ostream &BicubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      OSTREAM << m_Border[ coord::GetAxes().e_X]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Y]           << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Values                   << '\n';
      OSTREAM << m_Dsecox                   << '\n';
      OSTREAM << m_Dsecoy                   << '\n';
      OSTREAM << m_Dsecoxy                  << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_X]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Y]          << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].Second() << '\n';

      // end
      return OSTREAM;
    }

    //! read BicubicSpline from std::istream
    std::istream &BicubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_X] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Y] = SplineBorderType( border);
      ISTREAM >> m_Start[ coord::GetAxes().e_X];
      ISTREAM >> m_Delta[ coord::GetAxes().e_X];
      ISTREAM >> m_Start[ coord::GetAxes().e_Y];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Y];
      ISTREAM >> m_Values;
      ISTREAM >> m_Dsecox;
      ISTREAM >> m_Dsecoy;
      ISTREAM >> m_Dsecoxy;
      ISTREAM >> m_LinCont[ coord::GetAxes().e_X];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Y];
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].Second();
      m_DeltaNorm[ coord::GetAxes().e_X] = Sqr( m_Delta[ coord::GetAxes().e_X]) / 6.0;
      m_DeltaNorm[ coord::GetAxes().e_Y] = Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6.0;

      // end
      return ISTREAM;
    }

  ////////////////
  // operations //
  ////////////////

  } // namespace math
} // namespace bcl
