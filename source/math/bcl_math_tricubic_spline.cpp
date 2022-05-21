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
#include "math/bcl_math_tricubic_spline.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_cubic_spline.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> TricubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new TricubicSpline())
    );

  ////////////////
  // operations //
  ////////////////

    //! train TricubicSpline
    void TricubicSpline::Train
    (
      const SplineBorderType BORDER[ 3], const double START[ 3], const double DELTA[ 3], const Tensor< double> &RESULTS,
      const bool LINCONT[ 3], const storage::Pair< double, double> FIRSTBE[ 3]
    )
    {
      //check, if the points are given in positive direction
      BCL_Assert( DELTA[ coord::GetAxes().e_X]>0 , "deltax <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Y]>0 , "deltay <= 0 not supported");
      BCL_Assert( DELTA[ coord::GetAxes().e_Z]>0 , "deltaz <= 0 not supported");

      //determine values for all three dimensions
      const int dimz( RESULTS.GetNumberCols());
      const int dimy( RESULTS.GetNumberRows());
      const int dimx( RESULTS.NumberLayers());

      //assigning values
      m_Border[ coord::GetAxes().e_X]  = BORDER[ coord::GetAxes().e_X];
      m_Border[ coord::GetAxes().e_Y]  = BORDER[ coord::GetAxes().e_Y];
      m_Border[ coord::GetAxes().e_Z]  = BORDER[ coord::GetAxes().e_Z];
      m_Start[ coord::GetAxes().e_X]   = START[ coord::GetAxes().e_X];
      m_Start[ coord::GetAxes().e_Y]   = START[ coord::GetAxes().e_Y];
      m_Start[ coord::GetAxes().e_Z]   = START[ coord::GetAxes().e_Z];
      m_Delta[ coord::GetAxes().e_X]   = DELTA[ coord::GetAxes().e_X];
      m_Delta[ coord::GetAxes().e_Y]   = DELTA[ coord::GetAxes().e_Y];
      m_Delta[ coord::GetAxes().e_Z]   = DELTA[ coord::GetAxes().e_Z];

      BCL_Assert
      (
        m_Border[0] != e_NotAKnot &&
        m_Border[1] != e_NotAKnot &&
        m_Border[2] != e_NotAKnot ,
        "The specified boundary condition is not yet supported"
      );

      m_Values   = RESULTS;
      m_Dsecox   = RESULTS;
      m_Dsecoy   = RESULTS;
      m_Dsecoz   = RESULTS;
      m_Dsecoxy  = RESULTS;
      m_Dsecoxz  = RESULTS;
      m_Dsecoyz  = RESULTS;
      m_Dsecoxyz = RESULTS;

      m_LinCont[ coord::GetAxes().e_X] = LINCONT[ coord::GetAxes().e_X];
      m_LinCont[ coord::GetAxes().e_Y] = LINCONT[ coord::GetAxes().e_Y];
      m_LinCont[ coord::GetAxes().e_Z] = LINCONT[ coord::GetAxes().e_Z];
      m_FirstBe[ coord::GetAxes().e_X] = FIRSTBE[ coord::GetAxes().e_X];
      m_FirstBe[ coord::GetAxes().e_Y] = FIRSTBE[ coord::GetAxes().e_Y];
      m_FirstBe[ coord::GetAxes().e_Z] = FIRSTBE[ coord::GetAxes().e_Z];

      //train seven times for fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz
      //reduction to Spline2D by training only 2D-layers at the same time
      for( int layer( 0); layer < dimx; ++layer)
      {
        linal::Matrix< double> values( dimy, dimz);
        for( int row( 0); row < dimy; ++row)
        {
          for( int col( 0); col < dimz; ++col)
          {
            values( row, col) = m_Values( layer, row, col);
          }
        }
        BicubicSpline bs;
        SplineBorderType border[2]            = { BORDER[ coord::GetAxes().e_Y],   BORDER[ coord::GetAxes().e_Z]};
        double start[2]                       = { START[ coord::GetAxes().e_Y],    START[ coord::GetAxes().e_Z]};
        double delta[2]                       = { DELTA[ coord::GetAxes().e_Y],    DELTA[ coord::GetAxes().e_Z]};
        bool lin_cont[2]                      = { LINCONT[ coord::GetAxes().e_Y],  LINCONT[ coord::GetAxes().e_Z]};
        storage::Pair< double, double> firstbe[ 2] = { FIRSTBE[ coord::GetAxes().e_Y], FIRSTBE[ coord::GetAxes().e_Z]};

        bs.Train( border, start, delta, values, lin_cont, firstbe);
        m_Dsecoy.ReplaceLayer(  layer, bs.GetDsecox());
        m_Dsecoz.ReplaceLayer(  layer, bs.GetDsecoy());
        m_Dsecoyz.ReplaceLayer( layer, bs.GetDsecoxy());
      }

      for( int row( 0); row < dimy; ++row)
      {
        for( int col( 0); col < dimz; ++col)
        {
          linal::Vector< double> values( dimx), dsecoz( dimx), dsecoy( dimx), dsecoyz( dimx);
          for( int layer( 0); layer < dimx; ++layer)
          {
            values(   layer) = m_Values(  layer, row, col);
            dsecoy(   layer) = m_Dsecoy(  layer, row, col);
            dsecoz(   layer) = m_Dsecoz(  layer, row, col);
            dsecoyz(  layer) = m_Dsecoyz( layer, row, col);
          }
          CubicSpline cs, csz, csy, csyz;
          cs.Train(   BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], values , FIRSTBE[ coord::GetAxes().e_X]);
          csy.Train(  BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoy , FIRSTBE[ coord::GetAxes().e_X]);
          csz.Train(  BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoz , FIRSTBE[ coord::GetAxes().e_X]);
          csyz.Train( BORDER[ coord::GetAxes().e_X], START[ coord::GetAxes().e_X], DELTA[ coord::GetAxes().e_X], dsecoyz, FIRSTBE[ coord::GetAxes().e_X]);
          for( int layer( 0); layer < dimx; ++layer)
          {
            m_Dsecox(   layer, row, col) = cs.GetDsecox()(    layer);
            m_Dsecoxy(  layer, row, col) = csy.GetDsecox()(   layer);
            m_Dsecoxz(  layer, row, col) = csz.GetDsecox()(   layer);
            m_Dsecoxyz( layer, row, col) = csyz.GetDsecox()(  layer);
          }
        }
      }
      return;
    }

    //! return value at certain (x, y, z)
    double TricubicSpline::F( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return F( m_Start[ coord::GetAxes().e_X], y, z) + ( x - m_Start[ coord::GetAxes().e_X]) * dFdx( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return F( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z)
                + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
                * dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y]) || ( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return F( x, m_Start[ coord::GetAxes().e_Y], z) + ( y - m_Start[ coord::GetAxes().e_Y]) * dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return F( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z)
                + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
                * dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return F( x, y, m_Start[ coord::GetAxes().e_Z]) + ( z - m_Start[ coord::GetAxes().e_Z]) * dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return F( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                + ( z - m_Start[ coord::GetAxes().e_Z] - ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                * dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X] for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //determine j with m_Start[ coord::GetAxes().e_Y]+(i-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+i*m_Delta[ coord::GetAxes().e_Y] for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      // determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Deltaz < z < m_Start[ coord::GetAxes().e_Z]+k*m_Deltaz for the correct supporting points
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      // this formula was derived from the Spline2D formula
      // the idea is to 'combine every part of the 2D formula
      // with dzm, dzp, dz3m, dz3p'

      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - (i-1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - (j-1) * m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - (k-1) * m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp*dxp*dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm*dxm*dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp*dyp*dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym*dym*dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp*dzp*dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm*dzm*dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indices
      while( i<1) i += dimx;
      while( j<1) j += dimy;
      while( k<1) k += dimz;

      return
        dzm
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))

        +dzp
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))

        +dz3m
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +dz3p
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
      ;
    }

    //! return partial derivative at certain (x, y, z) for x
    double TricubicSpline::dFdx( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdx( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdx( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdx( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z] || m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdx( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdx( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start+(i-1)*m_Delta < x < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //the same for j
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //the same for k
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);
      double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - ( k - 1) * m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indizes
      while( i<1) i += dimx;
      while( j<1) j += dimy;
      while( k<1) k += dimz;

      return
        dzm
        * (
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))
        )
        +dzp
        *(
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))
        )
        +dz3m
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))
        )
        +dz3p
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))
        )
      ;
    }

    //! return partial derivative at certain (x, y, z) for y
    double TricubicSpline::dFdy( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdy( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdy( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z] || m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdy( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdy( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X]
      //for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y]
      //for the correct supporting points
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Delta[ coord::GetAxes().e_Z] < z < m_Start[ coord::GetAxes().e_Z]+k*m_Delta[ coord::GetAxes().e_Z]
      //for the correct supporting points
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      //generate some auxiliary variables
      const double delta_aktx( x - m_Start[ coord::GetAxes().e_X] - ( i - 1) * m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y - m_Start[ coord::GetAxes().e_Y] - ( j - 1) * m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z - m_Start[ coord::GetAxes().e_Z] - ( k - 1) * m_Delta[ coord::GetAxes().e_Z]);

      //see F(x, y, z) for a short explanation of the values
      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);
      const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
      const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

      //generate positive values to prevent some problems with the indizes
      while( i <= 0) i += dimx;
      while( j <= 0) j += dimy;
      while( k <= 0) k += dimz;

      return
        dzm
        * (
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dzp
        *(
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, k%dimz)+m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecox(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dz3m
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        +dz3p
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoxz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
      ;
    }

    //! return partial derivative at certain (x, y, z) for z
    double TricubicSpline::dFdz( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());

      //check if argument x is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
        && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x)))
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          return dFdz( m_Start[ coord::GetAxes().e_X], y, z);
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          return dFdz( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X], y, z);
        }
      }

      //check if argument y is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
        && ( ( y < m_Start[ coord::GetAxes().e_Y] || m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          return dFdz( x, m_Start[ coord::GetAxes().e_Y], z);
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          return dFdz( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
        }
      }

      //check if argument z is in range for non-periodic splines
      if
      (
            ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
        && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
      )
      {
        BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          return dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          return dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
        }
      }

      //determine i with m_Start+(i-1)*m_Delta < x < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
      while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

      //the same for j and y
      int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
      while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

      //the same for k and z
      int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
      while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

      //generate some auxiliary variables
      const double delta_aktx( x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X]);
      const double delta_akty( y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y]);
      const double delta_aktz( z-m_Start[ coord::GetAxes().e_Z]-(k-1)*m_Delta[ coord::GetAxes().e_Z]);

      const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
      const double dxm( 1 - dxp);
      const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

      const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
      const double dym( 1 - dyp);
      const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
      const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

      const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
      const double dzm( 1 - dzp);

      //generate positive values to prevent some problems with the indizes
      while( i <= 0) i += dimx;
      while( j <= 0) j += dimy;
      while( k <= 0) k += dimz;

      return
        -(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        +(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        -(3*dzm*dzm-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +(3*dzp*dzp-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
        + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
        +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
        + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
        +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
      ;
    }

    //! return value and partial derivatives at certain (x, y, z)
    storage::Pair< double, linal::Vector< double> > TricubicSpline::FdF( const double x, const double y, const double z) const
    {
      const int dimx( m_Values.NumberLayers());
      const int dimy( m_Values.GetNumberRows());
      const int dimz( m_Values.GetNumberCols());
      double fvalue( 0), dfdxvalue( 0), dfdyvalue( 0), dfdzvalue( 0);

      //check if argument is in range for non-periodic splines
      if
      (
        (
              ( m_Border[ coord::GetAxes().e_X] != e_Periodic)
          && ( ( x < m_Start[ coord::GetAxes().e_X]) || ( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] < x))
        )
        ||
        (
              ( m_Border[ coord::GetAxes().e_Y] != e_Periodic)
          && ( ( y < m_Start[ coord::GetAxes().e_Y]) || ( m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y] < y))
        )
        ||
        (
              ( m_Border[ coord::GetAxes().e_Z] != e_Periodic)
          && ( ( z < m_Start[ coord::GetAxes().e_Z]) || ( m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z] < z))
        )
      )
      {
        BCL_MessageDbg( "argument out of range, using linear continuation");
        if( x < m_Start[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( m_Start[ coord::GetAxes().e_X], y, z);
          dfdyvalue = dFdy( m_Start[ coord::GetAxes().e_X], y, z);
          dfdzvalue = dFdz( m_Start[ coord::GetAxes().e_X], y, z);
          fvalue    = F(    m_Start[ coord::GetAxes().e_X], y, z) + ( x - m_Start[ coord::GetAxes().e_X]) * dfdxvalue;
        }
        if( x > m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_X], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          dfdyvalue = dFdy( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          dfdzvalue = dFdz( m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z);
          fvalue    = F(    m_Start[ coord::GetAxes().e_X] + ( dimx - 1) * m_Delta[ coord::GetAxes().e_X] , y, z)
                      + ( x - m_Start[ coord::GetAxes().e_X] - ( dimx - 1) * m_Delta[ coord::GetAxes().e_X]) * dfdxvalue;
        }
        if( y < m_Start[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, m_Start[ coord::GetAxes().e_Y], z);
          dfdyvalue = dFdy( x, m_Start[ coord::GetAxes().e_Y], z);
          dfdzvalue = dFdz( x, m_Start[ coord::GetAxes().e_Y], z);
          fvalue    = F(    x, m_Start[ coord::GetAxes().e_Y], z) + ( y - m_Start[ coord::GetAxes().e_Y]) * dfdyvalue;
        }
        if( y > m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Y], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          dfdyvalue = dFdy( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          dfdzvalue = dFdz( x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z);
          fvalue    = F(    x, m_Start[ coord::GetAxes().e_Y] + ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y], z)
                      + ( y - m_Start[ coord::GetAxes().e_Y] - ( dimy - 1) * m_Delta[ coord::GetAxes().e_Y]) * dfdyvalue;
        }
        if( z < m_Start[ coord::GetAxes().e_Z])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, y, m_Start[ coord::GetAxes().e_Z]);
          dfdyvalue = dFdy( x, y, m_Start[ coord::GetAxes().e_Z]);
          dfdzvalue = dFdz( x, y, m_Start[ coord::GetAxes().e_Z]);
          fvalue    = F(    x, y, m_Start[ coord::GetAxes().e_Z]) + ( z - m_Start[ coord::GetAxes().e_Z]) * dfdzvalue;
        }
        if( z > m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
        {
          BCL_Assert( m_LinCont[ coord::GetAxes().e_Z], "argument out of range for non-periodic spline!");
          dfdxvalue = dFdx( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          dfdyvalue = dFdy( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          dfdzvalue = dFdz( x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]);
          fvalue    = F(    x, y, m_Start[ coord::GetAxes().e_Z] + ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z])
                      + ( z - m_Start[ coord::GetAxes().e_Z] - ( dimz - 1) * m_Delta[ coord::GetAxes().e_Z]) * dfdzvalue;
        }
      }
      else
      {
        //determine i with m_Start[ coord::GetAxes().e_X]+(i-1)*m_Delta[ coord::GetAxes().e_X] < x < m_Start[ coord::GetAxes().e_X]+i*m_Delta[ coord::GetAxes().e_X]
        //for the correct supporting points
        int i( int( floor( ( x - m_Start[ coord::GetAxes().e_X]) / m_Delta[ coord::GetAxes().e_X])));
        while( m_Start[ coord::GetAxes().e_X] + i * m_Delta[ coord::GetAxes().e_X] < x) i++;

        //determine j with m_Start[ coord::GetAxes().e_Y]+(j-1)*m_Delta[ coord::GetAxes().e_Y] < y < m_Start[ coord::GetAxes().e_Y]+j*m_Delta[ coord::GetAxes().e_Y]
        //for the correct supporting points
        int j( int( floor( ( y - m_Start[ coord::GetAxes().e_Y]) / m_Delta[ coord::GetAxes().e_Y])));
        while( m_Start[ coord::GetAxes().e_Y] + j * m_Delta[ coord::GetAxes().e_Y] < y) j++;

        //determine k with m_Start[ coord::GetAxes().e_Z]+(k-1)*m_Delta[ coord::GetAxes().e_Z] < z < m_Start[ coord::GetAxes().e_Z]+k*m_Delta[ coord::GetAxes().e_Z]
        //for the correct supporting points
        int k( int( floor( ( z - m_Start[ coord::GetAxes().e_Z]) / m_Delta[ coord::GetAxes().e_Z])));
        while( m_Start[ coord::GetAxes().e_Z] + k * m_Delta[ coord::GetAxes().e_Z] < z) k++;

        //see method F(x,y,z) for detailed formula

        const double delta_aktx( x-m_Start[ coord::GetAxes().e_X]-(i-1)*m_Delta[ coord::GetAxes().e_X]);
        const double delta_akty( y-m_Start[ coord::GetAxes().e_Y]-(j-1)*m_Delta[ coord::GetAxes().e_Y]);
        const double delta_aktz( z-m_Start[ coord::GetAxes().e_Z]-(k-1)*m_Delta[ coord::GetAxes().e_Z]);

        const double dxp( delta_aktx / m_Delta[ coord::GetAxes().e_X]);
        const double dxm( 1 - dxp);
        const double dx3p( ( dxp * dxp * dxp - dxp) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);
        const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta[ coord::GetAxes().e_X]) / 6);

        const double dyp( delta_akty / m_Delta[ coord::GetAxes().e_Y]);
        const double dym( 1 - dyp);
        const double dy3p( ( dyp * dyp * dyp - dyp) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);
        const double dy3m( ( dym * dym * dym - dym) * Sqr( m_Delta[ coord::GetAxes().e_Y]) / 6);

        const double dzp( delta_aktz / m_Delta[ coord::GetAxes().e_Z]);
        const double dzm( 1 - dzp);
        const double dz3p( ( dzp * dzp * dzp - dzp) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);
        const double dz3m( ( dzm * dzm * dzm - dzm) * Sqr( m_Delta[ coord::GetAxes().e_Z]) / 6);

        //generate positive values to prevent some problems with the indizes
        while( i < 1) i += dimx;
        while( j < 1) j += dimy;
        while( k < 1) k += dimz;

        //compute values
        fvalue =
        dzm
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))

        +dzp
        *(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))

        +dz3m
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +dz3p
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
        ;

        dfdxvalue =
        dzm
        * (
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))
        )

        +dzp
        *(
          -(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))
        )

        +dz3m
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))
        )

        +dz3p
        *(
          -(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
          -(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          +(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_X]
          - (3 * dxm*dxm - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          + (3 * dxp*dxp - 1) * m_Delta[ coord::GetAxes().e_X] / 6*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))
        )
        ;

        dfdyvalue =
        dzm
        * (
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Values(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dzp
        *(
            dxm*(-m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Values( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Values(i%dimx      , (j-1)%dimy, k%dimz)+m_Values(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecox(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoy(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxy(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dz3m
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )

        +dz3p
        *(
            dxm*(-m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxp*(-m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3m*(-m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          +dx3p*(-m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+m_Dsecoxz(i%dimx      , j%dimy, k%dimz))/m_Delta[ coord::GetAxes().e_Y]
          + dxm*(-(3*dym*dym-1)*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          + dxp*(-(3*dym*dym-1)*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3m*(-(3*dym*dym-1)*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
          +dx3p*(-(3*dym*dym-1)*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz))* m_Delta[ coord::GetAxes().e_Y]/ 6
        )
        ;

        dfdzvalue =
        -(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Values(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, (k-1)%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        +(dxm*(dym*m_Values( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Values( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Values(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Values(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecox( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecox( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecox(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecox(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoy( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoy(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoy(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxy( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxy(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxy(i%dimx    , j%dimy, k%dimz)))
        / m_Delta[ coord::GetAxes().e_Z]

        -(3*dzm*dzm-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, (k-1)%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, (k-1)%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, (k-1)%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, (k-1)%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, (k-1)%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, (k-1)%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, (k-1)%dimz)))

        +(3*dzp*dzp-1)*m_Delta[ coord::GetAxes().e_Z]/ 6
        *(dxm*(dym*m_Dsecoz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz( (i-1)%dimx  , j%dimy, k%dimz))
          + dxp*(dym*m_Dsecoz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoz(i%dimx      , j%dimy, k%dimz))
        +dx3m*(dym*m_Dsecoxz( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz( (i-1)%dimx  , j%dimy, k%dimz))
          +dx3p*(dym*m_Dsecoxz(i%dimx      , (j-1)%dimy, k%dimz)+dyp*m_Dsecoxz(i%dimx      , j%dimy, k%dimz))
        + dxm*(dy3m*m_Dsecoyz( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz( (i-1)%dimx , j%dimy, k%dimz))
          + dxp*(dy3m*m_Dsecoyz(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoyz(i%dimx     , j%dimy, k%dimz))
        +dx3m*(dy3m*m_Dsecoxyz( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz( (i-1)%dimx, j%dimy, k%dimz))
          +dx3p*(dy3m*m_Dsecoxyz(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*m_Dsecoxyz(i%dimx    , j%dimy, k%dimz)))
        ;
      }

      return storage::Pair< double, linal::Vector< double> >( fvalue, linal::MakeVector< double>( dfdxvalue, dfdyvalue, dfdzvalue));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read TricubicSpline from std::istream
    std::istream &TricubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_X] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Y] = SplineBorderType( border);
      ISTREAM >> border;
      m_Border[ coord::GetAxes().e_Z] = SplineBorderType( border);
      ISTREAM >> m_Start[ coord::GetAxes().e_X];
      ISTREAM >> m_Start[ coord::GetAxes().e_Y];
      ISTREAM >> m_Start[ coord::GetAxes().e_Z];
      ISTREAM >> m_Delta[ coord::GetAxes().e_X];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Y];
      ISTREAM >> m_Delta[ coord::GetAxes().e_Z];
      ISTREAM >> m_Values;
      ISTREAM >> m_Dsecox;
      ISTREAM >> m_Dsecoy;
      ISTREAM >> m_Dsecoxy;
      ISTREAM >> m_Dsecoz;
      ISTREAM >> m_Dsecoxz;
      ISTREAM >> m_Dsecoyz;
      ISTREAM >> m_Dsecoxyz;
      ISTREAM >> m_LinCont[ coord::GetAxes().e_X];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Y];
      ISTREAM >> m_LinCont[ coord::GetAxes().e_Z];
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_X].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Y].Second();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Z].First();
      ISTREAM >> m_FirstBe[ coord::GetAxes().e_Z].Second();

      // end
      return ISTREAM;
    }

    //! write TricubicSpline into std::ostream
    std::ostream &TricubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      OSTREAM << m_Border[ coord::GetAxes().e_X]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Y]           << '\n';
      OSTREAM << m_Border[ coord::GetAxes().e_Z]           << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Start[ coord::GetAxes().e_Z]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_X]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Y]            << '\n';
      OSTREAM << m_Delta[ coord::GetAxes().e_Z]            << '\n';
      OSTREAM << m_Values                   << '\n';
      OSTREAM << m_Dsecox                   << '\n';
      OSTREAM << m_Dsecoy                   << '\n';
      OSTREAM << m_Dsecoxy                  << '\n';
      OSTREAM << m_Dsecoz                   << '\n';
      OSTREAM << m_Dsecoxz                  << '\n';
      OSTREAM << m_Dsecoyz                  << '\n';
      OSTREAM << m_Dsecoxyz                 << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_X]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Y]          << '\n';
      OSTREAM << m_LinCont[ coord::GetAxes().e_Z]          << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_X].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Y].Second() << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Z].First()  << '\n';
      OSTREAM << m_FirstBe[ coord::GetAxes().e_Z].Second() << '\n';

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
