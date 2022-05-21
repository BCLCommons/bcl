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
#include "math/bcl_math_cubic_spline.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_smooth_data.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CubicSpline::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSpline())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSpline
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @return trained spline
    CubicSpline &CubicSpline::Train
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::VectorConstInterface< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      // check, if the points are given in positive direction
      BCL_Assert( DELTA > 0, "deltax <= 0 not supported");

      // determine value for dimension of x
      const int dim( RESULTS.GetSize());

      // assigning values
      m_Border = BORDER;
      m_Start  = START;
      m_Delta  = DELTA;

      BCL_Assert( BORDER != e_NotAKnot, "The specified boundary condition is not yet supported");

      // auxiliary variables
      const double delta1( DELTA / 6), delta2( 2 * DELTA / 3);

      m_Values = RESULTS;

      linal::Matrix< double> coeffs( dim, dim);
      linal::Vector< double> derivs( dim, 0.0);

      // train once for the values of fxx
      // those values are equivalent for every type of spline considered here
      for( int i( 1); i < dim - 1; ++i)
      {
        coeffs( i, i - 1) = delta1;
        coeffs( i, i    ) = delta2;
        coeffs( i, i + 1) = delta1;

        derivs( i) = ( m_Values( i + 1) - 2 * m_Values( i) + m_Values( i - 1)) / DELTA;
      }

      // setting the second order derivative on start and end to 0, "natural" cubic spline
      if( m_Border == e_Natural)
      {
        coeffs(     0,     0) = 1;
        coeffs( dim-1, dim-1) = 1;
        // natural borders yields a simple tridiagonal matrix that can be readily solved using this specialty
        // function
        m_Dsecox = linal::MatrixInversionInterface< double>::SolveTridiagonalMatrix( coeffs, derivs);
        return *this;
      }

      // periodic, the function starts over after reaching the ends, continuously differentiable everywhere
      if( m_Border == e_Periodic)
      {
        coeffs(     0, dim - 1) = delta1;
        coeffs(     0,       0) = delta2;
        coeffs(     0,       1) = delta1;
        derivs(     0)          = ( m_Values( 1) - 2 * m_Values( 0) + m_Values( dim-1)) / DELTA;

        coeffs( dim - 1, dim - 2) = delta1;
        coeffs( dim - 1, dim - 1) = delta2;
        coeffs( dim - 1,       0) = delta1;
        derivs( dim - 1)          = ( m_Values( 0) - 2 * m_Values( dim-1) + m_Values( dim-2)) / DELTA;
      }

      // set the first order derivative at x_0 to first_start and at x_dim-1 to first_end
      if( m_Border == e_FirstDer)
      {
        coeffs(       0,       0) = -delta2/2;
        coeffs(       0,       1) = -delta1;
        derivs(       0)          = FIRSTBE.First() - ( m_Values( 1) - m_Values( 0)) / DELTA;

        coeffs( dim - 1, dim - 1) = delta2 / 2;
        coeffs( dim - 1, dim - 2) = delta1;
        derivs( dim - 1)          = FIRSTBE.Second() - ( m_Values( dim - 1) - m_Values( dim - 2)) / DELTA;
      }

      // computation of the second order derivatives in every given point of the spline
      // The spline matrix is always very well-conditioned, so pivoting is a waste of time, so skip pivoting
      // in the GJ solver by passing it false
      m_Dsecox = linal::MatrixInversionGaussJordan< double>( coeffs, false).Solve( derivs);

      return *this;
    }

    //! @brief train CubicSpline with pre processed data
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @param PREPROC determines the degree of smoothing, should be between 0 and 1 to give the part/influence
    //! of the actual value to the pre processed value
    //! @return trained spline
    CubicSpline &CubicSpline::TrainWithPreprocessing
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::VectorConstInterface< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE, const double PREPROC
    )
    {
      // passing a smoothed vector to the Train function
      return Train( BORDER, START, DELTA, SmoothData::SmoothVector( RESULTS, PREPROC, true), FIRSTBE);
    }

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    double CubicSpline::dF( const double &ARGUMENT) const
    {

      // number of grid points
      const int dim( m_Values.GetSize());

      int i( int( floor( ( ARGUMENT - m_Start) / m_Delta)) + 1);

      // not within supporting points - left
      if( i < 1)
      {
        // if the spline is periodic, adjust i to be positive ( > 0) and within range
        if( m_Border == e_Periodic)
        {
          // bring close to actual range
          i %= dim;

          // if between end and start
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval
            if( dxp < 0.0)
            {
              dxp += 1.0;
            }

            // return
            return Derivative( dim - 1, 0, dxp);
          }
          // generate positive index value ( > 0)
          while( i < 1)
          {
            i += dim;
          }
        }
        else
        {
          return Derivative( 0, 1, double( 0));
        }
      }
      // not within supporting points - right
      else if( i >= dim)
      {
        const int end( dim - 1);

        // if the spline is periodic, adjust i to be positive ( > 0)
        if( m_Border == e_Periodic)
        {
          // generate index value within range
          i %= dim;

          // special case, where interpolation happens between end and beginning
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            const double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval

            // return
            return Derivative( end, 0, dxp);
          }
        }
        else
        {
          // derivative at last supporting points
          return Derivative( end - 1, end, 1.0);
        }
      }

      double fn_mod( fmod( ARGUMENT - m_Start, m_Delta));
      double dxp( fn_mod / m_Delta);

      // see Numerical recipes in C++, pages 116-118
      //double dxp( math::Fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // delta_akt is an offset from the beginning of the actual interval\n";

      // for numerical reason the fmod function sometimes give a value less than m_Delta
      // divison will have worked correctly, so need to reset dxp to 0.0.
      if( Absolute( fn_mod - m_Delta) < std::numeric_limits< double>::epsilon())
      {
        dxp = 0.0;
      }
      else if( dxp < 0.0)
      {
        dxp += 1.0;
      }

      return Derivative( i - 1, i, dxp);
    }

    //! @brief return derivative and value at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return value and derivative at ARGUMENT
    storage::Pair< double, double> CubicSpline::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSpline::operator()( const double &ARGUMENT) const
    {
      // number of grid points
      const int dim( m_Values.GetSize());

      // determine i with m_Start+(i-1)*m_Delta < ARGUMENT < m_Start+i*m_Delta for the correct supporting points
      int i( int( floor( ( ARGUMENT - m_Start) / m_Delta)) + 1);

      // not within supporting points - left
      if( i < 1)
      {
        // if the spline is periodic, adjust i to be positive ( > 0) and within range
        if( m_Border == e_Periodic)
        {
          // bring close to actual range
          i %= dim;

          // if between end and start
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval
            if( dxp < 0.0)
            {
              dxp += 1.0;
            }

            // return
            return Function( dim - 1, 0, dxp);
          }
          // generate positive index value ( > 0)
          while( i < 1)
          {
            i += dim;
          }
        }
        else
        {
          return Function( 0, 1, double( 0)) + ( ARGUMENT - m_Start) * Derivative( 0, 1, double( 0));
        }
      }
      // not within supporting points - right
      else if( i >= dim)
      {
        const int end( dim - 1);

        // if the spline is periodic, adjust i to be positive ( > 0)
        if( m_Border == e_Periodic)
        {
          // generate index value within range
          i %= dim;

          // special case, where interpolation happens between end and beginning
          if( i == 0)
          {
            // see Numerical recipes in C++, pages 116-118
            const double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // relative offset from the beginning of the actual interval

            // return
            return Function( end, 0, dxp);
          }
        }
        else
        {
          // derivative at last supporting points
          return Function( end - 1, end, 1.0) + ( ARGUMENT - m_Start - ( dim - 1) * m_Delta) * Derivative( end - 1, end, 1.0);
        }
      }

      // see Numerical recipes in C++, pages 116-118
      double dxp( fmod( ARGUMENT - m_Start, m_Delta) / m_Delta); // delta_akt is an offset from the beginning of the actual interval\n";
      if( dxp < 0.0)
      {
        dxp += 1.0;
      }

      return Function( i - 1, i, dxp);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSpline from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSpline::Read( std::istream &ISTREAM)
    {
      // read parameters
      int border;
      ISTREAM >> border;
      m_Border = SplineBorderType( border);
      io::Serialize::Read( m_Start , ISTREAM);
      io::Serialize::Read( m_Delta , ISTREAM);
      io::Serialize::Read( m_Values, ISTREAM);
      io::Serialize::Read( m_Dsecox, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSpline into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSpline::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Border, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Start , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Delta , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Values, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Dsecox, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSpline::Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      // see Numerical recipes in C++, pages 116-118
      // relative distance from right grid point
      const double dxm( 1 - DXP);

      return
          ( m_Values( INDEX_RIGHT) - m_Values( INDEX_LEFT)) / m_Delta
        - ( 3 * dxm * dxm - 1) / 6 * m_Delta * m_Dsecox( INDEX_LEFT)
        + ( 3 * DXP * DXP - 1) / 6 * m_Delta * m_Dsecox( INDEX_RIGHT);
    }

    //! @brief calculate function between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return function depending on relative distance DXP
    double CubicSpline::Function( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      // see Numerical recipes in C++, pages 116-118
      // relative distance from right grid point
      const double dxm( 1 - DXP);
      const double dx3p( ( DXP * DXP * DXP - DXP) * Sqr( m_Delta) / 6); // =0 at the gridpoints, adds cubic part of the spline
      const double dx3m( ( dxm * dxm * dxm - dxm) * Sqr( m_Delta) / 6); // =0 at the gridpoints, adds cubic part of the spline

      return
          dxm * m_Values( INDEX_LEFT) + DXP * m_Values( INDEX_RIGHT)
        + dx3m * m_Dsecox( INDEX_LEFT) + dx3p * m_Dsecox( INDEX_RIGHT);
    }

  } // namespace math
} // namespace bcl
