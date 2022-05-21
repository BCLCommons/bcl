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
#include "math/bcl_math_cubic_spline_damped.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CubicSplineDamped::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSplineDamped())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSplineDamped
    //! @param X abscissa Values of dataset
    //! @param Y ordinate values of dataset
    //! @return trained spline
    CubicSplineDamped &CubicSplineDamped::Train
    (
      const linal::VectorConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y,
      const double &DY_START,
      const double &DY_END
    )
    {
      m_ConstantDelta = false;
      m_Delta = 0.0;
      BCL_Assert( X.GetSize() == Y.GetSize(), "Train called with non-equal sized vectors!");
      // test whether the x-range is unique.
      size_t n_duplicates( 0);
      if( X.GetSize() > size_t( 1))
      {
        double prev( X( 0));
        for( size_t i( 1), n_values( X.GetSize()); i < n_values; ++i)
        {
          if( X( i) == prev)
          {
            ++n_duplicates;
          }
          BCL_Assert( X( i) >= prev, "cubic spline requires sorted X vector! Given vector was: " + util::Format()( X));
          prev = X( i);
        }
      }
      if( !n_duplicates)
      {
        m_X = X;
        m_Y = Y;
      }
      else
      {
        // duplicates were found. Copy unique m_X, average Y across duplicate m_X
        const size_t n_uniq( X.GetSize() - n_duplicates);
        m_X = linal::Vector< double>( n_uniq);
        m_Y = linal::Vector< double>( n_uniq);
        double prev( std::numeric_limits< double>::quiet_NaN());
        for( size_t i( 0), n_values( X.GetSize()), uniq_i( 0); i < n_values;)
        {
          if( X( i) == prev)
          {
            // average Y across all shared X values
            RunningAverage< double> ave;
            ave += Y( uniq_i - 1);
            ave += Y( i);
            while( ++i < n_values && X( i) == prev)
            {
              ave += Y( i);
            }
            m_Y( uniq_i - 1) = ave.GetAverage();
          }
          else
          {
            prev = X( i);
            m_X( uniq_i) = X( i);
            m_Y( uniq_i) = Y( i);
            ++i;
            ++uniq_i;
          }
        }
      }
      m_Delta = ( m_X.Last() - m_X.First()) / double( std::max( m_X.GetSize(), size_t( 2)) - 1);

      // determine derivatives
      const size_t size( m_X.GetSize()), sizem1( size - 1);
      m_ConstantDelta = EqualWithinMachineTolerance( m_X, linal::FillVector( size, m_X.First(), m_Delta));
      m_dY = linal::Vector< double>( size, double( 0.0));

      if( size <= size_t( 1))
      {
        m_A = m_B = linal::Vector< double>();
        return *this;
      }
      m_A = m_B = linal::Vector< double>( sizem1, double( 0.0));

      linal::Vector< double> secants( m_dY), intervals( m_dY);
      for( size_t i( 0); i < sizem1; ++i)
      {
        intervals( i) = m_X( i + 1) - m_X( i);
        secants( i) = ( m_Y( i + 1) - m_Y( i)) / intervals( i);
      }

      // boundary condition; specified derivative at RHS
      if( size > 2)
      {
        intervals( size - 1) = intervals( size - 2);
        secants( size - 1) = util::IsDefined( DY_END) ? DY_END : secants( size - 2);
      }
      else
      {
        secants( size - 1) = util::IsDefined( DY_END) ? DY_END : 0.0;
      }
      for( size_t i( 1); i < size; ++i)
      {
        const bool going_down( secants( i - 1) <= 0.0);
        if( going_down != ( secants( i) <= 0.0))
        {
          // local maxima or minima in the input data / derivative must be 0
          continue;
        }
        // compute the secant through i-1 - i+1
        const double secant_ddy
        (
          ( secants( i - 1) * intervals( i) + secants( i) * intervals( i - 1))
          /
          ( intervals( i) + intervals( i - 1))
        );
        // derivative magnitude is the smaller absolute of the twice the adjacent secants, or the average secant
        m_dY( i) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * std::max( secants( i - 1), secants( i)))
            : std::min( secant_ddy, 2.0 * std::min( secants( i - 1), secants( i)))
          );
      }
      if( !util::IsDefined( DY_START))
      {
        // compute the secant through 0 - 1
        // equation 24
        const double secant_ddy
        (
          ( secants( 0) * intervals( 0) * ( 2.0 + intervals( 1) / intervals( 0)) - secants( 1))
          /
          ( intervals( 0) + intervals( 1))
        );
        const bool going_down( secants( 0) <= 0.0);
        if( going_down != ( secant_ddy < 0.0))
        {
          m_dY( 0) = 0.0;
        }
        else
        {
          m_dY( 0) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * secants( 0))
            : std::min( secant_ddy, 2.0 * secants( 0))
          );
        }
      }
      else
      {
        m_dY( 0) = DY_START;
      }
      if( !util::IsDefined( DY_END) && size > 2)
      {
        // equation 25
        const double interval_left( intervals( size - 3)), interval_right( intervals( size - 2));
        const double secants_final( secants( size - 2));
        const double secant_ddy
        (
          ( secants_final * interval_right * ( 2.0 + interval_left / interval_right) - secants( size - 3))
          /
          ( interval_right + interval_left)
        );
        const bool going_down( secants_final <= 0.0);
        if( going_down != ( secant_ddy < 0.0))
        {
          m_dY( size - 1) = 0.0;
        }
        else
        {
          m_dY( size - 1) =
          (
            going_down
            ? std::max( secant_ddy, 2.0 * secants_final)
            : std::min( secant_ddy, 2.0 * secants_final)
          );
        }
      }
      for( size_t i( 0); i < sizem1; ++i)
      {
        m_B( i) = ( 3.0 * secants( i) - 2.0 * m_dY( i) - m_dY( i + 1)) / intervals( i);
        m_A( i) = ( m_dY( i) + m_dY( i + 1) - 2.0 * secants( i)) / ( intervals( i) * intervals( i));
      }
      return *this;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an value in the range of the trained cubic spline and returning what the function value
    //! @brief would be if there were a function
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSplineDamped::operator()( const double &ARGUMENT) const
    {
      const size_t length( m_X.GetSize());
      if( length == size_t( 1))
      {
        // 1 point; no spline, assume constant value
        return m_Y( 0);
      }

      if( ARGUMENT <= m_X( 0))
      {
        return m_Y( 0) + m_dY( 0) * ( m_X( 0) - ARGUMENT);
      }
      else if( ARGUMENT >= m_X( length - 1))
      {
        return m_Y( length - 1) + m_dY( length - 1) * ( ARGUMENT - m_X( length - 1));
      }
      else if( !util::IsDefined( ARGUMENT))
      {
        return util::GetUndefined< double>();
      }

      size_t index( 0);

      if( !m_ConstantDelta)
      {
        // find the interval which contains ARGUMENT ( the sampled point)
        while( index < length && ARGUMENT > m_X( index))
        {
          ++index;
        }

        // handle case point as at a gridpoint
        if( index < length && ARGUMENT == m_X( index))
        {
          return m_Y( index);
        }

        --index;
        if( index >= m_X.GetSize())
        {
          BCL_MessageStd( "BAD INDEX A: X: " + util::Format()( ARGUMENT) + " " + util::Format()( m_X( 0)));
        }
      }
      else
      {
        // determine i with m_Start+(i-1)*m_Delta < ARGUMENT < m_Start+i*m_Delta for the correct supporting points
        index = int( floor( ( ARGUMENT - m_X( 0)) / m_Delta));
        if( index >= m_X.GetSize())
        {
          BCL_MessageStd( "BAD INDEX B: X: " + util::Format()( ARGUMENT) + " " + util::Format()( m_X( 0)));
        }
      }

      // Horner's rule to calculate cubic polynomial
      const double x( ARGUMENT - m_X( index));
      return ( ( m_A( index) * x + m_B( index)) * x + m_dY( index)) * x + m_Y( index);
    }

    //! @brief train CubicSpline, requiring a maximally smooth, piecewise monotonic solution
    //! @param START x-value for the start
    //! @param DELTA difference between consecutive X-points
    //! @param Y ordinate values of dataset
    //! @param DY_START first order derivative at X(0). if nan/undefined, will be the same as first order derivative at X(1)
    //! @param DY_END   first order derivative at X_last. if nan/undefined, will be the same as first order derivative at X(last-1)
    //! @return *this
    //! Reference: M. Steffen. "A simple method for monotonic interpolation in one dimension" Astronomy and
    //! Astrophysics. 1990
    //! Note that it is not required that the Y be monotonic. All this function does is ensure that each cubic piece
    //! of the function is monotonic over the given interval and keeps the y' continuous. y'' will not, in general, be
    //! continuous if the spline is trained with this function
    CubicSplineDamped &CubicSplineDamped::Train
    (
      const double &START,
      const double &DELTA,
      const linal::VectorConstInterface< double> &Y,
      const double &DY_START,
      const double &DY_END
    )
    {
      return this->Train( linal::FillVector( Y.GetSize(), START, DELTA), Y, DY_START, DY_END);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSplineDamped from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSplineDamped::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_X , ISTREAM);
      io::Serialize::Read( m_Y , ISTREAM);
      io::Serialize::Read( m_dY, ISTREAM);
      io::Serialize::Read( m_A, ISTREAM);
      io::Serialize::Read( m_B, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSplineDamped into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSplineDamped::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_dY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_A, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_B, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    const double CubicSplineDamped::dF( const double &ARGUMENT) const
    {
      size_t length( m_X.GetSize());

      // if outside left boundary
      if( ARGUMENT < m_X( 0))
      {
        return m_dY( 0);
      }
      // if outside right boundary
      else if( ARGUMENT >= m_X( length - 1))
      {
        return m_dY( length - 1);
      }

      // locate left index
      size_t left_index( 0);
      if( !m_ConstantDelta)
      {
        for( size_t i( 0); m_X( i) <= ARGUMENT; ++i)
        {
          left_index = i;
        }
      }
      else
      {
        left_index = int( floor( ( ARGUMENT - m_X( 0)) / m_Delta));
      }

      // Set the offset from the beginning of the actual interval
      return Derivative( left_index, ( ARGUMENT - m_X( left_index)));
    }

    storage::Pair< double, double> CubicSplineDamped::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSplineDamped::Derivative( const int INDEX_LEFT, const double DXP) const
    {
      return m_dY( INDEX_LEFT) + DXP * ( 2.0 * m_B( INDEX_LEFT) + 3.0 * m_A( INDEX_LEFT) * DXP);
    }

  } // namespace math
} // namespace bcl
