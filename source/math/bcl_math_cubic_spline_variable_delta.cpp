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
#include "math/bcl_math_cubic_spline_variable_delta.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "math/bcl_math_running_average.h"
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
    const util::SiPtr< const util::ObjectInterface> CubicSplineVariableDelta::s_Instance
    (
      GetObjectInstances().AddInstance( new CubicSplineVariableDelta())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief train CubicSplineVariableDelta
    //! @param X abscissa Values of dataset
    //! @param Y ordinate values of dataset
    //! @return trained spline
    CubicSplineVariableDelta &CubicSplineVariableDelta::Train
    (
      const SplineBorderType BORDER,
      const linal::VectorConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      m_Border = BORDER;

      // test whether the x-range is unique.
      size_t n_duplicates( 0);
      if( X.GetSize())
      {
        double prev( X( 0));
        for( size_t i( 1), n_values( X.GetSize()); i < n_values; ++i)
        {
          if( X( i) == prev)
          {
            ++n_duplicates;
          }
          BCL_Assert( X( i) >= prev, "cubic spline requires sorted X vector!");
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

      BCL_Assert( BORDER == e_NotAKnot || BORDER == e_Natural, "The specified boundary condition is not yet supported");

      // Compute 2nd derivatives at each point
      GenerateSecondDerivative();

      return *this;
    }

    //! @brief train CubicSpline for Constant Delta.  This overloads the train function with the same input from the
    //         original implementation of cublic spline
    //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
    //! @param START the start of the interval the spline is defined on
    //! @param DELTA the distance between two support points of the spline
    //! @param RESULTS the function values at the support points of the spline
    //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
    //! @return trained spline
    CubicSplineVariableDelta &CubicSplineVariableDelta::Train
    (
      const SplineBorderType BORDER,
      const double START,
      const double DELTA,
      const linal::Vector< double> &RESULTS,
      const storage::Pair< double, double> &FIRSTBE
    )
    {
      m_X = linal::FillVector( RESULTS.GetSize(), START, DELTA);
      m_Y = RESULTS;
      m_Border = BORDER;

      // Compute 2nd derivatives at each point
      GenerateSecondDerivative();

      return *this;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an value in the range of the trained cubic spline and returning what the function value
    //! @brief would be if there were a function
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double CubicSplineVariableDelta::operator()( const double &ARGUMENT) const
    {
      const size_t length( m_X.GetSize());
      if( length == size_t( 1))
      {
        // 1 point; no spline, assume constant value
        return m_Y( 0);
      }

      // !Todo should we extend the spline to extrapolation in addition to interpolation
      BCL_Assert( ARGUMENT >= m_X( 0) && ARGUMENT <= m_X( length - 1), "The input is outside the trained spline");

      double A0, A1, A2, A3;
      size_t index( 0);

      // find the interval which contains ARGUMENT ( the sampled point)
      while( index < length && ARGUMENT > m_X( index))
      {
        ++index;
      }

      // handle case point as at a gridpoint. This is necessary only due to the bounds assertion above
      if( index < length && ARGUMENT == m_X( index))
      {
        return m_Y( index);
      }

      --index;

      const double dx( m_X( index + 1) - m_X( index));
      const double dy( m_Y( index + 1) - m_Y( index));

      A0 = m_Y( index);
      A1 = dy / dx - ( dx / 6.0) * ( m_SecondDerivative( index + 1) + 2.0 * m_SecondDerivative( index));
      A2 = m_SecondDerivative( index) / 2.0;
      A3 = ( m_SecondDerivative( index + 1) - m_SecondDerivative( index)) / ( 6.0 * dx);

      // Horner's rule to calculate cubic polynomial
      const double x( ARGUMENT - m_X( index));
      return ( ( A3 * x + A2) * x + A1) * x + A0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read CubicSplineVariableDelta from std::istream
    //! @param ISTREAM the stream to read the spline from
    //! @return istream after reading
    std::istream &CubicSplineVariableDelta::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_X , ISTREAM);
      io::Serialize::Read( m_Y , ISTREAM);
      io::Serialize::Read( m_SecondDerivative, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write CubicSplineVariableDelta into std::ostream
    //! @param OSTREAM the stream to write the spline to
    //! @param INDENT indentation of the spline
    //! @return ostream after writing
    std::ostream &CubicSplineVariableDelta::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Y, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SecondDerivative, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ///////////////////////
  // helper  functions //
  ///////////////////////

    //! @brief use m_x and m_y and initinal yd vectors to set m_SecondDerivative
    void CubicSplineVariableDelta::GenerateSecondDerivative()
    {
      size_t vector_size( m_X.GetSize());
      if( vector_size < size_t( 3))
      {
        // handle 0-2 points
        m_SecondDerivative = linal::Vector< double>( vector_size, double( 0.0));
        return;
      }

      linal::Vector< double> yd( vector_size);
      double h0, h1, r0, r1;

      // initialize first row data

      // delta x0
      h0 = m_X( 1) - m_X( 0);

      // delta x1
      h1 = m_X( 2) - m_X( 1);

      // delta y0
      r0 = ( m_Y( 1) - m_Y( 0)) / h0;

      // delta y1
      r1 = ( m_Y( 2) - m_Y( 1)) / h1;

      // This matrix will be filled with the coefficients of the derivative equations
      linal::Matrix< double> coeffs( vector_size, vector_size);

      if( m_Border == e_NotAKnot)
      {
        coeffs( 0, 0) = -h1;
        coeffs( 0, 1) = h0 + h1;
        coeffs( 0, 2) = -h0;
        yd( 0) = ( 0.0);
      }
      else // ( m_Border == e_Natural)
      {
        // Set B Value for initial point
         coeffs( 0, 0) = 1.0;
         yd( 0) = ( 0.0);
      }

      // initialize interior data
      for( size_t i( 1); i < vector_size - 1; i++)
      {
        h0 = m_X( i) - m_X( i - 1);
        h1 = m_X( i + 1) - m_X( i);
        r0 = ( m_Y( i) - m_Y( i - 1)) / h0;
        r1 = ( m_Y( i + 1) - m_Y( i)) / h1;

        // Set A Value
        coeffs( i, i - 1) = h0;

        // Set B Value
        coeffs( i, i) = 2 * ( h0 + h1);

        // Set C Value
        coeffs( i, i + 1) = h1;

        yd( i) = 6 * ( r1 - r0);
      }

      // initialize last row data

      if( m_Border == e_NotAKnot)
      {
        // Set A Value
        coeffs( vector_size - 1, vector_size - 3) = -h1;

        // Set B Value
        coeffs( vector_size - 1, vector_size - 2) = h0 + h1;

        // Set C Value
        coeffs( vector_size - 1, vector_size - 1) = -h0;

        yd( vector_size - 1) = 0.0;
        m_SecondDerivative = linal::MatrixInversionGaussJordan< double>( coeffs, false).Solve( yd);
      }
      else // if ( m_Border == e_Natural)
      {
        // Set C Value
        coeffs( vector_size - 1, vector_size - 1) = 1.0;
        yd( vector_size - 1) = 0.0;

        // computation of the second order derivatives in every given point of the spline
        // natural borders yields a simple tridiagonal matrix that can be readily solved using this specialty
        // function
        m_SecondDerivative = linal::MatrixInversionInterface< double>::SolveTridiagonalMatrix( coeffs, yd);
      }
    }

    //! @brief return derivative at certain ARGUMENT
    //! @param ARGUMENT x value
    //! @return derivative at ARGUMENT
    const double CubicSplineVariableDelta::dF( const double &ARGUMENT) const
    {

      size_t length( m_X.GetSize());

      // if outside left boundary
      if( ARGUMENT < m_X( 0))
      {
        return Derivative( 0, 1, double( 0.0));
      }
      // if outside right boundary
      else if( ARGUMENT >= m_X( length - 1))
      {
        const int end( length - 1);
        return Derivative( end - 1, end, 1.0);
      }

      // locate left index
      size_t i( 0);
      int left_index( 0);

      while( m_X( i) <= ARGUMENT)
      {
        left_index = i;
        ++i;
      }

      // Assign right index
      int right_index( left_index + 1);

      double left_border( m_X( left_index));
      double right_border( m_X( right_index));
      double delta( right_border - left_border);

      // Set the offset from the beginning of the actual interval
      double dxp( ( ARGUMENT - m_X( left_index)) / delta);

      return Derivative( left_index, right_index, dxp);
    }

    storage::Pair< double, double> CubicSplineVariableDelta::FdF( const double &ARGUMENT) const
    {
      return storage::Pair< double, double>( operator()( ARGUMENT), dF( ARGUMENT));
    }

    //! @brief calculate derivative between two cells
    //! @param INDEX_LEFT index of left grid point
    //! @param INDEX_RIGHT index of right grid point
    //! @param DXP relative distance from left grid point, must be element [0, 1]
    //! @return derivative depending on relative distance DXP
    double CubicSplineVariableDelta::Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const
    {
      const double dxm( 1 - DXP);
      const double delta( m_X( INDEX_RIGHT) - m_X( INDEX_LEFT));

      return
          ( m_Y( INDEX_RIGHT) - m_Y( INDEX_LEFT)) / delta
        - ( 3.0 * dxm * dxm - 1) / 6.0 * delta * m_SecondDerivative( INDEX_LEFT)
        + ( 3.0 * DXP * DXP - 1) / 6.0 * delta * m_SecondDerivative( INDEX_RIGHT);
    }

  } // namespace math
} // namespace bcl
