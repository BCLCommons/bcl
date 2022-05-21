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

#ifndef BCL_MATH_CUBIC_SPLINE_VARIABLE_DELTA_H_
#define BCL_MATH_CUBIC_SPLINE_VARIABLE_DELTA_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_spline_border_type.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CubicSplineVariableDelta
    //! @brief This class is a cubic spline function for irregulary-spaced points
    //! See Cubic Spline Interpolation: A Review Geroge Wolbert Septemeber 1988
    //! This uses the not-a-knot boundary condition for the endpoints
    //!
    //! @see @link example_math_cubic_spline_variable_delta.cpp @endlink
    //! @author putnamdk
    //! @date 6/12/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CubicSplineVariableDelta :
      public FunctionInterfaceSerializable< double, double>
    {
    private:

    //////////
    // data //
    //////////

      SplineBorderType         m_Border;           //!< controls the behavior of the function at the borders
      linal::Vector< double>   m_X;                //!< list of Abscissa points
      linal::Vector< double>   m_Y;                //!< list of Ordinate values at given x point
      linal::Vector< double>   m_SecondDerivative; //!< second derivative at each knot in the spline

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

       //! @brief copy constructor
      CubicSplineVariableDelta *Clone() const
      {
        return new CubicSplineVariableDelta( *this);
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

      //! @brief access to the values
      //! @return the x function values at the support points of the spline
      const linal::Vector< double> &GetXValues() const
      {
        return m_X;
      }

      //! @brief access to the values
      //! @return the y function values at the support points of the spline
      const linal::Vector< double> &GetYValues() const
      {
        return m_Y;
      }

      //! @brief access to the values
      //! @return the derivative function values at the support points of the spline
      const linal::Vector< double> &GetSecondDerivativeValues() const
      {
        return m_SecondDerivative;
      }

      //! @brief return derivative at certain ARGUMENT
      //! @param ARGUMENT x value
      //! @return derivative at ARGUMENT
      const double dF( const double &ARGUMENT) const;

      //! @brief return value and derivative at ARGUMENT
      //! @param ARGUMENT x value
      //! @return value and derivative at ARGUMENT
      storage::Pair< double, double> FdF( const double &ARGUMENT) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief train CubicSplineVariableDelta
      //! @param X abscissa Values of dataset
      //! @param Y ordinate values of dataset
      //! @param BORDER the variable to control the boundary condition.  The four conditions are:
      //  1) e_Natural
      //  2) e_Periodic
      //  3) e_FirstDer
      //  4) e_NotAKnot
      //! @return trained spline
      CubicSplineVariableDelta &Train
      (
        const SplineBorderType BORDER,
        const linal::VectorConstInterface< double> &X,
        const linal::VectorConstInterface< double> &Y,
        const storage::Pair< double, double> &FIRSTBE = ( storage::Pair< double, double>( 0.0, 0.0))
      );

      //! @brief train CubicSpline for Constant Delta.  This overloads the train function with the same input from the
      //         original implementation of cublic spline
      //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
      //! @param START the start of the interval the spline is defined on
      //! @param DELTA the distance between two support points of the spline
      //! @param RESULTS the function values at the support points of the spline
      //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
      //! @return trained spline
      CubicSplineVariableDelta &Train
      (
        const SplineBorderType BORDER,
        const double START,
        const double DELTA,
        const linal::Vector< double> &RESULTS,
        const storage::Pair< double, double> &FIRSTBE = ( storage::Pair< double, double>( 0.0, 0.0))
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      double operator()( const double &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read CubicSplineVariableDelta from std::istream
      //! @param ISTREAM the stream to read the spline from
      //! @return istream after reading
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write CubicSplineVariableDelta to std::ostream
      //! @param OSTREAM the stream to write the spline to
      //! @param INDENT indentation of the spline
      //! @return ostream after writing
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    ///////////////////////
    // helper  functions //
    ///////////////////////

      //! @brief use m_x and m_y and initinal YD vectors to set m_SecondDerivative
      void GenerateSecondDerivative();

      //! @brief calculate derivative between two cells
      //! @param INDEX_LEFT index of left grid point
      //! @param INDEX_RIGHT index of right grid point
      //! @param DXP relative distance from left grid point, must be element [0, 1]
      //! @return derivative depending on relative distance DXP
      double Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const;

    }; // class CubicSplineVariableDelta

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_CUBIC_SPLINE_VARIABLE_DELTA_H_
