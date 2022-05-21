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

#ifndef BCL_MATH_CUBIC_SPLINE_DAMPED_H_
#define BCL_MATH_CUBIC_SPLINE_DAMPED_H_

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
    //! @class CubicSplineDamped
    //! @brief This class is a cubic spline function for irregulary-spaced points
    //! The normal cubic spline with variable delta is subject to chaotic behavior with irregularly spaced points, for
    //! which may exhibit non-local effects and introduce novel extrema (minima or maxima) that were not in the in the source
    //! data. This class yields a spline that interpolates smoothly (at 1st derivative) between points, with a guarantee
    //! of not introducing extrema that are not in the data, and will not shift extrema that are in the data
    //! That is, the actual maxima in the data will be the actual maxima in the trained spline.
    //!
    //! @see @link example_math_cubic_spline_damped.cpp @endlink
    //! @author mendenjl
    //! @date Feb 11, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CubicSplineDamped :
      public FunctionInterfaceSerializable< double, double>
    {
    private:

    //////////
    // data //
    //////////

      linal::Vector< double>   m_X;                //!< list of Abscissa points
      linal::Vector< double>   m_Y;                //!< list of Ordinate values at given x point
      linal::Vector< double>   m_dY;               //!< derivative at each y-point
      linal::Vector< double>   m_A;                //!< a (cubic) term in the cubic polynomial
      linal::Vector< double>   m_B;                //!< b (quadratic) term in the cubic polynomial
      double                   m_Delta;            //!< Delta - average delta if points are not equidistant
      bool                     m_ConstantDelta;    //!< Whether delta is constant

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
      CubicSplineDamped *Clone() const
      {
        return new CubicSplineDamped( *this);
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

      //! @brief access to delta
      //! @return the delta value (average delta in case the deltas are not constant)
      const double &GetDelta() const
      {
        return m_Delta;
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

      //! @brief train CubicSpline, requiring a maximally smooth, piecewise monotonic solution
      //! @param X abscissa of dataset
      //! @param Y ordinate values of dataset
      //! @param DY_START first order derivative at X(0). if nan/undefined, will be the same as first order derivative at X(1)
      //! @param DY_END   first order derivative at X_last. if nan/undefined, will be the same as first order derivative at X(last-1)
      //! @return *this
      //! Reference: M. Steffen. "A simple method for monotonic interpolation in one dimension" Astronomy and
      //! Astrophysics. 1990
      //! Note that it is not required that the Y be monotonic. All this function does is ensure that each cubic piece
      //! of the function is monotonic over the given interval and keeps the y' continuous. y'' will not, in general, be
      //! continuous if the spline is trained with this function
      CubicSplineDamped &Train
      (
        const linal::VectorConstInterface< double> &X,
        const linal::VectorConstInterface< double> &Y,
        const double &DY_START = std::numeric_limits< double>::quiet_NaN(),
        const double &DY_END = std::numeric_limits< double>::quiet_NaN()
      );

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
      CubicSplineDamped &Train
      (
        const double &START,
        const double &DELTA,
        const linal::VectorConstInterface< double> &Y,
        const double &DY_START = std::numeric_limits< double>::quiet_NaN(),
        const double &DY_END = std::numeric_limits< double>::quiet_NaN()
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

      //! @brief read CubicSplineDamped from std::istream
      //! @param ISTREAM the stream to read the spline from
      //! @return istream after reading
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write CubicSplineDamped to std::ostream
      //! @param OSTREAM the stream to write the spline to
      //! @param INDENT indentation of the spline
      //! @return ostream after writing
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    ///////////////////////
    // helper  functions //
    ///////////////////////

      //! @brief calculate derivative between two cells
      //! @param INDEX_LEFT index of left grid point
      //! @param DXP relative distance from left grid point, must be element [0, 1]
      //! @return derivative depending on relative distance DXP
      double Derivative( const int INDEX_LEFT, const double DXP) const;

    }; // class CubicSplineDamped

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_CUBIC_SPLINE_DAMPED_H_
