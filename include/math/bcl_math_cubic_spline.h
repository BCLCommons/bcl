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

#ifndef BCL_MATH_CUBIC_SPLINE_H_
#define BCL_MATH_CUBIC_SPLINE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_spline_border_type.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CubicSpline
    //! @brief This class is a cubic spline interpolation function which is used to approximate data with a
    //! function with continuous derivative. Numerical recipes in C++ pp. 116-119 contains the theory behind this class
    //! The code in this class is original and no where near as buggy as the numerical recipes versions
    //!
    //! @see @link example_math_cubic_spline.cpp @endlink
    //! @author mueller
    //! @date 11/04/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CubicSpline :
      public FunctionInterfaceSerializable< double, double>
    {
    private:

    //////////
    // data //
    //////////

      SplineBorderType       m_Border;     //!< controls the behavior at x_0 and x_dim-1
      double         m_Start, m_Delta;     //!< gives the arguments as a sequence of equidistant points
      linal::Vector< double> m_Values;     //!< f(x)
      linal::Vector< double> m_Dsecox;     //!< second order derivatives

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct generic CubicSpline
      CubicSpline()
      {
      }

      //! @brief copy constructor
      CubicSpline *Clone() const
      {
        return new CubicSpline( *this);
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

      //! @brief get the second order derivatives of the spline
      //! @return the second order derivatives at the support points of the spline
      linal::Vector< double> const &GetDsecox() const
      {
        return m_Dsecox;
      }

      //! @brief access to the start value
      //! @return the start of the interval the spline is defined on
      double GetStart() const
      {
        return m_Start;
      }

      //! @brief access to the delta value
      //! @return the distance between two support points of the spline
      double GetDelta() const
      {
        return m_Delta;
      }

      //! @brief access to the values
      //! @return the function values at the support points of the spline
      const linal::Vector< double> &GetValues() const
      {
        return m_Values;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief return derivative at ARGUMENT
      //! @param ARGUMENT x value
      //! @return derivative at ARGUMENT
      double dF( const double &ARGUMENT) const;

      //! @brief return value and derivative at ARGUMENT
      //! @param ARGUMENT x value
      //! @return value and derivative at ARGUMENT
      storage::Pair< double, double> FdF( const double &ARGUMENT) const;

      //! @brief train CubicSpline
      //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
      //! @param START the start of the interval the spline is defined on
      //! @param DELTA the distance between two support points of the spline
      //! @param RESULTS the function values at the support points of the spline
      //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
      //! @return trained spline
      CubicSpline &Train
      (
        const SplineBorderType BORDER,
        const double START,
        const double DELTA,
        const linal::VectorConstInterface< double> &RESULTS,
        const storage::Pair< double, double> &FIRSTBE
      );

      //! @brief train CubicSpline with preprocessing of data to smooth the function values
      //! @param BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
      //! @param START the start of the interval the spline is defined on
      //! @param DELTA the distance between two support points of the spline
      //! @param RESULTS the function values at the support points of the spline
      //! @param FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
      //! @param PREPROC determines the degree of smoothing, should be between 0 and 1 to give the part/influence
      //! of the actual value to the pre processed value
      //! @return trained spline
      CubicSpline &TrainWithPreprocessing
      (
        const SplineBorderType BORDER,
        const double START,
        const double DELTA,
        const linal::VectorConstInterface< double> &RESULTS,
        const storage::Pair< double, double> &FIRSTBE, const double PREPOC
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

      //! @brief read CubicSpline from std::istream
      //! @param ISTREAM the stream to read the spline from
      //! @return istream after reading
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write CubicSpline to std::ostream
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
      //! @param INDEX_RIGHT index of right grid point
      //! @param DXP relative distance from left grid point, must be element [0, 1]
      //! @return derivative depending on relative distance DXP
      double Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const;

      //! @brief calculate function between two cells
      //! @param INDEX_LEFT index of left grid point
      //! @param INDEX_RIGHT index of right grid point
      //! @param DXP relative distance from left grid point, must be element [0, 1]
      //! @return function depending on relative distance DXP
      double Function( const int INDEX_LEFT, const int INDEX_RIGHT, const double DXP) const;

    }; // class CubicSpline

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_CUBIC_SPLINE_H_
