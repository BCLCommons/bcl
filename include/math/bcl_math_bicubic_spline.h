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

#ifndef BCL_MATH_BICUBIC_SPLINE_H_
#define BCL_MATH_BICUBIC_SPLINE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_binary_function_interface_serializable.h"
#include "bcl_math_spline_border_type.h"
#include "linal/bcl_linal_matrix.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BicubicSpline
    //! @brief This class is a bicubic spline interpolation function which is used to approximate
    //! data with a function with continous derivative. Numerical recipes in C++ pp. 116-119 contains the theory behind
    //! this class. The code in this class is original and no where near as buggy as the numerical recipes version
    //!
    //! @see @link example_math_bicubic_spline.cpp @endlink
    //! @author mueller
    //! @date 11/04/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BicubicSpline :
      public BinaryFunctionInterfaceSerializable< double, double, double>
    {
    private:

    //////////
    // data //
    //////////

      SplineBorderType m_Border[ 2]; //!< controls the behavior at x/y_0 and x/y_dim-1
      double           m_Start[  2];
      double           m_Delta[  2]; //!< gives the arguments as a sequence of equidistant points
      double           m_DeltaNorm[  2]; //!< gives the arguments as a sequence of equidistant points

      linal::Matrix< double>  m_Values;  //!< f(x)
      linal::Matrix< double>  m_Dsecox;  //!< second order derivatives for x
      linal::Matrix< double>  m_Dsecoy;  //!< second order derivatives for y
      linal::Matrix< double>  m_Dsecoxy; //!< second order derivatives for x and y

      storage::Pair< double, double> m_FirstBe[ 2]; //!< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for SplineBorderType FIRSTDER

      bool             m_LinCont[2]; //!< if the argument x is outside the range decide if the spline should be continued linearly

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct generic BicubicSpline
      BicubicSpline()
      {
      }

      //! copy constructor
      BicubicSpline *Clone() const
      {
        return new BicubicSpline( *this);
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

      //! get the values of the spline
      const linal::Matrix< double> &GetValues() const
      {
        return m_Values;
      }

      //! get the start for x
      double GetStartX() const
      {
        return m_Start[ 0];
      }

      //! get the start of y
      double GetStartY() const
      {
        return m_Start[ 1];
      }

      //! get the delta along x
      double GetDeltaX() const
      {
        return m_Delta[ 0];
      }

      //! get the delta along y
      double GetDeltaY() const
      {
        return m_Delta[ 1];
      }

      //! get the second order derivatives of the spline
      linal::Matrix< double> const &GetDsecox() const
      {
        return m_Dsecox;
      }

      linal::Matrix< double> const &GetDsecoy() const
      {
        return m_Dsecoy;
      }

      linal::Matrix< double> const &GetDsecoxy() const
      {
        return m_Dsecoxy;
      }

    //////////////
    // operator //
    //////////////

      //! @brief operator
      //! @param X argument x to spline
      //! @param Y argument y to spline
      //! @return the interpolate value at (x,y)
      double operator()( const double &X, const double &Y) const;

    ////////////////
    // operations //
    ////////////////

      //! return value at (x, y)
      double F( const linal::Vector< double> &ARGUMENTS) const;

      //! return partial derivative at (x, y) for x
      double dFdx( const linal::Vector< double> &ARGUMENTS) const;

      //! return partial derivative at (x, y) for y
      double dFdy( const linal::Vector< double> &ARGUMENTS) const;

      //! return value and derivative at (x, y)
      storage::Pair< double, linal::Vector< double> > FdF( const linal::Vector< double> &ARGUMENTS) const;

      //! train BicubicSpline
      void Train
      (
        const SplineBorderType BORDER[ 2],
        const double START[ 2],
        const double DELTA[ 2],
        const linal::Matrix< double> &RESULTS,
        const bool LINCONT[ 2],
        const storage::Pair< double, double> FIRSTBE[ 2]
      );

      //! train BicubicSpline with smoothing of the data
      void TrainWithPreprocessing
      (
        const SplineBorderType BORDER[ 2],
        const double START[ 2],
        const double DELTA[ 2],
        const linal::Matrix< double> &RESULTS,
        const bool LINCONT[ 2],
        const storage::Pair< double, double> FIRSTBE[ 2],
        const double PREPROC
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read BicubicSpline from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write BicubicSpline to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class BicubicSpline

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_BICUBIC_SPLINE_H_
