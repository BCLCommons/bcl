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

#ifndef BCL_MATH_TRICUBIC_SPLINE_H_
#define BCL_MATH_TRICUBIC_SPLINE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_spline_border_type.h"
#include "bcl_math_tensor.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TricubicSpline
    //! @brief This class is a tricubic spline interpolation function which is used to approximate data with a
    //! function with continuous derivative. Numerical recipes in C++ pp. 116-119 contains the theory behind
    //! this class. The code in this class is original
    //!
    //! @see @link example_math_tricubic_spline.cpp @endlink
    //! @author mueller
    //! @date 11/04/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TricubicSpline :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      SplineBorderType  m_Border[3];   //!< controls the behaviour at (x|y|z)_0 and (x|y|z)_dim-1, can be (NATURAL, PERIODIC, FIRSTDER)

      double m_Start[3], m_Delta[3];   //!< gives the arguments as a sequence of equidistant points

      Tensor< double> m_Values;        //!< the function values f(x) which will be splined
      Tensor< double> m_Dsecox;        //!< second order derivatives for x
      Tensor< double> m_Dsecoy;        //!< second order derivatives for y
      Tensor< double> m_Dsecoxy;       //!< second order derivatives for x and y
      Tensor< double> m_Dsecoz;        //!< second order derivatives for z
      Tensor< double> m_Dsecoxz;       //!< second order derivatives for x and z
      Tensor< double> m_Dsecoyz;       //!< second order derivatives for y and z
      Tensor< double> m_Dsecoxyz;      //!< second order derivatives for x, y and z

      storage::Pair< double, double> m_FirstBe[ 3]; //!< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for SplineBorderType FIRSTDER

      bool m_LinCont[3];    //!< if the argument is outside the range decide if the spline should be continued linearly

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! copy constructor
      TricubicSpline *Clone() const
      {
        return new TricubicSpline( *this);
      }

    //////////////
    // operator //
    //////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! return value at (x, y, z)
      virtual double F( const double x, const double y, const double z) const;

      //! return partial derivative at (x, y, z) for x
      virtual double dFdx( const double x, const double y, const double z) const;

      //! return partial derivative at (x, y, z) for y
      virtual double dFdy( const double x, const double y, const double z) const;

      //! return partial derivative at (x, y, z) for y
      virtual double dFdz( const double x, const double y, const double z) const;

      //! return value and derivative at (x, y, z)
      virtual storage::Pair< double, linal::Vector< double> > FdF( const double x, const double y, const double z) const;

      //! train TricubicSpline
      virtual void Train
      (
        const SplineBorderType BORDER[ 3], const double START[ 3], const double DELTA[ 3],
        const Tensor< double> &RESULTS, const bool LINCONT[ 3], const storage::Pair< double, double> FIRSTBE[ 3]
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read TricubicSpline from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write TricubicSpline to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_TRICUBIC_SPLINE_H_
