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

#ifndef BCL_MATH_QUATERNION_H_
#define BCL_MATH_QUATERNION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_rotation_matrix_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    BCL_API Quaternion Conjugate( const Quaternion &QUATERNION);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Quaternion
    //! @brief Quaternions are a four-dimensional generalization of complex numbers.
    //! @details They allow a simple representation of rotations, that is in many (not all!!!) cases faster than matrix-rotations.
    //! Main advantage is that there is (nearly) no numerical drift, that appears in matrix-rotations, due to the 9dim
    //! representation of a 3dim problem.
    //! Also, the small drift can be easily corrected through very simple normalization.
    //! It is a standalone class, as such not dependent on other classes, just the translation of quaternions into
    //! vectors or matrices require according classes to be included.
    //!
    //! @see @link example_math_quaternion.cpp @endlink
    //! @author staritrd
    //! @date 27. September 2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! The Quaternion Class
    class BCL_API Quaternion :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      double *m_Values; //!< data

    public:

      static const size_t s_Size = 4; //!< Number of Elements

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! The default constructor
      Quaternion() :
        m_Values( new double[ s_Size])
      {
        std::fill( m_Values, m_Values + s_Size, double( 0));
      }

      //! The copy constructor
      Quaternion( const Quaternion &QUATERNION) :
        m_Values( new double[ s_Size])
      {
        std::copy( QUATERNION.m_Values, QUATERNION.m_Values + s_Size, m_Values);
      }

      //! The value based constructor
      Quaternion( const double &R, const double &i_X, const double &i_Y, const double &i_Z) :
        m_Values( new double[ s_Size])
      {
        m_Values[ 0] =   R;
        m_Values[ 1] = i_X;
        m_Values[ 2] = i_Y;
        m_Values[ 3] = i_Z;
      }

      //! constructs quaternion from rotation axes (linal::Vector3D) and angle
      Quaternion( const linal::Vector3D &AXIS, const double ANGLE)
      {
        linal::Vector3D axis( AXIS);
        axis.Normalize();
        Equal( cos( 0.5 * ANGLE), sin( 0.5 * ANGLE) * axis);
      }

      //! copy constructor
      Quaternion *Clone() const
      {
        return new Quaternion( *this);
      }

      //! The destructor
      ~Quaternion()
      {
        delete [] m_Values;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator q1 = q2
      Quaternion &operator =( const Quaternion &QUATERNION);

      //! operator *= multiplies each value in quaternion with same factor
      Quaternion &operator *=( const double T);

      //! operator /= divides each value in quaternion with same factor
      Quaternion &operator /=( const double T);

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    //////////////
    // operator //
    //////////////

      //! access of the data
      double &operator[]( const size_t INDEX)
      {
        BCL_Assert( INDEX < s_Size, "access INDEX > 3");
        return m_Values[ INDEX];
      }

      //! read-only access of the data
      const double &operator[]( const size_t INDEX) const
      {
        BCL_Assert( INDEX < s_Size, "access INDEX > 3");
        return m_Values[ INDEX];
      }

    ////////////////
    // operations //
    ////////////////

      // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
      //! set random rotation
      Quaternion &SetRand();

      //! resets the values of the quaternion to the given four values
      Quaternion &Equal( const double R, const double X, const double Y, const double Z);

      //! resets the real part of the quaternion with single value and the imaginary part with a 3D vector
      Quaternion &Equal( const double R, const linal::Vector3D &VECTOR3D);

      //! normalizes the quaternion and returns it
      Quaternion &Normalize()
      {
        return operator /=( Norm());
      }

      //! returns the normalized quaternion, while the quaternion itself remains unchanged
      Quaternion Normalized() const
      {
        return Quaternion( *this).operator /=( Norm());
      }

      //! returns the norm of the quaternion
      double Norm() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! writes a linal::Vector3D into quaternion
      Quaternion &VectorToQuaternion( const linal::Vector3D &);

    protected:

      std::istream &Read( std::istream &ISTREAM);

      //! output data, default cout
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

    BCL_API Quaternion operator *( const Quaternion &QUAT_1, const Quaternion &QUAT_2);

    BCL_API Quaternion operator *( const linal::Vector3D &VECTOR3D, const Quaternion &QUATERNION);

    inline Quaternion operator *( const Quaternion &QUATERNION, const linal::Vector3D &VECTOR3D)
    {
      return VECTOR3D * QUATERNION;
    }

    inline Quaternion operator *( const Quaternion &QUATERNION, const double &T)
    {
      return Quaternion( QUATERNION).operator *=( T);
    }

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_QUATERNION_H_
