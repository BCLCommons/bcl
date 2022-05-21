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

#ifndef BCL_MATH_ROTATION_MATRIX_3D_H_
#define BCL_MATH_ROTATION_MATRIX_3D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix3x3.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotationMatrix3D
    //! @brief RotationMatrix3D for rotating a Vector in 3D space
    //! @details a 3x3 matrix contains all the euler angle derived coefficient, to rotate a coordinate in space around
    //!          the origin by a simple matrix vector multiplication
    //!
    //! @see @link example_math_rotation_matrix_3d.cpp @endlink
    //! @author karakam, woetzen
    //! @date Nov 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotationMatrix3D :
      public util::ObjectInterface
    {
    /////////////
    // friends //
    /////////////

      friend class linal::Vector3D;
      friend class TransformationMatrix3D;

    //////////
    // data //
    //////////

    private:

      linal::Matrix3x3< double> m_RotationMatrix3D; //!< RotationMatrix3D

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct empty RotationMatrix3D
      RotationMatrix3D();

      //! construct RotationMatrix3D from RotationAxis and RotationAngle
      RotationMatrix3D( const coord::Axis &AXIS, const double &ANGLE);

      //! construct RotationMatrix3D from any linal::Vector3D around which shall be rotated and the Angle.
      RotationMatrix3D( const linal::Vector3D &AXIS, const double ANGLE);

      //! construct RotationMatrix3D from three RotationAxis and three RotationAngles
      RotationMatrix3D( const coord::Axis AXIS[ 3], const double ANGLE[ 3]);

      //! construct RotationMatrix3D from three Euler angles
      RotationMatrix3D( const double PHI, const double PSI, const double THETA);

      //! construct RotationMatrix3D from a 3x3 matrix
      RotationMatrix3D( const linal::Matrix3x3< double> &ROTATION);

      //! copy constructor
      RotationMatrix3D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! return the matrix
      const linal::Matrix3x3< double> &GetMatrix() const;

      //! @brief return the X Y or Z axis in respect to their reference coordinate system x, y and z (euclidean unit vectors)
      //! @param AXIS axis that will be returned
      //! @return axis in coordinates of reference coordinate system
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! @brief return the axis of rotation; result is a unit vector
      //! @return axis of rotation; if this rotation matrix is applied to it, the vector will be unchanged
      linal::Vector3D GetRotationAxis() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return the euler angles for that rotation matrix
      //! http://www.gamedev.net/reference/articles/article1691.asp
      //! @return linal::Vector3D with 1st to 3rd element being alpha, beta, gamma - rotations around x y and z
      linal::Vector3D EulerAnglesXYZ() const;

      //! @brief return the Euler angles according to z-x'-z'' convention
      //! @return Euler angles according to z-x'-z'' convention
      storage::VectorND< 3, double> EulerAnglesZXZ() const;

      //http://journals.iucr.org/j/issues/2002/05/00/vi0166/index.html
      //! calculate effective rotation angle
      double EffectiveRotationAngle() const;

      // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
      //! set random rotation using a MAX_ROTATIONANGLE
      RotationMatrix3D &SetRand( const double MAX_ROTATIONANGLE = g_Pi);

    ///////////////
    // operators //
    ///////////////

      //! operator *= RotationMatrix3D
      RotationMatrix3D &operator *=( const RotationMatrix3D &MATRIX);

    //////////////////////
    // input and output //
    //////////////////////

      //! write RotationMatrix3D to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const;

    protected:

      //! write RotationMatrix3D to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read RotationMatrix3D from io::IFStream
      std::istream &Read( std::istream &ISTREAM);

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief construct from alpha, beta and gamma angles according to z-x-z convention
      //! @param ALPHA rotation angle around LCS z-axis
      //! @param BETA rotation angle around LCS x'-axis
      //! @param GAMMA rotation angle around LCS z''-axis
      static RotationMatrix3D CreateZXZ( const double ALPHA, const double BETA, const double GAMMA);

    private:

    ////////////////
    // operations //
    ////////////////

      //! private helper function to set rotation matrix data from RotationAxis and RotationAngles
      void Init( const coord::Axis &AXIS, const double &ANGLE);

    }; //class RotationMatrix3D

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_ROTATION_MATRIX_3D_H_
