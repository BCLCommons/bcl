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

#ifndef BCL_MATH_TRANSFORMATION_MATRIX_3D_H_
#define BCL_MATH_TRANSFORMATION_MATRIX_3D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_line_segment_3d.h"
#include "linal/bcl_linal_matrix.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransformationMatrix3D
    //! @brief This class is for 3-dimensional TransformationMatrix3D
    //!
    //! @see @link example_math_transformation_matrix_3d.cpp @endlink
    //! @author meilerj
    //! @date 21.08.2004
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TransformationMatrix3D :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

      friend class linal::Vector3D;

    private:

    //////////
    // data //
    //////////

      linal::Matrix< double> m_TransformationMatrix3D; //!< TransformationMatrix3D

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct TransformationMatrix3D and set to unitmatrix (empty constructor)
      TransformationMatrix3D();

      //! copy constructor
      TransformationMatrix3D( const TransformationMatrix3D &MATRIX);

      //! @brief move constructor
      TransformationMatrix3D( TransformationMatrix3D && MATRIX);

      //! construct from translation
      //! cannot be implemented here, because of cross dependency of TransformationMatrix3D and linal::Vector3D
      //! is indirectly implemented in linal::Vector3D as operator:
      //! linal::Vector3D::operator TransformationMatrix3D() const
      //TransformationMatrix3D( const linal::Vector3D &TRANSLATION);

      //! construct TransformationMatrix3D from 4linal::Vector3D ( x, y and z axis and the origin)
      TransformationMatrix3D
      (
        const linal::Vector3D &XAXIS,
        const linal::Vector3D &YAXIS,
        const linal::Vector3D &ZAXIS,
        const linal::Vector3D &ORIGIN
      );

      //! construct TransformationMatrix3D from MatrixND
      TransformationMatrix3D( const RotationMatrix3D &MATRIX);

      //! construct TransformationMatrix3D from Matrix< double> 4x4
      explicit TransformationMatrix3D( const linal::Matrix< double> &MATRIX);

      //! @brief construct TransformationMatrix3D from translation
      //! @param TRANSLATION the translation to perform
      explicit TransformationMatrix3D( const linal::Vector3D &TRANSLATION);

      //! construct TransformationMatrix3D from translation
      TransformationMatrix3D( const double &X, const double &Y, const double &Z);

      //! construct TransformationMatrix3D from RotationAxis and RotationAngle
      TransformationMatrix3D( const coord::Axis &AXIS, const double &ANGLE);

      //! construct TransformationMatrix3D from three RotationAxis and three RotationAngles
      TransformationMatrix3D( const coord::Axis AXES[ 3], const double ANGLE[ 3]);

      //! construct from vector with six components - three rotations and three translations
      //! @param VECTOR elements 1-3 rotation x, y, z; elements 4-6 translations x, y, z
      explicit TransformationMatrix3D( const linal::Vector< double> &VECTOR);

      //! @brief constructor from a definition state
      //! @param DEFINITION_STATE defined or undefined
      explicit TransformationMatrix3D( const util::UndefinedObject DEFINITION_STATE);

      //! @brief constructor taking the member variable parameters
      //! @param ROTATION_AXIS the axis around which the rotation will take place
      //! @param ROTATION_ORIGIN the origin point around which the rotation will occur
      //! @param ROTATION the amount of rotation
      TransformationMatrix3D
      (
        const coord::LineSegment3D &ROTATION_AXIS, const linal::Vector3D &ROTATION_ORIGIN, const double ROTATION
      );

      //! @brief constructor taking two lines, computes the transformation A->B
      //! @param A the starting line
      //! @param B the ending line
      TransformationMatrix3D( const coord::LineSegment3D &A, const coord::LineSegment3D &B);

      //! copy constructor
      TransformationMatrix3D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! access and modify the matrix
      void SetMatrix( const linal::Matrix< double> &MATRIX);

      //! read only access
      const linal::Matrix< double> &GetMatrix() const;

    ///////////////
    // operators //
    ///////////////

      //! operator = TransformationsMatrix3D
      virtual TransformationMatrix3D &operator =( const TransformationMatrix3D &MATRIX);

      //! operator = TransformationsMatrix3D
      TransformationMatrix3D &operator =( TransformationMatrix3D && MATRIX);

      //! operator *= TransformationMatrix3D
      TransformationMatrix3D &operator *=( const TransformationMatrix3D &MATRIX);

      //! operator += TRANSFORMATIONMATRIX3D
      TransformationMatrix3D &operator +=( const TransformationMatrix3D &MATRIX);

      //! operator /= SCALAR
      TransformationMatrix3D &operator /=( const double &SCALAR);

      //! operator *= SCALAR
      TransformationMatrix3D &operator *=( const double &SCALAR);

      //! apply translation
      TransformationMatrix3D &operator()( const TransformationMatrix3D &MATRIX);

      //!  apply translation
      TransformationMatrix3D &operator()( const double &X, const double &Y, const double &Z);

      //!  apply translation
      TransformationMatrix3D &operator()( const linal::Vector3D &TRANSLATION);

      //! apply rotation
      TransformationMatrix3D &operator()( const RotationMatrix3D &ROTATION);

      //!  apply RotationAxis and RotationAngle
      TransformationMatrix3D &operator()( const coord::Axis &AXIS, const double &ANGLE);

      //!  apply RotationAxis and three RotationAngles
      TransformationMatrix3D &operator()( const coord::Axis AXES[ 3], const double ANGLE[ 3]);

    ////////////////
    // operations //
    ////////////////

      //! Set Unit Matrix
      TransformationMatrix3D &SetUnit();

      //! Invert the transformation matrix
      TransformationMatrix3D &Invert();

      //! return rotationmatrix
      RotationMatrix3D GetRotation() const;

      //! return Translation
      linal::Vector3D GetTranslation() const;

      //! returns the axis information (see class Body)
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! returns the origin (see class Body)
      linal::Vector3D GetOrigin() const;

      //! @brief set the translation for this transformation matrix
      //! @param TRANSLATION translation vector
      void SetTranslation( const linal::Vector3D &TRANSLATION);

      //! @brief set the rotation matrix for this class
      //! @param ROTATION the rotation matrix to use
      void SetRotation( const RotationMatrix3D &ROTATION);

      //! modification of the axis (see class Body)
      void SetAxis( const coord::Axis &AXIS, const linal::Vector3D &NEW_AXIS);

      //! modification of the axis (see class Body)
      void SetOrigin( const linal::Vector3D &ORIGIN);

      //! @brief returns whether this transformation is defined or not by looking at the scaling factor
      //! @return whether this transformation is defined or not
      bool IsDefined() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! write TransformationMatrix3D to std::ostream using the given util::Format
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const;

    protected:

      //! write TransformationMatrix3D to std::ostream using the given util::Format
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read TransformationMatrix3D from io::IFStream
      virtual std::istream &Read( std::istream &ISTREAM);

      //! Set Unitary matrix, assuming that the matrix is currently composed of zeros
      void SetUnitFromZero();

    }; // class TransformationMatrix3D

    //! boolean operator TRANSFORMATIONMATRIX3D_A == TRANSFORMATIONMATRIX3D_B
    BCL_API bool operator ==
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    );

    //! boolean operator TRANSFORMATIONMATRIX3D_A != TRANSFORMATIONMATRIX3D_B
    BCL_API bool operator !=
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    );

    //! return Inverse of TransformationMatrix3D
    BCL_API TransformationMatrix3D Inverse( const TransformationMatrix3D &MATRIX);

    //! check if two transformation matrices are similar within a given rotational and translational tolerance
    BCL_API bool SimilarWithinTolerance
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B,
      const double TRANSLATION_TOLERANCE,
      const double ROTATION_TOLERANCE_RAD
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_TRANSFORMATION_MATRIX_3D_H_
