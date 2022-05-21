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
#include "math/bcl_math_transformation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> TransformationMatrix3D::s_Instance
    (
      GetObjectInstances().AddInstance( new TransformationMatrix3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct TransformationMatrix3D and set to unitmatrix (empty constructor)
    TransformationMatrix3D::TransformationMatrix3D() :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
    }

    //! copy constructor
    TransformationMatrix3D::TransformationMatrix3D( const TransformationMatrix3D &MATRIX) :
      m_TransformationMatrix3D( MATRIX.m_TransformationMatrix3D)
    {
    }

    //! @brief move constructor
    TransformationMatrix3D::TransformationMatrix3D( TransformationMatrix3D && MATRIX) :
      m_TransformationMatrix3D( std::move( MATRIX.m_TransformationMatrix3D))
    {
    }

    //! construct TransformationMatrix3D from 4linal::Vector3D ( x, y and z axis and the origin)
    TransformationMatrix3D::TransformationMatrix3D
    (
      const linal::Vector3D &XAXIS,
      const linal::Vector3D &YAXIS,
      const linal::Vector3D &ZAXIS,
      const linal::Vector3D &ORIGIN
    ) :
      m_TransformationMatrix3D( 4, 4)
    {
      double *trans( m_TransformationMatrix3D.Begin());
      for( const double *x_ptr = XAXIS.Begin(); x_ptr != XAXIS.End(); ++x_ptr)
      {
        ( *( trans++)) = ( *x_ptr);
      }

      trans++;
      for( const double *y_ptr = YAXIS.Begin(); y_ptr != YAXIS.End(); ++y_ptr)
      {
        ( *( trans++)) = ( *y_ptr);
      }

      trans++;
      for( const double *z_ptr = ZAXIS.Begin(); z_ptr != ZAXIS.End(); ++z_ptr)
      {
        ( *( trans++)) = ( *z_ptr);
      }

      trans++;
      for( const double *o_ptr = ORIGIN.Begin(); o_ptr != ORIGIN.End(); ++o_ptr)
      {
        ( *( trans++)) = ( *o_ptr);
      }

      //corner is 1
      ( *trans) = 1;
    }

    //! construct TransformationMatrix3D from MatrixND
    TransformationMatrix3D::TransformationMatrix3D( const RotationMatrix3D &MATRIX) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      SetRotation( MATRIX);
    }

    //! construct TransformationMatrix3D from Matrix< double> 4x4
    TransformationMatrix3D::TransformationMatrix3D( const linal::Matrix< double> &MATRIX) :
      m_TransformationMatrix3D( MATRIX)
    {
      BCL_Assert( MATRIX.IsSquare() && MATRIX.GetNumberRows() == size_t( 4), "Argument has to be 4x4 matrix");
    }

    //! @brief construct TransformationMatrix3D from translation
    //! @param TRANSLATION the translation to perform
    TransformationMatrix3D::TransformationMatrix3D( const linal::Vector3D &TRANSLATION) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      m_TransformationMatrix3D.ReplaceRow( 3, TRANSLATION);
    }

    //! construct TransformationMatrix3D from translation
    TransformationMatrix3D::TransformationMatrix3D( const double &X, const double &Y, const double &Z) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      m_TransformationMatrix3D( 3, 0) = X;
      m_TransformationMatrix3D( 3, 1) = Y;
      m_TransformationMatrix3D( 3, 2) = Z;
    }

    //! construct TransformationMatrix3D from coord::Axis and RotationAngle
    TransformationMatrix3D::TransformationMatrix3D( const coord::Axis &AXIS, const double &ANGLE) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      RotationMatrix3D a( AXIS, ANGLE);
      SetRotation( a.m_RotationMatrix3D);
    }

    //! construct TransformationMatrix3D from three coord::Axis and three RotationAngles
    TransformationMatrix3D::TransformationMatrix3D( const coord::Axis AXES[ 3], const double ANGLE[ 3]) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      SetRotation( RotationMatrix3D( AXES, ANGLE));
    }

    //! construct from vector with six components - three rotations and three translations
    //! @param VECTOR elements 1-3 rotation x, y, z; elements 4-6 translations x, y, z
    TransformationMatrix3D::TransformationMatrix3D( const linal::Vector< double> &VECTOR) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      BCL_Assert( VECTOR.GetSize() == 6, "argument vector needs to have 6 elements");
      SetUnitFromZero();
      SetRotation( RotationMatrix3D( &coord::GetAxes().e_X, VECTOR.Begin()));
      std::copy( VECTOR.Begin() + 3, VECTOR.End(), m_TransformationMatrix3D.Begin() + 12);
    }

    //! @brief constructor from a definition state
    //! @param DEFINITION_STATE defined or undefined
    TransformationMatrix3D::TransformationMatrix3D( const util::UndefinedObject DEFINITION_STATE) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnitFromZero();
      // set the scaling factor to undefined
      m_TransformationMatrix3D( 3, 3) = util::GetUndefinedDouble();
    }

    //! @brief constructor taking the member variable parameters
    //! @param ROTATION_AXIS the axis around which the rotation will take place
    //! @param ROTATION_ORIGIN the origin point around which the rotation will occur
    //! @param ROTATION the amount of rotation
    TransformationMatrix3D::TransformationMatrix3D
    (
      const coord::LineSegment3D &ROTATION_AXIS, const linal::Vector3D &ROTATION_ORIGIN, const double ROTATION
    ) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnit();
      m_TransformationMatrix3D.ReplaceRow( 3, -ROTATION_ORIGIN);

      // add the desired rotation according to "m_RotationAxis" and "m_Rotation" to "transform"
      this->operator ()( RotationMatrix3D( ROTATION_AXIS.GetDirection(), ROTATION));

      // add the translation of moving back to the original position
      this->operator ()( ROTATION_ORIGIN);
    }

    //! @brief constructor taking two lines, computes the transformation A->B
    //! @param A the starting line
    //! @param B the ending line
    TransformationMatrix3D::TransformationMatrix3D( const coord::LineSegment3D &A, const coord::LineSegment3D &B) :
      m_TransformationMatrix3D( 4, 4, 0.0)
    {
      SetUnit();
      // get the angle and axis between the lines

      // get projected angle between lines A_A-A_B and B_A-B_B
      const double angle( linal::ProjAngle( A.GetStartPoint(), A.GetEndPoint(), B.GetStartPoint(), B.GetEndPoint()));

      // get rotation axis perpendicular to the plane containing two lines
      linal::Vector3D rotation_axis
      (
        linal::CrossProduct( A.GetDirection(), B.GetDirection())
      );

      // translate the midpoint on the line between the two points to the origin
      const linal::Vector3D mid_point_b( ( B.GetStartPoint() + B.GetEndPoint()) / 2.0);
      ( *this)( -mid_point_b);

      // handle collinear case
      // Choose a random axis that is orthogonal to them
      while( rotation_axis.SquareNorm() == 0.0)
      {
        rotation_axis.SetRandomTranslation( 1.0);
        rotation_axis = linal::CrossProduct( rotation_axis, A.GetDirection());
      }

      RotationMatrix3D rotation( rotation_axis, angle);
      ( *this)( rotation);

      // translate the midpoint back to the midpoint of a
      const linal::Vector3D mid_point_a( ( A.GetStartPoint() + A.GetEndPoint()) / 2.0);
      ( *this)( mid_point_a);
    }

    //! copy constructor
    TransformationMatrix3D *TransformationMatrix3D::Clone() const
    {
      return new TransformationMatrix3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TransformationMatrix3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! access and modify the matrix
    void TransformationMatrix3D::SetMatrix( const linal::Matrix< double> &MATRIX)
    {
      BCL_Assert( MATRIX.IsSquare() && MATRIX.GetNumberRows() == size_t( 4), "Argument has to be 4x4 matrix");
      m_TransformationMatrix3D = MATRIX;
    }

    //! read only access
    const linal::Matrix< double> &TransformationMatrix3D::GetMatrix() const
    {
      return m_TransformationMatrix3D;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator = TransformationsMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator =( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D = MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator = TransformationsMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator =( TransformationMatrix3D && MATRIX)
    {
      m_TransformationMatrix3D = std::move( MATRIX.m_TransformationMatrix3D);
      return *this;
    }

    //! operator *= TransformationMatrix3D
    TransformationMatrix3D &TransformationMatrix3D::operator *=( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D = m_TransformationMatrix3D * MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator += TRANSFORMATIONMATRIX3D
    TransformationMatrix3D &TransformationMatrix3D::operator +=( const TransformationMatrix3D &MATRIX)
    {
      m_TransformationMatrix3D += MATRIX.m_TransformationMatrix3D;
      return *this;
    }

    //! operator /= SCALAR
    TransformationMatrix3D &TransformationMatrix3D::operator /=( const double &SCALAR)
    {
      m_TransformationMatrix3D /= SCALAR;
      return *this;
    }

    //! operator *= SCALAR
    TransformationMatrix3D &TransformationMatrix3D::operator *=( const double &SCALAR)
    {
      m_TransformationMatrix3D *= SCALAR;
      return *this;
    }

    //! apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const TransformationMatrix3D &MATRIX)
    {
      return operator *=( TransformationMatrix3D( MATRIX));
    }

    //!  apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const double &X, const double &Y, const double &Z)
    {
      return operator *=( TransformationMatrix3D( X, Y, Z));
    }

    //!  apply translation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const linal::Vector3D &TRANSLATION)
    {
      return operator *=( TransformationMatrix3D( TRANSLATION));
    }

    //! apply rotation
    TransformationMatrix3D &TransformationMatrix3D::operator()( const RotationMatrix3D &ROTATION)
    {
      return operator *=( TransformationMatrix3D( ROTATION));
    }

    //!  apply coord::Axis and RotationAngle
    TransformationMatrix3D &TransformationMatrix3D::operator()( const coord::Axis &AXIS, const double &ANGLE)
    {
      return operator *=( TransformationMatrix3D( AXIS, ANGLE));
    }

    //!  apply coord::Axis and three RotationAngles
    TransformationMatrix3D &TransformationMatrix3D::operator()( const coord::Axis AXES[ 3], const double ANGLE[ 3])
    {
      return operator *=( TransformationMatrix3D( AXES, ANGLE));
    }

  ////////////////
  // operations //
  ////////////////

    //! Set Unit Matrix
    TransformationMatrix3D &TransformationMatrix3D::SetUnit()
    {
      m_TransformationMatrix3D.SetZero();
      SetUnitFromZero();
      return *this;
    }

    //! Inverse
    TransformationMatrix3D &TransformationMatrix3D::Invert()
    {
      linal::MatrixInversionGaussJordan< double> inverter( m_TransformationMatrix3D);
      BCL_Assert( inverter.IsDefined(), "Could not compute the inverse of " + util::Format()( m_TransformationMatrix3D));
      m_TransformationMatrix3D = inverter.ComputeInverse();
      return *this;
    }

    //! return rotationmatrix
    RotationMatrix3D TransformationMatrix3D::GetRotation() const
    {
      linal::Matrix3x3< double> rotation;
      rotation.ReplaceRow( 0, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 0]));
      rotation.ReplaceRow( 1, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 1]));
      rotation.ReplaceRow( 2, linal::VectorConstReference< double>( 3, m_TransformationMatrix3D[ 2]));
      return RotationMatrix3D( rotation);
    }

    //! return Translation
    linal::Vector3D TransformationMatrix3D::GetTranslation() const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ 3]);
    }

    linal::Vector3D TransformationMatrix3D::GetAxis( const coord::Axis &AXIS) const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ AXIS]);
    }

    linal::Vector3D TransformationMatrix3D::GetOrigin() const
    {
      return linal::Vector3D( m_TransformationMatrix3D[ 3]);
    }

    //! @brief set the translation for this transformation matrix
    //! @param TRANSLATION translation vector
    void TransformationMatrix3D::SetTranslation( const linal::Vector3D &TRANSLATION)
    {
      m_TransformationMatrix3D.ReplaceRow( 3, TRANSLATION);
    }

    //! @brief set the rotation matrix for this class
    //! @param ROTATION the rotation matrix to use
    void TransformationMatrix3D::SetRotation( const RotationMatrix3D &ROTATION)
    {
      m_TransformationMatrix3D.ReplaceRow( 0, ROTATION.GetMatrix().GetRow( 0));
      m_TransformationMatrix3D.ReplaceRow( 1, ROTATION.GetMatrix().GetRow( 1));
      m_TransformationMatrix3D.ReplaceRow( 2, ROTATION.GetMatrix().GetRow( 2));
    }

    //! modification of the axis
    void TransformationMatrix3D::SetAxis( const coord::Axis &AXIS, const linal::Vector3D &NEW_AXIS)
    {
      for( size_t i = 0; i < 3; ++i)
      {
        m_TransformationMatrix3D( AXIS, i) = NEW_AXIS( i);
      }
    }

    void TransformationMatrix3D::SetOrigin( const linal::Vector3D &AXIS)
    {
      for( size_t i = 0; i < 3; ++i)
      {
        m_TransformationMatrix3D( 3, i) = AXIS( i);
      }
    }

    //! @brief returns whether this transformation is defined or not by looking at the scaling factor
    //! @return whether this transformation is defined or not
    bool TransformationMatrix3D::IsDefined() const
    {
      // return whether scaling factor is defined
      return util::IsDefined( m_TransformationMatrix3D( 3, 3));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write TransformationMatrix3D to std::ostream using the given util::Format
    std::ostream &TransformationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_TransformationMatrix3D, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

    //! write TransformationMatrix3D to std::ostream using the given util::Format
    std::ostream &TransformationMatrix3D::Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const
    {
      io::Serialize::Write( m_TransformationMatrix3D, OSTREAM, INDENT, FORMAT) << '\n';
      return OSTREAM;
    }

    //! read TransformationMatrix3D from io::IFStream
    std::istream &TransformationMatrix3D::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_TransformationMatrix3D, ISTREAM);
      return ISTREAM;
    }

    //! Set Unitary matrix, assuming that the matrix is currently composed of zeros
    void TransformationMatrix3D::SetUnitFromZero()
    {
      m_TransformationMatrix3D( 0, 0) = 1.0;
      m_TransformationMatrix3D( 1, 1) = 1.0;
      m_TransformationMatrix3D( 2, 2) = 1.0;
      m_TransformationMatrix3D( 3, 3) = 1.0;
    }

    //! boolean operator TRANSFORMATIONMATRIX3D_A == TRANSFORMATIONMATRIX3D_B
    bool operator ==
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    )
    {
      return
        std::equal
        (
          TRANSFORMATIONMATRIX3D_A.GetMatrix().Begin(),
          TRANSFORMATIONMATRIX3D_A.GetMatrix().End(),
          TRANSFORMATIONMATRIX3D_B.GetMatrix().Begin()
        );
    }

    //! boolean operator TRANSFORMATIONMATRIX3D_A != TRANSFORMATIONMATRIX3D_B
    bool operator !=
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B
    )
    {
      return !( TRANSFORMATIONMATRIX3D_A == TRANSFORMATIONMATRIX3D_B);
    }

    //! return Inverse of TransformationMatrix3D
    TransformationMatrix3D Inverse( const TransformationMatrix3D &MATRIX)
    {
      return TransformationMatrix3D( MATRIX).Invert();
    }

    //! check if two transformation matrices are similar within a given rotational and translational tolerance
    bool SimilarWithinTolerance
    (
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_A,
      const TransformationMatrix3D &TRANSFORMATIONMATRIX3D_B,
      const double TRANSLATION_TOLERANCE,
      const double ROTATION_TOLERANCE_RAD
    )
    {
      TransformationMatrix3D diff( TRANSFORMATIONMATRIX3D_A);
      diff( Inverse( TRANSFORMATIONMATRIX3D_B));

      // check for similar Translation and similar rotation
      return (
                 ( diff.GetTranslation().Norm() < TRANSLATION_TOLERANCE)
              && ( diff.GetRotation().EffectiveRotationAngle() < ROTATION_TOLERANCE_RAD)
             );
    }

  } // namespace math
} // namespace bcl
