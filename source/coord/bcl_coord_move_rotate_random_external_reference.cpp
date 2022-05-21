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
#include "coord/bcl_coord_move_rotate_random_external_reference.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //! @brief MethodType as string
    //! @param METHOD_TYPE the MethodType
    //! @return the string for the MethodType
    const std::string &MoveRotateRandomExternalReference::GetMethodDescriptor( const MethodType &METHOD_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "INTERNAL", "INTERNAL_ROTATE", "INTERNAL_TRANSLATE", "EXTERNAL", GetStaticClassName< MethodType>()
      };

      return s_descriptors[ METHOD_TYPE];
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveRotateRandomExternalReference::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateRandomExternalReference())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference() :
      m_MinRotationAngles(),
      m_MaxRotationAngles(),
      m_ReferenceOrientation(),
      m_Method( e_Internal)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const double MAX_ROTATION_ANGLE_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLES_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoveRotateRandomExternalReference
    MoveRotateRandomExternalReference *MoveRotateRandomExternalReference::Clone() const
    {
      return new MoveRotateRandomExternalReference( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveRotateRandomExternalReference::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateRandomExternalReference::GetAlias() const
    {
      static const std::string s_name( "MoveRotateRandomExternalReference");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateRandomExternalReference::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random rotation relative to an external coordinate system.");
      serializer.AddInitializer
      (
        "min rotation",
        "minimum rotation angle",
        io::Serialization::GetAgent( &m_MinRotationAngles)
      );
      serializer.AddInitializer
      (
        "max rotation",
        "maximum rotation angles",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      // serializer.AddInitializer
      // (
      //   "reference orientation",
      //   "external reference orientation",
      //   io::Serialization::GetAgent( &m_ReferenceOrientation)
      // );
      serializer.AddInitializer
      (
        "method",
        "rotation method",
        io::Serialization::GetAgent( &m_Method)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveRotateRandomExternalReference::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // construct a transformation matrix to be applied to the object
      math::TransformationMatrix3D transform;

      // apply the move using the method defined by "m_Method"
      switch( m_Method)
      {
        case e_Internal:
          // move to the external reference
          transform( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));
          transform( m_ReferenceOrientation);

          // generate random rotation
          transform( GenerateRandomRotation());

          // move back to the original position
          transform( MOVEABLE_OBJECT.GetOrientation());

          break;

        case e_InternalRotate:
          // rotate to match the external reference
          transform( math::Inverse( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
          transform( m_ReferenceOrientation.GetRotation());

          // generate random rotation
          transform( GenerateRandomRotation());

          // rotate back to the original position
          transform( MOVEABLE_OBJECT.GetOrientation().GetRotation());

          break;

        case e_InternalTranslate:
          // translate to the external reference
          transform( -MOVEABLE_OBJECT.GetCenter());
          transform( m_ReferenceOrientation.GetOrigin());

          // generate random rotation
          transform( GenerateRandomRotation());

          // move back to the original position
          transform( MOVEABLE_OBJECT.GetCenter());

          break;

        case e_External:
          // generate random rotation
          transform( GenerateRandomRotation());

          break;
        default: break;
      }

      // apply the transformation
      MOVEABLE_OBJECT.Transform( transform);

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a random rotation
    //! @return transformation matrix for the rotation
    math::TransformationMatrix3D MoveRotateRandomExternalReference::GenerateRandomRotation() const
    {
      // initialize transformation matrix by moving reference to origin
      math::TransformationMatrix3D transform( math::Inverse( m_ReferenceOrientation));

      // apply a random rotation
      transform( MoveRotateRandom::GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles));

      // move back to the reference position
      transform( m_ReferenceOrientation);

      // end
      return transform;
    }

  } // namespace coord
} // namespace bcl
