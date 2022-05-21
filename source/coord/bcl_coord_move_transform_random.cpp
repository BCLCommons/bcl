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
#include "coord/bcl_coord_move_transform_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveTransformRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTransformRandom())
    );

    //! @brief default constructor
    MoveTransformRandom::MoveTransformRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const double MAX_TRANSLATION,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation for each direction
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radians for rotation for each direction
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const linal::Vector3D &MAX_TRANSLATION,
      const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION        minimal distance for translation
    //! @param MAX_TRANSLATION        maximal distance for translation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION        minimal distance for translation for each direction
    //! @param MAX_TRANSLATION        maximal distance for translation for each direction
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation for each direction
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation for each direction
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const linal::Vector3D &MIN_TRANSLATION,
      const linal::Vector3D &MAX_TRANSLATION,
      const linal::Vector3D &MIN_ROTATION_ANGLE_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTransformRandom *MoveTransformRandom::Clone() const
    {
      return new MoveTransformRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTransformRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTransformRandom::GetAlias() const
    {
      static const std::string s_name( "MoveTransformRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTransformRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random transformation.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "min rotation",
        "minimum rotation",
        io::Serialization::GetAgent( &m_MinRotationAngles)
      );
      serializer.AddInitializer
      (
        "max rotation",
        "maximum rotation",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      serializer.AddInitializer
      (
        "internal transformation",
        "move objects to origin before applying transformation",
        io::Serialization::GetAgent( &m_TransformInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTransformRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // transform with generated transformationmatrix
      MOVEABLE_OBJECT.Transform( GenerateRandomTransformation( MOVEABLE_OBJECT));

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix
    //! @param MOVEABLE_OBJECT
    //! @return Transformationmatrix
    const math::TransformationMatrix3D MoveTransformRandom::GenerateRandomTransformation( const MovableInterface &MOVEABLE_OBJECT) const
    {
      // new transformation matrix
      math::TransformationMatrix3D transformation;
      math::TransformationMatrix3D orientation;

      // random rotation and translation
      const math::RotationMatrix3D rand_rot( MoveRotateRandom::GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles));
      linal::Vector3D rand_trans( MoveTranslateRandom::GenerateRandomTranslation( m_MaxTranslation, m_MinTranslation));

      // move Movable to origin
      if( m_TransformInternal)
      {
        // get the start orientation
        orientation = MOVEABLE_OBJECT.GetOrientation();
        // move to origin
        transformation( math::Inverse( orientation));
      }

      // apply random rotation
      transformation( rand_rot);

      // move Movable back to original position and orientation
      if( m_TransformInternal)
      {
        // move back to original position
        transformation( orientation);
        // correct the translation
        rand_trans.Rotate( rand_rot);
      }

      // apply random transformation
      transformation( rand_trans);

      // end
      return transformation;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
