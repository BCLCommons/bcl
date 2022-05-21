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
#include "coord/bcl_coord_move_rotate_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
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
    const util::SiPtr< const util::ObjectInterface> MoveRotateRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateRandom::MoveRotateRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_RotateInternal( INTERNAL)
    {

    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLES_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_RotateInternal( INTERNAL)
    {

    }

    //! @brief copy constructor
    //! @param MOVE_ROTATE_RANDOM MoveRotateRandom to be copied
    MoveRotateRandom::MoveRotateRandom( const MoveRotateRandom &MOVE_ROTATE_RANDOM) :
      m_MinRotationAngles( MOVE_ROTATE_RANDOM.m_MinRotationAngles),
      m_MaxRotationAngles( MOVE_ROTATE_RANDOM.m_MaxRotationAngles),
      m_RotateInternal( MOVE_ROTATE_RANDOM.m_RotateInternal)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveRotateRandom *MoveRotateRandom::Clone() const
    {
      return new MoveRotateRandom( *this);
    }

    //! @brief destructor
    MoveRotateRandom::~MoveRotateRandom()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveRotateRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateRandom::GetAlias() const
    {
      static const std::string s_name( "MoveRotateRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Random rotation of an object.");
      serializer.AddInitializer
      (
        "min angles",
        "minimum angles for rotation",
        io::Serialization::GetAgent( &m_MinRotationAngles),
        "(0,0,0)"
      );
      serializer.AddInitializer
      (
        "max angles",
        "maximum angles for rotation",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      serializer.AddInitializer
      (
        "internal",
        "move object to origin before rotation",
        io::Serialization::GetAgent( &m_RotateInternal),
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
    MovableInterface &MoveRotateRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_RotateInternal)
      {
        // transform with generated transformationmatrix
        MOVEABLE_OBJECT.Transform( GenerateRandomTransformationInternal( MOVEABLE_OBJECT));
      }
      else
      {
        // rotate with generate rotationmatrix
        MOVEABLE_OBJECT.Rotate( GenerateRandomRotation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix3D for internal coordinates (center of MOVEABLE_OBJECT as origin)
    //! @param MOVEABLE_OBJECT
    //! @return a translation vector
    const math::TransformationMatrix3D
    MoveRotateRandom::GenerateRandomTransformationInternal
    (
      const MovableInterface &MOVEABLE_OBJECT
    ) const
    {
      // apply transformation to origin
      math::TransformationMatrix3D new_matrix( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));

      // generate random rotation
      new_matrix( GenerateRandomRotation());

      // move back to original position in space
      new_matrix( MOVEABLE_OBJECT.GetOrientation());

      // end
      return new_matrix;
    }

    //! @brief generate a RotationMatrix3D
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateRandom::GenerateRandomRotation() const
    {
      return GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a RotationMatrix3D given the max rotation angles for each axis
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateRandom::GenerateRandomRotation
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD
    )
    {
      // calculate the difference
      const linal::Vector3D diff_angles( MAX_ROTATION_ANGLES_RAD - MIN_ROTATION_ANGLES_RAD);

      double rotation_around_x
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.X()
      );
      double rotation_around_y
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.Y()
      );
      double rotation_around_z
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.Z()
      );

      // add the minimum rotations
      rotation_around_x += rotation_around_x < 0.0 ? -MIN_ROTATION_ANGLES_RAD.X() : MIN_ROTATION_ANGLES_RAD.X();
      rotation_around_y += rotation_around_y < 0.0 ? -MIN_ROTATION_ANGLES_RAD.Y() : MIN_ROTATION_ANGLES_RAD.Y();
      rotation_around_z += rotation_around_z < 0.0 ? -MIN_ROTATION_ANGLES_RAD.Z() : MIN_ROTATION_ANGLES_RAD.Z();

      // construct the rotatio matrix from these three rotations
      math::RotationMatrix3D rotation( GetAxes().e_X, rotation_around_x);
      rotation *= math::RotationMatrix3D( GetAxes().e_Y, rotation_around_y);
      rotation *= math::RotationMatrix3D( GetAxes().e_Z, rotation_around_z);

      // end
      return rotation;
    }

  } // namespace coord
} // namespace bcl
