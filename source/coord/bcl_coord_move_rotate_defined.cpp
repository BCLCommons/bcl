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
#include "coord/bcl_coord_move_rotate_defined.h"

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
    const util::SiPtr< const util::ObjectInterface> MoveRotateDefined::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateDefined())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateDefined::MoveRotateDefined() :
      m_RotationAngle(),
      m_RotationAxis(),
      m_RotateInternal()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param ROTATION_AXIS axis to rotate around
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateDefined::MoveRotateDefined
    (
      const double ROTATION_ANGLE_RAD,
      const Axis &ROTATION_AXIS,
      const bool INTERNAL
    ) :
      m_RotationAngle( ROTATION_ANGLE_RAD),
      m_RotationAxis( ROTATION_AXIS),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief copy constructor
    //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
    MoveRotateDefined::MoveRotateDefined( const MoveRotateDefined &MOVE_ROTATE_DEFINED) :
      m_RotationAngle( MOVE_ROTATE_DEFINED.m_RotationAngle),
      m_RotationAxis( MOVE_ROTATE_DEFINED.m_RotationAxis),
      m_RotateInternal( MOVE_ROTATE_DEFINED.m_RotateInternal)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveRotateDefined *MoveRotateDefined::Clone() const
    {
      return new MoveRotateDefined( *this);
    }

    //! @brief destructor
    MoveRotateDefined::~MoveRotateDefined()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveRotateDefined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateDefined::GetAlias() const
    {
      static const std::string s_name( "MoveRotateDefined");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateDefined::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Rotate object around given axis and angle.");
      serializer.AddInitializer
      (
        "angle",
        "rotation angle",
        io::Serialization::GetAgent( &m_RotationAngle)
      );
      serializer.AddInitializer
      (
        "axis",
        "axis of rotation",
        io::Serialization::GetAgent( &m_RotationAxis)
      );
      serializer.AddInitializer
      (
        "internal",
        "move objects to origin before rotating them",
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
    MovableInterface &MoveRotateDefined::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_RotateInternal)
      {
        // transform with generated transformationmatrix
        MOVEABLE_OBJECT.Transform( GenerateTransformationInternal( MOVEABLE_OBJECT));
      }
      else
      {
        // rotate with generate rotationmatrix
        MOVEABLE_OBJECT.Rotate( GenerateRotation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix3D for internal coordinates (center of MOVEABLE_OBJECT as origin)
    //! @param MOVEABLE_OBJECT
    //! @return a translation vector
    const math::TransformationMatrix3D MoveRotateDefined::GenerateTransformationInternal( const MovableInterface &MOVEABLE_OBJECT) const
    {
      // apply transformation to origin
      math::TransformationMatrix3D new_matrix( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));

      // generate random rotation
      new_matrix( GenerateRotation());

      // move back to original position in space
      new_matrix( MOVEABLE_OBJECT.GetOrientation());

      // end
      return new_matrix;
    }

    //! @brief generate a RotationMatrix3D
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateDefined::GenerateRotation() const
    {
      return math::RotationMatrix3D( m_RotationAxis, m_RotationAngle);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
    MoveRotateDefined &MoveRotateDefined::operator =( const MoveRotateDefined &MOVE_ROTATE_DEFINED)
    {
      // update members
      m_RotationAngle  = MOVE_ROTATE_DEFINED.m_RotationAngle;
      m_RotationAxis   = MOVE_ROTATE_DEFINED.m_RotationAxis;
      m_RotateInternal = MOVE_ROTATE_DEFINED.m_RotateInternal;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
    //! @param AXIS Axis around which the flip is going to be applied
    //! @return an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
    MoveRotateDefined MoveRotateDefined::GetFlipMove( const Axis &AXIS)
    {
      // return MoveRotateDefined that rotates 180 degrees (flips) around provided AXIS
      return MoveRotateDefined( math::g_Pi, AXIS);
    }

  } // namespace coord
} // namespace bcl
