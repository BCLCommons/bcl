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

#ifndef BCL_COORD_MOVE_ROTATE_RANDOM_H_
#define BCL_COORD_MOVE_ROTATE_RANDOM_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_move_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveRotateRandom
    //! @brief MoveRotateRandom allows a random rotation of a given object
    //! @details This class rotates the given object randomly making sure the rotations are between the given
    //! min and max rotations allowed for each axis
    //!
    //! @see @link example_coord_move_rotate_random.cpp @endlink
    //! @author woetzen
    //! @date Feb 1, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveRotateRandom :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_MinRotationAngles; //!< minimal angle for rotation (around x-, y-, z-axis)
      linal::Vector3D m_MaxRotationAngles; //!< maximal angle for rotation (around x-, y-, z-axis)
      //! flag - if it is set to true, will move objects center to origin before applying rotation, and move it back
      bool           m_RotateInternal;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoveRotateRandom();

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
      MoveRotateRandom( const double MAX_ROTATION_ANGLE_RAD, const bool INTERNAL = true);

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
      //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
      MoveRotateRandom( const linal::Vector3D &MAX_ROTATION_ANGLES_RAD, const bool INTERNAL = true);

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
      MoveRotateRandom
      (
        const double MIN_ROTATION_ANGLE_RAD,
        const double MAX_ROTATION_ANGLE_RAD,
        const bool INTERNAL = true
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
      //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
      //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
      MoveRotateRandom
      (
        const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
        const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
        const bool INTERNAL = true
      );

      //! @brief copy constructor
      //! @param MOVE_ROTATE_RANDOM MoveRotateRandom to be copied
      MoveRotateRandom( const MoveRotateRandom &MOVE_ROTATE_RANDOM);

      //! @brief Clone function
      //! @return pointer to new Move
      MoveRotateRandom *Clone() const;

      //! @brief destructor
      ~MoveRotateRandom();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief move the MovableInterface derived object
      //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

      //! @brief generate a TransformationMatrix3D for internal coordinates (center of MOVEABLE_OBJECT as origin)
      //! @param MOVEABLE_OBJECT
      //! @return a translation vector
      const math::TransformationMatrix3D GenerateRandomTransformationInternal( const MovableInterface &MOVEABLE_OBJECT) const;

      //! @brief generate a RotationMatrix3D
      //! @return a translation vector
      const math::RotationMatrix3D GenerateRandomRotation() const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief generate a RotationMatrix3D given the max rotation angles for each axis
      //! @return a translation vector
      static const math::RotationMatrix3D GenerateRandomRotation
      (
        const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
        const linal::Vector3D &MIN_ROTATION_ANGLES_RAD = linal::Vector3D( 0.0)
      );

    }; // class MoveRotateRandom

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_ROTATE_RANDOM_H_
