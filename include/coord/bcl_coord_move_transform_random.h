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

#ifndef BCL_COORD_MOVE_TRANSFORM_RANDOM_H_
#define BCL_COORD_MOVE_TRANSFORM_RANDOM_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_move_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveTransformRandom
    //! @brief MoveInterface derived class, that transforms a MovableInterface randomly
    //! @details the applied transformation is randomized between a min and a max translation and rotation
    //!
    //! @see @link example_coord_move_transform_random.cpp @endlink
    //! @author woetzen
    //! @date Feb 1, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveTransformRandom :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_MinTranslation;    //!< minimal distance for translation in each direction (x, y, z)
      linal::Vector3D m_MaxTranslation;    //!< maximal distance for translation in each direction (x, y, z)
      linal::Vector3D m_MinRotationAngles; //!< minimal angle for rotation (around x-, y-, z-axis)
      linal::Vector3D m_MaxRotationAngles; //!< maximal angle for rotation (around x-, y-, z-axis)

      //! flag - if it is set to true, will move objects center to origin before applying translation, and move it back
      bool m_TransformInternal;

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
      MoveTransformRandom();

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_TRANSLATION        maximal distance for translation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
      MoveTransformRandom
      (
        const double MAX_TRANSLATION,
        const double MAX_ROTATION_ANGLE_RAD,
        const bool INTERNAL = true
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_TRANSLATION        maximal distance for translation for each direction
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation for each direction
      //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
      MoveTransformRandom
      (
        const linal::Vector3D &MAX_TRANSLATION,
        const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
        const bool INTERNAL = true
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_TRANSLATION        minimal distance for translation
      //! @param MAX_TRANSLATION        maximal distance for translation
      //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
      MoveTransformRandom
      (
        const double MIN_TRANSLATION,
        const double MAX_TRANSLATION,
        const double MIN_ROTATION_ANGLE_RAD,
        const double MAX_ROTATION_ANGLE_RAD,
        const bool INTERNAL = true
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_TRANSLATION        minimal distance for translation for each direction
      //! @param MAX_TRANSLATION        maximal distance for translation for each direction
      //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation for each direction
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation for each direction
      //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
      MoveTransformRandom
      (
        const linal::Vector3D &MIN_TRANSLATION,
        const linal::Vector3D &MAX_TRANSLATION,
        const linal::Vector3D &MIN_ROTATION_ANGLE_RAD,
        const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
        const bool INTERNAL = true
      );

      //! @brief Clone function
      //! @return pointer to new Move
      MoveTransformRandom *Clone() const;

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

      //! @brief generate a TransformationMatrix
      //! @param MOVEABLE_OBJECT
      //! @return Transformationmatrix
      const math::TransformationMatrix3D GenerateRandomTransformation( const MovableInterface &MOVEABLE_OBJECT) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MoveTransformRandom

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_TRANSFORM_RANDOM_H_
