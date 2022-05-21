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

#ifndef BCL_COORD_MOVE_ROTATE_DEFINED_H_
#define BCL_COORD_MOVE_ROTATE_DEFINED_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_axes.h"
#include "bcl_coord_move_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveRotateDefined
    //! @brief MoveInterface derived class that allows rotation around a given axis of given angles for a given object
    //! @details MoveRotateDefined allows a defined angles of rotation around a given axis. This rotation can be internal
    //! or external depending.
    //!
    //! @see @link example_coord_move_rotate_defined.cpp @endlink
    //! @author woetzen
    //! @date Feb 1, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveRotateDefined :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      double  m_RotationAngle; //!< angle for rotation (around x-, y-, z-axis)
      Axis    m_RotationAxis;  //!< axis around which rotation is done

      //! flag - if it is set to true, will move objects center to origin before applying rotation, and move it back
      bool    m_RotateInternal;

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
      MoveRotateDefined();

      //! @brief construct from undirected maximal translation and rotation
      //! @param ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param ROTATION_AXIS axis to be rotated around
      //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
      MoveRotateDefined
      (
        const double ROTATION_ANGLE_RAD,
        const Axis &ROTATION_AXIS,
        const bool INTERNAL = true
      );

      //! @brief copy constructor
      //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
      MoveRotateDefined( const MoveRotateDefined &MOVE_ROTATE_DEFINED);

      //! @brief Clone function
      //! @return pointer to new Move
      MoveRotateDefined *Clone() const;

      //! @brief destructor
      ~MoveRotateDefined();

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
      const math::TransformationMatrix3D GenerateTransformationInternal( const MovableInterface &MOVEABLE_OBJECT) const;

      //! @brief generate a RotationMatrix3D
      //! @return a translation vector
      const math::RotationMatrix3D GenerateRotation() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
      MoveRotateDefined &operator =( const MoveRotateDefined &MOVE_ROTATE_DEFINED);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief returns an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
      //! @param AXIS Axis around which the flip is going to be applied
      //! @return an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
      static MoveRotateDefined GetFlipMove( const Axis &AXIS);

    }; // class MoveRotateDefined

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_ROTATE_DEFINED_H_
