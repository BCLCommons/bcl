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

#ifndef BCL_COORD_MOVE_ROTATE_RANDOM_EXTERNAL_REFERENCE_H_
#define BCL_COORD_MOVE_ROTATE_RANDOM_EXTERNAL_REFERENCE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_move_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveRotateRandomExternalReference
    //! @brief performs a random rotation on an object using an external coordinate system
    //! @details This class performs a random rotation on an object using an external coordinate system that is defined
    //!          by the user.  An enum member allows for the specification of how the rotation is to be applied.
    //!
    //! @see @link example_coord_move_rotate_random_external_reference.cpp @endlink
    //! @author weinerbe
    //! @date Sep 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveRotateRandomExternalReference :
      public MoveInterface
    {

    public:

    //////////
    // data //
    //////////

      //! enumerator for the rotation method
      enum MethodType
      {
        e_Internal, //!< move the object to the reference orientation before applying rotation, this is the default
        e_InternalRotate, //!< rotate the object to match the reference orientation before applying rotation
        e_InternalTranslate, //!< translate the object to the reference orientation before applying rotation
        e_External, //!< rotate the object with respect to the reference orientation
        s_NumberMethodTypes //!< number of method types
      };

      //! @brief MethodType as string
      //! @param METHOD_TYPE the MethodType
      //! @return the string for the MethodType
      static const std::string &GetMethodDescriptor( const MethodType &METHOD_TYPE);

      //! @brief MethodTypeEnum is used for I/O of MethodType
      typedef util::WrapperEnum< MethodType, &GetMethodDescriptor, s_NumberMethodTypes> MethodTypeEnum;

    private:

      //! minimal angle for rotation (around x-, y-, z-axis)
      linal::Vector3D m_MinRotationAngles;

      //! maximal angle for rotation (around x-, y-, z-axis)
      linal::Vector3D m_MaxRotationAngles;

      //! external reference orientation
      math::TransformationMatrix3D m_ReferenceOrientation;

      //! rotation method to be used
      MethodTypeEnum m_Method;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoveRotateRandomExternalReference();

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param REFERENCE_ORIENTATION external reference orientation
      //! @param ROTATION_METHOD rotation method to be used
      MoveRotateRandomExternalReference
      (
        const double MAX_ROTATION_ANGLE_RAD,
        const math::TransformationMatrix3D &REFERENCE_ORIENTATION = math::TransformationMatrix3D(),
        const MethodType ROTATION_METHOD = e_Internal
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
      //! @param REFERENCE_ORIENTATION external reference orientation
      //! @param ROTATION_METHOD rotation method to be used
      MoveRotateRandomExternalReference
      (
        const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
        const math::TransformationMatrix3D &REFERENCE_ORIENTATION = math::TransformationMatrix3D(),
        const MethodType ROTATION_METHOD = e_Internal
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
      //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
      //! @param REFERENCE_ORIENTATION external reference orientation
      //! @param ROTATION_METHOD rotation method to be used
      MoveRotateRandomExternalReference
      (
        const double MIN_ROTATION_ANGLE_RAD,
        const double MAX_ROTATION_ANGLE_RAD,
        const math::TransformationMatrix3D &REFERENCE_ORIENTATION = math::TransformationMatrix3D(),
        const MethodType ROTATION_METHOD = e_Internal
      );

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
      //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
      //! @param REFERENCE_ORIENTATION external reference orientation
      //! @param ROTATION_METHOD rotation method to be used
      MoveRotateRandomExternalReference
      (
        const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
        const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
        const math::TransformationMatrix3D &REFERENCE_ORIENTATION = math::TransformationMatrix3D(),
        const MethodType ROTATION_METHOD = e_Internal
      );

      //! @brief Clone function
      //! @return pointer to new MoveRotateRandomExternalReference
      MoveRotateRandomExternalReference *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the rotation method
      //! @return the rotation method
      const MethodTypeEnum &GetMethod() const
      {
        return m_Method;
      }

      //! @brief sets the rotation method
      //! @param ROTATION_METHOD rotation method to be set
      void SetMethod( const MethodType ROTATION_METHOD)
      {
        m_Method = ROTATION_METHOD;
      }

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
      //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief generate a random rotation
      //! @return transformation matrix for the rotation
      math::TransformationMatrix3D GenerateRandomRotation() const;

    }; // class MoveRotateRandomExternalReference

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_MOVE_ROTATE_RANDOM_EXTERNAL_REFERENCE_H_
