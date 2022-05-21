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

#ifndef BCL_COORD_MOVABLE_INTERFACE_H_
#define BCL_COORD_MOVABLE_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_orientation_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MovableInterface
    //! @brief Interface that defines functions related to all movable classes
    //! This class provides the interface to functions that can applied to any class that is defined as Movable
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Nov 16, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MovableInterface :
      public OrientationInterface
    {

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Translation to be applied
      virtual void Translate( const linal::Vector3D &TRANSLATION) = 0;

      //! @brief transform the object by a given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      virtual void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D) = 0;

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      virtual void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D) = 0;

      //! @brief transform randomly by inner coordinates
      //! @param MAX_TRANSLATION maximal move distance
      //! @param MAX_ROTATION    maximal rotation angle
      void RandomTransformation( const double MAX_TRANSLATION, const double MAX_ROTATION);

    }; // class MovableInterface

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVABLE_INTERFACE_H_
