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

#ifndef BCL_COORD_ORIENTATION_INTERFACE_H_
#define BCL_COORD_ORIENTATION_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OrientationInterface
    //! @brief Interface that defines functions related to all classes with an orientation
    //! @details This class provides the interface to functions that can applied to any class that defines an Orientation
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date @date Aug 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API OrientationInterface :
      public virtual util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      virtual linal::Vector3D GetCenter() const = 0;

      //! @brief return the orientation of the object
      //! @return orientation
      virtual linal::Vector3D GetAxis( const Axis &AXIS) const
      {
        return *AXIS;
      }

      //! @brief return the orientation and Position as TransformationMatrix3D
      //! @return TransformationMatrix3D that defines orientation and position
      virtual const math::TransformationMatrix3D GetOrientation() const
      {
        // return
        return math::TransformationMatrix3D( GetCenter());
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief generate a random transformation around a given center
      //! @param MAX_TRANSLATION maximal translation
      //! @param MAX_ROTATION    maximal rotation angle (in rad)
      //! @param CENTER          origin around which the transformation will occur
      //! @return generated TransformationMatrix
      static
      math::TransformationMatrix3D
      GenerateRandomTransformationAroundCenter
      (
        const double MAX_TRANSLATION,
        const double MAX_ROTATION,
        const linal::Vector3D &CENTER = linal::Vector3D()
      );

    }; //class MovableInterface

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_ORIENTATION_INTERFACE_H_
