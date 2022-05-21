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

#ifndef BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_ROTATE_H_
#define BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_ROTATE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateTransformationMatrix3DRotate
    //! @brief MutateTransformationMatrix3DRotate is a class for mutating a transformation matrix
    //! by internally rotating it about a given axis
    //!
    //! @see @link example_restraint_mutate_transformation_matrix_3d_rotate.cpp @endlink
    //! @author alexanns
    //! @date 07/06/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateTransformationMatrix3DRotate :
      public math::MutateInterface< math::TransformationMatrix3D>
    {
    private:

    //////////
    // data //
    //////////

      coord::Axis m_Axis; //!< axis

      double m_MaxRotation; //!< max rotation
      double m_MinRotation; //!< min rotation

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateTransformationMatrix3DRotate();

      //! @brief construct from an axis and an amount of rotation
      //! @param AXIS is the math::RotationAxis about which the rotation will be performed
      //! @param MAX_ROTATION is a double denoting the maximum amount of of rotation
      //! @param MIN_ROTATION is a double denoting the minimum amount of of rotation
      MutateTransformationMatrix3DRotate( const coord::Axis &AXIS, const double MAX_ROTATION, const double MIN_ROTATION);

      //! virtual copy constructor
      MutateTransformationMatrix3DRotate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! virtual opartor taking an ARGUMENT and returning a mutate object of t_ArgumentType
      math::MutateResult< math::TransformationMatrix3D> operator()
      (
        const math::TransformationMatrix3D &MATRIX
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateTransformationMatrix3DRotate

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_ROTATE_H_
