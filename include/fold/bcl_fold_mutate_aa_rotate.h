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

#ifndef BCL_FOLD_MUTATE_AA_ROTATE_H_
#define BCL_FOLD_MUTATE_AA_ROTATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_line_segment_3d.h"
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateAARotate
    //! @brief This class is for rotating an amino acid around an arbitrary axis and around an arbitrary point of origin.
    //! @details The axis and origin of rotation are defined as member variables. The class also provides functionality to
    //! get the transformation matrix that results from the combination of the rotation axis, origin of rotation, and
    //! rotation amount. Rotation is in radians.
    //!
    //! @see @link example_fold_mutate_aa_rotate.cpp @endlink
    //! @author alexanns
    //! @date Sep 3, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateAARotate :
      public math::MutateInterface< biol::AABase>
    {

    private:

    //////////
    // data //
    //////////

      //! the axis around which the rotation will take place
      coord::LineSegment3D m_RotationAxis;

      //! the origin point around which the rotation will occur
      linal::Vector3D m_RotationOrigin;

      //! the amount of rotation
      double m_Rotation;

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
      MutateAARotate();

      //! @brief constructor taking the member variable parameters
      //! @param ROTATION_AXIS the axis around which the rotation will take place
      //! @param ROTATION_ORIGIN the origin point around which the rotation will occur
      //! @param ROTATION the amount of rotation
      MutateAARotate
      (
        const coord::LineSegment3D &ROTATION_AXIS, const linal::Vector3D &ROTATION_ORIGIN, const double ROTATION
      );

      //! @brief Clone function
      //! @return pointer to new MutateAARotate
      MutateAARotate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to rotate a residue
      //! @return the transformation matrix that is necessary to rotate according to member variables
      math::TransformationMatrix3D GetTransformationMatrix() const;

      //! @brief GetRotationAxis provides the axis around which the rotation will take place
      //! @return the axis around which the rotation will take place
      const coord::LineSegment3D &GetRotationAxis() const;

      //! @brief GetRotationOrigin gives the origin point around which the rotation will occur
      //! @return the origin point around which the rotation will occur
      const linal::Vector3D &GetRotationOrigin() const;

      //! @brief GetRotation gives the amount of rotation
      //! @return the amount of rotation
      double GetRotation() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param RESIDUE Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< biol::AABase> operator()( const biol::AABase &RESIDUE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateAARotate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_AA_ROTATE_H_ 
