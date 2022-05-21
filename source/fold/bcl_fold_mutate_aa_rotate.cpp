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
#include "fold/bcl_fold_mutate_aa_rotate.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateAARotate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAARotate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAARotate::MutateAARotate() :
      m_RotationAxis(),
      m_RotationOrigin(),
      m_Rotation()
    {
    }

    //! @brief constructor taking the member variable parameters
    //! @param ROTATION_AXIS the axis around which the rotation will take place
    //! @param ROTATION_ORIGIN the origin point around which the rotation will occur
    //! @param ROTATION the amount of rotation
    MutateAARotate::MutateAARotate
    (
      const coord::LineSegment3D &ROTATION_AXIS, const linal::Vector3D &ROTATION_ORIGIN, const double ROTATION
    ) :
      m_RotationAxis( ROTATION_AXIS),
      m_RotationOrigin( ROTATION_ORIGIN),
      m_Rotation( ROTATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAARotate
    MutateAARotate *MutateAARotate::Clone() const
    {
      return new MutateAARotate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAARotate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to rotate a residue
    //! @return the transformation matrix that is necessary to rotate "RESIDUE" according to member variables
    math::TransformationMatrix3D MutateAARotate::GetTransformationMatrix() const
    {
      // create TransformationMatrix3D "transform" and initialize with translation to the origin of rotation
      math::TransformationMatrix3D transform( -m_RotationOrigin);

      // add the desired rotation according to "m_RotationAxis" and "m_Rotation" to "transform"
      transform( math::RotationMatrix3D( m_RotationAxis.GetDirection(), m_Rotation));

      // add the translation of moving back to the original position
      transform( m_RotationOrigin);

      // return the transformation matrix "transform"
      return transform;
    }

    //! @brief GetRotationAxis provides the axis around which the rotation will take place
    //! @return the axis around which the rotation will take place
    const coord::LineSegment3D &MutateAARotate::GetRotationAxis() const
    {
      return m_RotationAxis;
    }

    //! @brief GetRotationOrigin gives the origin point around which the rotation will occur
    //! @return the origin point around which the rotation will occur
    const linal::Vector3D &MutateAARotate::GetRotationOrigin() const
    {
      return m_RotationOrigin;
    }

    //! @brief GetRotation gives the amount of rotation
    //! @return the amount of rotation
    double MutateAARotate::GetRotation() const
    {
      return m_Rotation;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param RESIDUE Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AABase> MutateAARotate::operator()( const biol::AABase &RESIDUE) const
    {
      // create ShPtr to new residue copy of "RESIDUE" which will be rotated
      util::ShPtr< biol::AABase> new_aa( RESIDUE.Clone());

      // transform "new_aa" according to the transformation matrix provided by GetTransformationMatrix
      new_aa->Transform( GetTransformationMatrix());

      // return the transformed "new_aa" in a MutateResult object
      return math::MutateResult< biol::AABase>( new_aa, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAARotate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RotationAxis, ISTREAM);
      io::Serialize::Read( m_RotationOrigin, ISTREAM);
      io::Serialize::Read( m_Rotation, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateAARotate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RotationAxis, OSTREAM, INDENT);
      io::Serialize::Write( m_RotationOrigin, OSTREAM, INDENT);
      io::Serialize::Write( m_Rotation, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
