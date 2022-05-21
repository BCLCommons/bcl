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
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_rotate.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateTransformationMatrix3DRotate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTransformationMatrix3DRotate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateTransformationMatrix3DRotate::MutateTransformationMatrix3DRotate() :
      m_Axis(),
      m_MaxRotation(),
      m_MinRotation()
    {
    }

    //! @brief construct from an axis and an amount of rotation
    //! @param AXIS is the math::RotationAxis about which the rotation will be performed
    //! @param MAX_ROTATION is a double denoting the amount of of rotation
    MutateTransformationMatrix3DRotate::MutateTransformationMatrix3DRotate( const coord::Axis &AXIS, const double MAX_ROTATION, const double MIN_ROTATION) :
      m_Axis( AXIS),
      m_MaxRotation( MAX_ROTATION),
      m_MinRotation( MIN_ROTATION)
    {
    }

    //! virtual copy constructor
    MutateTransformationMatrix3DRotate *MutateTransformationMatrix3DRotate::Clone() const
    {
      return new MutateTransformationMatrix3DRotate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateTransformationMatrix3DRotate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    math::MutateResult< math::TransformationMatrix3D> MutateTransformationMatrix3DRotate::operator()
    (
      const math::TransformationMatrix3D &MATRIX
    ) const
    {
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate::operator()");
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_Axis() " + util::Format()( m_Axis));
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_MaxRotation() " + util::Format()( m_MaxRotation));
      BCL_MessageDbg( "MutateTransformationMatrix3DRotate:: m_MinRotation() " + util::Format()( m_MinRotation));

      // get the transformation matrix of "MATRIX" at the origin
      util::ShPtr< math::TransformationMatrix3D> sp_origin( new math::TransformationMatrix3D());

      // rotate "origin" as defined by "m_Axis" and "m_MaxRotation", "m_MinRotation"
      if( m_MaxRotation == m_MinRotation)
      {
        sp_origin->operator()( m_Axis, m_MaxRotation);
      }
      else
      {
        sp_origin->operator()( m_Axis, random::GetGlobalRandom().Random< double>( m_MinRotation, m_MaxRotation));
      }

      // translate "origin" back to starting position
      sp_origin->operator()( MATRIX);

      // return unchanged "MATRIX"
      return math::MutateResult< math::TransformationMatrix3D>( sp_origin, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateTransformationMatrix3DRotate::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Axis, ISTREAM);
      io::Serialize::Read( m_MaxRotation, ISTREAM);
      io::Serialize::Read( m_MinRotation, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &MutateTransformationMatrix3DRotate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Axis, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxRotation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinRotation, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
