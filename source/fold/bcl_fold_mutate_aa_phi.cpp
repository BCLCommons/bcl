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
#include "fold/bcl_fold_mutate_aa_phi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "fold/bcl_fold_mutate_aa_rotate.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateAAPhi::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAAPhi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAAPhi::MutateAAPhi() :
      m_Rotation()
    {
    }

    //! @brief constructor taking a desired rotation amount in radians
    //! @param ROTATION the desired phi rotation amount in radians
    MutateAAPhi::MutateAAPhi( const double ROTATION) :
      m_Rotation( ROTATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAAPhi
    MutateAAPhi *MutateAAPhi::Clone() const
    {
      return new MutateAAPhi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAAPhi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to rotate phi of a residue
    //! @param RESIDUE the residue for which the transformation matrix will be calculated
    //! @return the transformation matrix that is necessary to rotate phi of "RESIDUE" by "m_Rotation"
    math::TransformationMatrix3D MutateAAPhi::GetTransformationMatrix( const biol::AABase &RESIDUE) const
    {
      // return the transformation matrix from the MutateAARotate created by GetMutateAARotate from RESIDUE
      return GetMutateAARotate( RESIDUE).GetTransformationMatrix();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param RESIDUE Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AABase> MutateAAPhi::operator()( const biol::AABase &RESIDUE) const
    {
      // create a new residue based on the desired phi change
      const math::MutateResult< biol::AABase> new_aa( GetMutateAARotate( RESIDUE)( RESIDUE));

      // return the "new_aa"
      return new_aa;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAAPhi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Rotation, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateAAPhi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Rotation, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetMutateAARotate provides the MutateAARotate object necesssary to do the desired phi change
    //! @param RESIDUE is the residue whose phi is going to be changed
    //! @return the MutateAARotate object necesssary to do the desired phi change
    MutateAARotate MutateAAPhi::GetMutateAARotate( const biol::AABase &RESIDUE) const
    {
      // get the coordinates of the nitrogen of "RESIDUE"
      const linal::Vector3D &n_atom_coordinates( RESIDUE.GetAtom( biol::GetAtomTypes().N).GetCoordinates());

      // get the ca coordinates from "RESIDUE"
      const linal::Vector3D &ca_atom_coordinates( RESIDUE.GetCA().GetCoordinates());

      // create a LineSegment3D "rotation_axis" which will be used as the phi rotation axis with N as origin
      const coord::LineSegment3D rotation_axis( n_atom_coordinates, ca_atom_coordinates);

      // create a MutateAARotate object out of "rotation_axis", "ca_atom_coordinates", and "m_Rotation"
      const MutateAARotate rotate( rotation_axis, ca_atom_coordinates, m_Rotation);

      // return the MutateAARotate which can be used to perform the desired phi rotation
      return rotate;
    }

  } // namespace fold
} // namespace bcl
