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
#include "fold/bcl_fold_mutate_aa_psi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateAAPsi::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAAPsi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAAPsi::MutateAAPsi() :
      m_Rotation()
    {
    }

    //! @brief constructor taking a desired rotation amount in radians
    //! @param ROTATION the desired psi rotation amount in radians
    MutateAAPsi::MutateAAPsi( const double ROTATION) :
      m_Rotation( ROTATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAAPsi
    MutateAAPsi *MutateAAPsi::Clone() const
    {
      return new MutateAAPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAAPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to rotate psi of a residue
    //! @param RESIDUE the residue for which the transformation matrix will be calculated
    //! @return the transformation matrix that is necessary to rotate psi of "RESIDUE" by "m_Rotation"
    math::TransformationMatrix3D MutateAAPsi::GetTransformationMatrix( const biol::AABase &RESIDUE) const
    {
      // get the C coordinates of "RESIDUE"
      const linal::Vector3D &c_atom_coordinates( RESIDUE.GetAtom( biol::GetAtomTypes().C).GetCoordinates());

      // get the CA coordinates of "RESIDUE"
      const linal::Vector3D &ca_atom_coordinates( RESIDUE.GetCA().GetCoordinates());

      // create a LineSegment3D "rotation_axis" which will be used as the psi rotation axis with CA as origin
      const coord::LineSegment3D rotation_axis( ca_atom_coordinates, c_atom_coordinates);

      // create TransformationMatrix3D and add the necessary transformations for rotating psi by "m_Rotation"
      math::TransformationMatrix3D transform( -ca_atom_coordinates);
      transform( math::RotationMatrix3D( rotation_axis.GetDirection(), m_Rotation));
      transform( ca_atom_coordinates);

      // return "transform"
      return transform;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param RESIDUE Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AABase> MutateAAPsi::operator()( const biol::AABase &RESIDUE) const
    {
      // create a ShPtr to a new residue
      util::ShPtr< biol::AABase> new_aa( RESIDUE.Clone());

      // just rotate the oxygen atom of the residue if it exists
      // get the Oxygen from this amino acid
      const biol::Atom &oxygen( RESIDUE.GetAtom( biol::GetAtomTypes().O));

      // if the type is defined, thus the oxygen atom exists
      if( oxygen.GetType().IsDefined())
      {
        // make a copy of the atom
        biol::Atom new_oxygen( oxygen);

        // apply the psi transformation
        new_oxygen.Transform( GetTransformationMatrix( RESIDUE));

        // set the oxygen of "new_aa" to "new_oxygen"
        new_aa->SetAtom( new_oxygen);
      }

      // create mutate result with "new_aa"
      const math::MutateResult< biol::AABase> mutate_result( new_aa, *this);

      // return "mutate_result"
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAAPsi::Read( std::istream &ISTREAM)
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
    std::ostream &MutateAAPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Rotation, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
