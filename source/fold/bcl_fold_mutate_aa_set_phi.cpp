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
#include "fold/bcl_fold_mutate_aa_set_phi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "fold/bcl_fold_mutate_aa_phi.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateAASetPhi::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAASetPhi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAASetPhi::MutateAASetPhi() :
      m_PreviousC(),
      m_DesiredPhi()
    {
    }

    //! @brief constructor taking member variables
    //! @param PREVIOUS_C the C atom of the residue previous to the residue whose phi will be changed
    //! @param DESIRED_PHI the desired phi angle value in radians
    MutateAASetPhi::MutateAASetPhi( const biol::Atom &PREVIOUS_C, const double DESIRED_PHI) :
      m_PreviousC( PREVIOUS_C),
      m_DesiredPhi( DESIRED_PHI)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAASetPhi
    MutateAASetPhi *MutateAASetPhi::Clone() const
    {
      return new MutateAASetPhi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAASetPhi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to set phi of a residue
    //! @param RESIDUE the residue for which the transformation matrix will be calculated
    //! @return the transformation matrix that is necessary to set phi of "RESIDUE" to "m_DesiredPhi"
    math::TransformationMatrix3D MutateAASetPhi::GetTransformationMatrix( const biol::AABase &RESIDUE) const
    {
      // return the transformation matrix that is necessary to set phi of "RESIDUE" to "m_DesiredPhi"
      return GetMutateAAPhi( RESIDUE).GetTransformationMatrix( RESIDUE);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param RESIDUE Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AABase> MutateAASetPhi::operator()( const biol::AABase &RESIDUE) const
    {
      // get the residue with phi mutated as desired in a MutateResult object
      const math::MutateResult< biol::AABase> mutated_phi( GetMutateAAPhi( RESIDUE)( RESIDUE));

      // return "mutated_phi"
      return mutated_phi;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAASetPhi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PreviousC, ISTREAM);
      io::Serialize::Read( m_DesiredPhi, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateAASetPhi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PreviousC, OSTREAM, INDENT);
      io::Serialize::Write( m_DesiredPhi, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetMutateAAPhi provides the MutateAAPhi object necesssary to set phi as desired
    //! @param RESIDUE is the residue whose phi is going to be set
    //! @return the MutateAAPhi object necesssary to set phi as desired
    MutateAAPhi MutateAASetPhi::GetMutateAAPhi( const biol::AABase &RESIDUE) const
    {
      // calculate the rotation that is needed to set phi to "m_DesiredPhi"
      double current_phi( RESIDUE.Phi());
      if( !util::IsDefined( current_phi))
      {
        current_phi = RESIDUE.CalculatePhi( m_PreviousC);
      }
      const double rotation( m_DesiredPhi - current_phi);

      // message the needed rotation to set phi to "m_DesiredPhi"
      BCL_MessageDbg
      (
        "current phi is " + util::Format()( current_phi)
        + " desired phi is " + util::Format()( m_DesiredPhi)
        + " needed rotation is " + util::Format()( rotation)
      );

      // create a MutateAAPhi object using "rotation"
      const MutateAAPhi phi_mutate( -rotation);

      // return MutateAAPhi which can be used to set phi of "RESIDUE" as desired
      return phi_mutate;
    }

  } // namespace fold
} // namespace bcl
