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
#include "fold/bcl_fold_mutate_aa_set_psi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "fold/bcl_fold_mutate_aa_psi.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateAASetPsi::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAASetPsi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAASetPsi::MutateAASetPsi() :
      m_FollowingN(),
      m_DesiredPsi()
    {
    }

    //! @brief constructor taking member variables
    //! @param FOLLOWING_N the N atom of the residue following the residue whose psi will be changed
    //! @param DESIRED_PSI the desired psi angle value in radians
    MutateAASetPsi::MutateAASetPsi( const biol::Atom &FOLLOWING_N, const double DESIRED_PSI) :
      m_FollowingN( FOLLOWING_N),
      m_DesiredPsi( DESIRED_PSI)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAASetPsi
    MutateAASetPsi *MutateAASetPsi::Clone() const
    {
      return new MutateAASetPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAASetPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to set psi of a residue
    //! @param RESIDUE the residue for which the transformation matrix will be calculated
    //! @return the transformation matrix that is necessary to set psi of "RESIDUE" to "m_DesiredPsi"
    math::TransformationMatrix3D MutateAASetPsi::GetTransformationMatrix( const biol::AABase &RESIDUE) const
    {
      // return the transformation matrix that is necessary to set psi of "RESIDUE" to "m_DesiredPsi"
      return GetMutateAAPsi( RESIDUE).GetTransformationMatrix( RESIDUE);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param RESIDUE Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AABase> MutateAASetPsi::operator()( const biol::AABase &RESIDUE) const
    {
      // get the residue with psi mutated as desired in a MutateResult object
      const math::MutateResult< biol::AABase> mutated_psi( GetMutateAAPsi( RESIDUE).operator()( RESIDUE));

      // return "mutated_psi"
      return mutated_psi;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAASetPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FollowingN, ISTREAM);
      io::Serialize::Read( m_DesiredPsi, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateAASetPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FollowingN, OSTREAM, INDENT);
      io::Serialize::Write( m_DesiredPsi, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetMutateAAPsi provides the MutateAAPsi object necesssary to set psi as desired
    //! @param RESIDUE is the residue whose psi is going to be set
    //! @return the MutateAAPsi object necesssary to set psi as desired
    MutateAAPsi MutateAASetPsi::GetMutateAAPsi( const biol::AABase &RESIDUE) const
    {
      // calculate the rotation that is needed to set psi to "m_DesiredPsi"
      double current_psi( RESIDUE.Psi());
      if( !util::IsDefined( current_psi))
      {
        current_psi = RESIDUE.CalculatePsi( m_FollowingN);
      }
      const double rotation( m_DesiredPsi - current_psi);

      // message the needed rotation to set phi to "m_DesiredPhi"
      BCL_MessageDbg
      (
        "current Psi is " + util::Format()( RESIDUE.CalculatePsi( m_FollowingN))
        + " desired Psi is " + util::Format()( m_DesiredPsi)
        + " needed rotation is " + util::Format()( rotation)
      );

      // create a MutateAAPsi object using "rotation"
      const MutateAAPsi psi_mutate( -rotation);

      // return MutateAAPsi which can be used to set psi of "RESIDUE" as desired
      return psi_mutate;
    }

  } // namespace fold
} // namespace bcl
