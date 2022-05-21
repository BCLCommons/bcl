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
#include "fold/bcl_fold_mutate_sse_bend_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateSSEBendRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSSEBendRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSSEBendRandom::MutateSSEBendRandom() :
      m_PhiChangeRange( math::Angle::Radian( -20.0), math::Angle::Radian( 20.0)),
      m_PsiChangeRange( math::Angle::Radian( -20.0), math::Angle::Radian( 20.0)),
      m_BendingDirection( biol::AASequenceFlexibility::e_Bidirectional),
      m_Scheme( GetStaticClassName< MutateSSEBendRandom>())
    {
    }

    //! @brief constructor phi/psi change ranges and a scheme
    //! @param PHI_CHANGE_RANGE range for phi changes
    //! @param PSI_CHANGE_RANGE range for psi changes
    //! @param BENDING_DIRECTION direction the phi/psi changes should be propagated towards
    //! @param SCHEME Scheme to be used
    MutateSSEBendRandom::MutateSSEBendRandom
    (
      const math::Range< double> &PHI_CHANGE_RANGE,
      const math::Range< double> &PSI_CHANGE_RANGE,
      const biol::AASequenceFlexibility::SequenceDirection BENDING_DIRECTION,
      const std::string &SCHEME
    ) :
      m_PhiChangeRange( PHI_CHANGE_RANGE),
      m_PsiChangeRange( PSI_CHANGE_RANGE),
      m_BendingDirection( BENDING_DIRECTION),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateSSEBendRandom *MutateSSEBendRandom::Clone() const
    {
      return new MutateSSEBendRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSSEBendRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get phi change range
    //! @return phi change range
    const math::Range< double> &MutateSSEBendRandom::GetPhiChangeRange() const
    {
      return m_PhiChangeRange;
    }

    //! @brief get phi change range
    //! @return phi change range
    const math::Range< double> &MutateSSEBendRandom::GetPsiChangeRange() const
    {
      return m_PsiChangeRange;
    }

    //! @brief return bending direction
    //! @return bending direction
    biol::AASequenceFlexibility::SequenceDirection MutateSSEBendRandom::GetBendingDirection() const
    {
      return m_BendingDirection;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that changes phi/psi angles of a random amino acid in the given SSE
    //! @param THIS_SSE SSE interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::SSE> MutateSSEBendRandom::operator()
    (
      const assemble::SSE &THIS_SSE
    ) const
    {
      // make a copy of the selected SSE
      util::ShPtr< assemble::SSE> new_sse( THIS_SSE.Clone());

      // decide on a random amino acid to pick
      const int seq_id
      (
        random::GetGlobalRandom().Integer
        (
          math::Range< int>( new_sse->GetFirstAA()->GetSeqID(), new_sse->GetLastAA()->GetSeqID())
        )
      );

      // calculate the phi psi change
      const storage::VectorND< 2, double> phi_psi_change
      (
        random::GetGlobalRandom().Double( m_PhiChangeRange), random::GetGlobalRandom().Double( m_PsiChangeRange)
      );

      // use the AASequence flexibility to bend the sequence
      biol::AASequenceFlexibility::ChangePhiPsi( *new_sse, seq_id, phi_psi_change, m_BendingDirection);

      // set the geometries
      new_sse->SetGeometry();

      // end
      return math::MutateResult< assemble::SSE>( new_sse, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSSEBendRandom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PhiChangeRange, ISTREAM);
      io::Serialize::Read( m_PsiChangeRange, ISTREAM);
      io::Serialize::Read( m_BendingDirection, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateSSEBendRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PhiChangeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PsiChangeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BendingDirection, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
