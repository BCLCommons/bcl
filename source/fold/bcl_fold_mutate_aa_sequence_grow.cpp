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
#include "fold/bcl_fold_mutate_aa_sequence_grow.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "fold/bcl_fold_mutation_residue.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateAASequenceGrow::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateAASequenceGrow())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateAASequenceGrow::MutateAASequenceGrow() :
      m_PhiPsiGenerator(),
      m_AnchorSequence()
    {
    }

    //! @brief constructor from member variables
    //! @param PHI_PSI_GENERATOR the method that will produce phi and psi angles as the sequence grows
    //! @param N_TERMINAL_AA_SEQUENCE the starting sequence which will be grown onto
    MutateAASequenceGrow::MutateAASequenceGrow
    (
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GENERATOR,
      const util::SiPtr< const biol::AASequence> &N_TERMINAL_AA_SEQUENCE
    ) :
      m_PhiPsiGenerator( PHI_PSI_GENERATOR),
      m_AnchorSequence( N_TERMINAL_AA_SEQUENCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateAASequenceInitializeCoordinates
    MutateAASequenceGrow *MutateAASequenceGrow::Clone() const
    {
      return new MutateAASequenceGrow( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateAASequenceGrow::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param SEQUENCE_TO_GROW Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< biol::AASequence> MutateAASequenceGrow::operator()
    (
      const biol::AASequence &SEQUENCE_TO_GROW
    ) const
    {
      const math::MutateResult< biol::AASequence> empty_result( util::ShPtr< biol::AASequence>(), *this);

      // make sure "m_AnchorSequence" is defined
      if( !m_AnchorSequence.IsDefined())
      {
        BCL_MessageCrt( "anchor sequence pointer is not defined");
        return empty_result;
      }

      // make sure that both of the sequences have at least one residue in them
      if( m_AnchorSequence->GetData().IsEmpty())
      {
        BCL_MessageCrt( "anchor sequence is empty");
        return empty_result;
      }
      if( SEQUENCE_TO_GROW.GetData().IsEmpty())
      {
        BCL_MessageCrt( "given sequence to grow is empty");
        return empty_result;
      }

      // create an sse of ideal strand sequence for the sequence to grow
      assemble::SSE new_sequence( SEQUENCE_TO_GROW, biol::GetSSTypes().STRAND);
      new_sequence.SetToIdealConformationAtOrigin();

      util::ShPtr< biol::AABase> previous_aa( m_AnchorSequence->GetLastAA());

      // iterate over new sequence and sample conformation
      for
      (
        biol::AASequence::iterator itr( new_sequence.Begin()), itr_end( new_sequence.End());
        itr != itr_end;
        previous_aa = *itr, ++itr
      )
      {
        util::ShPtr< biol::AABase> &current_aa( *itr);
        util::ShPtr< biol::AABase> next_aa;
        {
          biol::AASequence::iterator next_itr( itr + 1);
          if( next_itr != itr_end)
          {
            next_aa = *next_itr;
          }
        }

        MutationResidue current_mutation_residue( current_aa, previous_aa, next_aa);

        // generate phi psi
        const storage::VectorND< 2, double> new_phi_psi
        (
          m_PhiPsiGenerator->operator()( current_mutation_residue)
        );

        const storage::VectorND< 2, double> new_phi_psi_change( biol::AASequenceFlexibility::CalculatePhiPsiChange( new_sequence, current_aa->GetSeqID(), new_phi_psi));
        biol::AASequenceFlexibility::ChangePhiPsi( new_sequence, current_aa->GetSeqID(), new_phi_psi_change, biol::AASequenceFlexibility::e_CTerminal);

        // for first amino acid, use phi to connect to anchor sse
        if( itr == new_sequence.Begin())
        {
          BCL_MessageDbg( "attaching grow sequence on cterm of " + m_AnchorSequence->GetSequenceIdentification() + " with phi: " + util::Format()( new_phi_psi.First()));
          const math::TransformationMatrix3D transform( biol::AASequenceFactory::TransformationAppend( *m_AnchorSequence, *current_aa, new_phi_psi.First()));
          new_sequence.Transform( transform);
        }
      }

      util::ShPtr< biol::AASequence> sp_grown_seq( m_AnchorSequence->Clone());
      sp_grown_seq->AppendSequence( new_sequence);

      // return result
      return math::MutateResult< biol::AASequence>( sp_grown_seq, *this);
    }

    //! @brief grows the given sequence by connecting the cterminus to the anchor residue then growing nterminally
    //! @param SEQUENCE_TO_GROW sequence for which coordinates will be added
    //! @return MutateResult contains the sequence with newly assigned coordinates grown onto the anchor residue
    math::MutateResult< biol::AASequence>
    MutateAASequenceGrow::GrowTowardsNTerminus( const biol::AASequence &SEQUENCE_TO_GROW) const
    {
      const math::MutateResult< biol::AASequence> empty_result( util::ShPtr< biol::AASequence>(), *this);

      // make sure "m_AnchorSequence" is defined
      if( !m_AnchorSequence.IsDefined())
      {
        BCL_MessageCrt( "anchor sequence pointer is not defined");
        return empty_result;
      }

      // make sure that both of the sequences have at least one residue in them
      if( m_AnchorSequence->GetData().IsEmpty())
      {
        BCL_MessageCrt( "anchor sequence is empty");
        return empty_result;
      }
      if( SEQUENCE_TO_GROW.GetData().IsEmpty())
      {
        BCL_MessageCrt( "given sequence to grow is empty");
        return empty_result;
      }

      // create an sse of ideal strand sequence for the sequence to grow
      assemble::SSE new_sequence( SEQUENCE_TO_GROW, biol::GetSSTypes().STRAND);
      new_sequence.SetToIdealConformationAtOrigin();

      util::ShPtr< biol::AABase> previous_aa;

      // iterate over new sequence and sample conformation
      for
      (
        biol::AASequence::iterator itr( new_sequence.Begin()), itr_end( new_sequence.End());
        itr != itr_end;
        previous_aa = *itr, ++itr
      )
      {
        util::ShPtr< biol::AABase> &current_aa( *itr);
        util::ShPtr< biol::AABase> next_aa;
        {
          biol::AASequence::iterator next_itr( itr + 1);
          if( next_itr != itr_end)
          {
            next_aa = *next_itr;
          }
        }

        MutationResidue current_mutation_residue( current_aa, previous_aa, next_aa);

        // generate phi psi
        const storage::VectorND< 2, double> new_phi_psi
        (
          m_PhiPsiGenerator->operator()( current_mutation_residue)
        );

        const storage::VectorND< 2, double> new_phi_psi_change( biol::AASequenceFlexibility::CalculatePhiPsiChange( new_sequence, current_aa->GetSeqID(), new_phi_psi));
        biol::AASequenceFlexibility::ChangePhiPsi( new_sequence, current_aa->GetSeqID(), new_phi_psi_change, biol::AASequenceFlexibility::e_CTerminal);

        // for last amino acid, use phi to connect to anchor sse
        if( !next_aa.IsDefined())
        {
          BCL_MessageDbg( "attaching grow sequence on nterm of " + m_AnchorSequence->GetSequenceIdentification() + " with phi: " + util::Format()( new_phi_psi.First()));
          const math::TransformationMatrix3D transform( biol::AASequenceFactory::TransformationPrepend( *current_aa, *m_AnchorSequence, new_phi_psi.First()));
          new_sequence.Transform( transform);
        }
      }

      util::ShPtr< biol::AASequence> sp_grown_seq( new biol::AASequence( new_sequence));
      sp_grown_seq->AppendSequence( *m_AnchorSequence);

      // return result
      return math::MutateResult< biol::AASequence>( sp_grown_seq, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateAASequenceGrow::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PhiPsiGenerator, ISTREAM);
      io::Serialize::Read( m_AnchorSequence, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateAASequenceGrow::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PhiPsiGenerator, OSTREAM, INDENT);
      io::Serialize::Write( m_AnchorSequence, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
