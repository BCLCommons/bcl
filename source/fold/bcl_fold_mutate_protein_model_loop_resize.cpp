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
#include "fold/bcl_fold_mutate_protein_model_loop_resize.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelLoopResize::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelLoopResize())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelLoopResize::MutateProteinModelLoopResize() :
      MutateProteinModelSSEResize(),
      m_MaxSSESizes(),
      m_PhiPsiGenerator()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param LOCATOR locator that decides which sse in the protein model to mutate
    //! @param EXTEND_PROBABILITY probability for extending (1.0 - EXTEND_PROBABILITY for shrinking)
    //! @param LENGTH_CHANGE_RANGE range of number of residues that are to be added or removed in one mutate to one end
    //! @param MIN_SSE_SIZES minimum SSE sizes to be allowed when shrinking
    //! @param MAX_SSE_SIZES minimum SSE sizes to be allowed when extending
    //! @param GROWING_DIRECTION the side of the sse that changes will occur on
    //! @param PHI_PSI_GENERATOR the method that will be used in order to generate phi and psi angles if extending
    //! @param GROWING_DIRECTION the side of the sse that changes will occur on
    //! @param DISALLOW_OVERLAP true if extension of sse into neighboring sse is not allowed-else (replaces both sses)
    //! @param SCHEME the scheme of this mutate
    MutateProteinModelLoopResize::MutateProteinModelLoopResize
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
      const double &EXTEND_PROBABILITY,
      const math::Range< size_t> &LENGTH_CHANGE_RANGE,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
      const storage::Map< biol::SSType, size_t> &MAX_SSE_SIZES,
      const biol::AASequenceFlexibility::SequenceDirection &GROWING_DIRECTION,
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GENERATOR,
      const bool DISALLOW_OVERLAP,
      const std::string &SCHEME
    ) :
      MutateProteinModelSSEResize
      (
        LOCATOR, EXTEND_PROBABILITY, LENGTH_CHANGE_RANGE, GROWING_DIRECTION, false, MIN_SSE_SIZES, SCHEME
      ),
      m_MaxSSESizes( MAX_SSE_SIZES),
      m_PhiPsiGenerator( PHI_PSI_GENERATOR),
      m_DisallowOverlap( DISALLOW_OVERLAP)
    {
      BCL_Assert
      (
        m_Side == biol::AASequenceFlexibility::e_NTerminal ||
        m_Side == biol::AASequenceFlexibility::e_CTerminal,
        "Can only grow loops in either "
        + biol::AASequenceFlexibility::DirectionEnum( biol::AASequenceFlexibility::e_NTerminal).GetString()
        + " or "
        + biol::AASequenceFlexibility::DirectionEnum( biol::AASequenceFlexibility::e_CTerminal).GetString()
        + " direction but specified as "
        + biol::AASequenceFlexibility::DirectionEnum( m_Side).GetString()
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelLoopResize
    MutateProteinModelLoopResize *MutateProteinModelLoopResize::Clone() const
    {
      return new MutateProteinModelLoopResize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelLoopResize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelLoopResize::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model pointer
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // locate sse to grow
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator->Locate( PROTEIN_MODEL));

      // if sse cannot be located
      if( !located_sse.IsDefined())
      {
        // return empty result
        return empty_result;
      }

      // report selected sse to be extended/shrunk
      BCL_MessageDbg( "selected sse to be extended/shrunk: " + located_sse->GetIdentification());

      // determine whether to extend or shrink
      const bool extend( random::GetGlobalRandom().Double() <= m_ExtendProbability);

      // determine step size
      const size_t step_size( random::GetGlobalRandom().Random< size_t>( m_LengthChangeRange));

      // make a ShPtr to an SSE
      util::ShPtr< assemble::SSE> new_sse;

      // true if sse will shrink
      if( !extend)
      {
        // find the min sse size
        const storage::Map< biol::SSType, size_t>::const_iterator itr( m_MinSSESizes.Find( located_sse->GetType()));
        const size_t min_sse_size( itr == m_MinSSESizes.End() ? 0 : itr->second);

        // is sse long enough
        if( step_size > located_sse->GetSize())
        {
          return empty_result;
        }

        // calculate the size of SSE after the shrinking
        const size_t length_after_shrink( located_sse->GetSize() - step_size);

        // is the sse long enough after shrinking
        if( length_after_shrink < min_sse_size)
        {
          return empty_result;
        }

        // true if the n terminus of sse is being shrunk
        switch( m_Side)
        {
          case biol::AASequenceFlexibility::e_NTerminal:
          {
            new_sse =
              util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( located_sse->SubSequence( step_size, length_after_shrink), located_sse->GetType())
              );
            break;
          }
          case biol::AASequenceFlexibility::e_CTerminal:
          {
            new_sse =
              util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( located_sse->SubSequence( 0, length_after_shrink), located_sse->GetType())
              );
            break;
          }
          default:
            return empty_result;
        }
      }

      // sse will be extended
      else
      {
        // find the max sse size
        const storage::Map< biol::SSType, size_t>::const_iterator itr( m_MaxSSESizes.Find( located_sse->GetType()));
        const size_t max_sse_size( itr == m_MaxSSESizes.End() ? util::GetUndefined< size_t>() : itr->second);

        // determine how many residues to extend
        const size_t extend_size( random::GetGlobalRandom().Random< size_t>( m_LengthChangeRange));

        // sse too long
        if( located_sse->GetSize() + extend_size > max_sse_size)
        {
          return empty_result;
        }

        // check for overlap
        if( m_DisallowOverlap && WillOverlap( PROTEIN_MODEL, extend_size, *located_sse, m_Side))
        {
          // return empty result
          return empty_result;
        }

        // create a copy of the located SSE and store it in new_sse
        new_sse = util::ShPtr< assemble::SSE>( located_sse->Clone());

        // create a reference to the sequence
        const biol::AASequence &full_sequence( *PROTEIN_MODEL.GetChain( located_sse->GetChainID())->GetSequence());

        // true if extend on n terminal side
        if( m_Side == biol::AASequenceFlexibility::e_NTerminal)
        {
          const biol::AASequence subsequence( full_sequence.SubSequence( located_sse->GetFirstAA()->GetSeqID() - extend_size - 1, extend_size));
          BCL_MessageDbg( "subsequence size for prepending " + util::Format()( subsequence.GetSize()));
          new_sse = PrependSequence( *located_sse, subsequence);
        }

        // true if extend on c terminal side
        else if( m_Side == biol::AASequenceFlexibility::e_CTerminal)
        {
          const biol::AASequence subsequence( full_sequence.SubSequence( located_sse->GetLastAA()->GetSeqID(), extend_size));
          BCL_MessageDbg( "subsequence size for appending " + util::Format()( subsequence.GetSize()));
          BCL_MessageDbg( "subsequence appending " + util::Format()( subsequence.GetSequenceIdentification()));
          new_sse = AppendSequence( *located_sse, subsequence);
        }
      }

      // print out sse after extend/shrink
      BCL_MessageDbg( "sse after extended/shrunk: " + new_sse->GetIdentification());

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // insert the new SSE by replacing the old one
      if( new_sse->GetSize() > 0)
      {
        new_model->ReplaceWithOverlapping( new_sse);
      }
      else
      {
        new_model->Remove( *located_sse);
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelLoopResize::Read( std::istream &ISTREAM)
    {
      // read members
      MutateProteinModelSSEResize::Read( ISTREAM);
      io::Serialize::Read( m_MaxSSESizes    , ISTREAM);
      io::Serialize::Read( m_PhiPsiGenerator, ISTREAM);
      io::Serialize::Read( m_DisallowOverlap, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelLoopResize::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      MutateProteinModelSSEResize::Write(      OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxSSESizes    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PhiPsiGenerator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DisallowOverlap, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief prepends an aa sequence to the nterminal side of an sse
    //! @param SSE the sse to which a sequence will be prepended
    //! @param SEQUENCE the sequence that will be prepended to the sse
    //! @return sse which has been created by prepending the provided sequence onto n-terminus of the sse
    util::ShPtr< assemble::SSE>
    MutateProteinModelLoopResize::PrependSequence( const assemble::SSE &SSE, const biol::AASequence &SEQUENCE) const
    {
      util::ShPtr< assemble::SSE> new_sse( SSE.Clone());

      // iterate over the sequence to prepend it
      for
      (
        biol::AASequence::const_reverse_iterator aa_itr( SEQUENCE.ReverseBegin()), aa_itr_end( SEQUENCE.ReverseEnd());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // copy aa and set to strand ideal conformation
        util::ShPtr< biol::AABase> new_aa( ( *aa_itr)->Clone());
        new_aa->SetToIdealConformation( biol::GetSSTypes().STRAND, math::TransformationMatrix3D());

        // get the phi psi angles
        const MutationResidue mutation_resi( new_aa, util::ShPtr< biol::AABase>(), util::ShPtr< biol::AABase>());
        const storage::VectorND< 2, double> phi_psi( m_PhiPsiGenerator->operator ()( mutation_resi));
        BCL_MessageDbg( "PrependSequence phi psi is " + util::Format()( phi_psi));

        // prepend the current aa to the sse
        BCL_MessageDbg( "Prepending residue " + new_aa->GetIdentification() + " onto sse " + new_sse->GetIdentification());
        biol::AASequenceFactory::PrependAA( *new_aa, *new_sse, phi_psi.First(), phi_psi.Second());
        BCL_MessageDbg( "prepended " + new_aa->GetIdentification());
      }

      return new_sse;
    }

    //! @brief appends an aa sequence to the cterminal side of an sse
    //! @param SSE the sse to which a sequence will be appended
    //! @param SEQUENCE the sequence that will be appended to the sse
    //! @return sse which has been created by appending the provided sequence onto c-terminus of the sse
    util::ShPtr< assemble::SSE>
    MutateProteinModelLoopResize::AppendSequence( const assemble::SSE &SSE, const biol::AASequence &SEQUENCE) const
    {
      util::ShPtr< assemble::SSE> new_sse( SSE.Clone());

      // iterate over the sequence to append it
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // copy aa and set to strand ideal conformation
        util::ShPtr< biol::AABase> new_aa( ( *aa_itr)->Clone());
        new_aa->SetToIdealConformation( biol::GetSSTypes().STRAND, math::TransformationMatrix3D());

        // get the phi psi angles
        const MutationResidue mutation_resi( new_aa, util::ShPtr< biol::AABase>(), util::ShPtr< biol::AABase>());
        const storage::VectorND< 2, double> phi_psi( m_PhiPsiGenerator->operator ()( mutation_resi));
        BCL_MessageDbg( "AppendSequence phi psi is " + util::Format()( phi_psi));
        BCL_MessageDbg( "Appending residue " + new_aa->GetIdentification() + " onto sse " + new_sse->GetIdentification());

        // prepend the current aa to the sse
        biol::AASequenceFactory::AppendAA( *new_sse, *new_aa, phi_psi.First(), phi_psi.Second());
        BCL_MessageDbg( "appended " + new_aa->GetIdentification());
      }

      return new_sse;
    }

    bool MutateProteinModelLoopResize::WillOverlap
    (
      const assemble::ProteinModel &MODEL, const size_t EXTENSION_AMOUNT, const assemble::SSE &CURRENT_SSE,
      const biol::AASequenceFlexibility::SequenceDirection &EXTENSION_DIRECTION
    )
    {
      // make sure that extending won't overlap with neighboring sse
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > neighbor_sses
      (
        MODEL.GetNeighborSSEs( CURRENT_SSE)
      );

      // true if the left neighbor is defined and extension will occur on the nterminal side
      if( neighbor_sses.First().IsDefined() && EXTENSION_DIRECTION == biol::AASequenceFlexibility::e_NTerminal)
      {
        const assemble::SSE &left_sse( *neighbor_sses.First());

        const int current_sse_seq_id( CURRENT_SSE.GetFirstAA()->GetSeqID());
        const int left_sse_seq_id( left_sse.GetLastAA()->GetSeqID());

        // true if the sses will overlap
        if( left_sse_seq_id + int( EXTENSION_AMOUNT) >= current_sse_seq_id)
        {
          BCL_MessageDbg
          (
            "extending " + CURRENT_SSE.GetIdentification() + " by "
            + util::Format()( EXTENSION_AMOUNT) + " will overlap with sse " + left_sse.GetIdentification()
          );
          return true;
        }
      }
      else if( EXTENSION_DIRECTION == biol::AASequenceFlexibility::e_NTerminal)//< ntermini
      {
        // not enough residues prior in sequence to add
        if( int( CURRENT_SSE.GetFirstAA()->GetSeqID() - int( EXTENSION_AMOUNT)) < int( 1))
        {
          BCL_MessageDbg
          (
            "extending " + CURRENT_SSE.GetIdentification() + " by "
            + util::Format()( EXTENSION_AMOUNT) + " will run out of sequence to add"
          );
          return true;
        }
      }

      // true if the right neighbor is defined and extension will occur on the cterminal side
      if( neighbor_sses.Second().IsDefined() && EXTENSION_DIRECTION == biol::AASequenceFlexibility::e_CTerminal)
      {
        const assemble::SSE &right_sse( *neighbor_sses.Second());

        const int current_sse_seq_id( CURRENT_SSE.GetLastAA()->GetSeqID());
        const int right_sse_seq_id( right_sse.GetFirstAA()->GetSeqID());

        // true if the sses will overlap
        if( current_sse_seq_id + int( EXTENSION_AMOUNT) >= right_sse_seq_id)
        {
          BCL_MessageDbg
          (
            "extending " + CURRENT_SSE.GetIdentification() + " by "
            + util::Format()( EXTENSION_AMOUNT) + " will overlap with sse " + right_sse.GetIdentification()
          );
          return true;
        }
      }
      else if( EXTENSION_DIRECTION == biol::AASequenceFlexibility::e_CTerminal) //< ctermini
      {
        // not enough residues at end of sequence to add
        const int cterm_seq_id( MODEL.GetChain( CURRENT_SSE.GetChainID())->GetSequence()->GetLastAA()->GetSeqID());
        if
        (
          int( cterm_seq_id - CURRENT_SSE.GetLastAA()->GetSeqID()) < int( EXTENSION_AMOUNT)
        )
        {
          BCL_MessageDbg
          (
            "extending " + CURRENT_SSE.GetIdentification() + " by "
            + util::Format()( EXTENSION_AMOUNT) + " will run out of sequence to add"
          );
          return true;
        }
      }

      return false;
    }

  } // namespace fold
} // namespace bcl
