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
#include "fold/bcl_fold_mutate_protein_model_move_aa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelMoveAA::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new MutateProteinModelMoveAA
        (
          util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> >(),
          math::Range< size_t>( 1, 1),
          biol::AASequenceFlexibility::s_NumberSequenceDirections,
          storage::Map< biol::SSType, size_t>(),
          util::ShPtr< math::MutateInterface< assemble::SSE> >(),
          ""
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from max number of residues to move
    //! @param SSE_LOCATOR locator that decides which SSE in the protein model to find an adjacent one to
    //! @param RESIDUES_TO_MOVE_RANGE range of number of residues to move
    //! @param SEQUENCE_SIDE side on sequence to move aas
    //! @param MIN_SSE_SIZE minimum SSE sizes to be allowed as result after moving aas
    //! @param SP_MUTATE_LOCATED_SSE optional mutate that is applied to the located and modified sse before it is added to the model
    //! @param SCHEME the scheme
    MutateProteinModelMoveAA::MutateProteinModelMoveAA
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
      const math::Range< size_t> &RESIDUES_TO_MOVE_RANGE,
      const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_SIDE,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE,
      const util::ShPtr< math::MutateInterface< assemble::SSE> > &SP_MUTATE_LOCATED_SSE,
      const std::string &SCHEME
    ) :
      m_SSELocator( SSE_LOCATOR),
      m_ResdiuesToMoveRange( RESIDUES_TO_MOVE_RANGE),
      m_SSESide( SEQUENCE_SIDE),
      m_MinSSESize( MIN_SSE_SIZE),
      m_MutateLocatedSSE( SP_MUTATE_LOCATED_SSE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelMoveAA
    MutateProteinModelMoveAA *MutateProteinModelMoveAA::Clone() const
    {
      return new MutateProteinModelMoveAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelMoveAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelMoveAA::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a protein model by shifting AAs between SSEs
    //! @param PROTEIN_MODEL
    //! @return MutateResult that results from mutating to the PROTEIN_MODEL
    math::MutateResult< assemble::ProteinModel>
    MutateProteinModelMoveAA::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // empty result
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // locate the sse in the protein model
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator->Locate( PROTEIN_MODEL));

      // was sse located
      if( !located_sse.IsDefined())
      {
        return empty_result;
      }

      // find adjacent sse
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > adjacent_sses( PROTEIN_MODEL.GetAdjacentSSEs( *located_sse));

      util::SiPtr< const assemble::SSE> sse1;
      util::SiPtr< const assemble::SSE> sse2;

      size_t nr_res( random::GetGlobalRandom().Random( m_ResdiuesToMoveRange));

      // update sses according to sequence direction
      if( m_SSESide == biol::AASequenceFlexibility::e_NTerminal && adjacent_sses.First().IsDefined())
      {
        sse1 = adjacent_sses.First();
        sse2 = located_sse;
        // check if length after move is still acceptable
        const storage::Map< biol::SSType, size_t>::const_iterator min_sse_size_itr( m_MinSSESize.Find( sse1->GetType()));
        if
        (
             min_sse_size_itr != m_MinSSESize.End()
          && min_sse_size_itr->second + nr_res > sse1->GetSize()
        )
        {
          BCL_MessageVrb
          (
            "cannot move " + util::Format()( nr_res) + " aas, since the resulting sse would be too short: " +
            sse1->GetIdentification() + "\t" + " min size: " + util::Format()( min_sse_size_itr->second)
          )
          return empty_result;
        }
      }
      else if( m_SSESide == biol::AASequenceFlexibility::e_CTerminal && adjacent_sses.Second().IsDefined())
      {
        sse1 = located_sse;
        sse2 = adjacent_sses.Second();
        const storage::Map< biol::SSType, size_t>::const_iterator min_sse_size_itr( m_MinSSESize.Find( sse2->GetType()));
        if
        (
             min_sse_size_itr != m_MinSSESize.End()
          && min_sse_size_itr->second + nr_res > sse2->GetSize()
        )
        {
          BCL_MessageVrb
          (
            "cannot move " + util::Format()( nr_res) + " aas, since the resulting sse would be too short: " +
            sse2->GetIdentification() + "\t" + " min size: " + util::Format()( min_sse_size_itr->second)
          )
          return empty_result;
        }
      }
      // no adjacent sse for the sequence side was found
      else
      {
        return empty_result;
      }

      // create two new sses
      const biol::AASequenceFlexibility::SequenceDirection move_direction
      (
        m_SSESide == biol::AASequenceFlexibility::e_NTerminal ? biol::AASequenceFlexibility::e_CTerminal : biol::AASequenceFlexibility::e_NTerminal
      );
      storage::VectorND< 2, util::ShPtr< assemble::SSE> > moved_sses( MoveAAs( *sse1, *sse2, nr_res, move_direction));

      // unsuccessful move
      if( !moved_sses.First().IsDefined() && !moved_sses.Second().IsDefined())
      {
        return empty_result;
      }

      // mutate sse if desired
      if( m_MutateLocatedSSE.IsDefined())
      {
        // mutate the second sse, which was the previously located
        if( m_SSESide == biol::AASequenceFlexibility::e_NTerminal && moved_sses.Second().IsDefined())
        {
          const math::MutateResult< assemble::SSE> result_sse_mutate( m_MutateLocatedSSE->operator ()( *moved_sses.Second()));
          if( result_sse_mutate.GetArgument().IsDefined())
          {
            moved_sses.Second() = result_sse_mutate.GetArgument();
          }
        }
        // mutate the first sse, which was the previously located
        else if( m_SSESide == biol::AASequenceFlexibility::e_CTerminal && moved_sses.First().IsDefined())
        {
          const math::MutateResult< assemble::SSE> result_sse_mutate( m_MutateLocatedSSE->operator ()( *moved_sses.First()));
          if( result_sse_mutate.GetArgument().IsDefined())
          {
            moved_sses.First() = result_sse_mutate.GetArgument();
          }
        }
      }

      // create new protein model
      util::ShPtr< assemble::ProteinModel> sp_new_model( PROTEIN_MODEL.Clone());
      sp_new_model->Remove( *sse1);
      sp_new_model->Remove( *sse2);

      // replace the old sses with the new sses
      if( moved_sses.First().IsDefined())
      {
        sp_new_model->Insert( moved_sses.First());
      }
      if( moved_sses.Second().IsDefined())
      {
        sp_new_model->Insert( moved_sses.Second());
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( sp_new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelMoveAA::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSELocator         , ISTREAM);
      io::Serialize::Read( m_ResdiuesToMoveRange, ISTREAM);
      io::Serialize::Read( m_Scheme             , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelMoveAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSELocator         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ResdiuesToMoveRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme             , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief move a given number of amino acids in a given direction between two adjacent sses
    //! @param SSE_LEFT n-terminal SSE
    //! @param SSE_RIGHT c-terminal SSE
    //! @param NUMBER_AAS the number of amino acids to move from one sse to another
    //! @param DIRECTION the sequence direstion - c-terminal move from left to right, n-terminal right to left
    //! @return 2 new SSE with the same type as the given ones - if the number of AAs was larger than the size, the
    //!         ShPtr will be undefined
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateProteinModelMoveAA::MoveAAs
    (
      const assemble::SSE &SSE_LEFT,
      const assemble::SSE &SSE_RIGHT,
      const size_t NUMBER_AAS,
      const biol::AASequenceFlexibility::SequenceDirection &DIRECTION
    )
    {
      storage::VectorND< 2, util::ShPtr< assemble::SSE> > result;

      BCL_MessageDbg( "sses before moving aa:\n" + SSE_LEFT.GetIdentification() + '\n' + SSE_RIGHT.GetIdentification());

      // true if the distance is not defined
      if( !biol::AABase::AreAminoAcidsPeptideBonded( *SSE_LEFT.GetLastAA(), *SSE_RIGHT.GetFirstAA(), true))
      {
        BCL_MessageCrt( "cannot move amino acids between SSEs that are not peptide bonded!");
        return result;
      }

      // copy of nr aas, in case they need to be corrected
      size_t nr_res( NUMBER_AAS);

      // move residues up in sequence to next SSE
      if( DIRECTION == biol::AASequenceFlexibility::e_CTerminal)
      {
        BCL_MessageDbg( "sses after moving aa:");
        if( nr_res < SSE_LEFT.GetSize())
        {
          result.First() = util::ShPtr< assemble::SSE>( new assemble::SSE( SSE_LEFT.SubSequence( 0, SSE_LEFT.GetSize() - nr_res), SSE_LEFT.GetType()));
          BCL_MessageDbg( result.First()->GetIdentification());
        }

        // reduce the number of residues if necessary
        nr_res = std::min( nr_res, SSE_LEFT.GetSize());

        // sequence for other sse
        biol::AASequence seq2( SSE_RIGHT);
        seq2.PrependSequence( SSE_LEFT.SubSequence( SSE_LEFT.GetSize() - nr_res, nr_res));

        // new sse
        result.Second() = util::ShPtr< assemble::SSE>( new assemble::SSE( seq2, SSE_RIGHT.GetType()));
        BCL_MessageDbg( result.Second()->GetIdentification());
      }
      // move residues down in sequence to previous sse
      else if( DIRECTION == biol::AASequenceFlexibility::e_NTerminal)
      {
        BCL_MessageDbg( "sses after moving aa:");
        if( nr_res < SSE_RIGHT.GetSize())
        {
          result.Second() = util::ShPtr< assemble::SSE>( new assemble::SSE( SSE_RIGHT.SubSequence( nr_res, SSE_RIGHT.GetSize() - nr_res), SSE_RIGHT.GetType()));
          BCL_MessageDbg( result.Second()->GetIdentification());
        }

        // reduce the number of residues if necessary
        nr_res = std::min( nr_res, SSE_RIGHT.GetSize());
        biol::AASequence seq1( SSE_LEFT);
        seq1.AppendSequence( SSE_RIGHT.SubSequence( 0, nr_res));

        // new sse
        result.First() = util::ShPtr< assemble::SSE>( new assemble::SSE( seq1, SSE_LEFT.GetType()));
        BCL_MessageDbg( result.First()->GetIdentification());
      }

      // end
      return result;
    }

  } // namespace fold
} // namespace bcl
