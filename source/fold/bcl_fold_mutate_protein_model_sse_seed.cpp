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
#include "fold/bcl_fold_mutate_protein_model_sse_seed.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "fold/bcl_fold_mutate_protein_model_move_aa.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESeed::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSESeed())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSESeed::MutateProteinModelSSESeed()
    {
    }

    //! @brief constructor from a locator and other things
    //! @param SSE_LOCATOR locator that decides to which sse the seed sse is added
    //! @param SEED_LENGTH_RANGE number of residues to construct in seed sse, in addition to the cut
    //! @param DIRECTION the side of the located sse to which the seed sse is attached to
    //! @param SP_MUTATE_SEED optional mutate that is applied to the generated seed before it is added to the model
    //! @param CUT_IN_RANGE number of amino acids to cut into the located sse
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESeed::MutateProteinModelSSESeed
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
      const math::Range< size_t> &SEED_LENGTH_RANGE,
      const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
      const util::ShPtr< math::MutateInterface< assemble::SSE> > &SP_MUTATE_SEED,
      const math::Range< size_t> &CUT_IN_RANGE,
      const std::string &SCHEME
    ) :
      m_SSELocator( SSE_LOCATOR),
      m_SeedLengthRange( SEED_LENGTH_RANGE),
      m_Direction( DIRECTION),
      m_MutateSeed( SP_MUTATE_SEED),
      m_CutInRange( CUT_IN_RANGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSESeed
    MutateProteinModelSSESeed *MutateProteinModelSSESeed::Clone() const
    {
      return new MutateProteinModelSSESeed( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESeed::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelSSESeed::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an PROTEIN_MODEL and returning a protein model with an additional sse
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESeed::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate an sse and make copy
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator->Locate( PROTEIN_MODEL));

      // empty model for failed mutate
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // if sse cannot be located or the located sse is not helix or strand
      if( !located_sse.IsDefined() || located_sse->GetType() > biol::GetSSTypes().STRAND)
      {
        // warn user
        BCL_MessageVrb
        (
          "could not find sse to add seed sse to: " + util::Format()( *m_SSELocator)
        );

        // return empty result
        return empty_result;
      }

      // determine length of seed
      const size_t seed_length( random::GetGlobalRandom().Random< size_t>( m_SeedLengthRange));

      // create a reference to the sequence
      const biol::AASequence &full_sequence( *PROTEIN_MODEL.GetChain( located_sse->GetChainID())->GetSequence());

      biol::AASequence extended_sse( *located_sse);
      util::ShPtr< assemble::SSE> sp_new_sse;

      // seed length large 0
      if( seed_length > 0)
      {
        // grow on right side
        if( m_Direction == biol::AASequenceFlexibility::e_CTerminal)
        {
          // if there are not enough aas in the full sequence
          if( int( located_sse->GetLastAA()->GetSeqID() + seed_length) > full_sequence.GetLastAA()->GetSeqID())
          {
            BCL_MessageVrb( "full sequence does not have enough residues for seed sse on c term of: " + located_sse->GetIdentification());
            return empty_result;
          }

          // prepare seed
          assemble::SSE seed_seq( full_sequence.SubSequence( located_sse->GetLastAA()->GetSeqID(), seed_length), located_sse->GetType());
          seed_seq.SetToIdealConformationAtOrigin();

          // attach seed, so that it has correct peptide bond to previous SSE
          const math::TransformationMatrix3D trans( biol::AASequenceFactory::TransformationAppend( *located_sse, *seed_seq.GetFirstAA(), located_sse->GetType()->GetIdealPhi()));
          seed_seq.Transform( trans);
          sp_new_sse = util::ShPtr< assemble::SSE>( new assemble::SSE( seed_seq, biol::GetSSTypes().COIL));
        }
        // grow on left side
        else if( m_Direction == biol::AASequenceFlexibility::e_NTerminal)
        {
          // if there are not enough aas in the full sequence
          if( located_sse->GetFirstAA()->GetSeqID() < int( full_sequence.GetFirstAA()->GetSeqID() + seed_length))
          {
            BCL_MessageVrb( "full sequence does not have enough residues for seed sse on n term of: " + located_sse->GetIdentification());
            return empty_result;
          }

          // prepare seed
          assemble::SSE seed_seq( full_sequence.SubSequence( located_sse->GetFirstAA()->GetSeqID() - seed_length - 1, seed_length), located_sse->GetType());
          seed_seq.SetToIdealConformationAtOrigin();

          // calculate adequate tranformation for attachement
          const math::TransformationMatrix3D trans( biol::AASequenceFactory::TransformationPrepend( *seed_seq.GetLastAA(), *located_sse, located_sse->GetType()->GetIdealPhi()));
          seed_seq.Transform( trans);
          sp_new_sse = util::ShPtr< assemble::SSE>( new assemble::SSE( seed_seq, biol::GetSSTypes().COIL));
        }
        // unable to grow
        else
        {
          BCL_MessageVrb( "sequence direction incompatible");
          return empty_result;
        }
      }

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // cut in if desired
      const size_t cut_in_residues( random::GetGlobalRandom().Random( m_CutInRange));
      if( cut_in_residues >= located_sse->GetSize())
      {
        return empty_result;
      }

      if( cut_in_residues > 0)
      {
        BCL_MessageDbg
        (
          "cutting seed into located sse: " + located_sse->GetIdentification() + " by " + util::Format()( cut_in_residues) + " residues"
        );
        if( m_Direction == biol::AASequenceFlexibility::e_CTerminal)
        {
          if( sp_new_sse.IsDefined())
          {
            const storage::VectorND< 2, util::ShPtr< assemble::SSE> > moved_sses
            (
              MutateProteinModelMoveAA::MoveAAs( *located_sse, *sp_new_sse, cut_in_residues, m_Direction)
            );
            if( !moved_sses.First().IsDefined() && !moved_sses.Second().IsDefined())
            {
              return empty_result;
            }
            new_model->Remove( *located_sse);
            if( moved_sses.First().IsDefined())
            {
              new_model->Insert( moved_sses.First());
            }
            sp_new_sse = moved_sses.Second();
          }
          // split the sse
          else
          {
            const util::ShPtr< assemble::SSE> sp_sse_shorter( new assemble::SSE( located_sse->SubSequence( 0, located_sse->GetSize() - cut_in_residues), located_sse->GetType()));
            new_model->ReplaceWithOverlapping( sp_sse_shorter);
            sp_new_sse = util::ShPtr< assemble::SSE>( new assemble::SSE( located_sse->SubSequence( located_sse->GetSize() - cut_in_residues - 1, cut_in_residues), biol::GetSSTypes().COIL));
          }
        }
        else if( m_Direction == biol::AASequenceFlexibility::e_NTerminal)
        {
          if( sp_new_sse.IsDefined())
          {
            const storage::VectorND< 2, util::ShPtr< assemble::SSE> > moved_sses
            (
              MutateProteinModelMoveAA::MoveAAs( *sp_new_sse, *located_sse, cut_in_residues, m_Direction)
            );
            if( !moved_sses.First().IsDefined() && !moved_sses.Second().IsDefined())
            {
              return empty_result;
            }
            new_model->Remove( *located_sse);
            if( moved_sses.Second().IsDefined())
            {
              new_model->Insert( moved_sses.Second());
            }
            sp_new_sse = moved_sses.First();
          }
          else
          {
            const util::ShPtr< assemble::SSE> sp_sse_shorter( new assemble::SSE( located_sse->SubSequence( cut_in_residues, located_sse->GetSize() - cut_in_residues), located_sse->GetType()));
            new_model->ReplaceWithOverlapping( sp_sse_shorter);
            sp_new_sse = util::ShPtr< assemble::SSE>( new assemble::SSE( located_sse->SubSequence( 0, cut_in_residues), biol::GetSSTypes().COIL));
          }
        }
      }

      // mutate the seed if requested
      if( m_MutateSeed.IsDefined())
      {
        const math::MutateResult< assemble::SSE> result_seed_mutate( m_MutateSeed->operator ()( *sp_new_sse));
        if( result_seed_mutate.GetArgument().IsDefined())
        {
          sp_new_sse = result_seed_mutate.GetArgument();
        }
      }

      // insert the new SSE by replacing the old one
      if( !new_model->Insert( sp_new_sse))
      {
        BCL_MessageVrb( "seed sse overlaps with existing sse: " + sp_new_sse->GetIdentification());
        return empty_result;
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
    std::istream &MutateProteinModelSSESeed::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSELocator     , ISTREAM);
      io::Serialize::Read( m_SeedLengthRange, ISTREAM);
      io::Serialize::Read( m_Direction      , ISTREAM);
      io::Serialize::Read( m_MutateSeed     , ISTREAM);
      io::Serialize::Read( m_CutInRange     , ISTREAM);
      io::Serialize::Read( m_Scheme         , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelSSESeed::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSELocator     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeedLengthRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Direction      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MutateSeed     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CutInRange     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme         , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
