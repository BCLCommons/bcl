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
#include "fold/bcl_fold_mutate_protein_model_grow_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_mutate_aa_sequence_grow.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelGrowSSE::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelGrowSSE())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelGrowSSE::MutateProteinModelGrowSSE() :
      m_LocatorSSE(),
      m_PhiPsiGenerator(),
      m_AnchorAA(),
      m_GrowingDirection()
    {
    }

    //! @brief constructor taking parameters
    //! @param LOCATOR_SSE the locator which will be used to find the sse which is going to be grown
    //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
    //! @param ANCHOR_AA locator to residue to which the n-terminus of the growing sse will be anchored (connected to)
    //! @param GROWING_DIRECTION growing direction of the sse
    MutateProteinModelGrowSSE::MutateProteinModelGrowSSE
    (
      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface>
      > &LOCATOR_SSE,
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
      > &PHI_PSI_GEN,
      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel>
      > &ANCHOR_AA,
      const biol::AASequenceFlexibility::SequenceDirection GROWING_DIRECTION
    ) :
      m_LocatorSSE( LOCATOR_SSE),
      m_PhiPsiGenerator( PHI_PSI_GEN),
      m_AnchorAA( ANCHOR_AA),
      m_GrowingDirection( GROWING_DIRECTION)
    {
      BCL_Assert
      (
        GROWING_DIRECTION == biol::AASequenceFlexibility::e_NTerminal ||
        GROWING_DIRECTION == biol::AASequenceFlexibility::e_CTerminal,
        "SSEs can be grown only towards NTerminal or CTerminal, direction provided was: " +
        biol::AASequenceFlexibility::GetSequenceDirectionName( m_GrowingDirection)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelGrowSSE
    MutateProteinModelGrowSSE *MutateProteinModelGrowSSE::Clone() const
    {
      return new MutateProteinModelGrowSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelGrowSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelGrowSSE::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // true if the sse should be grown in the c terminal direction
      if( m_GrowingDirection == biol::AASequenceFlexibility::e_CTerminal)
      {
        return GrowTowardsCTerminus( PROTEIN_MODEL);
      }

      // at this point, the sse should be grown in the n nterminal direction
      return GrowTowardsNTerminus( PROTEIN_MODEL);
    }

    //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelGrowSSE::GrowTowardsCTerminus
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      const util::SiPtr< const assemble::SSE> sse( m_LocatorSSE->Locate( PROTEIN_MODEL));
      BCL_Assert( sse.IsDefined(), "sse could not be found");

      // sequence will be built up as the sse is grown
      biol::AASequence new_sse_sequence( util::ShPtrVector< biol::AABase>(), sse->GetChainID());

      // get the anchor aa
      const util::SiPtr< const biol::AABase> anchor_aa( m_AnchorAA->Locate( PROTEIN_MODEL));

      // assert the anchor residue could be found
      BCL_Assert( anchor_aa.IsDefined(), "could not locate residue with locator " + util::Format()( m_AnchorAA));

      // create sequence out of the anchor aa
      const biol::AASequence anchor_seq( util::ShPtrVector< biol::AABase>( 1, util::ShPtr< biol::AABase>( anchor_aa->Clone())), anchor_aa->GetChainID());

      // create a MutateAASequenceGrow to be used for growing the located SSE
      const MutateAASequenceGrow grower( m_PhiPsiGenerator, util::SiPtr< const biol::AASequence>( anchor_seq));

      // grow the sequence of interest
      const math::MutateResult< biol::AASequence> growing_result( grower( *sse));
      const util::ShPtr< biol::AASequence> &grown_sequence( growing_result.GetArgument());

      BCL_Assert( grown_sequence.IsDefined(), "grown sequence is not defined");

      // new protein model
      util::ShPtr< assemble::ProteinModel> new_protein_model( PROTEIN_MODEL.Clone());

      // create a new final SSE which contains all residues except the first one
      // since the first residue is actually the last residue of the anchor SSE
      util::ShPtr< assemble::SSE> final_new_sse
      (
        new assemble::SSE( grown_sequence->SubSequence( 1, grown_sequence->GetSize() - 1), biol::GetSSTypes().COIL)
      );

      // locate the anchor SSE from the model
      util::SiPtr< const assemble::SSE> anchor_sse
      (
        assemble::LocatorAA( anchor_aa->GetChainID(), anchor_aa->GetSeqID()).LocateSSE( PROTEIN_MODEL)
      );
      BCL_Assert( anchor_sse.IsDefined(), "Anchor sse was not found which contained " + anchor_aa->GetIdentification());
      // make a new SSE and update its last member with the first amino acid from grown sequence
      util::ShPtr< assemble::SSE> new_anchor_sse( anchor_sse->Clone());
      new_anchor_sse->GetData().LastElement() = grown_sequence->GetFirstAA();
      // replace it
      new_protein_model->ReplaceWithOverlapping( new_anchor_sse);
      // now we can add the final_new_sse
      new_protein_model->Replace( final_new_sse);
      // create mutate result from "replaced_sse_model"
      math::MutateResult< assemble::ProteinModel> mutate_result( new_protein_model, *this);

      // return
      return mutate_result;
    }

    //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelGrowSSE::GrowTowardsNTerminus
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      const util::SiPtr< const assemble::SSE> sse( m_LocatorSSE->Locate( PROTEIN_MODEL));
      BCL_Assert( sse.IsDefined(), "sse could not be found");

      // sequence will be built up as the sse is grown
      biol::AASequence new_sse_sequence( util::ShPtrVector< biol::AABase>(), sse->GetChainID());

      // get the anchor aa
      const util::SiPtr< const biol::AABase> anchor_aa( m_AnchorAA->Locate( PROTEIN_MODEL));

      // assert the anchor residue could be found
      BCL_Assert( anchor_aa.IsDefined(), "could not locate residue with locator " + util::Format()( m_AnchorAA));

      // create sequence out of the anchor aa
      const biol::AASequence anchor_seq
      (
        util::ShPtrVector< biol::AABase>( 1, util::ShPtr< biol::AABase>( anchor_aa->Clone())), anchor_aa->GetChainID()
      );

      // create a MutateAASequenceGrow to be used for growing the located SSE
      const MutateAASequenceGrow grower( m_PhiPsiGenerator, util::SiPtr< const biol::AASequence>( anchor_seq));

      // grow the sequence of interest
      const math::MutateResult< biol::AASequence> growing_result( grower.GrowTowardsNTerminus( *sse));

      const util::ShPtr< biol::AASequence> &grown_sequence( growing_result.GetArgument());

      BCL_Assert( grown_sequence.IsDefined(), "grown sequence is not defined");

      // new protein model
      util::ShPtr< assemble::ProteinModel> new_protein_model( PROTEIN_MODEL.Clone());

      // create a new final SSE which contains all residues except the last one
      // since the last residue is actually the first residue of the anchor SSE
      util::ShPtr< assemble::SSE> final_new_sse
      (
        new assemble::SSE( grown_sequence->SubSequence( 0, grown_sequence->GetSize() - 1), biol::GetSSTypes().COIL)
      );

      // locate the anchor SSE from the model
      util::SiPtr< const assemble::SSE> anchor_sse
      (
        assemble::LocatorAA( anchor_aa->GetChainID(), anchor_aa->GetSeqID()).LocateSSE( PROTEIN_MODEL)
      );
      BCL_Assert( anchor_sse.IsDefined(), "Anchor sse was not found which contained " + anchor_aa->GetIdentification());
      // make a new SSE and update its first member with the last amino acid from grown sequence
      util::ShPtr< assemble::SSE> new_anchor_sse( anchor_sse->Clone());
      new_anchor_sse->GetData().FirstElement() = grown_sequence->GetLastAA();
      // replace it
      new_protein_model->ReplaceWithOverlapping( new_anchor_sse);
      // now we can add the final_new_sse
      new_protein_model->Replace( final_new_sse);
      // create mutate result from "replaced_sse_model"
      math::MutateResult< assemble::ProteinModel> mutate_result( new_protein_model, *this);

      // return
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelGrowSSE::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LocatorSSE,       ISTREAM);
      io::Serialize::Read( m_PhiPsiGenerator,  ISTREAM);
      io::Serialize::Read( m_AnchorAA,         ISTREAM);
      io::Serialize::Read( m_GrowingDirection, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelGrowSSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LocatorSSE,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PhiPsiGenerator,  OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AnchorAA,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GrowingDirection, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief static function that grows all the missing coil regions in the given ProteinModel
    //!        splits the coil and grows first part N to C and the second part C to N
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
    //! @return ShPtr to ProteinModel with the loops grown
    util::ShPtr< assemble::ProteinModel> MutateProteinModelGrowSSE::GrowAllMissingCoilsBidirectional
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN
    )
    {
      // create a model to keep track as we grow loops
      util::ShPtr< assemble::ProteinModel> sp_model( PROTEIN_MODEL.Clone());
      sp_model->AddLoops( true, false);

      // iterate over chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the SSEs for this chain
        const util::SiPtrVector< const assemble::SSE> all_sses( ( *chain_itr)->GetSSEs());

        // iterate over all coils
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // create reference
          const assemble::SSE &this_sse( **sse_itr);

          // check if this coil has undefined coordinates
          if( this_sse.GetType() != biol::GetSSTypes().COIL || this_sse.HasDefinedCoordinates())
          {
            continue;
          }

          // true if at nterminal sse i.e. very first sse in chain
          if( sse_itr == all_sses.Begin())
          {
            BCL_MessageDbg( "growing nterminal sse " + this_sse.GetIdentification());
            // grow sse
            GrowSSE( **sse_itr, sp_model, biol::AASequenceFlexibility::e_NTerminal, PHI_PSI_GEN);
            continue;
          }
          // true if at cterminal sse i.e. last sse in the chain
          else if( sse_itr == --all_sses.End())
          {
            BCL_MessageDbg( "growing cterminal sse " + this_sse.GetIdentification());
            // grow sse
            GrowSSE( **sse_itr, sp_model, biol::AASequenceFlexibility::e_CTerminal, PHI_PSI_GEN);
            continue;
          }

          BCL_Assert( sp_model->Remove( this_sse), "could not remove " + this_sse.GetIdentification());

          // create n terminal segment
          const size_t nterm_segment_size( this_sse.GetSize() - this_sse.GetSize() / 2);
          if( nterm_segment_size > 0)
          {
            const biol::AASequence new_nterm_seq( this_sse.SubSequence( 0, nterm_segment_size));
            util::ShPtr< assemble::SSE> nterm_sse( new assemble::SSE( new_nterm_seq, this_sse.GetType()));
            BCL_Assert( sp_model->Insert( nterm_sse), "could not insert " + nterm_sse->GetIdentification());
            GrowSSE( *nterm_sse, sp_model, biol::AASequenceFlexibility::e_CTerminal, PHI_PSI_GEN);
          }

          // create c terminal segment
          const size_t cterm_segment_size( this_sse.GetSize() - nterm_segment_size);
          if( cterm_segment_size > 0)
          {
            const biol::AASequence new_cterm_seq( this_sse.SubSequence( nterm_segment_size, cterm_segment_size));
            util::ShPtr< assemble::SSE> cterm_sse( new assemble::SSE( new_cterm_seq, this_sse.GetType()));
            BCL_Assert( sp_model->Insert( cterm_sse), "could not insert " + cterm_sse->GetIdentification());
            GrowSSE( *cterm_sse, sp_model, biol::AASequenceFlexibility::e_NTerminal, PHI_PSI_GEN);
          }
        }
      }

      // end
      return sp_model;
    }

    //! @brief grows coordinates for specified sses
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
    //! @param LOCATOR_SSES the specific sses that will be grown
    //! @return ptr to protein model that has the grown sses as specifed
    util::ShPtr< assemble::ProteinModel> MutateProteinModelGrowSSE::GrowSpecifiedCoilsBidirectional
    (
      const assemble::ProteinModel &START_MODEL,
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN,
      const util::SiPtrVector< const assemble::LocatorSSE> &LOCATOR_SSES
    )
    {
      // create a model to keep track as we grow loops
      util::ShPtr< assemble::ProteinModel> sp_model( START_MODEL.Clone());
      sp_model->AddLoops( true, false);

      // iterate over the sse locators
      for
      (
        util::SiPtrVector< const assemble::LocatorSSE>::const_iterator
        locator_itr( LOCATOR_SSES.Begin()), locator_itr_end( LOCATOR_SSES.End());
        locator_itr != locator_itr_end; ++locator_itr
      )
      {
        util::SiPtr< const assemble::SSE> located_sse( ( *locator_itr)->Locate( START_MODEL));
        BCL_Assert( located_sse.IsDefined(), "could not locate sse " + ( *locator_itr)->GetIdentification());

        // create reference
        const assemble::SSE &this_sse( *located_sse);

        // true if at nterminal sse i.e. very first sse in chain
        const util::SiPtrVector< const assemble::SSE> chain_sses
        (
          START_MODEL.GetChain( ( *locator_itr)->GetChainID())->GetSSEs()
        );
        if( located_sse == chain_sses.FirstElement())
        {
          BCL_MessageDbg( "growing nterminal sse " + this_sse.GetIdentification());
          // grow sse and get locator
          GrowSSE( this_sse, sp_model, biol::AASequenceFlexibility::e_NTerminal, PHI_PSI_GEN);
          continue;
        }
        // true if at cterminal sse i.e. last sse in the chain
        else if( located_sse == chain_sses.LastElement())
        {
          BCL_MessageDbg( "growing cterminal sse " + this_sse.GetIdentification());
          // grow sse and get locator
          GrowSSE( this_sse, sp_model, biol::AASequenceFlexibility::e_CTerminal, PHI_PSI_GEN);
          continue;
        }

        // size of n and c terminal segments
        const size_t nterm_segment_size( this_sse.GetSize() - this_sse.GetSize() / 2);
        const size_t cterm_segment_size( this_sse.GetSize() - nterm_segment_size);

        // construct the n and c terminal sub sequences
        const biol::AASequence new_nterm_seq( this_sse.SubSequence( 0, nterm_segment_size));
        const biol::AASequence new_cterm_seq( this_sse.SubSequence( nterm_segment_size, cterm_segment_size));
        BCL_Assert( ( new_nterm_seq.GetSize() + new_cterm_seq.GetSize()) == this_sse.GetSize(), "size wrong");

        // construct n and c terminal sses and put them into the model after removing the sse they were created from
        util::ShPtr< assemble::SSE> nterm_sse( new assemble::SSE( new_nterm_seq, this_sse.GetType()));
        util::ShPtr< assemble::SSE> cterm_sse( new assemble::SSE( new_cterm_seq, this_sse.GetType()));
        BCL_Assert( sp_model->Remove( this_sse), "could not remove " + this_sse.GetIdentification());
        BCL_Assert( sp_model->Insert( nterm_sse), "could not insert " + nterm_sse->GetIdentification());
        BCL_Assert( sp_model->Insert( cterm_sse), "could not insert " + cterm_sse->GetIdentification());

        GrowSSE( *nterm_sse, sp_model, biol::AASequenceFlexibility::e_CTerminal, PHI_PSI_GEN);
        GrowSSE( *cterm_sse, sp_model, biol::AASequenceFlexibility::e_NTerminal, PHI_PSI_GEN);
      }

      // end
      return sp_model;
    }

    //! @brief grows a given sse in a protein model
    //! @param SSE the sse in the protein model that will be grown
    //! @param MODEL the protein model in which the sse will be grown
    //! @param DIRECTION the direction the sse needs to be grown
    //! @param PHI_PSI_GEN the method for generating the phi and psi angles to be assigned to the sse
    void MutateProteinModelGrowSSE::GrowSSE
    (
      const assemble::SSE &SSE,
      util::ShPtr< assemble::ProteinModel> &MODEL,
      const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN
    )
    {
      // construct locator for SSE
      const util::ShPtr< assemble::LocatorSSE> sp_sse_locator
      (
        new assemble::LocatorSSE
        (
          SSE.GetChainID(), SSE.GetFirstAA()->GetSeqID(), SSE.GetLastAA()->GetSeqID()
        )
      );

      // construct locator for anchor res
      util::ShPtr< assemble::LocatorAA> sp_anchor_aa_locator;

      // true if C to N direction - anchor aa is next aa in sequence
      if( DIRECTION == biol::AASequenceFlexibility::e_NTerminal)
      {
        sp_anchor_aa_locator =
          util::ShPtr< assemble::LocatorAA>
         (
           new assemble::LocatorAA( SSE.GetChainID(), SSE.GetLastAA()->GetSeqID() + 1)
         );
      }
      else //< N to C direction - anchor aa is previous aa in sequence
      {
        sp_anchor_aa_locator =
          util::ShPtr< assemble::LocatorAA>
         (
           new assemble::LocatorAA( SSE.GetChainID(), SSE.GetFirstAA()->GetSeqID() - 1)
         );
      }

      BCL_MessageDbg
      (
        "\tanchor aa locator => chain: " + util::Format()( sp_anchor_aa_locator->GetLocatorChain().GetChainID()) +
        " id: " + util::Format()( sp_anchor_aa_locator->GetAAID())
      );
      BCL_MessageDbg
      (
        "\tdirection: " + biol::AASequenceFlexibility::GetSequenceDirectionName( DIRECTION)
      );

      // construct a mutate and apply it
      MODEL =
        MutateProteinModelGrowSSE( sp_sse_locator, PHI_PSI_GEN, sp_anchor_aa_locator, DIRECTION)( *MODEL).GetArgument();
    }

  } // namespace fold
} // namespace bcl
