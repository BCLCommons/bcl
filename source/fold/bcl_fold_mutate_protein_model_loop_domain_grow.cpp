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
#include "fold/bcl_fold_mutate_protein_model_loop_domain_grow.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutate_aa_sequence_grow.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelLoopDomainGrow::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelLoopDomainGrow())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelLoopDomainGrow::MutateProteinModelLoopDomainGrow() :
      m_LocatorLoopDomain(),
      m_PhiPsiGenerator()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param LOOP_DOMAIN_LOCATOR the locator which will be used to find a Loop Domain in a protein model
    //! @param PHI_PSI_GENERATOR method to be used in order to generate phi and psi angles as the loop domain is grown
    MutateProteinModelLoopDomainGrow::MutateProteinModelLoopDomainGrow
    (
      const util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > &LOOP_DOMAIN_LOCATOR,
      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GENERATOR
    ) :
      m_LocatorLoopDomain( LOOP_DOMAIN_LOCATOR),
      m_PhiPsiGenerator( PHI_PSI_GENERATOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelLoopDomainGrow
    MutateProteinModelLoopDomainGrow *MutateProteinModelLoopDomainGrow::Clone() const
    {
      return new MutateProteinModelLoopDomainGrow( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelLoopDomainGrow::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelLoopDomainGrow::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // create LoopDomain from "PROTEIN_MODEL" with "m_LocatorLoopDomain"
      util::ShPtr< LoopDomain> loop_domain( m_LocatorLoopDomain->Locate( PROTEIN_MODEL));

      if( loop_domain->GetSequenceDirection() != biol::AASequenceFlexibility::e_CTerminal)
      {
        BCL_MessageCrt( "cannot work with loop domains for nterminal flexibility.");
        return empty_result;
      }

      // get a const reference to the loop segments in "loop_domain"
      const storage::Set< LoopSegment, LoopSegmentSequenceOrder> &loop_segments( loop_domain->GetSegments());

      // make sure the first loop segment is not rigid
      if( loop_segments.Begin()->IsRigid())
      {
        BCL_MessageCrt( "first loop segment should not be rigid");
        return empty_result;
      }

      // find adjacent sse
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > adjecent_sses( PROTEIN_MODEL.GetAdjacentSSEs( *loop_segments.Begin()->GetSSE()));

      if( !adjecent_sses.First().IsDefined())
      {
        BCL_MessageCrt( "cannot grow domain without an adjacent sse in protein model");
        return empty_result;
      }

      // connect the first segment to the n-terminal anchor sse
      const util::ShPtr< assemble::SSE> connected_n_c_sse
      (
        ConnectNonRigidSSE( *adjecent_sses.First(), *loop_segments.Begin()->GetSSE())
      );

      loop_domain->ReplaceSegment( LoopSegment( connected_n_c_sse, loop_segments.Begin()->IsRigid()));

      // now iterate through the rest of the loop segments and grow or connect them
      for
      (
        storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator
          segment_itr_a( loop_segments.Begin()), segment_itr_b( ++loop_segments.Begin()),
          segment_itr_end( loop_segments.End());
        segment_itr_a != segment_itr_end && segment_itr_b != segment_itr_end;
      )
      {
        // get pointer to the n terminarl sse which will be used as the anchor point
        const util::ShPtr< assemble::SSE> &n_terminal_sse( segment_itr_a->GetSSE());

        const bool n_terminal_sse_is_rigid( segment_itr_a->IsRigid());

        // Consecutive rigid loop segments following the loop segment denoted by "segment_itr_b" are moved
        // in conjunction with "segment_itr_b" since they must be a part of the same rigid domain as "segment_itr_b"
        // is in. Rigid segments should have coordinates, so they can all be connected as a rigid body.
        // If "segment_itr_b" is not rigid, it is assumed to not have coordinates and will be grown from scratch
        // from the n terminal anchor sse denoted by "n_terminal_sse". Even if it has coordinates, if it is
        // not rigid, it will be grown from scratch.

        // true if the segment denoted by "segment_itr_b" is not rigid
        if( !segment_itr_b->IsRigid())
        {
          BCL_MessageStd( "segment_itr_b is not rigid " + segment_itr_b->GetSSE()->GetIdentification());

          const storage::VectorND< 2, util::ShPtr< assemble::SSE> > connected_n_c_sses
          (
            ConnectNonRigidSSE( *n_terminal_sse, *segment_itr_b->GetSSE())
          );
          loop_domain->ReplaceSegment( LoopSegment( connected_n_c_sses.First(), segment_itr_a->IsRigid()));

            std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool> status
            (
              loop_domain->ReplaceSegment( LoopSegment( connected_n_c_sses.Second(), segment_itr_b->IsRigid()))
            );

          // set a to position of last replaced segment and b to next segment
          segment_itr_a = status.first;
          segment_itr_b = segment_itr_a;
          ++segment_itr_b;
        }
        else //< segment pointed by to "segment_itr_b" is rigid
        {
          // create ShPtr list of sses which will hold the sses that need to be connect along with "segment_itr_b"
          storage::List< LoopSegment> c_terminal_sses( 1, *segment_itr_b);

          // create c-terminal sequence which will have concatenated sequences of sses connected with "segment_itr_b"
          biol::AASequence c_terminal_sequence( *segment_itr_b->GetSSE());

          BCL_MessageStd( "segment_itr_b is rigid " + segment_itr_b->GetSSE()->GetIdentification());

          // increase the iterators
          ++segment_itr_a;
          ++segment_itr_b;

          // stop if segment_itr_b is at the end
          if( segment_itr_b == segment_itr_end)
          {
            // connect sequences
            const storage::Pair< util::ShPtr< assemble::SSE>, storage::List< LoopSegment> > connected_n_c_sses
            (
              ConnectRigidSequence( *n_terminal_sse, c_terminal_sequence, c_terminal_sses)
            );

            loop_domain->ReplaceSegment( LoopSegment( connected_n_c_sses.First(), segment_itr_a->IsRigid()));
            loop_domain->ReplaceSegment( connected_n_c_sses.Second());

            // break out of for loop
            break;
          }

          // collect consecutive segments that are also rigid so that they can all be moved together
          while( segment_itr_b->IsRigid())
          {
            // add sse denoted by "segment_itr_b" to "c_terminal_sses"
            c_terminal_sses.PushBack( *segment_itr_b);

            // append the sequence of "segment_itr_b" to "c_terminal_sequence"
            c_terminal_sequence.AppendSequence( *segment_itr_b->GetSSE());

            // increase iterators
            ++segment_itr_a;
            ++segment_itr_b;

            // stop if segment_itr_b is at the end
            if( segment_itr_b == segment_itr_end)
            {
              // break out of while loop
              break;
            }
          }

          // connect to "n_terminal_sse" all the sequences in "c_terminal_sequence" as a rigid body
          const storage::Pair< util::ShPtr< assemble::SSE>, storage::List< LoopSegment> > connected_n_c_sses
          (
            ConnectRigidSequence( *n_terminal_sse, c_terminal_sequence, c_terminal_sses)
          );

          loop_domain->ReplaceSegment( LoopSegment( connected_n_c_sses.First(), n_terminal_sse_is_rigid));

          std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool> status
          (
            loop_domain->ReplaceSegment( connected_n_c_sses.Second())
          );

          // set a to position of last replaced segment and b to next segment
          segment_itr_a = status.first;
          segment_itr_b = segment_itr_a;
          ++segment_itr_b;
        }
      }

      // new protein model
      util::ShPtr< assemble::ProteinModel> new_protein_model( PROTEIN_MODEL.Clone());

      // iterate through the loop domain and replace the sses in the protein model
      util::ShPtr< assemble::ProteinModel> replaced_sse_model
      (
        loop_domain->UpdateProteinModel( *new_protein_model)
      );

      // create mutate result from "replaced_sse_model"
      math::MutateResult< assemble::ProteinModel> mutate_result( replaced_sse_model, *this);

      // return "mutate_result"
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelLoopDomainGrow::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LocatorLoopDomain, ISTREAM);
      io::Serialize::Read( m_PhiPsiGenerator, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelLoopDomainGrow::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LocatorLoopDomain, OSTREAM, INDENT);
      io::Serialize::Write( m_PhiPsiGenerator, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief ConnectNonRigidSSE connects a non rigid c-terminal sse to an n-terminal anchor sse
    //! @param N_TERMINAL_SSE the sse which will have a c-terminal portion connected to it
    //! @param C_TERMINAL_SSE the non-rigid c-terminal sse which will be connected to "N_TERMINAL_SSE"
    util::ShPtr< assemble::SSE> MutateProteinModelLoopDomainGrow::ConnectNonRigidSSE
    (
      const assemble::SSE &N_TERMINAL_SSE, const assemble::SSE &C_TERMINAL_SSE
    ) const
    {
      // create a MutateAASequenceGrow to be used for growing the c-terminal sequence
      const MutateAASequenceGrow grower( m_PhiPsiGenerator, util::SiPtr< const biol::AASequence>( N_TERMINAL_SSE));

      // grow the c-terminal sequence onto the n-terminal sequence and get the resulting aasequence
      const util::ShPtr< biol::AASequence> connected_sequence( grower( C_TERMINAL_SSE).GetArgument());

      // make sure the new aa sequence is defined
      BCL_Assert( connected_sequence.IsDefined(), "connected sequence is not defined");

      // make sure the sizes add up for the connected sequence and the two unconnected sequences
      BCL_Assert
      (
        connected_sequence->GetSize() == ( N_TERMINAL_SSE.GetSize() + C_TERMINAL_SSE.GetSize()),
        "size is " + util::Format()( connected_sequence->GetSize()) + " but should be " +
        util::Format()( N_TERMINAL_SSE.GetSize() + C_TERMINAL_SSE.GetSize())
      );

      // make new sse for c-terminal loop segment sequence
      const util::ShPtr< assemble::SSE> new_cterminal_loop_segment
      (
        new assemble::SSE
        (
          connected_sequence->SubSequence( N_TERMINAL_SSE.GetSize(), C_TERMINAL_SSE.GetSize()),
          C_TERMINAL_SSE.GetType()
        )
      );

      return new_cterminal_loop_segment;
    }

    //! @brief ConnectRigidSequence connects a non rigid c-terminal sse to an n-terminal anchor sse
    //! @param N_TERMINAL_SSE the sse which will have a c-terminal portion connected to it
    //! @param C_TERMINAL_SEQUENCE
    //! @param C_TERMINAL_SSES the non-rigid c-terminal sse which will be connected to "N_TERMINAL_SSE"
    storage::Pair< util::ShPtr< assemble::SSE>, storage::List< LoopSegment> >
    MutateProteinModelLoopDomainGrow::ConnectRigidSequence
    (
      const assemble::SSE &N_TERMINAL_SSE,
      const biol::AASequence &C_TERMINAL_SEQUENCE,
      const storage::List< LoopSegment> &C_TERMINAL_SSES
    ) const
    {
      // create a MutationResidue object for the first residue in "C_TERMINAL_SEQUENCE"
      MutationResidue c_terminal_start_mutation_residue
      (
        C_TERMINAL_SEQUENCE.GetFirstAA(),
        N_TERMINAL_SSE.GetLastAA(),
        util::ShPtr< biol::AABase>()
      );

      // create a MutationResidue object for the last residue in "n_terminal_sse"
      MutationResidue n_terminal_end_mutation_residue
      (
        N_TERMINAL_SSE.GetLastAA(),
        util::ShPtr< biol::AABase>(),
        C_TERMINAL_SEQUENCE.GetFirstAA()
      );

      // get the phi that the first residue in "C_TERMINAL_SEQUENCE" will have
      const double phi( m_PhiPsiGenerator->operator()( c_terminal_start_mutation_residue).First());

      // make sure phi is defined
      BCL_Assert( util::IsDefined( phi), "phi is not defined");

      // get the psi that the last residue in "n_terminal_sse" will have
      const double psi( m_PhiPsiGenerator->operator()( n_terminal_end_mutation_residue).Second());

      // make sure psi is defined
      BCL_Assert( util::IsDefined( psi), "psi is not defined");

      // connect "C_TERMINAL_SEQUENCE" and "n_terminal_sse" and get the ShPtr to the new connected sequence
      util::ShPtr< biol::AASequence> connected_sequence( new biol::AASequence( N_TERMINAL_SSE));
      biol::AASequenceFactory::AppendSequence( *connected_sequence, C_TERMINAL_SEQUENCE, phi, psi);

      // set a size_t to zero to which will keep track of the starting positions of sses in "connected_sequence"
      size_t sse_start_index( 0);

      // update the nterminal sse
      // make new sse for nterminal sse
      util::ShPtr< assemble::SSE> new_nterminal_sse
      (
        new assemble::SSE
        (
          connected_sequence->SubSequence( 0, N_TERMINAL_SSE.GetSize()), N_TERMINAL_SSE.GetType()
        )
      );

      // increase sse start index by the size of "n_terminal_sse"
      sse_start_index += N_TERMINAL_SSE.GetSize();

      storage::List< LoopSegment> new_sses;

      // update the c-terminal sses
      for
      (
        storage::List< LoopSegment>::const_iterator
          c_terminal_sses_itr( C_TERMINAL_SSES.Begin()), c_terminal_sses_itr_end( C_TERMINAL_SSES.End());
        c_terminal_sses_itr != c_terminal_sses_itr_end;
        ++c_terminal_sses_itr
      )
      {
        // get shptr to current sse denoted by "c_terminal_sses_itr"
        const assemble::SSE &current_sse( *c_terminal_sses_itr->GetSSE());

        // get identifier of "current_sse_identifier" which will be used as a correctness check
        const std::string current_sse_identifier( current_sse.GetIdentification());

        // make new sse for current sse
        util::ShPtr< assemble::SSE> new_current_sse
        (
          new assemble::SSE
          (
            connected_sequence->SubSequence( sse_start_index, current_sse.GetSize()), current_sse.GetType()
          )
        );

        // increase "sse_start_index" to the position of the start of the next sse
        sse_start_index += current_sse.GetSize();

        new_sses.PushBack( LoopSegment( new_current_sse, c_terminal_sses_itr->IsRigid()));

        // get the identification of the changed "current_sse"
        const std::string new_current_sse_identifier( current_sse.GetIdentification());

        // make sure the two identifications are the same, otherwise something went wrong
        BCL_Assert
        (
          current_sse_identifier == new_current_sse_identifier, "after setting sse the identifier has changed from " +
          current_sse_identifier + " to " + new_current_sse_identifier
        );
      }

      return storage::Pair< util::ShPtr< assemble::SSE>, storage::List< LoopSegment> >( new_nterminal_sse, new_sses);
    }

  } // namespace fold
} // namespace bcl
