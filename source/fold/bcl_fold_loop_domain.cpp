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
#include "fold/bcl_fold_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_locator_loop_domain.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopDomain::s_Instance
    (
      GetObjectInstances().AddInstance( new LoopDomain())
    );

    //! @brief return command line flag for specifying loop domains via a file
    //! @return command line flag for specifying loop domains via a file
    util::ShPtr< command::FlagInterface> &LoopDomain::GetFlagLoopDomainFilename()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "loop_domain",
          "\tSpecific loop domain(s) should be built. If this flag is not given, then all loops will be built."
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter
        (
          "loop_domain_filename", "\tfilename for loop domain file", "loop_domain.txt"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
      }

      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopDomain::LoopDomain() :
      m_LoopSegments()
    {
    }

    //! @brief constructor taking parameters to create member variables
    //! @param SEGMENTS list of LoopSegments which will be used to create "m_LoopSegments"
    LoopDomain::LoopDomain
    (
      const storage::List< LoopSegment> &SEGMENTS
    ) :
      m_LoopSegments()
    {
      BCL_Assert( !SEGMENTS.IsEmpty(), "need at least one segment in loop domain");

      // build up m_LoopSegments
      for
      (
        storage::List< LoopSegment>::const_iterator
          segment_itr( SEGMENTS.Begin()), segment_itr_end( SEGMENTS.End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // create const reference to current loop segment
        const LoopSegment &loop_segment( *segment_itr);

        // insert "loop_segment" into "loop_segment"
        m_LoopSegments.InsertElement( loop_segment);
      }
    }

    //! @brief copy constructor
    //! @param OTHER loop domain to be copied
    LoopDomain::LoopDomain( const LoopDomain &OTHER) :
      m_LoopSegments( OTHER.m_LoopSegments)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LoopSegment
    LoopDomain *LoopDomain::Clone() const
    {
      return new LoopDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LoopDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the directionality of the loop domain
    //! @return sequence direction in this case e_CTerminal since changes are propagated cterminally
    biol::AASequenceFlexibility::SequenceDirection LoopDomain::GetSequenceDirection() const
    {
      return biol::AASequenceFlexibility::e_CTerminal;
    }

    //! @brief GetSegments gives the set LoopSegments ordered by sequence that make up the loop domain
    //! @return the set LoopSegments ordered by sequence that make up the loop domain
    const storage::Set< LoopSegment, LoopSegmentSequenceOrder> &LoopDomain::GetSegments() const
    {
      return m_LoopSegments;
    }

    //! @brief finds the segment and replaces it with the given segment
    //! @param LOOP_SEGMENT segment that will replace the currently matching segment
    //! @return pair of iterator pointing to iterator where the replaced segment is and bool indicating success or not
    std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool>
    LoopDomain::ReplaceSegment( const LoopSegment &LOOP_SEGMENT)
    {
      // get iterator to segment where "LOOP_SEGMENT" is
      const storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator segment_itr
      (
        m_LoopSegments.Find( LOOP_SEGMENT)
      );

      // true if could not find loop segment containign "RESIDUE"
      if( segment_itr == m_LoopSegments.End())
      {
        BCL_MessageCrt( "could not find segment with sse " + LOOP_SEGMENT.GetSSE()->GetIdentification());
      }

      // remove the old segment and insert the new loop segment into m_LoopSegments
      m_LoopSegments.RemoveElement( segment_itr);
      std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool>
        insert_status( m_LoopSegments.Insert( LOOP_SEGMENT));
      BCL_Assert
      (
        insert_status.second, "failed to insert segment with sse " +
        LOOP_SEGMENT.GetSSE()->GetIdentification()
      );

      return insert_status;
    }

    //! @brief finds the segments and replaces them with the given segments
    //! @param LOOP_SEGMENTS segments that will replace the currently matching segments
      //! @return returns iter and bool to last replaced element
    std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool>
    LoopDomain::ReplaceSegment( const storage::List< LoopSegment> &LOOP_SEGMENTS)
    {
      std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool> status;
      for
      (
        storage::List< LoopSegment>::const_iterator
          segments_itr( LOOP_SEGMENTS.Begin()), segments_itr_end( LOOP_SEGMENTS.End());
        segments_itr != segments_itr_end;
        ++segments_itr
      )
      {
        status = ReplaceSegment( *segments_itr);
        BCL_Assert( status.second, "could not replace segment " + segments_itr->GetSSE()->GetIdentification());
      }

      return status;
    }

    //! @brief FindSegment provides an iterator to the segment containing a residue of interest
    //! @param RESIDUE the residue of interest for which the segment it is in is desired to be known
    //! @return iterator to the LoopSegment that contains RESIDUE
    storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator LoopDomain::FindSegment
    (
      const biol::AABase &RESIDUE
    ) const
    {
      // get the seq id of "RESIDUE"
      const int resi_seq_id( RESIDUE.GetSeqID());

      // get the chain id of RESIDUE
      const char resi_chain_id( RESIDUE.GetChainID());

      // find the segment where "RESIDUE" is and return iterator to it
      return FindSegment( resi_chain_id, resi_seq_id);
    }

    //! @brief FindSegment finds a residue of interest in the LoopDomain
    //! @param RESIDUE the residue of interest which should be found
    //! @return ShPtr pointing to the residue of interest
    util::ShPtr< biol::AABase>
    LoopDomain::FindResidue( const biol::AABase &RESIDUE) const
    {
      // find "RESIDUE"
      return FindResidue( assemble::LocatorAA( RESIDUE.GetChainID(), RESIDUE.GetSeqID()));
    }

    //! @brief FindSegment finds a residue of interest in the LoopDomain as described by a LocatorAA
    //! @param LOCATOR_AA the LocatorAA which describes the residue of iterest
    //! @return ShPtr pointing to the residue of interest as described by "LOCATOR_AA"
    util::ShPtr< biol::AABase> LoopDomain::FindResidue( const assemble::LocatorAA &LOCATOR_AA) const
    {
      // get the chain and seq ids out of "LOCATOR_AA"
      const int resi_seq_id( LOCATOR_AA.GetAAID());
      const char resi_chain_id( LOCATOR_AA.GetLocatorChain().GetChainID());

      // get an iterator to the segment where "resi_chain_id" and "resi_seq_id" are within
      const storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator segment_itr
      (
        FindSegment( resi_chain_id, resi_seq_id)
      );

      // true if the segment where the residue indicated by "LOCATOR_AA" would be could not be found
      if( segment_itr == m_LoopSegments.End())
      {
        // return
        BCL_MessageCrt( "could not locate segment for " + LOCATOR_AA.GetIdentification());

        // return empty pointer to an aabase
        return util::ShPtr< biol::AABase>();
      }

      // get the sse in the loopsegment denoted by "segment_itr"
      const util::ShPtr< assemble::SSE> segment_sse( segment_itr->GetSSE());

      // get ShPtr to the residue described by "LOCATOR_AA"
      const util::ShPtr< biol::AABase> located_residue( FindResidue( *segment_sse, resi_chain_id, resi_seq_id));

      // return "located_residue"
      return located_residue;
    }

    //! @brief SetResidue finds the pointer to residue of interest in the LoopDomain and sets the pointer to given ptr
    //! @param RESIDUE the ptr to residue of interest which should be found
    void LoopDomain::SetPtrToResidue( const util::ShPtr< biol::AABase> &RESIDUE)
    {
      BCL_Assert( RESIDUE.IsDefined(), " RESIDUE not defined");

      // get the chain and seq ids out of "LOCATOR_AA"
      const int resi_seq_id( RESIDUE->GetSeqID());
      const char resi_chain_id( RESIDUE->GetChainID());

      // get an iterator to the segment where "resi_chain_id" and "resi_seq_id" are within
      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator segment_itr
      (
        FindSegment( resi_chain_id, resi_seq_id)
      );

      // true if the segment where the residue could not be found
      if( segment_itr == m_LoopSegments.End())
      {
        // return
        BCL_MessageCrt( "could not locate segment for " + RESIDUE->GetIdentification());
        return;
      }

      LoopSegment new_loop_segment
      (
        util::ShPtr< assemble::SSE>( segment_itr->GetSSE()->Clone()), segment_itr->IsRigid()
      );

      // get the sse in the loopsegment denoted by "segment_itr"
      assemble::SSE &segment_sse( new_loop_segment.GetSSEReference());

      // iterate through the sequence of "segment_sse"
      for
      (
        biol::AASequence::iterator sse_itr( segment_sse.Begin()), sse_itr_end( segment_sse.End());
          sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // true if chain id and seq id of "RESIDUE" and the amino acid denoted by "sse_itr" are the same
        if( ( *sse_itr)->GetChainID() == resi_chain_id && ( *sse_itr)->GetSeqID() == resi_seq_id)
        {
          // set shptr to aa denoted by "sse_itr"
          ( *sse_itr) = RESIDUE;
        }
      }

      // remove the old segment and insert the new loop segment into m_LoopSegments
      m_LoopSegments.RemoveElement( segment_itr);
      BCL_Assert
      (
        m_LoopSegments.Insert( new_loop_segment).second, "failed to insert segment with sse " +
        new_loop_segment.GetSSE()->GetIdentification()
      );
    }

    //! @brief finds the pointer to residue of interest in the LoopDomain and gives iterator to it
    //! @param RESIDUE the ptr to residue of interest which should be found
    //! @return iterator to residue
    biol::AASequence::const_iterator LoopDomain::FindResidueIterator( const biol::AABase &RESIDUE) const
    {
      // get the chain and seq ids out of "LOCATOR_AA"
      const int resi_seq_id( RESIDUE.GetSeqID());
      const char resi_chain_id( RESIDUE.GetChainID());

      // get an iterator to the segment where "resi_chain_id" and "resi_seq_id" are within
      const storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator segment_itr
      (
        FindSegment( resi_chain_id, resi_seq_id)
      );

      BCL_Assert( segment_itr != m_LoopSegments.End(), "could not locate segment for " + RESIDUE.GetIdentification());

      // get the sse in the loopsegment denoted by "segment_itr"
      const assemble::SSE &segment_sse( *segment_itr->GetSSE());

      // iterate through the sequence of "segment_sse"
      for
      (
        biol::AASequence::const_iterator sse_itr( segment_sse.Begin()), sse_itr_end( segment_sse.End());
          sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // true if chain id and seq id of "RESIDUE" and the amino acid denoted by "sse_itr" are the same
        if( ( *sse_itr)->GetChainID() == resi_chain_id && ( *sse_itr)->GetSeqID() == resi_seq_id)
        {
          // return shptr to aa denoted by "sse_itr"
          return sse_itr;
        }
      }

      return segment_sse.End();
    }

    //! @brief FindSegment gives iterator to the segment containing an AA of interest described by chain and seq id
    //! @param RESI_CHAIN_ID the chain id of the residue of interest
    //! @param RESI_SEQ_ID the seq id of the residue of interest
    //! @return iterator to the LoopSegment that contains the residue described by "RESI_CHAIN_ID" and "RESI_SEQ_ID"
    storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator LoopDomain::FindSegment
    (
      const char RESI_CHAIN_ID, const int RESI_SEQ_ID
    ) const
    {
      // loop through the segments until the one containing "RESIDUE" is found
      for
      (
        storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator
          segment_itr( m_LoopSegments.Begin()), segment_itr_end( m_LoopSegments.End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // get the start and ending seq ids and chain id of the sse denoted by "segment_itr"
        const int sse_start_seq_id( segment_itr->GetSSE()->GetFirstAA()->GetSeqID());
        const int sse_last_seq_id( segment_itr->GetSSE()->GetLastAA()->GetSeqID());
        const char sse_chain_id( segment_itr->GetSSE()->GetChainID());

        // true if chain ids match and "RESI_SEQ_ID" is between "sse_start_seq_id" and "sse_last_seq_id"
        if( sse_start_seq_id <= RESI_SEQ_ID && sse_last_seq_id >= RESI_SEQ_ID && RESI_CHAIN_ID == sse_chain_id)
        {
          // return the iterator
          return segment_itr;
        }
      }

      // message that the segment was not found
      BCL_MessageDbg
      (
        "did not find segment with residue with chain id " + util::Format()( RESI_CHAIN_ID) +
        " and seq id " + util::Format()( RESI_SEQ_ID)
      );

      // return iterator to the end of "m_LoopSegments"
      return m_LoopSegments.End();
    }

    //! @brief FindResidue finds a residue of interest in an SSE of the loop domain as described by a chain and seq id
    //! @param SEGMENT_SSE the sse which contains the residue of interest
    //! @param CHAIN_ID the chain id of the residue of interest
    //! @param SEQ_ID the seq id of the residue of interest
    //! @return ShPtr pointing to the residue of interest
    util::ShPtr< biol::AABase> LoopDomain::FindResidue
    (
      const assemble::SSE &SEGMENT_SSE, const char CHAIN_ID, const int SEQ_ID
    ) const
    {
      // iterate through the sequence of "segment_sse"
      for
      (
        biol::AASequence::const_iterator sse_itr( SEGMENT_SSE.Begin()), sse_itr_end( SEGMENT_SSE.End());
          sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // true if chain id and seq id of "RESIDUE" and the amino acid denoted by "sse_itr" are the same
        if( ( *sse_itr)->GetChainID() == CHAIN_ID && ( *sse_itr)->GetSeqID() == SEQ_ID)
        {
          // return shptr to aa denoted by "sse_itr"
          return *sse_itr;
        }
      }

      // message that the residue could not be found
      BCL_MessageDbg
      (
        "could not find residue with chain id " + util::Format()( CHAIN_ID) + " and seq id "
        + util::Format()( SEQ_ID) + " in sse " + SEGMENT_SSE.GetIdentification()
      );

      // return empty shptr
      return util::ShPtr< biol::AABase>();
    }

    //! @brief GetChainID provides the chain id of the LoopDomain. Since residues must be from the same chain, it
    //! can get the chain id from any residue.
    //! @return the chain id of the loop domain
    const char &LoopDomain::GetChainID() const
    {
      // return the chain id of any of the components of the loop domain since it is checked that they are all the
      // same
      return m_LoopSegments.Begin()->GetSSE()->GetChainID();
    }

    //! @brief returns the residues in the loop domain in sequence order including anchor sse anchor residue and the
    //!        created psuedo residue
    //! @return map with the amino acids and a bool indicating whether or not the aa is part of a rigid segment or not
    storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > LoopDomain::GetResidues() const
    {
      // map with residues and bool indicating if the residue is rigid or not
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > residues;

      // iterate through all the segments
      for
      (
        storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator
          segment_itr( GetSegments().Begin()), segment_itr_end( GetSegments().End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // get reference to the sse denoted by "segment_itr"
        const assemble::SSE &sse( *segment_itr->GetSSE());

        const bool is_rigid( segment_itr->IsRigid());

        // iterate through the amino acids of the current segment
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator aa_itr( sse.Begin()), aa_itr_end( sse.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // insert the residue and whether or not it is rigid
          residues.PushBack( storage::Pair< util::ShPtr< biol::AABase>, bool>( *aa_itr, is_rigid));
        } // iterate through sse aas
      } // iterate through segments

      return residues;
    }

    //! @brief gives the residue that is attached to the anchor sse
    //!        i.e. the residue that is closest to point of attachment
    //! @return ShPtr to residue that is attached to the anchor sse
    const util::ShPtr< biol::AABase> &LoopDomain::GetMostProximalLoopSegmentAA() const
    {
      return GetSegments().Begin()->GetSSE()->GetFirstAA();
    }

    //! @brief gives the loop segment that is most distant in sequence to attachment to the anchor sse
    //!        this is the sse that the pseudo residue attaches to
    //! @return sse that the pseudo residue attaches to and is the most distant in sequence from anchor sse
    const assemble::SSE &LoopDomain::GetMostDistalLoopSegment() const
    {
      return ( --GetSegments().End())->GetConstSSEReference();
    }

    //! @brief gives a ShPtr to the residue following the given residue in sequence
    //! @param AA the residue for which the following residue will be found
    //! @return ShPtr to residue that follows AA in sequence
    util::ShPtr< biol::AABase> LoopDomain::FindFollowingResidue( const biol::AABase &AA) const
    {
      const storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > residues( GetResidues());

      // need at least one residue, so that Last() is valid
      if( residues.IsEmpty())
      {
        return util::ShPtr< biol::AABase>();
      }

      // iterate through the list of residues until AA is found
      // go to one before end since if the last residue is reached, there is no following residue
      for
      (
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator
          resi_itr( residues.Begin()), resi_itr_end( residues.Last());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        const biol::AABase &resi( *resi_itr->First());

        if( ResiduesMatch( AA, resi))
        {
          ++resi_itr;
          // increment the iterator to go to next residue and return the following residue
          return resi_itr->First();
        }
      }

      return util::ShPtr< biol::AABase>();
    }

    //! @brief gives a ShPtr to the residue preceding the given residue in sequence
    //! @param AA the residue for which the preceding residue will be found
    //! @return ShPtr to residue that precedes AA in sequence
    util::ShPtr< biol::AABase> LoopDomain::FindPreviousResidue( const biol::AABase &AA) const
    {
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > residues( GetResidues());

      // need at least one residue, so that Last() is valid
      if( residues.IsEmpty())
      {
        return util::ShPtr< biol::AABase>();
      }

      // iterate through the list of residues until AA is found
      // start at one after begin since if there is no previous residue for the first residue
      for
      (
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator
          resi_itr( ++residues.Begin()), resi_itr_end( residues.End());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        const biol::AABase &resi( *resi_itr->First());

        if( ResiduesMatch( AA, resi))
        {
          // decrement the iterator to go to previous residue and return the previous residue
          return ( --resi_itr)->First();
        }
      }

      return util::ShPtr< biol::AABase>();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief sets the psi of a residue in the loop domain
    //! @param AA the amino acid that will be found and mutated
    //! @param PSI the psi angle that residue will be set to
    void LoopDomain::SetPsi( const biol::AABase &AA, const double PSI)
    {
      // get the current psi angle
      double current_psi( AA.Psi());

      // phi within aa is not defined
      if( !util::IsDefined( current_psi))
      {
        // use the neighbor for the phi
        const util::ShPtr< biol::AABase> following_residue( FindFollowingResidue( AA));

        // neighbor is defined
        if( following_residue.IsDefined())
        {
          current_psi = AA.CalculatePsi( following_residue->GetAtom( biol::GetAtomTypes().N));
        }
        // no neighbor defined, psi is undefined, no loop domain change required
        else
        {
          return;
        }
      }

      const storage::VectorND< 2, double> phi_psi_change( 0.0, PSI - current_psi);

      BCL_MessageDbg
      (
        "changing psi by: " + util::Format()( phi_psi_change.Second()) + " from AA: " + AA.GetIdentification() +
        " in cterm direction"
      );

      storage::Set< LoopSegment, LoopSegmentSequenceOrder> new_segments;

      // iterate over loop segments
      storage::VectorND< 2, math::TransformationMatrix3D> transformations;

      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::iterator itr( m_LoopSegments.Begin()), itr_end( m_LoopSegments.End());
      for( ; itr != itr_end; ++itr)
      {
        if( itr->IsRigid())
        {
          new_segments.Insert( *itr);
          continue;
        }

        if( itr->GetSSE()->FindAABySeqID( AA.GetSeqID()) == itr->GetSSE()->End())
        {
          new_segments.Insert( *itr);
          continue;
        }

        // amino acid is part of that non rigid loop segment
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        transformations =
          biol::AASequenceFlexibility::ChangePhiPsi
          (
            *new_sse,
            AA.GetSeqID(),
            phi_psi_change,
            biol::AASequenceFlexibility::e_CTerminal
          );
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
        ++itr;

        break;
      }

      for( ; itr != itr_end; ++itr)
      {
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        new_sse->Transform( transformations.Second());
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
      }

      m_LoopSegments = new_segments;
    }

    //! @brief sets the phi of a residue in the loop domain
    //! @param AA the amino acid that will be found and mutated
    //! @param PHI the phi angle that residue will be set to
    void LoopDomain::SetPhi( const biol::AABase &AA, const double PHI)
    {
      // get the current phi angle
      double current_phi( AA.Phi());

      // phi within aa is not defined
      if( !util::IsDefined( current_phi))
      {
        // use the neighbor for the phi
        const util::ShPtr< biol::AABase> previous_residue( FindPreviousResidue( AA));

        // neighbor is defined
        if( previous_residue.IsDefined())
        {
          current_phi = AA.CalculatePhi( previous_residue->GetAtom( biol::GetAtomTypes().C));
        }
        // no neighbor defined, phi is undefined, no loop domain change required
        else
        {
          return;
        }
      }

      const storage::VectorND< 2, double> phi_psi_change( PHI - current_phi, 0.0);

      BCL_MessageDbg
      (
        "changing phi by: " + util::Format()( phi_psi_change.First()) + " from AA: " + AA.GetIdentification() +
        " in cterm direction"
      );

      storage::Set< LoopSegment, LoopSegmentSequenceOrder> new_segments;

      // iterate over loop segments
      storage::VectorND< 2, math::TransformationMatrix3D> transformations;

      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::iterator itr( m_LoopSegments.Begin()), itr_end( m_LoopSegments.End());
      for( ; itr != itr_end; ++itr)
      {
        if( itr->IsRigid())
        {
          new_segments.Insert( *itr);
          continue;
        }

        if( itr->GetSSE()->FindAABySeqID( AA.GetSeqID()) == itr->GetSSE()->End())
        {
          new_segments.Insert( *itr);
          continue;
        }

        // amino acid is part of that non rigid loop segment
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        transformations =
          biol::AASequenceFlexibility::ChangePhiPsi
          (
            *new_sse,
            AA.GetSeqID(),
            phi_psi_change,
            biol::AASequenceFlexibility::e_CTerminal
          );
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
        ++itr;

        break;
      }

      for( ; itr != itr_end; ++itr)
      {
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        new_sse->Transform( transformations.Second());
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
      }

      m_LoopSegments = new_segments;
    }

    //! @brief return the points that can be used for CCD
    //! @param PROTEIN_MODEL the model, to find the adjacent sse in to connect to
    //! @return target and moving points
    storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> LoopDomain::TargetAndMovingPointsForCCD
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> points;

      util::SiPtr< const assemble::SSE> sp_left_sse;
      util::SiPtr< const assemble::SSE> sp_right_sse;

      // get the sses - c is flexible side
      if( GetSequenceDirection() == biol::AASequenceFlexibility::e_CTerminal)
      {
        sp_left_sse = ( --m_LoopSegments.End())->GetSSE();
        storage::VectorND< 2, util::SiPtr< const assemble::SSE> > adjacent_sses( PROTEIN_MODEL.GetAdjacentSSEs( *sp_left_sse));
        if( !adjacent_sses.Second().IsDefined())
        {
          return points;
        }
        sp_right_sse = adjacent_sses.Second();
      }
      // n is the flexible side
      else if( GetSequenceDirection() == biol::AASequenceFlexibility::e_NTerminal)
      {
        sp_right_sse = m_LoopSegments.Begin()->GetSSE();
        storage::VectorND< 2, util::SiPtr< const assemble::SSE> > adjacent_sses( PROTEIN_MODEL.GetAdjacentSSEs( *sp_right_sse));
        if( !adjacent_sses.First().IsDefined())
        {
          return points;
        }
        sp_left_sse = adjacent_sses.First();
      }
      else
      {
        return points;
      }

      // amino acids
      const biol::AABase &left_aa( *sp_left_sse->GetLastAA());
      const biol::AABase &right_aa( *sp_right_sse->GetFirstAA());

      // pseudo atoms to superimpose
      const storage::Map< biol::AtomType, biol::Atom> pseudo_atoms( biol::AABackBoneCompleter::GenerateHNCA( left_aa));
      const storage::Map< biol::AtomType, biol::Atom>::const_iterator itr_end( pseudo_atoms.End());

      // hydrogen
      {
        const storage::Map< biol::AtomType, biol::Atom>::const_iterator itr( pseudo_atoms.Find( biol::GetAtomTypes().H));
        const linal::Vector3D &actual_hydrogen_coord( right_aa.GetAtom( biol::GetAtomTypes().H).GetCoordinates());

        if( itr != itr_end && itr->second.GetCoordinates().IsDefined() && actual_hydrogen_coord.IsDefined())
        {
          if( GetSequenceDirection() == biol::AASequenceFlexibility::e_CTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                actual_hydrogen_coord,
                itr->second.GetCoordinates()
              )
            );
          }
          else if( GetSequenceDirection() == biol::AASequenceFlexibility::e_NTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                itr->second.GetCoordinates(),
                actual_hydrogen_coord
              )
            );
          }
        }
      }

      // nitrogen
      {
        const storage::Map< biol::AtomType, biol::Atom>::const_iterator itr( pseudo_atoms.Find( biol::GetAtomTypes().N));
        const linal::Vector3D &actual_nitrogen_coord( right_aa.GetAtom( biol::GetAtomTypes().N).GetCoordinates());

        if( itr != itr_end && itr->second.GetCoordinates().IsDefined() && actual_nitrogen_coord.IsDefined())
        {
          if( GetSequenceDirection() == biol::AASequenceFlexibility::e_CTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                actual_nitrogen_coord,
                itr->second.GetCoordinates()
              )
            );
          }
          else if( GetSequenceDirection() == biol::AASequenceFlexibility::e_NTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                itr->second.GetCoordinates(),
                actual_nitrogen_coord
              )
            );
          }
        }
      }

      // carbon alpha
      {
        const storage::Map< biol::AtomType, biol::Atom>::const_iterator itr( pseudo_atoms.Find( biol::GetAtomTypes().CA));
        const linal::Vector3D &actual_calpha_coord( right_aa.GetAtom( biol::GetAtomTypes().CA).GetCoordinates());

        if( itr != itr_end && itr->second.GetCoordinates().IsDefined() && actual_calpha_coord.IsDefined())
        {
          if( GetSequenceDirection() == biol::AASequenceFlexibility::e_CTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                actual_calpha_coord,
                itr->second.GetCoordinates()
              )
            );
          }
          else if( GetSequenceDirection() == biol::AASequenceFlexibility::e_NTerminal)
          {
            points.PushBack
            (
              coord::CyclicCoordinateDescent::TargetAndMovingPointPair
              (
                itr->second.GetCoordinates(),
                actual_calpha_coord
              )
            );
          }
        }
      }

      // end
      return points;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LoopDomain::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_LoopSegments, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LoopDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_LoopSegments, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief places the sses from a loop domain into a protein model
    //! @param PROTEIN_MODEL the protein model whose sses will be replaced
    //! @return ShPtr to ProteinModel which has had its sses replaced with the sses of LOOP_DOMAIN
    util::ShPtr< assemble::ProteinModel> LoopDomain::UpdateProteinModel
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create a ShPtr copy to "PROTEIN_MODEL"
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // iterate over the loop segments of "LOOP_DOMAIN" in order to replace them in "PROTEIN_MODEL"
      for
      (
        storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator
          segment_itr( GetSegments().Begin()), segment_itr_end( GetSegments().End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // create const reference to the sse in the loop segment currently denoted by "segment_itr"
        const util::ShPtr< assemble::SSE> &current_sse( segment_itr->GetSSE());

        // replace "current_sse" in "temp_model" and make sure the replace was successful
        BCL_Assert( new_model->Replace( current_sse), "failed to replace " + current_sse->GetIdentification());
      }

      // return "temp_model" which has had its sses replaced by the sses of "LOOP_DOMAIN"
      return new_model;
    }

    //! @brief determines if the seq id and chain id of two residues are the same
    //! @param RESI_A first residue
    //! @param RESI_B second residue
    //! @return boolean true if the seqid and chain ids match - false otherwise
    bool LoopDomain::ResiduesMatch( const biol::AABase &RESI_A, const biol::AABase &RESI_B)
    {
      return RESI_A.GetChainID() == RESI_B.GetChainID() && RESI_A.GetSeqID() == RESI_B.GetSeqID();
    }

  } // namespace fold
} // namespace bcl
