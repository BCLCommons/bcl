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
#include "fold/bcl_fold_loop_domain_c_to_n.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopDomainCToN::s_Instance
    (
      GetObjectInstances().AddInstance( new LoopDomainCToN())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopDomainCToN::LoopDomainCToN() :
      LoopDomain()
    {
    }

    //! @brief constructor taking parameters to create member variables
    //! @param SEGMENTS list of LoopSegments which will be used to create "m_LoopSegments"
    LoopDomainCToN::LoopDomainCToN
    (
      const storage::List< LoopSegment> &SEGMENTS
    ) :
      LoopDomain( SEGMENTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LoopDomainCToN
    LoopDomainCToN *LoopDomainCToN::Clone() const
    {
      return new LoopDomainCToN( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopDomainCToN::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the directionality of the loop domain
    //! @return sequence direction in this case e_NTerminal since changes are propagated nterminally
    biol::AASequenceFlexibility::SequenceDirection LoopDomainCToN::GetSequenceDirection() const
    {
      return biol::AASequenceFlexibility::e_NTerminal;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the residues in the loop domain in sequence order including anchor sse anchor residue and the
    //!        created psuedo residue
    //! @return map with the amino acids and a bool indicating whether or not the aa is part of a rigid segment or not
    storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > LoopDomainCToN::GetResidues() const
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
    const util::ShPtr< biol::AABase> &LoopDomainCToN::GetMostProximalLoopSegmentAA() const
    {
      return ( --GetSegments().End())->GetSSE()->GetLastAA();
    }

    //! @brief gives the loop segment that is most distant in sequence to attachment to the anchor sse
    //!        this is the sse that the pseudo residue attaches to
    //! @return sse that the pseudo residue attaches to and is the most distant in sequence from anchor sse
    const assemble::SSE &LoopDomainCToN::GetMostDistalLoopSegment() const
    {
      return GetSegments().Begin()->GetConstSSEReference();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief sets the phi of a residue in the loop domain
    //! @param AA the amino acid that will be found and mutated
    //! @param PHI the phi angle that residue will be set to
    void LoopDomainCToN::SetPhi( const biol::AABase &AA, const double PHI)
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
        " in nterm direction"
      );

      storage::Set< LoopSegment, LoopSegmentSequenceOrder> new_segments;

      // iterate over loop segments
      storage::VectorND< 2, math::TransformationMatrix3D> transformations;

      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::reverse_iterator itr( m_LoopSegments.ReverseBegin()), itr_end( m_LoopSegments.ReverseEnd());
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
            biol::AASequenceFlexibility::e_NTerminal
          );
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
        ++itr;

        break;
      }

      for( ; itr != itr_end; ++itr)
      {
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        new_sse->Transform( transformations.First());
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
      }

      m_LoopSegments = new_segments;
    }

    //! @brief sets the psi of a residue in the loop domain
    //! @param AA the amino acid that will be found and mutated
    //! @param PSI the psi angle that residue will be set to
    void LoopDomainCToN::SetPsi( const biol::AABase &AA, const double PSI)
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
        " in nterm direction"
      );

      storage::Set< LoopSegment, LoopSegmentSequenceOrder> new_segments;

      // iterate over loop segments
      storage::VectorND< 2, math::TransformationMatrix3D> transformations;

      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::reverse_iterator itr( m_LoopSegments.ReverseBegin()), itr_end( m_LoopSegments.ReverseEnd());
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
            biol::AASequenceFlexibility::e_NTerminal
          );
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
        ++itr;
        break;
      }

      for( ; itr != itr_end; ++itr)
      {
        util::ShPtr< assemble::SSE> new_sse( itr->GetSSE().HardCopy());
        new_sse->Transform( transformations.First());
        new_segments.Insert( LoopSegment( new_sse, itr->IsRigid()));
      }

      m_LoopSegments = new_segments;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LoopDomainCToN::Read( std::istream &ISTREAM)
    {
      // read members
      LoopDomain::Read( ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LoopDomainCToN::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      LoopDomain::Write( OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
