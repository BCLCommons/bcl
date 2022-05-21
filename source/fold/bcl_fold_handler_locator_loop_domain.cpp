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
#include "fold/bcl_fold_handler_locator_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "fold/bcl_fold_collector_unconnected_sses.h"
#include "fold/bcl_fold_handler_locator_loop_segment.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_loop_segment.h"
#include "fold/bcl_fold_loop_segment_sequence_order.h"
#include "fold/bcl_fold_mutation_residue.h"
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"
#include "fold/bcl_fold_protocol_loop_coordinate_add.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! the string used to indicate the end of a loop domain
      const std::string HandlerLocatorLoopDomain::s_DomainEndString( "DomainLocatorEnd");

      //! single instance of that class
      const util::SiPtr< const util::ObjectInterface> HandlerLocatorLoopDomain::s_Instance
      (
        GetObjectInstances().AddInstance( new HandlerLocatorLoopDomain())
      );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param RANDOMIZE_PHI_PSI bool true if the phi and psi values held in the pseudo residue should be randomized
    HandlerLocatorLoopDomain::HandlerLocatorLoopDomain( const bool RANDOMIZE_PHI_PSI) :
      m_RandomizePseudoResiduePhiPsi( RANDOMIZE_PHI_PSI)
    {
    }

    //! @brief Clone function
    //! @return pointer to new HandlerLocatorLoopDomain
    HandlerLocatorLoopDomain *HandlerLocatorLoopDomain::Clone() const
    {
      return new HandlerLocatorLoopDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &HandlerLocatorLoopDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief HandleRead creates a LocatorLoopDomain when from an istream and protein model
    //! @param ISTREAM the file stream which contains the loop domain locator information
    //! @param MODEL the protein model which will be used to help create the LocatorLoopDomain
    //! @return LocatorLoopDomain which has been created from ISTREAM and MODEL
    LocatorLoopDomain HandlerLocatorLoopDomain::HandleRead
    (
      std::istream &ISTREAM, const assemble::ProteinModel &MODEL
    ) const
    {
      // get the list of LocatorLoopSegments from "ISTREAM"
      const storage::List< LocatorLoopSegment> loop_segment_locators( CreateLocatorLoopSegments( ISTREAM));

      // true if "loop_segment_locators" has no loop segment locators in it
      if( loop_segment_locators.IsEmpty())
      {
        // return empty LocatorLoopDomain
        return LocatorLoopDomain();
      }

      // create loop domain locator from "loop_segment_locators"
      const LocatorLoopDomain domain_locator( loop_segment_locators);

      // return "domain_locator"
      return domain_locator;
    }

    //! @brief HandleReadMultiple can create a list of LocatorLoopDomain objects from a file and a protein model
    //! @param ISTREAM the file stream which contains the loop domain locators information
    //! @param MODEL the protein model which will be used to help create the LocatorLoopDomains
    //! @return list of LocatorLoopDomains which has been created from ISTREAM and MODEL
    util::ShPtrList< LocatorLoopDomain> HandlerLocatorLoopDomain::HandleReadMultiple
    (
      std::istream &ISTREAM, const assemble::ProteinModel &MODEL
    ) const
    {
      // create list to hold all of the LocatorLoopDomain objects
      util::ShPtrList< LocatorLoopDomain> locator_list;

      // read in the first LocatorLoopDomain from ISTREAM
      const util::ShPtr< LocatorLoopDomain> current_locator( HandleRead( ISTREAM, MODEL).Clone());

      // add "current_locator" to "locator_list"
      locator_list.PushBack( current_locator);

      // print the number of LocatorLoopDomains read in so far
      BCL_MessageDbg
      (
        "read in " + util::Format()( locator_list.GetSize()) + " locator domains so far"
      );

      // while the ISTREAM still has stuff in it
      while( !ISTREAM.eof() && ISTREAM.good() && !ISTREAM.fail())
      {
        // print the number of LocatorLoopDomains read in so far
        BCL_MessageDbg
        (
          "read in " + util::Format()( locator_list.GetSize()) + " locator domains so far"
        );

        // create the current LocatorLoopDomain from ISTREAM
        const util::ShPtr< LocatorLoopDomain> current_locator( HandleRead( ISTREAM, MODEL).Clone());

        // true if "current_locator" is not empty
        if( !current_locator->GetLoopSegments().IsEmpty())
        {
          // add "current_locator" to "locator_list"
          locator_list.PushBack( current_locator);
        }
      }

      // return the created list of LocatorLoopDomains
      return locator_list;
    }

    //! @brief CreateLocatorLoopDomainsForInteriorCoil creates a list of loop domain locators from a protein model
    //! It creates a loop domain locator for every coil sse that is in the protein model
    //! @param MODEL the model for which a list of loop domain locators will be created for all its coil sses
    //! @return a list of LocatorLoopDomains created from MODEL - one for each coild sse
    const util::ShPtrList< LocatorLoopDomain> HandlerLocatorLoopDomain::CreateLocatorLoopDomainsForInteriorCoil
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // create list which will hold the loop domain locators created from MODEL
      util::ShPtrList< LocatorLoopDomain> loop_domain_list;

      // get the list of coil sses in "MODEL". The sses will be ordered by chain and sequence
      const util::SiPtrVector< const assemble::SSE> coil_sses( MODEL.GetSSEs( biol::GetSSTypes().COIL));

      // iterate through "coil_sses"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( coil_sses.Begin()), sse_itr_end( coil_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // get the starting seq id of the sse
        const int sse_start_seq_id( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // get the ending seq id of the sse
        const int sse_end_seq_id( ( *sse_itr)->GetLastAA()->GetSeqID());

        // get the chain id of the sse
        const char chain_id( ( *sse_itr)->GetChainID());

        // true if the current sse is not the first or last sse in the sse's chain
        // this is important since only interior regions of the protein will be built
        if
        (
          sse_start_seq_id != MODEL.GetChain( chain_id)->GetSequence()->GetFirstAA()->GetSeqID() &&
          sse_end_seq_id != MODEL.GetChain( chain_id)->GetSequence()->GetLastAA()->GetSeqID()
        )
        {
          // create a LocatorLoopSegment for the sse currently denoted by "sse_itr"
          const LocatorLoopSegment current_loop_segment
          (
            assemble::LocatorSSE( chain_id, sse_start_seq_id, sse_end_seq_id),
            false //< loopsegment is not rigid
          );

          // put "current_loop_segment" into a list
          const storage::List< LocatorLoopSegment> loop_segment_list( 1, current_loop_segment);

          // create the loop domain locator
          const util::ShPtr< LocatorLoopDomain> domain_locator
          (
            new LocatorLoopDomain( loop_segment_list)
          );

          // add "domain_locator" to "loop_domain_list"
          loop_domain_list.PushBack( domain_locator);
        }
      }

      // return the list of loop domains where every coil sse is a loop domain
      return loop_domain_list;
    }

    //! @brief creates a list of loop domain locators from a protein model with bidirectionally grown loops
    //! every loop is assumed to be made up of two coils with on attached n and the other attached c terminally
    //! @param MODEL the model for which a list of loop domain locators will be created for all its coil sses
    //! @return a list of LocatorLoopDomains created from MODEL
    util::ShPtr< util::ShPtrList< LocatorLoopDomain> >
    HandlerLocatorLoopDomain::CreateBidirectionalLocatorsForInteriorCoil
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // create list which will hold the loop domain locators created from MODEL
      util::ShPtr< util::ShPtrList< LocatorLoopDomain> > loop_domain_list( new util::ShPtrList< LocatorLoopDomain>());

      const CollectorUnconnectedSSE collector_unconnected_sse_n_flex
      (
        biol::AASequenceFlexibility::e_NTerminal,
        true,
        storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
        true
      );
      const CollectorUnconnectedSSE collector_unconnected_sse_c_flex
      (
        biol::AASequenceFlexibility::e_CTerminal,
        true,
        storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
        true
      );

      const util::SiPtrList< const assemble::SSE> n_flex_coils( collector_unconnected_sse_n_flex.Collect( MODEL));
      const util::SiPtrList< const assemble::SSE> c_flex_coils( collector_unconnected_sse_c_flex.Collect( MODEL));

      // iterate over n flex coils
      for( util::SiPtrList< const assemble::SSE>::const_iterator itr( n_flex_coils.Begin()), itr_end( n_flex_coils.End()); itr != itr_end; ++itr)
      {
        // insert locator into list of loop domains
        BCL_MessageVrb( "CreateCToNLocator adding loop domain locator for sse " + ( *itr)->GetIdentification());
        const util::ShPtr< LocatorLoopDomain> cterm_nflex_locator( CreateCToNLocator( **itr, MODEL));
        loop_domain_list->PushBack( cterm_nflex_locator);
      }

      // iterate over c flex coils
      for( util::SiPtrList< const assemble::SSE>::const_iterator itr( c_flex_coils.Begin()), itr_end( c_flex_coils.End()); itr != itr_end; ++itr)
      {
        // insert locator into list of loop domains
        BCL_MessageVrb( "CreateNToCLocator adding loop domain locator for sse " + ( *itr)->GetIdentification());
        const util::ShPtr< LocatorLoopDomain> nterm_cflex_locator( CreateNToCLocator( **itr, MODEL));
        loop_domain_list->PushBack( nterm_cflex_locator);
      }

      BCL_MessageVrb( "loop domain list size is " + util::Format()( loop_domain_list->GetSize()));

      // return the list of loop domains where every coil sse is a loop domain
      return loop_domain_list;
    }

    //! @brief creates a locator loop domain that will give an N to C direction loop domain
    //! @param SSE the sse that will be used to make a loop domain locator
    //! @param MODEL the model that will be used to make LocatorLoopDomain
    //! @return ShPtr to a LocatorLoopDomain
    util::ShPtr< LocatorLoopDomain>
    HandlerLocatorLoopDomain::CreateNToCLocator( const assemble::SSE &SSE, const assemble::ProteinModel &MODEL) const
    {
      // get the starting seq id of the sse
      const int sse_start_seq_id( SSE.GetFirstAA()->GetSeqID());

      // get the ending seq id of the sse
      const int sse_end_seq_id( SSE.GetLastAA()->GetSeqID());

      // get the chain id of the sse
      const char chain_id( SSE.GetChainID());

      // create a LocatorLoopSegment for the sse currently denoted by "sse_itr"
      const LocatorLoopSegment current_loop_segment
      (
        assemble::LocatorSSE( chain_id, sse_start_seq_id, sse_end_seq_id),
        false //< loopsegment is not rigid
      );

      // put "current_loop_segment" into a list
      const storage::List< LocatorLoopSegment> loop_segment_list( 1, current_loop_segment);

      // create the loop domain locator
      const util::ShPtr< LocatorLoopDomain> domain_locator
      (
        new LocatorLoopDomain( loop_segment_list)
      );

      return domain_locator;
    }

    //! @brief creates a locator loop domain that will give an C to N direction loop domain
    //! @param SSE the sse that will be used to make a loop domain locator
    //! @param MODEL the model that will be used to make LocatorLoopDomain
    //! @return ShPtr to a LocatorLoopDomain
    util::ShPtr< LocatorLoopDomain>
    HandlerLocatorLoopDomain::CreateCToNLocator( const assemble::SSE &SSE, const assemble::ProteinModel &MODEL) const
    {
      // get the starting seq id of the sse
      const int sse_start_seq_id( SSE.GetFirstAA()->GetSeqID());

      // get the ending seq id of the sse
      const int sse_end_seq_id( SSE.GetLastAA()->GetSeqID());

      // get the chain id of the sse
      const char chain_id( SSE.GetChainID());

      // create a LocatorLoopSegment for the sse currently denoted by "sse_itr"
      const LocatorLoopSegment current_loop_segment
      (
        assemble::LocatorSSE( chain_id, sse_start_seq_id, sse_end_seq_id),
        false //< loopsegment is not rigid
      );

      // put "current_loop_segment" into a list
      const storage::List< LocatorLoopSegment> loop_segment_list( 1, current_loop_segment);

      // get anchor aa locator
      const assemble::LocatorAA aa_anchor_locator( SSE.GetChainID(), SSE.GetFirstAA()->GetSeqID());

    ///////////////////////////
    // create pseudo residue //
    ///////////////////////////

      const util::SiPtr< const biol::AABase> pseudo_residue_basis
      (
        assemble::LocatorAA( SSE.GetChainID(), SSE.GetFirstAA()->GetSeqID() - 1).Locate( MODEL)
      );

      BCL_Assert( pseudo_residue_basis.IsDefined(), "pseudo residue basis is not defined");

      // create a ShPtr to the residue that is the basis for the pseudo residue
      const util::ShPtr< biol::AABase> dummy_pseudo_residue_basis( pseudo_residue_basis->Clone());

      // define phi and psi
      MutationResidue temp( dummy_pseudo_residue_basis, util::ShPtr< biol::AABase>(), util::ShPtr< biol::AABase>());
      storage::VectorND< 2, double> phi_psi
      (
        PhiPsiGeneratorRamachandran::GetDefaultInstance().operator()( temp)
      );

      // create the loop domain locator
      const util::ShPtr< LocatorLoopDomain> domain_locator
      (
        new LocatorLoopDomain( loop_segment_list, false)
      );

      return domain_locator;
    }

    //! @brief GetFormat gives the format that this handler needs in order to work
    //! @return string which describes the format needed by this handler in order for it to work
    std::string HandlerLocatorLoopDomain::GetFormat() const
    {
      return "Every loop domain should be separated by \"DomainLocatorEnd\"\n"
        "example\n"
        "'A' 114 114 false\n"
        "DomainLocatorEnd\n"
        "'A' 124 125 false\n"
        "DomainLocatorEnd\n"
        "'A' 126 134 false\n"
        "'A' 135 136 false\n"
        "DomainLocatorEnd";
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerLocatorLoopDomain::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &HandlerLocatorLoopDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief CreateLocatorLoopSegments creates a list of loop segment locators from an istream
    //! @param ISTREAM the file stream which contains the information to create a list of loop segment locators
    //! @return a list of LocatorLoopSegments created from ISTREAM
    const storage::List< LocatorLoopSegment> HandlerLocatorLoopDomain::CreateLocatorLoopSegments
    (
      std::istream &ISTREAM
    ) const
    {
      // create list which will hold the loop segment locators created from ISTREAM
      storage::List< LocatorLoopSegment> loop_segment_list;

      // create string "current" line which will hold each of the current lines as read in from ISTREAM
      std::string current_line;

      // get the current line from ISTREAM and put it into "current_line"
      std::getline( ISTREAM, current_line);

      // trim "current_line"
      util::TrimString( current_line);

      // print "current_line"
      BCL_MessageDbg( "current_line |" + current_line + "|");

      // go until the end of file has been reached or "s_DomainEndString" has been read in
      while( current_line != s_DomainEndString && !ISTREAM.eof())
      {
        // create LoopSegment and add it to "loop_segment_list"
        loop_segment_list.PushBack( HandlerLocatorLoopSegment().HandleRead( current_line));

        // get the next line from ISREAM and put it into "current_line"
        std::getline( ISTREAM, current_line);

        // trim "current_line"
        util::TrimString( current_line);

        // print "current_line
        BCL_MessageDbg( "current_line |" + current_line + "|");
      }

      // return the list of created loop segment locators
      return loop_segment_list;
    }

    //! @brief CreateNTerminalSSELocator creates the locator to the nterminal anchor sse of a loop domain
    //! @param LOOP_SEGMENT_LOCATORS list of loop segment locators which make up a loop domain
    //! @param MODEL the protein model which the loop segment locators and the loop domain refer
    //! @return LocatorSSE which can locate the nterminal anchor sse of the loop domain
    const assemble::LocatorSSE HandlerLocatorLoopDomain::CreateNTerminalSSELocator
    (
      const storage::List< LocatorLoopSegment> &LOOP_SEGMENT_LOCATORS, const assemble::ProteinModel &MODEL
    ) const
    {
      // order the loop segments by sequence
      const storage::Set< LocatorLoopSegment, LoopSegmentSequenceOrder> ordered_segments
      (
        LOOP_SEGMENT_LOCATORS.Begin(), LOOP_SEGMENT_LOCATORS.End()
      );

      // make sure "ordered_segments" and "LOOP_SEGMENT_LOCATORS" have the same size
      BCL_Assert( ordered_segments.GetSize() == LOOP_SEGMENT_LOCATORS.GetSize(), "sizes differ");

      // get the first loop segment locator according to sequence
      const LocatorLoopSegment &locator_first_loop_segment( *ordered_segments.Begin());

      const char chain_id( locator_first_loop_segment.GetLocatorSSE().GetChainID());

      // get the sses from "MODEL" in a set sorted by sequence
      const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &model_sses( MODEL.GetChain( chain_id)->GetData());

      // locate the first segment in MODEL
      const LoopSegment first_loop_segment( locator_first_loop_segment.Locate( MODEL));

      // get an iterator to the first loop segment as it is in "MODEL"
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator itr_first_loop_sse
      (
        model_sses.Find( first_loop_segment.GetSSE())
      );

      // make sure that the sse could be found in model
      BCL_Assert
      (
        itr_first_loop_sse != model_sses.End(), "could not find " + first_loop_segment.GetSSE()->GetIdentification()
      );

      // make sure that the first loop segment is not the first sse in the model
      BCL_Assert
      (
        itr_first_loop_sse != model_sses.Begin(), "first loop segment is the last sse in the protein model"
       );

        // move "itr_first_loop_sse" to the previous sse
        --itr_first_loop_sse;

        // get the chain id and first and last seq ids of the sse denoted by "itr_first_loop_sse"
        const char nterminal_sse_chain_id( ( *itr_first_loop_sse)->GetChainID());
        const int  nterminal_sse_start_seq_id( ( *itr_first_loop_sse)->GetFirstAA()->GetSeqID());
        const int  nterminal_sse_last_seq_id( ( *itr_first_loop_sse)->GetLastAA()->GetSeqID());

        // create locator for locating the nterminal anchor sse of the loop domain
        const assemble::LocatorSSE nterminal_sse_locator
        (
          nterminal_sse_chain_id, nterminal_sse_start_seq_id, nterminal_sse_last_seq_id
        );

        // return the locator
        return nterminal_sse_locator;
      }

  } // namespace fold
} // namespace bcl
