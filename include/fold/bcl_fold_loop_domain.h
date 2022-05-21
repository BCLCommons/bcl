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

#ifndef BCL_FOLD_LOOP_DOMAIN_H_
#define BCL_FOLD_LOOP_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_segment.h"
#include "bcl_fold_loop_segment_sequence_order.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "coord/bcl_coord_cyclic_coordinate_descent.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopDomain
    //! @brief is for representing loop domains.
    //! @details A loop domain is a collection of LoopSegments which are
    //! secondary structure elements where the dihedral angles of each SSE are fixed (rigid) or allowed to be changed.
    //! The LoopDomain has a LoopPseudoResidue which represents the cterminal anchor residue of the LoopDomain. This
    //! is the residue the LoopDomain should connect to.
    //! The LoopDomain also has a locator to the sse just before it in sequence. This is the nterminal anchor sse to
    //! which the LoopDomain is connected on the n-terminal side.
    //!
    //!
    //! @see @link example_fold_loop_domain.cpp @endlink
    //! @author alexanns
    //! @date Aug 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopDomain :
      public util::ObjectInterface
    {

    protected:

    //////////
    // data //
    //////////

      //! set LoopSegments ordered by sequence that make up the loop domain
      storage::Set< LoopSegment, LoopSegmentSequenceOrder> m_LoopSegments;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief return command line flag for specifying loop domains via a file
      //! @return command line flag for specifying loop domains via a file
      static util::ShPtr< command::FlagInterface> &GetFlagLoopDomainFilename();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopDomain();

      //! @brief constructor taking parameters to create member variables
      //! @param SEGMENTS list of LoopSegments which will be used to create "m_LoopSegments"
      LoopDomain
      (
        const storage::List< LoopSegment> &SEGMENTS
      );

      //! @brief copy constructor
      //! @param OTHER loop domain to be copied
      LoopDomain( const LoopDomain &OTHER);

      //! @brief Clone function
      //! @return pointer to new LoopSegment
      LoopDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gives the directionality of the loop domain
      //! @return sequence direction in this case e_CTerminal since changes are propagated cterminally
      virtual biol::AASequenceFlexibility::SequenceDirection GetSequenceDirection() const;

      //! @brief GetSegments gives the set LoopSegments ordered by sequence that make up the loop domain
      //! @return the set LoopSegments ordered by sequence that make up the loop domain
      const storage::Set< LoopSegment, LoopSegmentSequenceOrder> &GetSegments() const;

      //! @brief finds the segment and replaces it with the given segment
      //! @param LOOP_SEGMENT segment that will replace the currently matching segment
      //! @return pair of iterator pointing to iterator where the replaced segment is and bool indicating success or not
      std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool>
      ReplaceSegment( const LoopSegment &LOOP_SEGMENT);

      //! @brief finds the segments and replaces them with the given segments
      //! @param LOOP_SEGMENTS segments that will replace the currently matching segments
      //! @return returns iter and bool to last replaced element
      std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool>
      ReplaceSegment( const storage::List< LoopSegment> &LOOP_SEGMENTS);

      //! @brief FindSegment provides an iterator to the segment containing a residue of interest
      //! @param RESIDUE the residue of interest for which the segment it is in is desired to be known
      //! @return iterator to the LoopSegment that contains RESIDUE
      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator FindSegment
      (
        const biol::AABase &RESIDUE
      ) const;

      //! @brief finds a residue of interest in the LoopDomain
      //! @param RESIDUE the residue of interest which should be found
      //! @return ShPtr pointing to the residue of interest
      util::ShPtr< biol::AABase> FindResidue( const biol::AABase &RESIDUE) const;

      //! @brief finds a residue of interest in the LoopDomain as described by a LocatorAA
      //! @param LOCATOR_AA the LocatorAA which describes the residue of iterest
      //! @return ShPtr pointing to the residue of interest as described by "LOCATOR_AA"
      util::ShPtr< biol::AABase> FindResidue( const assemble::LocatorAA &LOCATOR_AA) const;

      //! @brief SetResidue finds the pointer to residue of interest in the LoopDomain and sets the pointer to given ptr
      //! @param RESIDUE the ptr to residue of interest which should be found
      void SetPtrToResidue( const util::ShPtr< biol::AABase> &RESIDUE);

      //! @brief finds the pointer to residue of interest in the LoopDomain and gives iterator to it
      //! @param RESIDUE the ptr to residue of interest which should be found
      //! @return iterator to residue
      biol::AASequence::const_iterator FindResidueIterator( const biol::AABase &RESIDUE) const;

      //! @brief FindSegment gives iterator to the segment containing an AA of interest described by chain and seq id
      //! @param RESI_CHAIN_ID the chain id of the residue of interest
      //! @param RESI_SEQ_ID the seq id of the residue of interest
      //! @return iterator to the LoopSegment that contains the residue described by "RESI_CHAIN_ID" and "RESI_SEQ_ID"
      storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator FindSegment
      (
        const char RESI_CHAIN_ID, const int RESI_SEQ_ID
      ) const;

      //! @brief FindResidue finds a residue of interest in an SSE of the loop domain as described by a chain and seq id
      //! @param SEGMENT_SSE the sse which contains the residue of interest
      //! @param CHAIN_ID the chain id of the residue of interest
      //! @param SEQ_ID the seq id of the residue of interest
      //! @return ShPtr pointing to the residue of interest
      util::ShPtr< biol::AABase> FindResidue
      (
        const assemble::SSE &SEGMENT_SSE, const char CHAIN_ID, const int SEQ_ID
      ) const;

      //! @brief GetChainID provides the chain id of the LoopDomain. Since residues must be from the same chain, it
      //! can get the chain id from any residue.
      //! @return the chain id of the loop domain
      const char &GetChainID() const;

      //! @brief returns the residues in the loop domain in sequence order including anchor sse anchor residue and the
      //!        created psuedo residue
      //! @return map with the amino acids and a bool indicating whether or not the aa is part of a rigid segment or not
      virtual storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > GetResidues() const;

      //! @brief gives the residue that is attached to the anchor sse
      //!        i.e. the residue that is closest to point of attachment
      //! @return ShPtr to residue that is attached to the anchor sse
      virtual const util::ShPtr< biol::AABase> &GetMostProximalLoopSegmentAA() const;

      //! @brief gives the loop segment that is most distant in sequence to attachment to the anchor sse
      //!        this is the sse that the pseudo residue attaches to
      //! @return sse that the pseudo residue attaches to and is the most distant in sequence from anchor sse
      virtual const assemble::SSE &GetMostDistalLoopSegment() const;

      //! @brief gives a ShPtr to the residue following the given residue in sequence
      //! @param AA the residue for which the following residue will be found
      //! @return ShPtr to residue that follows AA in sequence
      util::ShPtr< biol::AABase> FindFollowingResidue( const biol::AABase &AA) const;

      //! @brief gives a ShPtr to the residue preceding the given residue in sequence
      //! @param AA the residue for which the preceding residue will be found
      //! @return ShPtr to residue that precedes AA in sequence
      util::ShPtr< biol::AABase> FindPreviousResidue( const biol::AABase &AA) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the psi of a residue in the loop domain
      //! @param AA the amino acid that will be found and mutated
      //! @param PSI the psi angle that residue will be set to
      virtual void SetPsi( const biol::AABase &AA, const double PSI);

      //! @brief sets the phi of a residue in the loop domain
      //! @param AA the amino acid that will be found and mutated
      //! @param PHI the phi angle that residue will be set to
      virtual void SetPhi( const biol::AABase &AA, const double PHI);

      //! @brief return the points that can be used for CCD
      //! @param PROTEIN_MODEL the model, to find the adjacent sse in to connect to
      //! @return target and moving points
      storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> TargetAndMovingPointsForCCD
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief places the sses from a loop domain into a protein model
      //! @param PROTEIN_MODEL the protein model whose sses will be replaced
      //! @return ShPtr to ProteinModel which has had its sses replaced with the sses of LOOP_DOMAIN
      util::ShPtr< assemble::ProteinModel> UpdateProteinModel
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief determines if the seq id and chain id of two residues are the same
      //! @param RESI_A first residue
      //! @param RESI_B second residue
      //! @return boolean true if the seqid and chain ids match - false otherwise
      static bool ResiduesMatch( const biol::AABase &RESI_A, const biol::AABase &RESI_B);

    }; // class LoopDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_DOMAIN_H_ 
