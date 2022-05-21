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
#include "fold/bcl_fold_mutate_domain_sse_pair_trim.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateDomainSSEPairTrim::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDomainSSEPairTrim())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDomainSSEPairTrim::MutateDomainSSEPairTrim() :
      m_Scheme(),
      m_MinSSESizes(),
      m_NumberResiduesToTrim()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SCHEME the scheme for this mutate
    //! @param MIN_SSE_SIZES minimum SSE sizes to be allowed when shrinking
    //! @param NUMBER_RESIS_TO_TRIM the total number of residues that will possibly be trimmed
    MutateDomainSSEPairTrim::MutateDomainSSEPairTrim
    (
      const size_t NUMBER_RESIS_TO_TRIM,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_MinSSESizes( MIN_SSE_SIZES),
      m_NumberResiduesToTrim( NUMBER_RESIS_TO_TRIM)
    {

    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainSSEPairTrim
    MutateDomainSSEPairTrim *MutateDomainSSEPairTrim::Clone() const
    {
      return new MutateDomainSSEPairTrim( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDomainSSEPairTrim::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDomainSSEPairTrim::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a domain
    //! @param DOMAIN domain which will be mutated
    //! @return MutateResult with the mutated domain
    math::MutateResult< assemble::Domain>
    MutateDomainSSEPairTrim::operator()( const assemble::Domain &THIS_DOMAIN) const
    {
      // static empty domain
      static util::ShPtr< assemble::Domain> s_empty_domain;

      // get the sses from this domain - these are ordered by sequence by the domain
      const util::SiPtrVector< const assemble::SSE> domain_sses( THIS_DOMAIN.GetSSEs());

      // true if there are not two sses in the domain
      if( domain_sses.GetSize() != 2)
      {
        return math::MutateResult< assemble::Domain>( s_empty_domain, *this);
      }

      // get the combination of trimming that gives the lowest ratio of euclidian distance to connecting residues
      storage::VectorND< 2, util::ShPtr< assemble::SSE> > best_trimmed_sses
      (
        GetBestTrimmedSSEs( *domain_sses.FirstElement(), *domain_sses.LastElement())
      );

      // true if either of the best sses are not defined
      if( !best_trimmed_sses.First().IsDefined() || !best_trimmed_sses.Second().IsDefined())
      {
        return math::MutateResult< assemble::Domain>( s_empty_domain, *this);
      }

      util::ShPtr< assemble::Domain> new_domain( THIS_DOMAIN.Clone());

      // remove the two old sses
      BCL_Assert
      (
        new_domain->Remove( *domain_sses.FirstElement()),
        "could not remove sse " + domain_sses.FirstElement()->GetIdentification()
      );
      BCL_Assert
      (
        new_domain->Remove( *domain_sses.LastElement()),
        "could not remove sse " + domain_sses.LastElement()->GetIdentification()
      );

      // insert both new sses
      BCL_Assert
      (
        new_domain->Insert( best_trimmed_sses.First()),
        "could not Insert sse " + best_trimmed_sses.First()->GetIdentification()
      );
      BCL_Assert
      (
        new_domain->Insert( best_trimmed_sses.Second()),
        "could not Insert sse " + best_trimmed_sses.Second()->GetIdentification()
      );

      return math::MutateResult< assemble::Domain>( new_domain, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainSSEPairTrim::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_MinSSESizes, ISTREAM);
      io::Serialize::Read( m_NumberResiduesToTrim, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainSSEPairTrim::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_MinSSESizes, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_NumberResiduesToTrim, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateDomainSSEPairTrim::GetBestTrimmedSSEs
    (
      const assemble::SSE &SSE_N_TERMINAL, const assemble::SSE &SSE_C_TERMINAL
    ) const
    {
      // to hold the best combination of trimmed sses
      util::ShPtr< assemble::SSE> best_sse_a;
      util::ShPtr< assemble::SSE> best_sse_b;

      // to keep track of the smallest ratio of distance per residue observed
      double best_dist_per_resi( std::numeric_limits< double>::max());

      // iterate through the possible number of residues to trim
      for( size_t aa_removed_a( 0); aa_removed_a <= m_NumberResiduesToTrim; ++aa_removed_a)
      {
        // trim the number of residues off of first sse. c terminus of sse coming first in sequence will be changed
        util::ShPtr< assemble::SSE> trim_sse_n_terminal
        (
          TrimSSE( SSE_N_TERMINAL, biol::AASequenceFlexibility::e_CTerminal, aa_removed_a)
        );

        if( !trim_sse_n_terminal.IsDefined())
        {
          continue;
        }

        for
        (
          size_t aa_removed_b( 0); aa_removed_b <= ( m_NumberResiduesToTrim - aa_removed_a);
          ++aa_removed_b
        )
        {
          // trim number of residues off of second sse. n terminus of sse coming second in sequence will be changed
          util::ShPtr< assemble::SSE> trim_sse_c_terminal
          (
            TrimSSE( SSE_C_TERMINAL, biol::AASequenceFlexibility::e_NTerminal, aa_removed_b)
          );

          if( !trim_sse_c_terminal.IsDefined())
          {
            continue;
          }

          // get the ratio of euclidian distance to number of residues
          const double dist_per_resi( CalculateDistancePerResidue( *trim_sse_n_terminal, *trim_sse_c_terminal));

          BCL_MessageDbg( "trim_sse_n_terminal " + trim_sse_n_terminal->GetIdentification());
          BCL_MessageDbg( "trim_sse_c_terminal " + trim_sse_c_terminal->GetIdentification());
          BCL_MessageDbg( "dist_per_resi " + util::Format()( dist_per_resi));

          // true if this is the best combination seen so far
          if( dist_per_resi < best_dist_per_resi)
          {
            best_dist_per_resi = dist_per_resi;
            best_sse_a = trim_sse_n_terminal;
            best_sse_b = trim_sse_c_terminal;
          }
        }
      }

      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > best_trimmed_sses( best_sse_a, best_sse_b);

      return best_trimmed_sses;
    }

    util::ShPtr< assemble::SSE> MutateDomainSSEPairTrim::TrimSSE
    (
      const assemble::SSE &SSE, const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION,
      const size_t NUM_RESI_TO_REMOVE
    ) const
    {
      // true if the trimmed sse will be too short according to the desired minimum sse sizes
      storage::Map< biol::SSType, size_t>::const_iterator type_size( m_MinSSESizes.Find( SSE.GetType()));
      BCL_Assert
      (
        type_size != m_MinSSESizes.End(), "could not find ss type " + SSE.GetType().GetName() + " in min sse sizes "
        + util::Format()( m_MinSSESizes)
      );

      if( int( SSE.GetSize()) - int( NUM_RESI_TO_REMOVE) < int( type_size->second))
      {
        return util::ShPtr< assemble::SSE>();
      }

      // index of starting subsequence for trimmed sse
      const size_t start_pos
      (
        SEQUENCE_DIRECTION == biol::AASequenceFlexibility::e_NTerminal ? NUM_RESI_TO_REMOVE : 0
      );

      // the length of the subsequence for the trimmed sse
      const size_t length( SSE.GetSize() - NUM_RESI_TO_REMOVE);

      // the sequence for the trimed sse
      const biol::AASequence sequence( SSE.SubSequence( start_pos, length));

      // the trimmed sse
      const util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( sequence, SSE.GetType()));

      // return the new sse
      return new_sse;
    }

    double MutateDomainSSEPairTrim::CalculateDistancePerResidue
    (
      const assemble::SSE &SSE_N_TERMINAL, const assemble::SSE &SSE_C_TERMINAL
    ) const
    {
      // get c terminal residue of n terminal sse
      const biol::AABase &c_terminal_aa( *SSE_N_TERMINAL.GetLastAA());

      // get n terminal residue of c terminal sse
      const biol::AABase &n_terminal_aa( *SSE_C_TERMINAL.GetFirstAA());

      const double euclidian_distance( biol::GetPeptideBondLength( c_terminal_aa, n_terminal_aa));
      const size_t sequence_distance( biol::CalculateSequenceDistance( SSE_N_TERMINAL, SSE_C_TERMINAL));
      BCL_MessageDbg( "euclidian_distance " + util::Format()( euclidian_distance));
      BCL_MessageDbg( "sequence_distance " + util::Format()( sequence_distance));

      return euclidian_distance / sequence_distance;
    }

  } // namespace fold
} // namespace bcl
