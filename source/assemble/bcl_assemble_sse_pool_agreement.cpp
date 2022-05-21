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
#include "assemble/bcl_assemble_sse_pool_agreement.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolAgreement::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPoolAgreement())
    );

    //! @brief default scheme as string
    //! @return string for default scheme
    const std::string &SSEPoolAgreement::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "sse_pool_agreement");

      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param SCHEME scheme for function
    //! @param LOG before adding differences, take log
    SSEPoolAgreement::SSEPoolAgreement( const bool LOG, const std::string &SCHEME) :
      m_Scheme( SCHEME),
      m_Log( LOG)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolAgreement
    SSEPoolAgreement *SSEPoolAgreement::Clone() const
    {
      return new SSEPoolAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPoolAgreement::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief virtual operator taking two ARGUMENTs and returning a t_ResultType object
    //! @param POOL the pool to be evaluated
    //! @param POOL_TEMPLATE the template to evaluate against
    //! @return double that expresses the agreement of the pool to the template
    double SSEPoolAgreement::AgreementToTemplate( const SSEPool &POOL, const SSEPool &POOL_TEMPLATE) const
    {
      const double missing_sse_penalty_multiplier( 3.0); // multiplier
      const int tolerate_shift_residues( 1); // number of residues difference in beginning and end, that is not penalized

      // initialize agreement
      double agreement( 0.0);

      // iterate over the sses in the template
      for
      (
        SSEPool::const_iterator temp_itr( POOL_TEMPLATE.Begin()), temp_itr_end( POOL_TEMPLATE.End());
        temp_itr != temp_itr_end; ++temp_itr
      )
      {
        // if the current SSE is empty
        if( ( *temp_itr)->GetSize() == 0)
        {
          // warn user and skip this SSE
          BCL_MessageCrt( "Empty SSE found " + ( *temp_itr)->GetIdentification());
          continue;
        }

        // ignore coils for agreement evaluation
        if( !( *temp_itr)->GetType()->IsStructured())
        {
          continue;
        }

        // find overlapping sses in the POOL in question include identical, exclude different type
        const util::SiPtrList< const SSE> overlap_sses( POOL.GetOverlappingSSEs( **temp_itr, false, true));

        // number of overlapping
        const size_t nr_overlapping( overlap_sses.GetSize());

        // no sse found - total miss
        if( nr_overlapping == 0)
        {
          BCL_MessageDbg
          (
            "no overlapping sse found for template: " + ( *temp_itr)->GetIdentification()
          );
          double penalty( ( *temp_itr)->GetSize());
          if( m_Log)
          {
            penalty = std::log( penalty + 1);
          }
          agreement += missing_sse_penalty_multiplier * penalty;
          continue;
        }

        // initialize this agreement
        double this_agreement( 0.0);

        // iterate over overlap
        for
        (
          util::SiPtrList< const SSE>::const_iterator over_itr( overlap_sses.Begin()), over_itr_end( overlap_sses.End());
          over_itr != over_itr_end;
          ++over_itr
        )
        {
          // calculate the overlap and store it
          storage::VectorND< 2, int> this_overlap( Overlap( **over_itr, **temp_itr));

          // calculate the score and sum it up
          double shift_left(   std::max( 0, math::Absolute( this_overlap.First()) - tolerate_shift_residues));
          double shift_right(  std::max( 0, math::Absolute( this_overlap.Second()) - tolerate_shift_residues));
          double shift_length( math::Absolute( this_overlap.First() + this_overlap.Second()));

          if( m_Log)
          {
            shift_left   = std::log( shift_left + 1);
            shift_right  = std::log( shift_right + 1);
            shift_length = std::log( shift_length + 1);
          }

          BCL_MessageDbg
          (
            "Agreement between " + ( *temp_itr)->GetIdentification() + " and " +
            ( *over_itr)->GetIdentification() + " ==> " +
            util::Format()( shift_left) + "\t" +
            util::Format()( shift_right) + "\t" +
            util::Format()( shift_length)
          );

          // add the agreement scores
          this_agreement += shift_left + shift_right + shift_length;
        }

        // normalize this agreement sum by number of overlap and add it to the global agreement value
        agreement += this_agreement / double( nr_overlapping);
      }

      return agreement;
    }

    //! @brief Q3 score for two SSEPools
    //! @param POOL_A the pool a
    //! @param POOL_B the pool b
    //! @return number_correct_sse/total_number_res * 100%
    double SSEPoolAgreement::Q3Score( const SSEPool &POOL_A, const SSEPool &POOL_B) const
    {
      storage::Map< util::SiPtr< const biol::AAData>, biol::SSType> per_res_sstype_a;
      for( SSEPool::const_iterator sse_itr( POOL_A.Begin()), sse_itr_end( POOL_A.End()); sse_itr != sse_itr_end; ++sse_itr)
      {
        const biol::SSType &current_type( ( *sse_itr)->GetType());
        for( SSE::const_iterator aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End()); aa_itr != aa_itr_end; ++aa_itr)
        {
          if( !per_res_sstype_a.Insert( std::pair< util::SiPtr< const biol::AAData>, biol::SSType>( ( *aa_itr)->GetData(), current_type)).second)
          {
            BCL_MessageCrt( "supplied pool with overlapping sse definitions at res: " + ( *aa_itr)->GetIdentification())
            return util::GetUndefined< double>();
          }
        }
      }

      int total_number( per_res_sstype_a.GetSize());
      int correct_number( 0);

      for( SSEPool::const_iterator sse_itr( POOL_B.Begin()), sse_itr_end( POOL_B.End()); sse_itr != sse_itr_end; ++sse_itr)
      {
        const biol::SSType &current_type( ( *sse_itr)->GetType());
        for( SSE::const_iterator aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End()); aa_itr != aa_itr_end; ++aa_itr)
        {
          const storage::Map< util::SiPtr< const biol::AAData>, biol::SSType>::const_iterator find_itr( per_res_sstype_a.Find( ( *aa_itr)->GetData()));

          // residue in pool a that is not present in the pool b
          if( find_itr == per_res_sstype_a.End())
          {
            ++total_number;
          }
          else if( find_itr->second == current_type)
          {
            ++correct_number;
          }
        }
      }

      return double( 100 * correct_number) / double( total_number);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking two pools and returning the symmetric agreement, by taking the average of the
    //!        agreement if each pool is treated as template to the other
    //! @param POOL_A the pool a
    //! @param POOL_B the pool b
    //! @return double that expresses the agreement between two pools
    double SSEPoolAgreement::operator()( const SSEPool &POOL_A, const SSEPool &POOL_B) const
    {
      const double agreement_a_b( AgreementToTemplate( POOL_A, POOL_B));
      const double agreement_b_a( AgreementToTemplate( POOL_B, POOL_A));
      return agreement_a_b + agreement_b_a;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolAgreement::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolAgreement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the overlap between two sequences
    //! @param SEQUENCE the sequence of interest
    //! @param SEQUENCE_TEMPLATE the reference sequence
    //! @return pair where first is the overlap on the left (first seq id difference) and the right (last seqid difference),
    //!         where the sign indicates if the SEQUENCE is shorter (negative) or longer (positive) then the template, sum
    //!         of those two is the overall length difference
    storage::VectorND< 2, int> SSEPoolAgreement::Overlap
    (
      const biol::AASequence &SEQUENCE,
      const biol::AASequence &SEQUENCE_TEMPLATE
    )
    {
      BCL_Assert( SEQUENCE.GetSize() > 0 && SEQUENCE_TEMPLATE.GetSize() > 0, "empty sequences supplied");
      return
        storage::VectorND< 2, int>
        (
          SEQUENCE_TEMPLATE.GetFirstAA()->GetSeqID() - SEQUENCE.GetFirstAA()->GetSeqID(),
          SEQUENCE.GetLastAA()->GetSeqID() - SEQUENCE_TEMPLATE.GetLastAA()->GetSeqID()
        );
    }

  } // namespace assemble
} // namespace bcl
