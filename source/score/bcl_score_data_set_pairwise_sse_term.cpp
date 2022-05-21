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
#include "score/bcl_score_data_set_pairwise_sse_term.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSSETerm::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSSETerm())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSSETerm::GetDefaultScheme()
    {
      static const std::string s_scheme( "sse_term");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSETerm::DataSetPairwiseSSETerm( const std::string &SCHEME) :
      m_SSEPool(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param SSE_POOL pool to use as sse definitions
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSETerm::DataSetPairwiseSSETerm
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL, const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSSETerm
    DataSetPairwiseSSETerm *DataSetPairwiseSSETerm::Clone() const
    {
      return new DataSetPairwiseSSETerm( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSSETerm::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSSETerm::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseSSETerm::operator()( const restraint::DataSetPairwise &DATA) const
    {
      if( DATA.IsEmpty())
      {
        return 0;
      }

      // get random non overlapping set of sses from the pool
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
      (
        m_SSEPool->GetRandomNonOverlappingSet()
      );

      // number of spin labels which is two times the dataset size
      const size_t l( DATA.GetSize() * 2);

      // number of sses
      const size_t s( sses.GetSize());

      // Q is the number of spin labels per SSE - these are the acceptable integer values of Q
      const size_t Q_prime( l / s);
      const size_t Q_prime_prime( Q_prime + 1);

      // the remainder of l / s
      const size_t R( l % s);

      // to keep track of how many labels are in each sse
      const storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> sse_count
      (
        GetSSECounts( sses, DATA)
      );

      // calculate Ssse(L) which is the secondary structure score component L
      const double s_sse_L( CalculateSSEScoreComponentL( sse_count, Q_prime_prime, l));

      // calculate Ssse(S) which is the secondary structure score component S
      const double s_sse_S( CalculateSSEScoreComponentS( sse_count, Q_prime, Q_prime_prime, R, s));

      const double composite_score( -( s_sse_L + s_sse_S) / 2.0);

      // return the score
      return composite_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSSETerm::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEPool, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseSSETerm::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief for each sse gives the number of data points that fall within that sse
    //! @param SSES the sses which will be used to count how many data points fall within them
    //! @param DATA the data set which will be used to count how many data points fall within the SSES
    //! @return map which has for each sse the number of data points that fall within that sse
    storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap>
    DataSetPairwiseSSETerm::GetSSECounts
    (
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> &SSES,
      const restraint::DataSetPairwise &DATA
    )
    {
      // to keep track of how many labels are in each sse
      storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> sse_count;

      // iterate through the sses to get the number of labels in each sse
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( SSES.Begin()), sse_itr_end( SSES.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const util::SiPtr< const assemble::SSE> &sse( *sse_itr);
        // add the sse to the counting map
        std::pair
        <
          storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap>::iterator, bool
        > insert_status( sse_count.Insert( std::pair< util::SiPtr< const assemble::SSE>, size_t>( sse, 0)));

        BCL_Assert( insert_status.second, "could not insert " + sse->GetIdentification());

        // iterate through the data set to see which sse each point of each data pair is in
        for
        (
          restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
          data_itr != data_itr_end;
          ++data_itr
        )
        {
          // true if the first data point is in the sse
          if( data_itr->First()->IsWithin( *sse))
          {
            // if so add a count for the current SSE
            ++insert_status.first->second;
          }
          // true if the second data point is in the sse
          if( data_itr->Second()->IsWithin( *sse))
          {
            // if so add a count for the current SSE
            ++insert_status.first->second;
          }
        }
      }

      // return the map of sses and how many data points from the data set are in each sse
      return sse_count;
    }

    //! @brief calculates the average percentage of points positioned in each SSE up to the ideal value, Qprimeprime
    //! @param SSE_COUNT map with a count for each sse of the number of points in each sse
    //! @param Q_PRIME_PRIME see reference
    //! @param NUMBER_OF_POINTS the number of points in the data set i.e. 2 * (dataset size) ( this is l in reference)
    //! @return double which is the L component of the sse term score
    double DataSetPairwiseSSETerm::CalculateSSEScoreComponentL
    (
      const storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> &SSE_COUNT,
      const size_t Q_PRIME_PRIME,
      const size_t NUMBER_OF_POINTS
    )
    {
      double score( 0);

      // iterate through the map
      for
      (
        storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( SSE_COUNT.Begin()), sse_itr_end( SSE_COUNT.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // number of labels in the current sse
        const size_t li( sse_itr->second);
        const size_t min( std::min( li, Q_PRIME_PRIME));
        score += min;
      }

      // divide by number of points
      score /= double( NUMBER_OF_POINTS);

      return score;
    }

    //! @brief calculates S component of SSE score, derived from frachtion of SSEs containing exactly
    //! @param SSE_COUNT map with a count for each sse of the number of points in each sse
    //! @param Q_PRIME_PRIME see reference
    //! @param REMAINDER "R" in reference
    //! @param NUMBER_SSES "s" in reference
    //! @return double which is component S of the sse score
    double DataSetPairwiseSSETerm::CalculateSSEScoreComponentS
    (
      const storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> &SSE_COUNT,
      const size_t Q_PRIME,
      const size_t Q_PRIME_PRIME,
      const size_t REMAINDER,
      const size_t NUMBER_SSES
    )
    {
      // the number of sses with Q_PRIME data points
      size_t E_prime( 0);

      // the number of sses with Q_PRIME_PRIME data points
      size_t E_prime_prime( 0);

      // iterate through the sses and their data point counts to determine E_prime and E_prime_prime
      for
      (
        storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( SSE_COUNT.Begin()), sse_itr_end( SSE_COUNT.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const size_t current_count( sse_itr->second);

        if( current_count == Q_PRIME)
        {
          ++E_prime;
        }

        if( current_count == Q_PRIME_PRIME)
        {
          ++E_prime_prime;
        }
      }

      // min( Eprime,( s - R))
      size_t min_E_prime( std::min( E_prime, ( NUMBER_SSES - REMAINDER)));

      // min( E_prime_prime, R)
      size_t min_E_prime_prime( std::min( E_prime_prime, REMAINDER));

      // S component of SSE score
      const double SsseS( double( min_E_prime + min_E_prime_prime) / double( NUMBER_SSES));

      return SsseS;
    }

  } // namespace score
  
} // namespace bcl
