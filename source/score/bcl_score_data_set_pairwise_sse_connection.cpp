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
#include "score/bcl_score_data_set_pairwise_sse_connection.h"

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
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSSEConnection::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSSEConnection())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSSEConnection::GetDefaultScheme()
    {
      static const std::string s_scheme( "sse_connection");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSEConnection::DataSetPairwiseSSEConnection( const std::string &SCHEME) :
      m_SSEPool(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param SSE_POOL pool to use as sse definitions
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSSEConnection::DataSetPairwiseSSEConnection
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL, const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSSEConnection
    DataSetPairwiseSSEConnection *DataSetPairwiseSSEConnection::Clone() const
    {
      return new DataSetPairwiseSSEConnection( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSSEConnection::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSSEConnection::GetScheme() const
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
    double DataSetPairwiseSSEConnection::operator()( const restraint::DataSetPairwise &DATA) const
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

      BCL_Assert( !sses.IsEmpty(), "no sses");

      const storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> sse_pair_counts
      (
        GetSSEPairCounts( sses, DATA)
      );

      // the number of sses
      const size_t s( sses.GetSize());

      // the number of restraints
      const size_t r( DATA.GetSize());

      // the number of sse pairs
      const size_t p( s * ( s - 1) / 2);

      // minimum acceptable integer value for C, the ideal number of connections for each SSE pair
      const size_t C_prime( r / p);

      // maximum acceptable integer value for C, the ideal number of connections for each SSE pair
      // if evenly divisible is it necessary to increment
      const size_t C_prime_prime( C_prime + 1);

      // remainder of r and p
      const size_t M( r % p);

      // sse connection score term R, average percent of restraints in each SSE pair up to ideal value C_prime_prime
      const double S_EC_R( CalculateSSEConnectionScoreComponentR( sse_pair_counts, C_prime_prime, r));

      // sse connection score term C, fraction of SSE pairs that contain exactly the ideal number of restraints
      const double S_EC_C( CalculateSSEConnectionScoreComponentC( sse_pair_counts, C_prime, C_prime_prime, p, M));

      const double composite_score( -( S_EC_R + S_EC_C) / 2.0);

      // return the score
      return composite_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSSEConnection::Read( std::istream &ISTREAM)
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
    std::ostream &DataSetPairwiseSSEConnection::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

    //! @brief Gives information on how many times each SSE pair is connected by a data pair
    //! @param SSES the sses which will be used to count how many datapoints fall in each
    //! @param DATA the dataset which will be used to see how many of each data point are in each sse
    //! @return map of sse pairs and how many times they are connected
    storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>
    DataSetPairwiseSSEConnection::GetSSEPairCounts
    (
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> &SSES,
      const restraint::DataSetPairwise &DATA
    )
    {
      // to keep track of  how many measurements connect the each sse pair
      storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> sse_pairs_counts;

      // iterate through the sses
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( SSES.Begin()), sse_itr_end( SSES.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        const util::SiPtr< const assemble::SSE> &sse_a( *sse_itr_a);

        // iterate through the sses again to make sse pairs - pairs of the same sse are not considered
        for
        (
          storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr_b
            (
              ++storage::Set
              <
                util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap
              >::const_iterator( sse_itr_a)
            );
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          const util::SiPtr< const assemble::SSE> &sse_b( *sse_itr_b);

          // make the pair of sses
          const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > sse_pair( sse_a, sse_b);

          // insert the current pair into sse_pairs_counts
          std::pair
          <
            storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>::iterator, bool
          > insert_status
          (
            sse_pairs_counts.Insert
            (
              std::pair< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>( sse_pair, 0)
            )
          );

          BCL_Assert
          (
            insert_status.second, "could not insert " + sse_a->GetIdentification() + " and "
            + sse_b->GetIdentification()
          );

          // iterate through the data set to see which sse each point of each data pair is in
          for
          (
            restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
            data_itr != data_itr_end;
            ++data_itr
          )
          {
            // true if the current data pair connect the sse pair
            // they should be sorted by sequence so can just compare first with first and second with second
            if( data_itr->First()->IsWithin( *sse_a) && data_itr->Second()->IsWithin( *sse_b))
            {
              // if so add a count for the current SSE pair
              ++insert_status.first->second;
            }

            // make sure sorting is working as expected
            BCL_Assert
            (
              !( data_itr->First()->IsWithin( *sse_b) && data_itr->Second()->IsWithin( *sse_a)), "sorting is broken"
            );
          } //< iterate through dataset
        } //< iterate through sses
      } //< iterate through sses

      // return the map of sse pairs and how many times they are connected
      return sse_pairs_counts;
    }

    //! @brief calculates the R component of the SSE connection score
    //! @param SSE_PAIRS_COUNTS map of sse pairs and how many times they are connected
    //! @param C_PRIME_PRIME maximum acceptable integer value for C, the ideal number of connections for each SSE pair
    //! @param NUMBER_RESTRAINTS this is "r" in the reference
    //! @return double which is the R component of the SSE connection score
    double DataSetPairwiseSSEConnection::CalculateSSEConnectionScoreComponentR
    (
      const storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> &SSE_PAIRS_COUNTS,
      const size_t C_PRIME_PRIME,
      const size_t NUMBER_RESTRAINTS
    )
    {
      double score( 0);

      // iterate through the map
      for
      (
        storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>::const_iterator
          sse_itr( SSE_PAIRS_COUNTS.Begin()), sse_itr_end( SSE_PAIRS_COUNTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // number of connections between the current sse pair
        const size_t ri( sse_itr->second);
        const size_t min( std::min( ri, C_PRIME_PRIME));
        score += min;
      }

      // divide by number of points
      score /= double( NUMBER_RESTRAINTS);

      return score;
    }

    //! @brief sse connection score term C, fraction of SSE pairs that contain exactly the ideal number of restraints
    //! @param SSE_PAIRS_COUNTS map of sse pairs and how many times they are connected
    //! @param C_PRIME minimum acceptable integer value for C, the ideal number of connections for each SSE pair
    //! @param C_PRIME_PRIME maximum acceptable integer value for C, the ideal number of connections for each SSE pair
    //! @param NUM_SSE_PAIRS this is "p" in the reference
    //! @param REMAINDER this is "M" in the reference
    //! @return score for sse connection score component C
    double DataSetPairwiseSSEConnection::CalculateSSEConnectionScoreComponentC
    (
      const storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> &SSE_PAIRS_COUNTS,
      const size_t C_PRIME,
      const size_t C_PRIME_PRIME,
      const size_t NUM_SSE_PAIRS,
      const size_t REMAINDER
    )
    {
      // the number of SSE pairs with C_PRIME restraints
      size_t F_prime( 0);

      // the number of SSE pairs with C_PRIME_PRIME restraints
      size_t F_prime_prime( 0);

      // iterate through the sses and their data point counts to determine F_prime and F_prime_prime
      for
      (
        storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>::const_iterator
          sse_itr( SSE_PAIRS_COUNTS.Begin()), sse_itr_end( SSE_PAIRS_COUNTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const size_t current_count( sse_itr->second);

        if( current_count == C_PRIME)
        {
          ++F_prime;
        }

        else if( current_count == C_PRIME_PRIME)
        {
          ++F_prime_prime;
        }
      }

      // min( Fprime,( p - M))
      size_t min_F_prime( std::min( F_prime, ( NUM_SSE_PAIRS - REMAINDER)));

      // min( F_prime_prime, M)
      size_t min_F_prime_prime( std::min( F_prime_prime, REMAINDER));

      // C component of SSE connection score
      const double component_C( double( min_F_prime + min_F_prime_prime) / double( NUM_SSE_PAIRS));

      return component_C;
    }

  } // namespace score
  
} // namespace bcl
