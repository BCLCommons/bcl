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
#include "assemble/bcl_assemble_sse_pool_split_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! @brief the default scheme for this class
    const std::string &SSEPoolSplitSSE::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "pool_split_sse");
      return s_default_scheme;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolSplitSSE::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new SSEPoolSplitSSE
        (
          sspred::GetMethods().e_Undefined,
          util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >(),
          util::ShPtr< math::MutateInterface< SSE> >()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from locator and mutate
    //! @param SS_METHOD the method to use to locate the largest drop in the ss prediction for the located sse
    //! @param SP_LOCATE_SSE locator of sse from domain
    //! @param SP_MUTATE_SSE mutate to apply to sse that has the lower of the predictions
    //! @param SCHEME the scheme
    SSEPoolSplitSSE::SSEPoolSplitSSE
    (
      const sspred::Method &SS_METHOD,
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
      const util::ShPtr< math::MutateInterface< SSE> > &SP_MUTATE_SSE,
      const std::string &SCHEME
    ) :
      m_Method( SS_METHOD),
      m_SSELocator( SP_LOCATE_SSE),
      m_SSEMutate( SP_MUTATE_SSE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolSplitSSE
    SSEPoolSplitSSE *SSEPoolSplitSSE::Clone() const
    {
      return new SSEPoolSplitSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolSplitSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &SSEPoolSplitSSE::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a pool by splitting a single sse
    //! @param SSE_POOL SSEPool to be mutated
    //! @return MutateResult that results from mutating to the SSE_POOL
    math::MutateResult< SSEPool> SSEPoolSplitSSE::operator()( const SSEPool &SSE_POOL) const
    {
      // locate a sse from the pool
      const util::SiPtr< const SSE> located_sse( m_SSELocator->Locate( SSE_POOL));

      // was locating successful
      if( !located_sse.IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // length of sse to split
      const size_t seq_length( located_sse->GetSize());

      // need at least two residues to split
      if( seq_length < 2)
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // find drop in prediction
      const storage::Pair< util::SiPtr< const biol::AABase>, double> aa_pos_left_diff( FindLargestSSPredDifferencePosition( *located_sse, m_Method));

      // no drop could be identified
      if( !aa_pos_left_diff.First().IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // either left or right sse to mutate (the one with the lower prediction)
      const bool left( aa_pos_left_diff.Second() > 0.0);
      const size_t split_pos( aa_pos_left_diff.First()->GetSeqID() - located_sse->GetFirstAA()->GetSeqID());

      util::ShPtr< SSE> sp_sse_a( new SSE( located_sse->SubSequence( 0, split_pos), located_sse->GetType()));
      util::ShPtr< SSE> sp_sse_b( new SSE( located_sse->SubSequence( split_pos, seq_length - split_pos), located_sse->GetType()));

      // mutate either sse
      if( m_SSEMutate.IsDefined())
      {
        if( left)
        {
          math::MutateResult< SSE> result( m_SSEMutate->operator ()( *sp_sse_a));
          if( result.GetArgument().IsDefined())
          {
            sp_sse_a = result.GetArgument();
          }
        }
        else
        {
          math::MutateResult< SSE> result( m_SSEMutate->operator ()( *sp_sse_b));
          if( result.GetArgument().IsDefined())
          {
            sp_sse_b = result.GetArgument();
          }
        }
      }

      BCL_MessageDbg
      (
        "split sse: " + located_sse->GetIdentification() + "\nat: " +
        aa_pos_left_diff.First()->GetIdentification() +
        "\ninto\n" +
        sp_sse_a->GetIdentification() + "\n" + sp_sse_b->GetIdentification()
      );

      // replace with new sses
      util::SiPtrList< const SSE> pool_sses( SSE_POOL.Begin(), SSE_POOL.End());
      util::SiPtrList< const SSE>::iterator new_end
      (
        std::remove_if
        (
          pool_sses.Begin(), pool_sses.End(),
          SSECompare( *located_sse)
        )
      );

      util::ShPtr< SSEPool> sp_pool( new SSEPool( util::SiPtrList< const SSE>( pool_sses.Begin(), new_end), false));
      sp_pool->Insert( sp_sse_a);
      sp_pool->Insert( sp_sse_b);

      // end
      return math::MutateResult< SSEPool>( sp_pool, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolSplitSSE::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method    , ISTREAM);
      io::Serialize::Read( m_SSELocator, ISTREAM);
      io::Serialize::Read( m_SSEMutate , ISTREAM);
      io::Serialize::Read( m_Scheme    , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolSplitSSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSELocator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSEMutate , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme    , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find the amino acid right to the highest difference for the given sspred::Method and SSType and report the difference
    //! @param ELEMENT the secondary structure element
    //! @param SSMETHOD the ssprediction method
    //! @return Pair of difference (right-left) and the amino acid to the left
    storage::Pair< util::SiPtr< const biol::AABase>, double> SSEPoolSplitSSE::FindLargestSSPredDifferencePosition
    (
      const SSE &ELEMENT, const sspred::Method &SSMETHOD
    )
    {
      // initialize the pair
      storage::Pair< util::SiPtr< const biol::AABase>, double> aa_diff( util::SiPtr< const biol::AABase>(), 0.0);

      // sequence too short
      if( ELEMENT.GetSize() < 2)
      {
        return aa_diff;
      }

      // sstype
      const biol::SSType sstype( ELEMENT.GetType());

      // iterator for amino acids
      SSE::const_iterator itr( ELEMENT.Begin()), itr_end( ELEMENT.End());

      // first aa
      aa_diff.First() = *itr;
      double left_prediction( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));

      // next aa
      ++itr;

      // iterate over all aas
      while( itr != itr_end)
      {
        const double prediction( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));
        const double difference( prediction - left_prediction);

        // update the highest difference if necessary
        if( math::Absolute( difference) > math::Absolute( aa_diff.Second()))
        {
          aa_diff.First() = *itr;
          aa_diff.Second() = difference;
        }

        // go to next
        left_prediction = prediction;
        ++itr;
      }

      // if largest difference is 0, do not return a position
      if( aa_diff.Second() == 0.0)
      {
        aa_diff.First() = util::SiPtr< const biol::AABase>();
      }

      // end
      return aa_diff;
    }

  } // namespace assemble
} // namespace bcl
