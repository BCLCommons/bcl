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
#include "score/bcl_score_sse_pool_sses.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolSSEs::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new SSEPoolSSEs
        (
          util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >(),
          false, false
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from sse score and normalization
    SSEPoolSSEs::SSEPoolSSEs
    (
      const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_SCORE_SSE,
      const bool NORMALIZE_SSE,
      const bool NORMALIZE_NUMBER_SSES
    ) :
      m_SSEScore( SP_SCORE_SSE),
      m_NormalizeSSE( NORMALIZE_SSE),
      m_NormalizeNumberSSEs( NORMALIZE_NUMBER_SSES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolSSEs
    SSEPoolSSEs *SSEPoolSSEs::Clone() const
    {
      return new SSEPoolSSEs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolSSEs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPoolSSEs::GetScheme() const
    {
      return m_SSEScore->GetScheme();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the pool
    //! @param POOL the sse pool to be scored
    //! @param MEMBRANE membrane object
    //! @return the sum of the normalized scores
    double SSEPoolSSEs::operator()( const assemble::SSEPool &POOL, const biol::Membrane &MEMBRANE) const
    {
      double score( 0);
      size_t number( 0);
      size_t entities( 0);

      // iterate over all sses in the pool
      for( assemble::SSEPool::const_iterator itr( POOL.Begin()), itr_end( POOL.End()); itr != itr_end; ++itr)
      {
        // score that sses
        storage::Pair< double, size_t> result( m_SSEScore->operator ()( **itr, MEMBRANE));

        // at least one entity needs to be scored
        if( result.Second() == 0)
        {
          continue;
        }

        if( !util::IsDefined( result.First()))
        {
          BCL_MessageDbg
          (
            "scoring SSE: " + ( *itr)->GetIdentification() +
            " with score: " + m_SSEScore->GetScheme() + " returned undefined"
          );
          continue;
        }
//        // normalize by the entities scored
//        if( m_NormalizeSSE)
//        {
//          result.First() /= double( result.Second());
//        }

        score    += result.First();
        entities += result.Second();
        ++number;
      }

      if( m_NormalizeSSE && entities > 0)
      {
        score /= double( entities);
      }
      // normalize by the number of sses
      else if( m_NormalizeNumberSSEs && number > 0)
      {
        score /= double( number);
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param POOL Argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @param MEMBRANE membrane object
    //! @return std::ostream which was written to
    std::ostream &SSEPoolSSEs::WriteDetailedSchemeAndValues
    (
      const assemble::SSEPool &POOL,
      const biol::Membrane &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      // iterate over all sses in the pool
      for( assemble::SSEPool::const_iterator itr( POOL.Begin()), itr_end( POOL.End()); itr != itr_end; ++itr)
      {
        // score that sses
        storage::Pair< double, size_t> result( m_SSEScore->operator ()( **itr, MEMBRANE));

        // normalize by the entities scored
        if( m_NormalizeSSE)
        {
          result.First() /= double( result.Second());
        }

        OSTREAM << ( *itr)->GetIdentification() << '\t' << result.First() << '\t' << result.Second() << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolSSEs::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEScore           , ISTREAM);
      io::Serialize::Read( m_NormalizeSSE       , ISTREAM);
      io::Serialize::Read( m_NormalizeNumberSSEs, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolSSEs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEScore           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NormalizeSSE       , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NormalizeNumberSSEs, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace score
} // namespace bcl
