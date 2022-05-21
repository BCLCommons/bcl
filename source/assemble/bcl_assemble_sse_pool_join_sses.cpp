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
#include "assemble/bcl_assemble_sse_pool_join_sses.h"

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
    const std::string &SSEPoolJoinSSEs::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "pool_join_sses");
      return s_default_scheme;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolJoinSSEs::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPoolJoinSSEs( util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >(), false, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from scheme
    //! @param SP_LOCATE_SSE locator of sse from domain
    //! @param JOIN_LEFT join the sse left to the located; default: true
    //! @param JOIN_RIGHT join the sse right to the located; default: true
    //! @param SCHEME the scheme
    SSEPoolJoinSSEs::SSEPoolJoinSSEs
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
      const bool JOIN_LEFT,
      const bool JOIN_RIGHT,
      const std::string &SCHEME
    ) :
      m_SSELocator( SP_LOCATE_SSE),
      m_JoinLeft(   JOIN_LEFT),
      m_JoinRight(  JOIN_RIGHT),
      m_Scheme(     SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolJoinSSEs
    SSEPoolJoinSSEs *SSEPoolJoinSSEs::Clone() const
    {
      return new SSEPoolJoinSSEs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolJoinSSEs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &SSEPoolJoinSSEs::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a pool by joining multiple sses
    //! @param SSE_POOL
    //! @return MutateResult that results from mutating to the SSE_POOL
    math::MutateResult< SSEPool> SSEPoolJoinSSEs::operator()( const SSEPool &SSE_POOL) const
    {
      // locate an sses
      const util::SiPtr< const SSE> located_sse( m_SSELocator->Locate( SSE_POOL));

      // was locating successful
      if( !located_sse.IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // create new sequence
      biol::AASequence new_sequence( *located_sse);

      // neighbor sses
      const storage::VectorND< 2, util::SiPtr< const SSE> > neighbors( SSE_POOL.GetAdjacentSSEs( *located_sse));

      if( !neighbors.First().IsDefined() || !neighbors.Second().IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      bool joined_left( false);
      bool joined_right( false);

      util::SiPtrList< const SSE> pool_sses( SSE_POOL.Begin(), SSE_POOL.End());
      util::SiPtrList< const SSE>::iterator new_end( std::remove_if( pool_sses.Begin(), pool_sses.End(), SSECompare( *located_sse)));

      storage::Set< biol::SSType> ss_types;
      ss_types.Insert( located_sse->GetType());

      // find left adjacent sse
      if( m_JoinLeft && neighbors.First().IsDefined())
      {
        new_sequence.PrependSequence( *neighbors.First());
        joined_left = true;
        new_end = std::remove_if( pool_sses.Begin(), new_end, SSECompare( *neighbors.First()));
        ss_types.Insert( neighbors.First()->GetType());
      }

      // find right adjacent sse
      if( m_JoinRight && neighbors.Second().IsDefined())
      {
        new_sequence.AppendSequence( *neighbors.Second());
        joined_right = true;
        new_end = std::remove_if( pool_sses.Begin(), new_end, SSECompare( *neighbors.Second()));
        ss_types.Insert( neighbors.Second()->GetType());
      }

      // if there was no join, skip
      if( !joined_left && !joined_right)
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // joined sse
      util::ShPtr< SSE> sp_sse( new SSE( new_sequence, *random::GetGlobalRandom().Iterator( ss_types.Begin(), ss_types.End(), ss_types.GetSize())));

      BCL_MessageDbg
      (
        "joined sses: " +
        std::string( joined_left ? "\nleft: " + neighbors.First()->GetIdentification() : "") +
        "\nmiddle: " + located_sse->GetIdentification() +
        std::string( joined_right ? "\nright: " + neighbors.Second()->GetIdentification() : "") +
        "\ninto: " + sp_sse->GetIdentification()
      )

      // new pool
      util::ShPtr< SSEPool> sp_pool( new SSEPool( util::SiPtrList< const SSE>( pool_sses.Begin(), new_end), false));
      sp_pool->Insert( sp_sse);

      // end
      return math::MutateResult< SSEPool>( sp_pool, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolJoinSSEs::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSELocator, ISTREAM);
      io::Serialize::Read( m_JoinLeft  , ISTREAM);
      io::Serialize::Read( m_JoinRight , ISTREAM);
      io::Serialize::Read( m_Scheme    , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolJoinSSEs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSELocator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_JoinLeft  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_JoinRight , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme    , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
