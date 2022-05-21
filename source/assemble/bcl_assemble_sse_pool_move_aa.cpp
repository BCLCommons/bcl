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
#include "assemble/bcl_assemble_sse_pool_move_aa.h"

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
    const std::string &SSEPoolMoveAA::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "pool_move_aa");
      return s_default_scheme;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolMoveAA::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPoolMoveAA( math::Range< size_t>( 1, 1)))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from max number of residues to move
    //! @param RESIDUES_TO_MOVE_RANGE range of number of residues to move
    //! @param SCHEME the scheme
    SSEPoolMoveAA::SSEPoolMoveAA
    (
      const math::Range< size_t> &RESIDUES_TO_MOVE_RANGE,
      const std::string &SCHEME
    ) :
      m_ResdiuesToMoveRange( RESIDUES_TO_MOVE_RANGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolMoveAA
    SSEPoolMoveAA *SSEPoolMoveAA::Clone() const
    {
      return new SSEPoolMoveAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolMoveAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &SSEPoolMoveAA::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a pool by shifting aas between sses
    //! @param SSE_POOL
    //! @return MutateResult that results from mutating to the SSE_POOL
    math::MutateResult< SSEPool> SSEPoolMoveAA::operator()( const SSEPool &SSE_POOL) const
    {
      // all sses from the pool
      util::SiPtrList< const SSE> pool_sses( SSE_POOL.Begin(), SSE_POOL.End());

      if( pool_sses.GetSize() < 2)
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // pick an sse randomly
      const size_t index( random::GetGlobalRandom().Random( pool_sses.GetSize() - 2));

      // identify two adjacent sses
      util::SiPtrList< const SSE>::const_iterator itr1( pool_sses.Begin());
      util::SiPtrList< const SSE>::const_iterator itr2( pool_sses.Begin());
      std::advance( itr1, index);
      std::advance( itr2, index + 1);

      const SSE &sse1( **itr1);
      const SSE &sse2( **itr2);

      // check that the sses are consecutive
      if( ( ( *itr2)->GetFirstAA()->GetSeqID() - ( *itr1)->GetLastAA()->GetSeqID()) != 1)
      {
        BCL_MessageDbg
        (
          "this mutate only works for pools of consecutive SSEs, but encountered:\n" +
          ( *itr1)->GetIdentification() + '\n' + ( *itr2)->GetIdentification()
        );
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // create new pool
      util::SiPtrList< const SSE>::iterator new_end
      (
        std::remove_if
        (
          pool_sses.Begin(), pool_sses.End(),
          SSECompare( sse1)
        )
      );

      new_end = std::remove_if
      (
        pool_sses.Begin(), new_end,
        SSECompare( sse2)
      );

      util::ShPtr< SSEPool> sp_pool( new SSEPool( util::SiPtrList< const SSE>( pool_sses.Begin(), new_end), false));

      // create two new sses
      size_t nr_res( random::GetGlobalRandom().Random( m_ResdiuesToMoveRange));
      const bool move_up( random::GetGlobalRandom().Boolean());

      BCL_MessageDbg( "sses before moving aa:\n" + sse1.GetIdentification() + '\n' + sse2.GetIdentification());

      // move residues up in sequence to next SSE
      if( move_up)
      {
        BCL_MessageDbg( "sses after moving aa:");
        if( nr_res < sse1.GetSize())
        {
          util::ShPtr< SSE> sp_sse( new SSE( sse1.SubSequence( 0, sse1.GetSize() - nr_res), sse1.GetType()));
          sp_pool->Insert( sp_sse);
          BCL_MessageDbg( sp_sse->GetIdentification());
        }

        // reduce the number of residues if necessary
        nr_res = std::min( nr_res, sse1.GetSize());

        // sequence for other sse
        biol::AASequence seq2( sse2);
        seq2.PrependSequence( sse1.SubSequence( sse1.GetSize() - nr_res, nr_res));

        // new sse
        util::ShPtr< SSE> sp_sse( new SSE( seq2, sse2.GetType()));
        sp_pool->Insert( sp_sse);
        BCL_MessageDbg( sp_sse->GetIdentification());
      }
      // move residues down in sequence to previous sse
      else
      {
        BCL_MessageDbg( "sses after moving aa:");
        if( nr_res < sse2.GetSize())
        {
          util::ShPtr< SSE> sp_sse( new SSE( sse2.SubSequence( nr_res, sse2.GetSize() - nr_res), sse2.GetType()));
          sp_pool->Insert( sp_sse);
          BCL_MessageDbg( sp_sse->GetIdentification());
        }

        // reduce the number of residues if necessary
        nr_res = std::min( nr_res, sse2.GetSize());
        biol::AASequence seq1( sse1);
        seq1.AppendSequence( sse2.SubSequence( 0, nr_res));

        // new sse
        util::ShPtr< SSE> sp_sse( new SSE( seq1, sse1.GetType()));
        sp_pool->Insert( sp_sse);
        BCL_MessageDbg( sp_sse->GetIdentification());
      }

      // end
      return math::MutateResult< SSEPool>( sp_pool, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolMoveAA::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ResdiuesToMoveRange, ISTREAM);
      io::Serialize::Read( m_Scheme             , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolMoveAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ResdiuesToMoveRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme             , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
