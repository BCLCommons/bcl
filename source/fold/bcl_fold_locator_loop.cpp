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
#include "fold/bcl_fold_locator_loop.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> LocatorLoop::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorLoop())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param IGNORE_TERMS ignore terminal loops
    //! @param LOCATE_UNDEFINED locate undefined loop regions
    LocatorLoop::LocatorLoop( bool IGNORE_TERMS, bool LOCATE_UNDEFINED) :
      m_IgnoreTerms( IGNORE_TERMS),
      m_LocateUndefined( LOCATE_UNDEFINED)
    {
    }

    //! @brief clone function
    //! @return pointer to a new LocatorLoop
    LocatorLoop *LocatorLoop::Clone() const
    {
      return new LocatorLoop( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LocatorLoop::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates loops in a protein model and computes a parameterization of the segment
    //! @param MODEL protein model for which to find loops
    //! @return parameterizations of loops in the given protein model
    util::ShPtrVector< LoopParameters> LocatorLoop::Locate
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // get all loops in the model and filter out the ones not consistent with given options
      const util::SiPtrVector< const assemble::SSE> loops( MODEL.GetSSEs( biol::GetSSTypes().COIL));
      util::SiPtrVector< const assemble::SSE> loops_filtered;
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        const assemble::SSE &sse( **loop_it);
        const int num_residues( MODEL.GetChain( sse.GetChainID())->GetSequence()->GetSize());
        if( !m_IgnoreTerms || ( sse.GetFirstAA()->GetSeqID() != 1 && sse.GetLastAA()->GetSeqID() != num_residues))
        {
          if( IsDefined( sse) != m_LocateUndefined)
          {
            loops_filtered.PushBack( *loop_it);
          }
        }
      }

      // compute the parameterizations for the filtered loops
      util::ShPtrVector< LoopParameters> segments;
      for
      (
        auto loop_it( loops_filtered.Begin()), loop_it_end( loops_filtered.End());
        loop_it != loop_it_end;
        ++loop_it
      )
      {
        // get the anchor residues for this loop
        const assemble::SSE &sse( **loop_it);
        const char chain_id( sse.GetChainID());
        const util::SiPtrVector< const biol::AABase> &residues
        (
          MODEL.GetChain( chain_id)->GetSequence()->GetMembers()
        );
        const int num_residues( residues.GetSize());
        const int anchor_1_seq_id( std::max( 1, ( **loop_it).GetFirstAA()->GetSeqID() - 1));
        const int anchor_2_seq_id( std::min( num_residues, ( **loop_it).GetLastAA()->GetSeqID() + 1));
        const biol::AABase &anchor_1( *residues( anchor_1_seq_id - 1));
        const biol::AABase &anchor_2( *residues( anchor_2_seq_id - 1));

        // compute the parameterization of the loop
        if( anchor_2.GetSeqID() - anchor_1.GetSeqID() > 0)
        {
          segments.PushBack( LoopParameters::Create( anchor_1, anchor_2));
        }
      }

      return segments;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &LocatorLoop::Read( std::istream &ISTREAM)
    {
      // read members from input stream
      io::Serialize::Read( m_IgnoreTerms, ISTREAM);
      io::Serialize::Read( m_LocateUndefined, ISTREAM);

      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &LocatorLoop::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members into output stream
      io::Serialize::Write( m_IgnoreTerms, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_LocateUndefined, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines if the given loop has fully defined coordinates
    //! @param LOOP loop for which to determine if it is fully defined
    //! @return true if the given loop has fully defined coordinates
    bool LocatorLoop::IsDefined( const assemble::SSE &LOOP)
    {
      const util::ShPtrVector< biol::AABase> &residues( LOOP.GetData());
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &current_res( **res_it);
        if( !current_res.HasDefinedCoordinates())
        {
          return false;
        }
      }

      return true;
    }

  } // namespace fold
} // namespace bcl
