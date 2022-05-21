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
#include "fold/bcl_fold_loop_segment.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopSegment::s_Instance
    (
      GetObjectInstances().AddInstance( new LoopSegment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopSegment::LoopSegment() :
      m_SSE(),
      m_IsRigid()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SSE the sse which the loop segment refers to
    //! @param IS_RIGID boolean indicating if the dihedral angles of the LoopSegment can be changed or not
    LoopSegment::LoopSegment( const util::ShPtr< assemble::SSE> &SSE, const bool IS_RIGID) :
      m_SSE( SSE),
      m_IsRigid( IS_RIGID)
    {
    }

    //! @brief copy constructor
    //! @param OTHER loop segment to be copied
    LoopSegment::LoopSegment( const LoopSegment &OTHER) :
      m_SSE(),
      m_IsRigid( OTHER.m_IsRigid)
    {
      if( OTHER.m_SSE.IsDefined())
      {
        m_SSE = util::ShPtr< assemble::SSE>( OTHER.m_SSE->Clone());
      }
    }

    //! @brief Clone function
    //! @return pointer to new LoopSegment
    LoopSegment *LoopSegment::Clone() const
    {
      return new LoopSegment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LoopSegment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetSSE gives the sse which the loop segment refers to
    //! @return "m_SSE" the sse which the loop segment refers to
    const util::ShPtr< assemble::SSE> &LoopSegment::GetSSE() const
    {
      return m_SSE;
    }

    //! @brief GetSSE gives the sse which the loop segment refers to
    //! @return "m_SSE" the sse which the loop segment refers to
    util::ShPtr< assemble::SSE> &LoopSegment::GetSSE()
    {
      return m_SSE;
    }

    //! @brief GetSSE gives the sse which the loop segment refers to
    //! @return "m_SSE" the sse which the loop segment refers to
    assemble::SSE &LoopSegment::GetSSEReference()
    {
      return *m_SSE;
    }

    //! @brief GetSSE gives the sse which the loop segment refers to
    //! @return "m_SSE" the sse which the loop segment refers to
    const assemble::SSE &LoopSegment::GetConstSSEReference() const
    {
      return *m_SSE;
    }

    //! @brief IsRigid gives boolean indicating if the dihedral angles of the LoopSegment can be changed or not
    //! @return "m_IsRigid" the boolean indicating if the dihedral angles of the LoopSegment can be changed or not
    bool LoopSegment::IsRigid() const
    {
      return m_IsRigid;
    }

    //! @brief gives iterator to residue provided
    //! @param AA_BASE residue for which an iterator will be given
    //! @return iterator to the provided residue
    biol::AASequence::iterator LoopSegment::GetAAIterator( const biol::AABase &AA_BASE)
    {
      // iterate over the sse
      for
      (
        biol::AASequence::iterator itr( m_SSE->Begin()), itr_end( m_SSE->End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetChainID() == AA_BASE.GetChainID() && ( *itr)->GetSeqID() == AA_BASE.GetSeqID())
        {
          return itr;
        }
      }

      return m_SSE->End();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LoopSegment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSE, ISTREAM);
      io::Serialize::Read( m_IsRigid, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LoopSegment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSE, OSTREAM, INDENT);
      io::Serialize::Write( m_IsRigid, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
