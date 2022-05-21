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
#include "assemble/bcl_assemble_sse_compare.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if SSE_A comes before SSE_B in sequence
    //! @param SSE_A first SSE
    //! @param SSE_B second SSE
    //! @return true if SSE_A comes before SSE_B in sequence
    bool SSELessThan::operator()( const SSE &SSE_A, const SSE &SSE_B) const
    {
      // if in a previous chain
      if( SSE_A.GetChainID() < SSE_B.GetChainID())
      {
        return true;
      }

      // if in the same chain
      if( SSE_A.GetChainID() == SSE_B.GetChainID())
      {
        // check that there are amino acids available - if they are nt, compare size
        if( SSE_A.GetData().IsEmpty() || SSE_B.GetData().IsEmpty())
        {
          return SSE_A.GetSize() < SSE_B.GetSize();
        }

        // if sse_a starts first return true
        if( SSE_A.GetFirstAA()->GetSeqID() < SSE_B.GetFirstAA()->GetSeqID())
        {
          return true;
        }
        // if sse_a starts at the same position
        if( SSE_A.GetFirstAA()->GetSeqID() == SSE_B.GetFirstAA()->GetSeqID())
        {
          // if sse_a ends earlier
          if( SSE_A.GetLastAA()->GetSeqID() < SSE_B.GetLastAA()->GetSeqID())
          {
            return true;
          }
          // if they end at the same seq id
          if( SSE_A.GetLastAA()->GetSeqID() == SSE_B.GetLastAA()->GetSeqID())
          {
            // check if ss_type is smaller
            if( SSE_A.GetType() < SSE_B.GetType())
            {
              return true;
            }
          }
        }
      }

      // for all other cases return false
      return false;
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
    bool SSELessThan::operator()
    (
      const util::PtrInterface< SSE> &PTR_SSE_A,
      const util::PtrInterface< SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
    bool SSELessThan::operator()
    (
      const util::PtrInterface< const SSE> &PTR_SSE_A,
      const util::PtrInterface< const SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
    bool SSELessThan::operator()
    (
      const util::PtrInterface< const SSE> &PTR_SSE_A,
      const util::PtrInterface< SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
    bool SSELessThan::operator()
    (
      const util::PtrInterface< SSE> &PTR_SSE_A,
      const util::PtrInterface< const SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if SSE_A comes before SSE_B in sequence, false if overlap
    //! @param SSE_A first SSE
    //! @param SSE_B second SSE
    //! @return true if SSE_A comes before SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()( const SSE &SSE_A, const SSE &SSE_B) const
    {
      if( SSE_A.GetChainID() < SSE_B.GetChainID())
      {
        return true;
      }
      else if( SSE_A.GetChainID() > SSE_B.GetChainID())
      {
        return false;
      }

      // an empty sse is smaller
      if( SSE_A.GetData().IsEmpty())
      {
        return true;
      }

      // an empty is smaller
      if( SSE_B.GetData().IsEmpty())
      {
        return false;
      }

      if
      (
        ( SSE_A.GetFirstAA()->GetSeqID() < SSE_B.GetFirstAA()->GetSeqID())
        &&
        ( SSE_A.GetLastAA()->GetSeqID() < SSE_B.GetFirstAA()->GetSeqID()))
      {
        return true;
      }

      // else return larger
      return false;
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()
    (
      const util::PtrInterface< SSE> &PTR_SSE_A,
      const util::PtrInterface< SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()
    (
      const util::PtrInterface< const SSE> &PTR_SSE_A,
      const util::PtrInterface< const SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
    //! @param SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if SSE_A comes before PTR_SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()( const SSE &SSE_A, const util::PtrInterface< SSE> &PTR_SSE_B) const
    {
      return operator()( SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
    //! @param SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if SSE_A comes before PTR_SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()( const SSE &SSE_A, const util::PtrInterface< const SSE> &PTR_SSE_B) const
    {
      return operator()( SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
    //! @param PTR_SSE_A first SSE
    //! @param SSE_B second SSE
    //! @return true if PTR_SSE_A comes before SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()( const util::PtrInterface< SSE> &PTR_SSE_A, const SSE &SSE_B) const
    {
      return operator()( *PTR_SSE_A, SSE_B);
    }

    //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
    //! @param PTR_SSE_A first SSE
    //! @param SSE_B second SSE
    //! @return true if PTR_SSE_A comes before SSE_B in sequence, false if overlap
    bool SSELessThanNoOverlap::operator()( const util::PtrInterface< const SSE> &PTR_SSE_A, const SSE &SSE_B) const
    {
      return operator()( *PTR_SSE_A, SSE_B);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if SSE_A is smaller than SSE_B in size
    //! @param SSE_A first SSE
    //! @param SSE_B second SSE
    //! @return true if SSE_A is smaller than SSE_B in size
    bool SSELessThanBySize::operator()( const SSE &SSE_A, const SSE &SSE_B) const
    {
      // if has a smaller size
      if( SSE_A.GetSize() < SSE_B.GetSize())
      {
        return true;
      }
      // if the sizes of both SSEs are same
      if( SSE_A.GetSize() == SSE_B.GetSize())
      {
        // if in a previous chain
        if( SSE_A.GetChainID() < SSE_B.GetChainID())
        {
          return true;
        }
        // if in the same chain
        if( SSE_A.GetChainID() == SSE_B.GetChainID())
        {
          // if sse_a starts first return true
          if( SSE_A.GetFirstAA()->GetSeqID() < SSE_B.GetFirstAA()->GetSeqID())
          {
            return true;
          }
          // if sse_a starts at the same position
          if( SSE_A.GetFirstAA()->GetSeqID() == SSE_B.GetFirstAA()->GetSeqID())
          {
            // if sse_a ends earlier
            if( SSE_A.GetLastAA()->GetSeqID() < SSE_B.GetLastAA()->GetSeqID())
            {
              return true;
            }
            // if they end at the same seq id
            if( SSE_A.GetLastAA()->GetSeqID() == SSE_B.GetLastAA()->GetSeqID())
            {
              // check if ss_type is smaller
              if( SSE_A.GetType() < SSE_B.GetType())
              {
                return true;
              }
            }
          }
        }
      }

      // for all other cases return false
      return false;
    }

    //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
    bool SSELessThanBySize::operator()
    (
      const util::PtrInterface< SSE> &PTR_SSE_A,
      const util::PtrInterface< SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
    bool SSELessThanBySize::operator()
    (
      const util::PtrInterface< const SSE> &PTR_SSE_A,
      const util::PtrInterface< const SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
    bool SSELessThanBySize::operator()
    (
      const util::PtrInterface< const SSE> &PTR_SSE_A,
      const util::PtrInterface< SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

    //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
    //! @param PTR_SSE_A first SSE
    //! @param PTR_SSE_B second SSE
    //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
    bool SSELessThanBySize::operator()
    (
      const util::PtrInterface< SSE> &PTR_SSE_A,
      const util::PtrInterface< const SSE> &PTR_SSE_B
    ) const
    {
      return operator()( *PTR_SSE_A, *PTR_SSE_B);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a SSE reference
    //! @param THIS_SSE reference to SSE to be compared against
    SSECompare::SSECompare( const SSE &THIS_SSE) :
      m_SSE( THIS_SSE)
    {
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compares the provided THIS_SSE with m_SSE to check it is the same sse
    //! @param THIS_SSE SSE to be compared
    //! @return true if THIS_SSE is equal to m_SSE
    bool SSECompare::operator()( const SSE &THIS_SSE) const
    {
      return THIS_SSE == m_SSE;
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
    //! @param PTR_SSE SSE to be compared
    //! @return true if PTR_SSE is equal to m_SSE
    bool SSECompare::operator()( const util::PtrInterface< SSE> &PTR_SSE) const
    {
      return *PTR_SSE == m_SSE;
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
    //! @param PTR_SSE SSE to be compared
    //! @return true if PTR_SSE is equal to m_SSE
    bool SSECompare::operator()( const util::PtrInterface< const SSE> &PTR_SSE) const
    {
      return *PTR_SSE == m_SSE;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a begin seq id, end seq id, chain id and SSType
    //! @param BEGIN_SEQ_ID seq id of the first amino acid  of the SSE searched
    //! @param END_SEQ_ID seq id of the last amino acid of the SSE searched
    //! @param CHAIN_ID chain id of the SSE searched
    //! @param SS_TYPE type of the SSE searched
    SSECompareByIdentity::SSECompareByIdentity
    (
      const size_t        BEGIN_SEQ_ID,
      const size_t        END_SEQ_ID,
      const char          CHAIN_ID,
      const biol::SSType &SS_TYPE
    ) :
      m_BeginSeqID( BEGIN_SEQ_ID),
      m_EndSeqID( END_SEQ_ID),
      m_ChainID( CHAIN_ID),
      m_Type( SS_TYPE)
    {
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compares the provided THIS_SSE with m_SSE to check it is the same sse
    //! @param THIS_SSE SSE to be compared
    //! @return true if THIS_SSE is equal to m_SSE
    bool SSECompareByIdentity::operator()( const SSE &THIS_SSE) const
    {
      return
        size_t( THIS_SSE.GetFirstAA()->GetSeqID()) == m_BeginSeqID &&
        size_t( THIS_SSE.GetLastAA()->GetSeqID()) == m_EndSeqID &&
        THIS_SSE.GetChainID() == m_ChainID &&
        THIS_SSE.GetType() == m_Type;
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
    //! @param PTR_SSE SSE to be compared
    //! @return true if PTR_SSE is equal to m_SSE
    bool SSECompareByIdentity::operator()( const util::PtrInterface< SSE> &PTR_SSE) const
    {
      return operator()( *PTR_SSE);
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
    //! @param PTR_SSE SSE to be compared
    //! @return true if PTR_SSE is equal to m_SSE
    bool SSECompareByIdentity::operator()( const util::PtrInterface< const SSE> &PTR_SSE) const
    {
      return operator()( *PTR_SSE);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a SSE reference
    //! @param THIS_SSE reference to SSE to be compared against
    SSECompareOverlap::SSECompareOverlap( const SSE &THIS_SSE) :
      m_SSE( THIS_SSE)
    {
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compares the provided THIS_SSE with m_SSE to check for overlap
    //! @param THIS_SSE SSE to be checked for overlap
    //! @return true if THIS_SSE overlaps with m_SSE
    bool SSECompareOverlap::operator()( const SSE &THIS_SSE) const
    {
      // return overlap
      return biol::DoOverlap( THIS_SSE, m_SSE);
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check for overlap
    //! @param PTR_SSE SSE to be checked for overlap
    //! @return true if PTR_SSE overlaps with m_SSE
    bool SSECompareOverlap::operator()( const util::PtrInterface< SSE> &PTR_SSE) const
    {
      // return overlap
      return biol::DoOverlap( *PTR_SSE, m_SSE);
    }

    //! @brief compares the provided PTR_SSE with m_SSE to check for overlap
    //! @param PTR_SSE SSE to be checked for overlap
    //! @return true if PTR_SSE overlaps with m_SSE
    bool SSECompareOverlap::operator()( const util::PtrInterface< const SSE> &PTR_SSE) const
    {
      // return overlap
      return biol::DoOverlap( *PTR_SSE, m_SSE);
    }

  } // namespace assemble
} // namespace bcl
