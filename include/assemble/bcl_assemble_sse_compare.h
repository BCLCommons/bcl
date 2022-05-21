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

#ifndef BCL_ASSEMBLE_SSE_COMPARE_H_
#define BCL_ASSEMBLE_SSE_COMPARE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSELessThan
    //! @brief This is a function class for comparing two SSEs in less-than fashion
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date 22.02.2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSELessThan
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_A comes before SSE_B in sequence
      //! @param SSE_A first SSE
      //! @param SSE_B second SSE
      //! @return true if SSE_A comes before SSE_B in sequence
      bool operator()( const SSE &SSE_A, const SSE &SSE_B) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
      bool operator()
      (
        const util::PtrInterface< SSE> &PTR_SSE_A,
        const util::PtrInterface< SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
      bool operator()
      (
        const util::PtrInterface< const SSE> &PTR_SSE_A,
        const util::PtrInterface< const SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
      bool operator()
      (
        const util::PtrInterface< const SSE> &PTR_SSE_A,
        const util::PtrInterface< SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence
      bool operator()
      (
        const util::PtrInterface< SSE> &PTR_SSE_A,
        const util::PtrInterface< const SSE> &PTR_SSE_B
      ) const;

    }; // class SSELessThan

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSELessThanNoOverlap
    //! @brief Function class for comparing two SSEs in less-than fashion, but returns false for overlapping SSEs
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date 18.05.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSELessThanNoOverlap
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_A comes before SSE_B in sequence, false if overlap
      //! @param SSE_A first SSE
      //! @param SSE_B second SSE
      //! @return true if SSE_A comes before SSE_B in sequence, false if overlap
      bool operator()( const SSE &SSE_A, const SSE &SSE_B) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
      bool operator()( const util::PtrInterface< SSE> &PTR_SSE_A, const util::PtrInterface< SSE> &PTR_SSE_B) const;

      //! @brief returns true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A comes before PTR_SSE_B in sequence, false if overlap
      bool operator()( const util::PtrInterface< const SSE> &PTR_SSE_A, const util::PtrInterface< const SSE> &PTR_SSE_B) const;

      //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
      //! @param SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if SSE_A comes before PTR_SSE_B in sequence, false if overlap
      bool operator()( const SSE &SSE_A, const util::PtrInterface< SSE> &PTR_SSE_B) const;

      //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
      //! @param SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if SSE_A comes before PTR_SSE_B in sequence, false if overlap
      bool operator()( const SSE &SSE_A, const util::PtrInterface< const SSE> &PTR_SSE_B) const;

      //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
      //! @param PTR_SSE_A first SSE
      //! @param SSE_B second SSE
      //! @return true if PTR_SSE_A comes before SSE_B in sequence, false if overlap
      bool operator()( const util::PtrInterface< SSE> &PTR_SSE_A, const SSE &SSE_B) const;

      //! @brief returns true if SSE_A comes before SSE_B in sequence, false if overlap
      //! @param PTR_SSE_A first SSE
      //! @param SSE_B second SSE
      //! @return true if PTR_SSE_A comes before SSE_B in sequence, false if overlap
      bool operator()( const util::PtrInterface< const SSE> &PTR_SSE_A, const SSE &SSE_B) const;

    }; // class SSELessThanNoOverlap

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSELessThanBySize
    //! @brief This is a function class for comparing two SSEs in less-than by size fashion
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date 22.02.2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSELessThanBySize
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_A is smaller than SSE_B in size
      //! @param SSE_A first SSE
      //! @param SSE_B second SSE
      //! @return true if SSE_A is smaller than SSE_B in size
      bool operator()( const SSE &SSE_A, const SSE &SSE_B) const;

      //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
      bool operator()
      (
        const util::PtrInterface< SSE> &PTR_SSE_A,
        const util::PtrInterface< SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
      bool operator()
      (
        const util::PtrInterface< const SSE> &PTR_SSE_A,
        const util::PtrInterface< const SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
      bool operator()
      (
        const util::PtrInterface< const SSE> &PTR_SSE_A,
        const util::PtrInterface< SSE> &PTR_SSE_B
      ) const;

      //! @brief returns true if PTR_SSE_A is smaller than PTR_SSE_B in size
      //! @param PTR_SSE_A first SSE
      //! @param PTR_SSE_B second SSE
      //! @return true if PTR_SSE_A is smaller than PTR_SSE_B in size
      bool operator()
      (
        const util::PtrInterface< SSE> &PTR_SSE_A,
        const util::PtrInterface< const SSE> &PTR_SSE_B
      ) const;

    }; // class SSELessThanBySize

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSECompare
    //! @brief Function to be used in std::find_if to check for finding SSEs equal to the provided SSE
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date 22.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSECompare
    {
    public:

    //////////
    // data //
    //////////

      //! SSE to compare to
      const SSE &m_SSE;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a SSE reference
      //! @param THIS_SSE reference to SSE to be compared against
      SSECompare( const SSE &THIS_SSE);

    ///////////////
    // operators //
    ///////////////

      //! @brief compares the provided THIS_SSE with m_SSE to check it is the same sse
      //! @param THIS_SSE SSE to be compared
      //! @return true if THIS_SSE is equal to m_SSE
      bool operator()( const SSE &THIS_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
      //! @param PTR_SSE SSE to be compared
      //! @return true if PTR_SSE is equal to m_SSE
      bool operator()( const util::PtrInterface< SSE> &PTR_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
      //! @param PTR_SSE SSE to be compared
      //! @return true if PTR_SSE is equal to m_SSE
      bool operator()( const util::PtrInterface< const SSE> &PTR_SSE) const;

    }; // class SSECompare

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSECompareByIdentity
    //! @brief Function for searching SSEs that has same seqids, chain id and SSType with the given SSE
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date Nov 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSECompareByIdentity
    {
    public:

    //////////
    // data //
    //////////

      const size_t      m_BeginSeqID; //!< seq id of the first amino acid  of the SSE searched
      const size_t      m_EndSeqID;   //!< seq id of the last amino acid of the SSE searched
      const char        m_ChainID;    //!< chain id of the SSE searched
      const biol::SSType m_Type;      //!< type of the SSE searched

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a begin seq id, end seq id, chain id and SSType
      //! @param BEGIN_SEQ_ID seq id of the first amino acid  of the SSE searched
      //! @param END_SEQ_ID seq id of the last amino acid of the SSE searched
      //! @param CHAIN_ID chain id of the SSE searched
      //! @param SS_TYPE type of the SSE searched
      SSECompareByIdentity
      (
        const size_t        BEGIN_SEQ_ID,
        const size_t        END_SEQ_ID,
        const char          CHAIN_ID,
        const biol::SSType &SS_TYPE
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief compares the provided THIS_SSE with m_SSE to check it is the same sse
      //! @param THIS_SSE SSE to be compared
      //! @return true if THIS_SSE is equal to m_SSE
      bool operator()( const SSE &THIS_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
      //! @param PTR_SSE SSE to be compared
      //! @return true if PTR_SSE is equal to m_SSE
      bool operator()( const util::PtrInterface< SSE> &PTR_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check it is the same sse
      //! @param PTR_SSE SSE to be compared
      //! @return true if PTR_SSE is equal to m_SSE
      bool operator()( const util::PtrInterface< const SSE> &PTR_SSE) const;

    }; // class SSECompareByIdentitye

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSECompareOverlap
    //! @brief Function to be used in std::find_if queries in sets and other containers where SSEs that overlap with
    //! with the given SSE is being searched
    //!
    //! @see @link example_assemble_sse_compare.cpp @endlink
    //! @author karakam
    //! @date 22.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSECompareOverlap
    {
    public:

    //////////
    // data //
    //////////

      //! SSE to compare to
      const SSE &m_SSE;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a SSE reference
      //! @param THIS_SSE reference to SSE to be compared against
      SSECompareOverlap( const SSE &THIS_SSE);

    ///////////////
    // operators //
    ///////////////

      //! @brief compares the provided THIS_SSE with m_SSE to check for overlap
      //! @param THIS_SSE SSE to be checked for overlap
      //! @return true if THIS_SSE overlaps with m_SSE
      bool operator()( const SSE &THIS_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check for overlap
      //! @param PTR_SSE SSE to be checked for overlap
      //! @return true if PTR_SSE overlaps with m_SSE
      bool operator()( const util::PtrInterface< SSE> &PTR_SSE) const;

      //! @brief compares the provided PTR_SSE with m_SSE to check for overlap
      //! @param PTR_SSE SSE to be checked for overlap
      //! @return true if PTR_SSE overlaps with m_SSE
      bool operator()( const util::PtrInterface< const SSE> &PTR_SSE) const;

    }; // class SSECompareOverlap

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_SSE_COMPARE_H_
