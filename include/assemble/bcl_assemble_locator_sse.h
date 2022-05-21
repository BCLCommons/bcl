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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSE
    //! @brief This class is used for locating a specified SSE from a given Protein Model
    //! @details This class uses the member chain locator to find the corresponding chain with the correct chain id
    //! and then iterates over SSEs in that chain to find the one with the correct seq ids for the first and the last
    //! amino acids in that SSE and returns it
    //!
    //! @see @link example_assemble_locator_sse.cpp @endlink
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSE :
      public find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      char m_ChainID;   //!< Chain id
      int  m_StartAAID; //!< the SeqID or PDBID of the first amino acid in the SSE
      int  m_EndAAID;   //!< the SeqID or PDBID of the last amino acid in the SSE
      bool m_UsePDBID;  //!< use the pdb id instead of the seq id to locate the amino acid

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSSE() :
        m_ChainID(),
        m_StartAAID(),
        m_EndAAID(),
        m_UsePDBID()
      {
      }

      //! @brief construct from a ChainID, and the SeqIDs of the first and last amino acids in the SSE
      //! @param CHAINID char which indicates the chain
      //! @param SSE_START starting residue identifier
      //! @param SSE_END ending residue identifier
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorSSE( const char CHAINID, const int SSE_START, const int SSE_END, const bool USE_PDB_ID = false) :
        m_ChainID( CHAINID),
        m_StartAAID( SSE_START),
        m_EndAAID( SSE_END),
        m_UsePDBID( USE_PDB_ID)
      {
      }

      //! @brief clone constructor
      LocatorSSE *Clone() const
      {
        return new LocatorSSE( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetIdentification
      //! @return the GetIdentification
      const std::string GetIdentification() const;

      //! GetSSEID gives the LocatorSSE
      //! @return returns storage::VectorND< 2, int> with the start and end SeqIDs of the amino acids
      storage::VectorND< 2, int> GetSSEID() const;

      //! @brief gives an identifying string for this locator
      //! @return returns an identifying string for this locator
      std::string GetSSEIDString() const;

      //! SetSSEID changes the SSE identifiers
      //! @param NEW_SSE_START size_t which indicates the new SSE start
      //! @param NEW_SSE_END size_t which indicates the new SSE end
      void SetSSEID( const int NEW_SSE_START, const int NEW_SSE_END);

      //! returns the chain locator
      //! @return returns the const reference to the chain locator
      const char &GetChainID() const;

      //! @brief get the aa locator for the start residue
      //! @return a locator to locate the start AA in the SSE
      const LocatorAA StartAALocator() const;

      //! @brief get the aa locator for the end residue
      //! @return a locator to locate the end AA in the SSE
      const LocatorAA EndAALocator() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate a sse in the domain
      //! @param SSE_DOMAIN domain which the LocatorSSE refers to
      //! @return returns SiPtr to the SSE denoted by the LocatorSSE
      util::SiPtr< const SSE> Locate( const DomainInterface &SSE_DOMAIN) const;

    }; // class LocatorSSE

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_SSE_H_
