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
#include "assemble/bcl_assemble_locator_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LocatorSSE::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSSE())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetIdentification
    //! @return the GetIdentification
    const std::string LocatorSSE::GetIdentification() const
    {
      return util::Format()( m_ChainID) + " " + util::Format()( m_StartAAID) + " " + util::Format()( m_EndAAID);
    }

    //! GetSSEID gives the LocatorSSE
    //! @return returns storage::Pair< size_t, size_t> with the start and end SeqIDs of the amino acids
    storage::VectorND< 2, int> LocatorSSE::GetSSEID() const
    {
      return storage::VectorND< 2, int>( m_StartAAID, m_EndAAID);
    }

    //! @brief gives an identifying string for this locator
    //! @return returns an identifying string for this locator
    std::string LocatorSSE::GetSSEIDString() const
    {
      // construct and return the id string
      std::string id
      (
        util::Format()(           m_ChainID) + " " +
        util::Format().W( 4).R()( m_StartAAID) + " " +
        util::Format().W( 4).R()( m_EndAAID)
      );
      if( m_UsePDBID)
      {
        id += " use_pdb_id";
      }

      return id;
    }

    //! SetSSEID changes the SSE identifiers
    //! @param NEW_SSE_START size_t which indicates the new SSE start
    //! @param NEW_SSE_END size_t which indicates the new SSE end
    void LocatorSSE::SetSSEID( const int NEW_SSE_START, const int NEW_SSE_END)
    {
      // set start and end ids of SSE
      m_StartAAID = NEW_SSE_START;
      m_EndAAID = NEW_SSE_END;
    }

    //! returns the chain locator
    //! @return returns the const reference to the chain locator
    const char &LocatorSSE::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief get the aa locator for the start residue
    //! @return a locator to locate the start AA in the SSE
    const LocatorAA LocatorSSE::StartAALocator() const
    {
      return LocatorAA( m_ChainID, m_StartAAID, m_UsePDBID);
    }

    //! @brief get the aa locator for the end residue
    //! @return a locator to locate the end AA in the SSE
    const LocatorAA LocatorSSE::EndAALocator() const
    {
      return LocatorAA( m_ChainID, m_EndAAID, m_UsePDBID);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &LocatorSSE::GetAlias() const
    {
      static const std::string s_name( "LocatorSSE");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Locates a specified SSE in a protein model.");
      serializer.AddInitializer
      (
        "chain id",
        "chain id of the SSE",
        io::Serialization::GetAgent( &m_ChainID)
      );
      serializer.AddInitializer
      (
        "start id",
        "sequence ID of the first residue in the SSE",
        io::Serialization::GetAgent( &m_StartAAID)
      );
      serializer.AddInitializer
      (
        "end id",
        "sequence ID of the last residue in the SSE",
        io::Serialization::GetAgent( &m_EndAAID)
      );
      serializer.AddInitializer
      (
        "use pdb id",
        "use PDB ID instead of sequence ID to locate the residues",
        io::Serialization::GetAgent( &m_UsePDBID)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate a sse in the domain
    //! @param SSE_DOMAIN domain which the LocatorSSE refers to
    //! @return returns SiPtr to the SSE denoted by the LocatorSSE
    util::SiPtr< const SSE> LocatorSSE::Locate( const DomainInterface &SSE_DOMAIN) const
    {
      // get all sses
      const util::SiPtrVector< const SSE> sses( SSE_DOMAIN.GetSSEs());

      // iterate over the SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator itr( sses.Begin()), itr_end( sses.End());
        itr != itr_end;
        ++itr
      )
      {
        // check to see if the SeqIDs of the first and last amino acids of current SSE match m_StartAAID and m_EndAAID
        int aa_id_first( ( *itr)->GetFirstAA()->GetSeqID());
        int aa_id_last( ( *itr)->GetLastAA()->GetSeqID());

        if( m_UsePDBID)
        {
          BCL_MessageDbg( "using pbdb id");
          aa_id_first = ( *itr)->GetFirstAA()->GetPdbID();
          aa_id_last = ( *itr)->GetLastAA()->GetPdbID();
        }
        BCL_MessageDbg
        (
          "trying to find " + GetSSEIDString() + " comparing against " +
          util::Format()( aa_id_first) + " " + util::Format()( aa_id_last) + " "
          + util::Format()( ( *itr)->GetChainID()) + " which is " + ( *itr)->GetIdentification()
        );
        if( ( *itr)->GetChainID() == m_ChainID && aa_id_first == m_StartAAID && aa_id_last == m_EndAAID)
        {
          // if so then return this SSE
          return *itr;
        }
      }

      // if not found, return empty SiPtr< const SSE>
      BCL_MessageDbg
      (
        "sse starting with seqID " + util::Format()( m_StartAAID) + " and ending with "
        + util::Format()( m_EndAAID) + " in chain " + m_ChainID + " does not exist in domain " + SSE_DOMAIN.GetIdentification()
      );

      return util::SiPtr< const SSE>();
    }

  } // namespace assemble
} // namespace bcl
