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
#include "assemble/bcl_assemble_locator_aa.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> LocatorAA::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorAA())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string description
    //! @return gives formatted string description
    std::string LocatorAA::GetIdentification() const
    {
      std::stringstream write;
      io::Serialize::Write( GetLocatorChain().GetChainID(), write, 1);
      io::Serialize::Write( GetAAID(), write, 1);
      io::Serialize::Write( std::string( "UsePDBID"), write, 1);
      io::Serialize::Write( m_UsePDBID, write, 1);
      return write.str();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorAA::ReadIdentification( std::istream &ISTREAM)
    {
      char chain_id;
      io::Serialize::Read( chain_id, ISTREAM);
      m_LocatorChain.SetChainID( chain_id);
      std::string use_pdb_id;
      io::Serialize::Read( m_AA_ID, ISTREAM);
      io::Serialize::Read( use_pdb_id, ISTREAM);
      io::Serialize::Read( m_UsePDBID, ISTREAM);
      return ISTREAM;
    }

    //! @brief returns a bool indicating if the pdb ids will be used instead of the seq ids
    //! @return boolean true if pdb ids should be used - false otherwise
    bool LocatorAA::GetUsePDBID() const
    {
      return m_UsePDBID;
    }

    //! @brief returns aa type, if known
    //! @return aa type, if known
    const biol::AAType &LocatorAA::GetAAType() const
    {
      return m_AAType;
    }

    //! @brief GetAAID gives the seq id of the amino acid to be located
    //! @return returns "m_AA_ID" int
    int const &LocatorAA::GetAAID() const
    {
      return m_AA_ID;
    }

    //! @brief SetAAID changes the amino acid id
    //! @param AA_ID int which indicates the new amino acid id
    void LocatorAA::SetAAID( const int AA_ID)
    {
      // set "m_AA_ID" to "AA_ID"
      m_AA_ID = AA_ID;
    }

    //! @brief returns the chain locator
    //! @return returns the const reference to the chain locator
    const LocatorChain &LocatorAA::GetLocatorChain() const
    {
      return m_LocatorChain;
    }

    //! @brief SetLocatorChain the chain locator to a new one
    //! @param LOCATOR_CHAIN the new locator chain which "m_LocatorChain" will be set to
    void LocatorAA::SetLocatorChain( const LocatorChain &LOCATOR_CHAIN)
    {
      // set "m_LocatorChain" to "LOCATOR_CHAIN"
      m_LocatorChain = LOCATOR_CHAIN;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorAA::GetAlias() const
    {
      static const std::string s_Name( "LocatorAA");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetAA finds the amino acid denoted by the LocatorAA
    //! @param SSE SSE which the LocatorAA refers to
    //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
    util::SiPtr< const biol::AABase> LocatorAA::Locate( const SSE &SSE) const
    {
      // iterate through the residues of "SSE"
      for
      (
        biol::AASequence::const_iterator
          itr_aa( SSE.Begin()), itr_aa_end( SSE.End());
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        // true if the seq id and the chain ids match
        if
        (
          ( m_UsePDBID ? ( *itr_aa)->GetPdbID() : ( *itr_aa)->GetSeqID()) == m_AA_ID
          && ( *itr_aa)->GetChainID() == m_LocatorChain.GetChainID()
        )
        {
          BCL_MessageDbg( ( *itr_aa)->GetIdentification());

          BCL_Assert
          (
            !m_AAType.IsDefined() || !m_AAType->IsNaturalAminoAcid() || m_AAType == ( *itr_aa)->GetType(),
            "AA types did not match " + m_AAType->GetName() + " != " + ( *itr_aa)->GetType().GetName()
            + " for residue: " + util::Format()( m_AA_ID)
            + " ; check input files and formats "
          );
          // return si ptr to the residue denoted by "itr_aa"
          return util::SiPtr< const biol::AABase>( *itr_aa);
        }
      }

      // return empty siptr since the residue was not found in "SSE"
      return util::SiPtr< const biol::AABase>();
    }

    //! @brief GetAA finds the amino acid denoted by the LocatorAA
    //! @param SSE SSE which the LocatorAA refers to
    //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
    util::SiPtr< biol::AABase> LocatorAA::Locate( SSE &SSE) const
    {
      // iterate through the residues of "SSE"
      for
      (
        biol::AASequence::iterator
          itr_aa( SSE.Begin()), itr_aa_end( SSE.End());
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        // true if the seq id and the chain ids match
        if
        (
          ( m_UsePDBID ? ( *itr_aa)->GetPdbID() : ( *itr_aa)->GetSeqID()) == m_AA_ID
          && ( *itr_aa)->GetChainID() == m_LocatorChain.GetChainID()
        )
        {
          BCL_MessageDbg( ( *itr_aa)->GetIdentification());

          BCL_Assert
          (
            !m_AAType.IsDefined() || !m_AAType->IsNaturalAminoAcid() || m_AAType == ( *itr_aa)->GetType(),
            "AA types did not match " + m_AAType->GetName() + " != " + ( *itr_aa)->GetType().GetName()
            + " for residue: " + util::Format()( m_AA_ID)
            + " ; check input files and formats "
          );
          // return si ptr to the residue denoted by "itr_aa"
          return util::SiPtr< biol::AABase>( *itr_aa);
        }
      }

      // return empty siptr since the residue was not found in "SSE"
      return util::SiPtr< biol::AABase>();
    }

    //! @brief Locate finds the amino acid denoted by the LocatorAA in a list of residues
    //! @param RESIDUE_LIST SSE the list of residues which will be searched for the residue of interest
    //! @return returns SiPtr to the amino acid denoted by the LocatorAA
    util::SiPtr< const biol::AABase> LocatorAA::Locate( const util::SiPtrVector< const biol::AABase> &RESIDUE_LIST) const
    {
      // iterate through the residues in "RESIDUE_LIST"
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          resi_itr( RESIDUE_LIST.Begin()), resi_itr_end( RESIDUE_LIST.End());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        // get seq id or pdb id
        int identifier;
        if( m_UsePDBID)
        {
          identifier = ( *resi_itr)->GetPdbID();
        }
        else
        {
          identifier = ( *resi_itr)->GetSeqID();
        }
        // true if the seq id and the chain id of the residue denoted by "resi_itr" matches "CHAIN_ID" and "SEQ_ID"
        if( identifier == m_AA_ID && ( *resi_itr)->GetChainID() == m_LocatorChain.GetChainID())
        {
          BCL_Assert
          (
            !m_AAType.IsDefined() || !m_AAType->IsNaturalAminoAcid() || m_AAType == ( *resi_itr)->GetType(),
            "AA types did not match " + m_AAType->GetName() + " != " + ( *resi_itr)->GetType().GetName()
            + " for residue: " + util::Format()( m_AA_ID)
            + " ; check input files and formats "
          );
          // return "resi_itr" indicating that the residue of interest exists in "RESI_LIST"
          return *resi_itr;
        }
      }

      // residue with chain id "SEQ_ID" and seq id "SEQ_ID" was not found in "RESI_LIST" so return empty SiPtr
      // indicating that it was not found
      return util::SiPtr< const biol::AABase>();
    }

    //! @brief Locate finds the amino acid denoted by the LocatorAA in a list of residues
    //! @param RESIDUE_LIST SSE the list of residues which will be searched for the residue of interest
    //! @return returns SiPtr to the amino acid denoted by the LocatorAA
    util::SiPtr< biol::AABase> LocatorAA::Locate( util::SiPtrVector< biol::AABase> &RESIDUE_LIST) const
    {
      // iterate through the residues in "RESIDUE_LIST"
      for
      (
        util::SiPtrVector< biol::AABase>::iterator
          resi_itr( RESIDUE_LIST.Begin()), resi_itr_end( RESIDUE_LIST.End());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        // get seq id or pdb id
        int identifier( m_UsePDBID ? ( *resi_itr)->GetPdbID() : ( *resi_itr)->GetSeqID());
        // true if the seq id and the chain id of the residue denoted by "resi_itr" matches "CHAIN_ID" and "SEQ_ID"
        if( identifier == m_AA_ID && ( *resi_itr)->GetChainID() == m_LocatorChain.GetChainID())
        {
          if( m_AAType.IsDefined() && m_AAType->IsNaturalAminoAcid() && m_AAType != ( *resi_itr)->GetType())
          {
            BCL_MessageVrb
            (
              "AA types did not match " + m_AAType->GetName() + " != " + ( *resi_itr)->GetType().GetName()
              + " for residue: " + util::Format()( m_AA_ID)
              + " ; check input files and formats "
            );
            return util::SiPtr< biol::AABase>();
          }
          // return "resi_itr" indicating that the residue of interest exists in "RESI_LIST"
          return *resi_itr;
        }
      }

      // residue with chain id "SEQ_ID" and seq id "SEQ_ID" was not found in "RESI_LIST" so return empty SiPtr
      // indicating that it was not found
      return util::SiPtr< biol::AABase>();
    }

    //! @brief GetSSE finds the SSE which contains the amino acid denoted by the LocatorAA
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
    //! @return returns SharedPointer to the SSE containing the amino acid denoted by the LocatorAA
    util::SiPtr< const SSE> LocatorAA::LocateSSE( const ProteinModel &PROTEIN_MODEL) const
    {
      // iterate over the SSEs of the chain to see which SSE contains the amino acid
      auto chain( m_LocatorChain.Locate( PROTEIN_MODEL));
      if( !chain.IsDefined())
      {
        return util::SiPtr< const SSE>();
      }
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          itr_begin( chain->GetData().Begin()),
          itr_end( chain->GetData().End());
        itr_begin != itr_end;
        ++itr_begin
      )
      {
        // get seq id or pdb id of the first and last residues of the sse
        int first_identifier, last_identifier;

        if( m_UsePDBID)
        {
          first_identifier = ( *itr_begin)->GetFirstAA()->GetPdbID();
          last_identifier = ( *itr_begin)->GetLastAA()->GetPdbID();
        }
        else
        {
          first_identifier = ( *itr_begin)->GetFirstAA()->GetSeqID();
          last_identifier = ( *itr_begin)->GetLastAA()->GetSeqID();
        }

        // check to see if the SeqIDs of the first and last amino acids of current SSE contain "m_AA_ID"
        if
        (
          first_identifier <= m_AA_ID &&
          last_identifier >= m_AA_ID
        )
        {
          // if so then return the SSE denoted by "itr_begin"
          return *itr_begin;
        }
      }

      // return empty shared pointer if SSE containing the amino acid is not found
      BCL_MessageDbg
      (
        "sse containing " + util::Format()( m_AA_ID) + " does not exist in protein model"
      );
      return util::SiPtr< const SSE>();
    }

    //! @brief GetAA finds the amino acid denoted by the LocatorAA
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
    //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
    util::SiPtr< const biol::AABase> LocatorAA::Locate( const ProteinModel &PROTEIN_MODEL) const
    {
      const auto located_sse( this->LocateSSE( PROTEIN_MODEL));
      if( located_sse.IsDefined())
      {
        return this->Locate( *located_sse);
      }

      // if "sse" is not defined
      // return empty simple pointer
      BCL_MessageDbg
      (
        "amino acid with seqid " + util::Format()( m_AA_ID) + " in chain " + m_LocatorChain.GetChainID()
        + " does not exist in protein model"
      );
      return util::SiPtr< const biol::AABase>();
    }

    //! @brief GetSSE finds the SSE which contains the amino acid denoted by the LocatorAA
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
    //! @return returns SharedPointer to the SSE containing the amino acid denoted by the LocatorAA
    util::SiPtr< SSE> LocatorAA::LocateSSE( ProteinModel &PROTEIN_MODEL) const
    {
      // iterate over the SSEs of the chain to see which SSE contains the amino acid
      auto chain( m_LocatorChain.Locate( PROTEIN_MODEL));
      if( !chain.IsDefined())
      {
        return util::SiPtr< SSE>();
      }
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator
          itr_begin( chain->GetData().Begin()),
          itr_end( chain->GetData().End());
        itr_begin != itr_end;
        ++itr_begin
      )
      {
        // get seq id or pdb id of the first and last residues of the sse
        int first_identifier, last_identifier;

        if( m_UsePDBID)
        {
          first_identifier = ( *itr_begin)->GetFirstAA()->GetPdbID();
          last_identifier = ( *itr_begin)->GetLastAA()->GetPdbID();
        }
        else
        {
          first_identifier = ( *itr_begin)->GetFirstAA()->GetSeqID();
          last_identifier = ( *itr_begin)->GetLastAA()->GetSeqID();
        }

        // check to see if the SeqIDs of the first and last amino acids of current SSE contain "m_AA_ID"
        if
        (
          first_identifier <= m_AA_ID &&
          last_identifier >= m_AA_ID
        )
        {
          // if so then return the SSE denoted by "itr_begin"
          util::ShPtr< SSE> sse( *itr_begin);
          return util::SiPtr< SSE>( sse);
        }
      }

      // return empty shared pointer if SSE containing the amino acid is not found
      BCL_MessageDbg
      (
        "sse containing " + util::Format()( m_AA_ID) + " does not exist in protein model"
      );
      return util::SiPtr< SSE>();
    }

    //! @brief GetAA finds the amino acid denoted by the LocatorAA
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
    //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
    util::SiPtr< biol::AABase> LocatorAA::Locate( ProteinModel &PROTEIN_MODEL) const
    {
      util::SiPtr< SSE> located_sse( this->LocateSSE( PROTEIN_MODEL));
      if( located_sse.IsDefined())
      {
        return this->Locate( *located_sse);
      }

      // if "sse" is not defined
      // return empty simple pointer
      BCL_MessageDbg
      (
        "amino acid with seqid " + util::Format()( m_AA_ID) + " in chain " + m_LocatorChain.GetChainID()
        + " does not exist in protein model"
      );
      return util::SiPtr< biol::AABase>();
    }

    //! @brief GetAA finds the amino acid denoted by the LocatorAA
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
    //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
    util::SiPtr< const biol::AABase> LocatorAA::Locate( const DomainInterface &PROTEIN_MODEL) const
    {
      auto sses( PROTEIN_MODEL.GetSSEs());
      for( auto itr_sse( sses.Begin()), itr_sse_end( sses.End()); itr_sse != itr_sse_end; ++itr_sse)
      {
        if( ( *itr_sse)->GetChainID() != m_LocatorChain.GetChainID())
        {
          continue;
        }
        const auto &last_aa( *( *itr_sse)->GetLastAA());
        const int last_id( m_UsePDBID ? last_aa.GetPdbID() : last_aa.GetSeqID());
        if( last_id >= m_AA_ID)
        {
          const auto &first_aa( *( *itr_sse)->GetFirstAA());
          const int first_id( m_UsePDBID ? first_aa.GetPdbID() : first_aa.GetSeqID());
          if( first_id <= m_AA_ID)
          {
            return this->Locate( **itr_sse);
          }
        }
      }

      // if "sse" is not defined
      // return empty simple pointer
      BCL_MessageDbg
      (
        "amino acid with seqid " + util::Format()( m_AA_ID) + " in chain " + m_LocatorChain.GetChainID()
        + " does not exist in protein model"
      );
      return util::SiPtr< const biol::AABase>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorAA::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Is used for locating a specified amino acid from a protein model. This class is used for locating "
        "a specified amino acid from a given protein model. It uses first the member chain locator to find the "
        "corresponding chain and then iterates over the amino acids of the SSEs to find the amino acid with the"
        "specified amino acid seq id."
      );

      parameters.AddInitializer
      (
        "locator_chain",
        "the chain locator that should be used to get the chain with the desired residue",
        io::Serialization::GetAgentWithRange( &m_LocatorChain.GetChainID(), 'A', 'Z')
      );

      parameters.AddInitializer
      (
        "seq_id",
        "the seq id of the desired residue",
        io::Serialization::GetAgent( &m_AA_ID),
        "1"
      );

      parameters.AddInitializer
      (
        "aa type",
        "the aa type; useful to verify that the same sequences/ids/formats are being used",
        io::Serialization::GetAgent( &m_AAType),
        "Undefined"
      );

      parameters.AddInitializer
      (
        "use_pdb_id",
        "boolean if true, the pdb id will be used for location (as opposed to the seq id). 1=true;0=false",
        io::Serialization::GetAgent( &m_UsePDBID),
        "0"
      );

      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorAA::Read( std::istream &ISTREAM)
    {
      return util::SerializableInterface::Read( ISTREAM);
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return util::SerializableInterface::Write( OSTREAM, INDENT);
    }

    //! @brief less than operator for comparing two locator aa objects by seq and chain id
    //! @param LOCATOR_LHS the first locator
    //! @param LOCATOR_RHS the second locator
    //! @return boolean true if LOCATOR_LHS is less than LOCATOR_RHS - false otherwise
    bool operator<( const LocatorAA &LOCATOR_LHS, const LocatorAA &LOCATOR_RHS)
    {
      // compare chain id
      if( LOCATOR_LHS.GetLocatorChain().GetChainID() < LOCATOR_RHS.GetLocatorChain().GetChainID())
      {
        return true;
      }

      if( LOCATOR_LHS.GetLocatorChain().GetChainID() > LOCATOR_RHS.GetLocatorChain().GetChainID())
      {
        return false;
      }

      // chain id match
      // compare seq id
      return LOCATOR_LHS.GetAAID() < LOCATOR_RHS.GetAAID();
    }

  } // namespace assemble
} // namespace bcl
