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
#include "pdb/bcl_pdb_residue_simple.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "pdb/bcl_pdb_entry_types.h"
#include "pdb/bcl_pdb_line_criterium.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ResidueSimple::s_Instance
    (
      GetObjectInstances().AddInstance( new ResidueSimple())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ResidueSimple::ResidueSimple() :
      m_ResidueName(),
      m_ChainID( ' '),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' ')
    {
    }

    //! @brief construct ResidueSimple from ResidueName and ChainID
    //! @param RESIDUE_NAME name
    //! @param CHAIN_ID chain id
    ResidueSimple::ResidueSimple
    (
      const std::string &RESIDUE_NAME,
      const char CHAIN_ID
    ) :
      m_ResidueName( RESIDUE_NAME),
      m_ChainID( CHAIN_ID),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' ')
    {
    }

    //! @brief construct ResidueSimple from ResidueName, ChainID, ID, ICode
    //! @param RESIDUE_NAME name
    //! @param CHAIN_ID chain id
    //! @param PDB_ID pdb sequence id
    //! @parmam I_CODE pdb residue insertion code
    ResidueSimple::ResidueSimple
    (
      const std::string &RESIDUE_NAME,
      const char CHAIN_ID,
      const int PDB_ID,
      const char I_CODE
    ) :
      m_ResidueName( RESIDUE_NAME),
      m_ChainID( CHAIN_ID),
      m_PDBID( PDB_ID),
      m_ICode( I_CODE)
    {
    }

    //! @brief copy constructor
    ResidueSimple::ResidueSimple( const ResidueSimple &RESIDUE) :
      m_ResidueName( RESIDUE.m_ResidueName),
      m_ChainID( RESIDUE.m_ChainID),
      m_PDBID( RESIDUE.m_PDBID),
      m_ICode( RESIDUE.m_ICode)
    {
    }

    //! copy constructor
    ResidueSimple *ResidueSimple::Clone() const
    {
      return new ResidueSimple( *this);
    }

    //! destructor
    ResidueSimple::~ResidueSimple()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ResidueSimple::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to pdb ID
    //! @return pdb sequence id
    int ResidueSimple::GetPDBID() const
    {
      return m_PDBID;
    }

    //! @brief set ID
    //! @param PDB_ID new pdb id for residue
    void ResidueSimple::SetPDBID( const int PDB_ID)
    {
      m_PDBID = PDB_ID;
    }

    //! @brief return insertion code
    //! @return pdb sequence insertion code, if multiple residues have the same pdb id
    char ResidueSimple::GetICode() const
    {
      return m_ICode;
    }

    //! @brief set insertion code
    //! @param I_CODE insertion code for that residue
    void ResidueSimple::SetICode( const char I_CODE)
    {
      m_ICode = I_CODE;
    }

    //! @brief return chain id
    //! @return the chain id this residue belongs to
    char ResidueSimple::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief set chain id
    //! @param CHAIN_ID the chain id
    void ResidueSimple::SetChainID( const char CHAIN_ID)
    {
      m_ChainID = CHAIN_ID;
    }

    //! @brief return residue name
    //! @return three letter code residue name
    const std::string &ResidueSimple::GetResidueName() const
    {
      return m_ResidueName;
    }

    //! @brief set residue name
    //! @param RESIDUE_NAME new residue name for that residue
    void ResidueSimple::SetResidueName( const std::string &RESIDUE_NAME)
    {
      m_ResidueName = RESIDUE_NAME;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> ResidueSimple::AtomSerials() const
    {
      return storage::Set< size_t>();
    }

    //! @brief line criterium to locate atom/hetatom lines for that residue
    //! @param LINE_TYPE the line type of interest
    //! @return criterium to check if a line corresponds to that residue
    LineCriterium ResidueSimple::GetCriterium( const LineType &LINE_TYPE) const
    {
      LineCriterium criterium;

      if( LINE_TYPE == GetLineTypes().ATOM)
      {
        criterium.AddCriterium( GetEntryTypes().ATOMResidueName      , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().ATOMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().ATOMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().ATOMInsertionCode    , m_ICode);
      }
      else if( LINE_TYPE == GetLineTypes().HETATM)
      {
        criterium.AddCriterium( GetEntryTypes().HETATMResidueName      , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().HETATMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().HETATMInsertionCode    , m_ICode);
      }

      // end
      return criterium;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write ResidueSimple to STREAM
    std::ostream &ResidueSimple::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ResidueName, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID    , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID      , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode      , OSTREAM) << '\n';

      // end
      return OSTREAM;
    }

    //! read ResidueSimple from std::istream
    std::istream &ResidueSimple::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_ResidueName, ISTREAM);
      io::Serialize::Read( m_ChainID    , ISTREAM);
      io::Serialize::Read( m_PDBID      , ISTREAM);
      io::Serialize::Read( m_ICode      , ISTREAM);

      // return
      return ISTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool ResidueSimple::operator ==( const ResidueInterface &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.GetChainID()
             && m_PDBID   == RESIDUE_RHS.GetPDBID()
             && m_ICode   == RESIDUE_RHS.GetICode();
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool ResidueSimple::operator !=( const ResidueInterface &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool ResidueSimple::operator <( const ResidueInterface &RESIDUE_RHS) const
    {
      // compare chainid
      const char rhs_chain_id( RESIDUE_RHS.GetChainID());
      if( m_ChainID < rhs_chain_id)
      {
        return true;
      }

      if( m_ChainID > rhs_chain_id)
      {
        return false;
      }

      // compare pdb id
      const int rhs_pdb_id( RESIDUE_RHS.GetPDBID());
      if( m_PDBID < rhs_pdb_id)
      {
        return true;
      }

      if( m_PDBID > rhs_pdb_id)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.GetICode())
      {
        return true;
      }

      return false;
    }

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_LHS lhs residue
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool ResidueSimple::operator ==( const ResidueSimple &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.m_ChainID
             && m_PDBID   == RESIDUE_RHS.m_PDBID
             && m_ICode   == RESIDUE_RHS.m_ICode;
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool ResidueSimple::operator !=( const ResidueSimple &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool ResidueSimple::operator <( const ResidueSimple &RESIDUE_RHS) const
    {
      // compare chain id
      if( m_ChainID < RESIDUE_RHS.m_ChainID)
      {
        return true;
      }

      if( m_ChainID > RESIDUE_RHS.m_ChainID)
      {
        return false;
      }

      // compare pdb id
      if( m_PDBID < RESIDUE_RHS.m_PDBID)
      {
        return true;
      }

      if( m_PDBID > RESIDUE_RHS.m_PDBID)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.m_ICode)
      {
        return true;
      }

      return false;
    }

    //! @brief compare with biol amino acid if they are equal
    //! @param AMINO_ACID right hand side amino acid to compare with
    //! @return true, if chain id, pdb sequence id and insertion code match
    bool ResidueSimple::operator ==( const biol::AABase &AMINO_ACID) const
    {
      return    m_ChainID == AMINO_ACID.GetChainID()
             && m_PDBID   == AMINO_ACID.GetPdbID()
             && m_ICode   == AMINO_ACID.GetPdbICode();
    }

    //! @brief compare with biol amino acid if they are not equal
    //! @param AMINO_ACID right hand side amino acid to compare with
    //! @return true, if either chain id, pdb sequence id and insertion code do not match
    bool ResidueSimple::operator !=( const biol::AABase &AMINO_ACID) const
    {
      return !operator==( AMINO_ACID);
    }

  } // namespace pdb
} // namespace bcl
