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
#include "pdb/bcl_pdb_residue.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> Residue::s_Instance
    (
      GetObjectInstances().AddInstance( new Residue())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Residue::Residue() :
      m_ResidueName(),
      m_ChainID( ' '),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' '),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName and ChainID
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' '),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName, ChainID, ID, ICode
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID,
      const int  &PDBID,
      const char &ICODE
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( PDBID),
      m_ICode( ICODE),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName, ChainID, ID, ICOde, and Lines
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID,
      const int  &PDBID,
      const char &ICODE,
      const util::ShPtrList< Line> &LINES
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( PDBID),
      m_ICode( ICODE),
      m_Lines( LINES)
    {
    }

    //! @brief copy constructor
    Residue::Residue( const Residue &RESIDUE) :
      m_ResidueName( RESIDUE.m_ResidueName),
      m_ChainID( RESIDUE.m_ChainID),
      m_PDBID( RESIDUE.m_PDBID),
      m_ICode( RESIDUE.m_ICode),
      m_Lines( RESIDUE.m_Lines)
    {
    }

    //! @brief construct form Residue
    //! @param RESIDUE the residue to copy
    Residue::Residue( const ResidueInterface &RESIDUE) :
      m_ResidueName( RESIDUE.GetResidueName()),
      m_ChainID( RESIDUE.GetChainID()),
      m_PDBID( RESIDUE.GetPDBID()),
      m_ICode( RESIDUE.GetICode()),
      m_Lines()
    {
    }

    //! copy constructor
    Residue *Residue::Clone() const
    {
      return new Residue( *this);
    }

    //! destructor
    Residue::~Residue()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Residue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return Lines
    util::ShPtrList< Line> const &Residue::GetLines() const
    {
      return m_Lines;
    }

    //! change Lines
    util::ShPtrList< Line> &Residue::ChangeLines()
    {
      return m_Lines;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> Residue::AtomSerials() const
    {
      storage::Set< size_t> serials;

      // iterate over all atom lines
      for( util::ShPtrList< Line>::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End()); itr != itr_end; ++itr)
      {
        serials.Insert( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().ATOMSerial));
      }

      // end
      return serials;
    }

    //! @brief line criterium to locate atom/hetatom lines for that residue
    //! @param LINE_TYPE the line type of interest
    //! @return criterium to check if a line corresponds to that residue
    LineCriterium Residue::GetCriterium( const LineType &LINE_TYPE) const
    {
      LineCriterium criterium;

      if( LINE_TYPE == GetLineTypes().ATOM)
      {
        criterium.AddCriterium( GetEntryTypes().ATOMName             , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().ATOMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().ATOMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().ATOMInsertionCode    , m_ICode);
      }
      else if( LINE_TYPE == GetLineTypes().HETATM)
      {
        criterium.AddCriterium( GetEntryTypes().HETATMName             , m_ResidueName);
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

    //! write Residue to STREAM
    std::ostream &Residue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ResidueName, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID    , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID      , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode      , OSTREAM) << '\n';
      io::Serialize::Write( m_Lines      , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Residue from std::istream
    std::istream &Residue::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_ResidueName, ISTREAM);
      io::Serialize::Read( m_ChainID    , ISTREAM);
      io::Serialize::Read( m_PDBID      , ISTREAM);
      io::Serialize::Read( m_ICode      , ISTREAM);
      io::Serialize::Read( m_Lines      , ISTREAM);

      // return
      return ISTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool Residue::operator ==( const ResidueInterface &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.GetChainID()
             && m_PDBID   == RESIDUE_RHS.GetPDBID()
             && m_ICode   == RESIDUE_RHS.GetICode();
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool Residue::operator !=( const ResidueInterface &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool Residue::operator <( const ResidueInterface &RESIDUE_RHS) const
    {
      // compar chainid
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
    bool Residue::operator ==( const Residue &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.m_ChainID
             && m_PDBID   == RESIDUE_RHS.m_PDBID
             && m_ICode   == RESIDUE_RHS.m_ICode;
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool Residue::operator !=( const Residue &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool Residue::operator <( const Residue &RESIDUE_RHS) const
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

  } // namespace pdb
} // namespace bcl
