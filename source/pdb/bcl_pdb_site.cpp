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
#include "pdb/bcl_pdb_site.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    //! @brief EvidenceCode as string
    //! @param EVIDENCE_CODE the EvidenceCode
    //! @return the string for the EvidenceCode
    const std::string &Site::GetEvidenceCodeDescriptor( const Site::EvidenceCode &EVIDENCE_CODE)
    {
      static const std::string s_descriptors[] =
      {
          "SOFTWARE",
          "AUTHOR",
          "UNKNOWN",
          GetStaticClassName< EvidenceCodeEnum>()
      };

      return s_descriptors[ EVIDENCE_CODE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Site::s_Instance
    (
      GetObjectInstances().AddInstance( new Site())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Site::Site()
    {
    }

    //! @brief construct from identifier, evidence code and description
    Site::Site( const std::string &NAME, const EvidenceCodeEnum &EVIDENCE_CODE, const std::string &DESCRIPTION) :
      m_Identifier( NAME),
      m_EvidenceCode( EVIDENCE_CODE),
      m_Description( DESCRIPTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Site
    Site *Site::Clone() const
    {
      return new Site( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Site::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief name/identifier of the site
    //! @return the name as string
    const std::string &Site::GetName() const
    {
      return m_Identifier;
    }

    //! @brief get evidence code
    //! @return enum EvicenceCodeEnum
    const Site::EvidenceCodeEnum &Site::GetEvidenceCode() const
    {
      return m_EvidenceCode;
    }

    //! @brief site description
    //! @return description of site
    const std::string &Site::GetDescription() const
    {
      return m_Description;
    }

    //! @brief access to the residues comprising the site that are within a chain
    //! @return ShPtrList to residues
    const storage::List< ResidueSimple> &Site::GetChainResidues() const
    {
      return m_ChainResidues;
    }

    //! @brief ligand that might be associated with that site
    //! @return ShPtr to Ligand; null if there is no ligand
    const util::ShPtr< Ligand> &Site::GetLigand() const
    {
      return m_Ligand;
    }

    //! @brief add a residue from a chain
    //! @param RESIDUE the residue to  add to the site
    void Site::AddChainResidue( const ResidueSimple &RESIDUE)
    {
      m_ChainResidues.PushBack( RESIDUE);
    }

    //! @brief add a residue from hetatm section
    //! @param RESIDUE the residue to  add to the site
    void Site::AddHetatmResidue( const ResidueSimple &RESIDUE)
    {
      m_HetResidues.PushBack( RESIDUE);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locators for the residues associated with the site definition
    //! @return the amino acid locators in a list
    storage::List< assemble::LocatorAA> Site::AALocators() const
    {
      storage::List< assemble::LocatorAA> aa_locators;

      // iterate through all residues
      for
      (
        storage::List< ResidueSimple>::const_iterator itr( m_ChainResidues.Begin()), itr_end( m_ChainResidues.End());
        itr != itr_end;
        ++itr
      )
      {
        aa_locators.PushBack( assemble::LocatorAA( itr->GetChainID(), itr->GetPDBID(), true));
      }

      // end
      return aa_locators;
    }

    //! @brief find ligand for site
    //! @param LIGANDS list of ligands, that might be associated with that site
    //! @return true, if a ligand was found, false otherwise
    bool Site::FindLigand( const util::ShPtrList< Ligand> &LIGANDS)
    {
      // check if this is a binding site
      static const std::string s_binding_site_string( "BINDING SITE FOR RESIDUE");

      // check if site is binding site for ligand
      if( m_Description.find( s_binding_site_string) == std::string::npos)
      {
        return false;
      }

      // ligand information
      const size_t pos_ligand_info( s_binding_site_string.length() + 1);
      const size_t ligand_info_length( 12);

      // name, chain and seq id
      std::stringstream ligand_info( util::TrimString( m_Description.substr( pos_ligand_info, ligand_info_length)));

      // ligand name
      std::string lig_name;
      ligand_info >> lig_name;
      lig_name = util::Format().R().W( 3).Fill( ' ')( lig_name);

      // chainid
      char chain_id( ' ');
      while( ligand_info.good() && chain_id == ' ')
      {
        ligand_info.get( chain_id);
      }

      // still expecting seq id
      if( !ligand_info.good())
      {
        return false;
      }

      std::string pdb_id_string;
      ligand_info >> pdb_id_string;

      if( util::TrimString( pdb_id_string).empty())
      {
        return false;
      }

      // insertion code could be the last character
      char pdb_icode( pdb_id_string[ pdb_id_string.length() - 1]);

      // if the insertion code is a digit, it is actually part of the seqid, and no insertion code is given
      if( std::isdigit( pdb_icode))
      {
        pdb_icode = ' ';
      }
      else
      {
        // found insertion code, the seq id string is one shorter
        pdb_id_string.erase( pdb_id_string.length() - 1, 1);
      }

      // the remaining string should be numerical
      if( !util::IsNumerical( pdb_id_string))
      {
        return false;
      }

      // the actual sequence id can be extracted
      const int seq_id( util::ConvertStringToNumericalValue< int>( pdb_id_string));

      BCL_MessageVrb
      (
        "site binds ligand: " + lig_name + ' ' + chain_id + ' ' + util::Format()( seq_id) + ' ' + pdb_icode
      );

      // iterate over ligands to find the correct one
      for
      (
        util::ShPtrList< Ligand>::const_iterator itr( LIGANDS.Begin()), itr_end( LIGANDS.End()); itr != itr_end; ++itr
      )
      {
        if( ( *itr)->GetChainID()     != chain_id)  continue;
        if( ( *itr)->GetResidueName() != lig_name)  continue;
        if( ( *itr)->GetPDBID()       != seq_id)    continue;
        if( ( *itr)->GetICode()       != pdb_icode) continue;

        // ligand found
        m_Ligand = *itr;
        return true;
      }

      // nothing found
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Site::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Identifier   , ISTREAM);
      io::Serialize::Read( m_EvidenceCode , ISTREAM);
      io::Serialize::Read( m_Description  , ISTREAM);
      io::Serialize::Read( m_ChainResidues, ISTREAM);
      io::Serialize::Read( m_HetResidues  , ISTREAM);
      io::Serialize::Read( m_Ligand       , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Site::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Identifier   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EvidenceCode , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainResidues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HetResidues  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ligand       , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
