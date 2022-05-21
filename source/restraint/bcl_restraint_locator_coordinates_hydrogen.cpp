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
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <cctype>
#include <functional>

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesHydrogen::s_Instance
    (
      util::Enumerated< find::LocatorCoordinatesInterface< assemble::ProteinModel> >::AddInstance
      (
        util::Enumerated< assemble::LocatorAtomCoordinatesInterface>::AddInstance
        (
          new LocatorCoordinatesHydrogen()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorCoordinatesHydrogen::LocatorCoordinatesHydrogen() :
      m_ChainID( 'A'),
      m_SeqID( util::GetUndefined< int>()),
      m_AtomTypeString( biol::GetAtomTypes().e_Undefined),
      m_AtomType( biol::GetAtomTypes().e_Undefined)
    {
    }

    //! @brief construct from chain id, sequence id, and atom type
    //! @param CHAIN_ID chain id
    //! @param SEQ_ID sequence id
    //! @param ATOM_TYPE atom type string
    LocatorCoordinatesHydrogen::LocatorCoordinatesHydrogen
    (
      const char CHAIN_ID,
      const int SEQ_ID,
      const std::string &ATOM_TYPE
    ) :
      m_ChainID( CHAIN_ID),
      m_SeqID( SEQ_ID),
      m_AtomTypeString( ATOM_TYPE),
      m_AtomType( biol::GetAtomTypes().e_Undefined)
    {
      // make string upper case
      std::transform
      (
        m_AtomTypeString.begin(),
        m_AtomTypeString.end(),
        m_AtomTypeString.begin(),
        []( unsigned char c) { return std::toupper( c);}
      );

      // try to get the real atom type
      m_AtomType = GetAtomTypeFromString();
    }

    //! @brief Clone function
    //! @return pointer to new LocatorCoordinatesHydrogen
    LocatorCoordinatesHydrogen *LocatorCoordinatesHydrogen::Clone() const
    {
      return new LocatorCoordinatesHydrogen( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorCoordinatesHydrogen::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorCoordinatesHydrogen::GetAlias() const
    {
      static std::string s_alias( "Hydrogen");
      return s_alias;
    }

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorCoordinatesHydrogen::GetIdentification() const
    {
      return util::Format()( m_ChainID) + " " + util::Format()( m_SeqID) + " " + m_AtomTypeString;
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorCoordinatesHydrogen::ReadIdentification( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_AtomTypeString, ISTREAM);

      return ISTREAM;
    }

    //! @brief returns reference to undefined AA TYPE, this function is required by the interface but not used for
    //!        this class
    //! @return reference to undefined AA TYPE
    const biol::AAType &LocatorCoordinatesHydrogen::GetAAType() const
    {
      return biol::GetAATypes().e_Undefined;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates the desired coordinates from a protein model
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @return the coordinates of the atom denoted by the LocatorAtom
    linal::Vector3D LocatorCoordinatesHydrogen::Locate( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      return LocateAtomCopy( PROTEIN_MODEL).GetCoordinates();
    }

    //! @brief locates an atom from a protein model, creating H or HA atoms as needed - other hydrogen atoms will have
    //!        m_AtomType as an atom type but with coords of the first side chain atom
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return copy of atom that has been located from PROTEIN_MODEL
    biol::Atom LocatorCoordinatesHydrogen::LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize static undefined atom
      static const biol::Atom s_undefined_atom;

      // try to locate the residue that contains the atom
      const util::SiPtr< const biol::AABase> sp_residue( LocateAA( PROTEIN_MODEL));

      // if the residue was not located
      if( !sp_residue.IsDefined())
      {
        // return an undefined atom
        return s_undefined_atom;
      }

      // convert the string to an atom type
      biol::AtomType atom_type( m_AtomType.IsDefined() ? m_AtomType : GetAtomTypeFromString( sp_residue->GetType()));

      // if the atom type is undefined
      if( atom_type == biol::GetAtomTypes().e_Undefined)
      {
        // return undefined atom
        return s_undefined_atom;
      }

      // try to locate the atom
      const util::SiPtr< const biol::Atom> sp_atom
      (
        assemble::LocatorAtom::LocateAtomFromModel( PROTEIN_MODEL, m_ChainID, m_SeqID, atom_type)
      );

      // if the atom was located
      if( sp_atom.IsDefined() && sp_atom->GetCoordinates().IsDefined())
      {
        // return the atom
        return *sp_atom;
      }

      // if the atom to be located is not hydrogen
      if
      (
        atom_type->GetElementType() != chemistry::GetElementTypes().e_Hydrogen &&
        atom_type != biol::GetAtomTypes().O1)
      {
        // return a default atom
        return s_undefined_atom;
      }

      // if the atom type is H
      if( atom_type == biol::GetAtomTypes().H)
      {
        // locate the previous residue
        const assemble::LocatorAA aa_locator( m_ChainID, m_SeqID - 1);
        const util::SiPtr< const biol::AABase> sp_prev_aa( aa_locator.Locate( PROTEIN_MODEL));

        // if the previous AA was not located
        if( !sp_prev_aa.IsDefined())
        {
          // return an undefined atom
          return s_undefined_atom;
        }

        // generate an H atom
        return biol::AABackBoneCompleter::GenerateHydrogen
        (
          *sp_residue, sp_prev_aa->GetAtom( biol::GetAtomTypes().C).GetCoordinates()
        );
      }
      // if the atom type is HA
      else if( atom_type == biol::GetAtomTypes().HA || atom_type == biol::GetAtomTypes().HA3)
      {
        // generate an HA atom
        return biol::AABackBoneCompleter::GenerateHA( *sp_residue);
      }

      // return an atom with the first side chain atom coords but with the requested atom type in ATOM_LOCATOR
      return biol::Atom( sp_residue->GetFirstSidechainAtom().GetCoordinates(), atom_type);
    }

    //! @brief takes a string of the Atom type as well as the amino acid type and converts them to the
    //!        correct IUPAC type
    //! @param AA_TYPE the AAType of the atom from the restraint file
    //! @return the IUPAC formatted atom type
    biol::AtomType LocatorCoordinatesHydrogen::GetAtomTypeFromString( const biol::AAType &AA_TYPE) const
    {
      // initialize return type
      biol::AtomType atom_type
      (
        AA_TYPE.IsDefined() ?
          AA_TYPE->GetAtomTypeFromAtomName( m_AtomTypeString) :
          biol::GetAtomTypes().TypeFromPDBAtomName( m_AtomTypeString)
      );

      // if the atom type is defined
      if( atom_type.IsDefined())
      {
        // return it
        return atom_type;
      }

      // find the AA in the map
      const storage::Map< std::string, biol::AtomType> &atom_strings_map( GetPseudoAtomMap().GetValue( AA_TYPE));

      // check to see if the atom_string is in the map
      if( atom_strings_map.Has( m_AtomTypeString))
      {
        // return the type
        return atom_strings_map.GetValue( m_AtomTypeString);
      }

      // if atom is not in the map but has the *, # or % symbol after the atom type
      if
      (
        !m_AtomTypeString.empty() &&
        (
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '*' ||
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '#' ||
          m_AtomTypeString[ m_AtomTypeString.length() - 1] == '%'
        )
      )
      {
        // try replacing the last char
        atom_type = ReplacePseudoAtomChar( m_AtomTypeString, AA_TYPE);

        // return the type if it is defined
        if( atom_type.IsDefined())
        {
          return atom_type;
        }

        // iterate from 1 to 3
        for( char i( '1'); i != '4'; ++i)
        {
          // set the new string HG* to HG11
          std::string atom_string( m_AtomTypeString.substr( 0, m_AtomTypeString.length() - 1) + i + i);

          // try replacing the last char again
          atom_type = ReplacePseudoAtomChar( atom_string, AA_TYPE);

          // return the type if it is defined
          if( atom_type.IsDefined())
          {
            return atom_type;
          }
        }
      }

      // end
      return atom_type;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorCoordinatesHydrogen::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_AtomTypeString, ISTREAM);
      io::Serialize::Read( m_AtomType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorCoordinatesHydrogen::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeqID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypeString, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief tries replacing the last character in an atom name string with 1, 2, or 3 and getting the atom type
    //! @param ATOM_NAME atom name string
    //! @param AA_TYPE the AAType of the atom from the restraint file
    //! @return atom type
    biol::AtomType LocatorCoordinatesHydrogen::ReplacePseudoAtomChar
    (
      const std::string &ATOM_NAME,
      const biol::AAType &AA_TYPE
    )
    {
      // set atom type to undefined
      biol::AtomType atom_type( biol::GetAtomTypes().e_Undefined);

      // copy the string
      std::string atom_string( ATOM_NAME);

      // if the AA type is not defined
      if( AA_TYPE == biol::GetAATypes().e_Undefined)
      {
        // change the last char to 1 and return
        atom_string[ atom_string.length() - 1] = '1';
        return biol::GetAtomTypes().TypeFromPDBAtomName( atom_string);
      }

      // iterate from 1 to 3
      for( char i( '1'); i != '4'; ++i)
      {
        // change the last character
        atom_string[ atom_string.length() - 1] = i;

        // try to get the type again
        atom_type = AA_TYPE->GetAtomTypeFromAtomName( atom_string);

        // if the atom_type is now defined
        if( atom_type.IsDefined())
        {
          // return it
          return atom_type;
        }
      }

      // return undefined
      return atom_type;
    }

    //! @brief gets a map of PsuedoAtoms that can be mapped back to an IUPAC type
    //! @return the map of these PsuedoAtoms
    const storage::Map
    <
      biol::AAType, storage::Map< std::string, biol::AtomType>
    > &LocatorCoordinatesHydrogen::GetPseudoAtomMap()
    {
      // create a map to store all of the data
      static storage::Map< biol::AAType, storage::Map< std::string, biol::AtomType> > s_map;
      if( s_map.IsEmpty())
      {
        // undefined (backbone or not residue-specific)
        storage::Map< std::string, biol::AtomType> undefined_data;
        undefined_data[ "HN"] = biol::GetAtomTypes().H;
        undefined_data[ "MTS"] = biol::GetAtomTypes().O1;
        undefined_data[ "MTSL"] = biol::GetAtomTypes().O1;
        s_map[ biol::GetAATypes().e_Undefined] = undefined_data;
        // GLY
        storage::Map< std::string, biol::AtomType> gly_data;
        gly_data[ "HN"] = biol::GetAtomTypes().H;
        gly_data[ "QPA"] = biol::GetAtomTypes().HA2;
        gly_data[ "QA"] = biol::GetAtomTypes().HA2;
        s_map[ biol::GetAATypes().GLY] = gly_data;
        // ALA
        storage::Map< std::string, biol::AtomType> ala_data;
        ala_data[ "HN"] = biol::GetAtomTypes().H;
        ala_data[ "MB"] = biol::GetAtomTypes().HB1;
        ala_data[ "QB"] = biol::GetAtomTypes().HB1;
        s_map[ biol::GetAATypes().ALA] = ala_data;
        // VAL
        storage::Map< std::string, biol::AtomType> val_data;
        val_data[ "HN"] = biol::GetAtomTypes().H;
        val_data[ "MG1"] = biol::GetAtomTypes().HG11;
        val_data[ "QG1"] = biol::GetAtomTypes().HG11;
        val_data[ "MG2"] = biol::GetAtomTypes().HG11;
        val_data[ "QG2"] = biol::GetAtomTypes().HG11;
        val_data[ "QQG"] = biol::GetAtomTypes().HG11;
        val_data[ "HG#"] = biol::GetAtomTypes().HG11;
        s_map[ biol::GetAATypes().VAL] = val_data;
        // ILE
        storage::Map< std::string, biol::AtomType> ile_data;
        ile_data[ "HN"] = biol::GetAtomTypes().H;
        ile_data[ "HG2"] = biol::GetAtomTypes().HG12;
        ile_data[ "HG3"] = biol::GetAtomTypes().HG13;
        ile_data[ "QPG"] = biol::GetAtomTypes().HG12;
        ile_data[ "QG1"] = biol::GetAtomTypes().HG12;
        ile_data[ "MG2"] = biol::GetAtomTypes().HG21;
        ile_data[ "QG2"] = biol::GetAtomTypes().HG21;
        ile_data[ "MG"] = biol::GetAtomTypes().HG21;
        ile_data[ "QG"] = biol::GetAtomTypes().HG21;
        ile_data[ "MD"] = biol::GetAtomTypes().HG11;
        ile_data[ "QD1"] = biol::GetAtomTypes().HG11;
        ile_data[ "HD#"] = biol::GetAtomTypes().HD11;
        s_map[ biol::GetAATypes().ILE] = ile_data;
        // LEU
        storage::Map< std::string, biol::AtomType> leu_data;
        leu_data[ "HN"] = biol::GetAtomTypes().H;
        leu_data[ "QPB"] = biol::GetAtomTypes().HB2;
        leu_data[ "QB"] = biol::GetAtomTypes().HB2;
        leu_data[ "MD1"] = biol::GetAtomTypes().HD11;
        leu_data[ "QD1"] = biol::GetAtomTypes().HD11;
        leu_data[ "MD2"] = biol::GetAtomTypes().HD21;
        leu_data[ "QD2"] = biol::GetAtomTypes().HD21;
        leu_data[ "CD1"] = biol::GetAtomTypes().CD1;
        leu_data[ "QQD"] = biol::GetAtomTypes().HD11;
        leu_data[ "MDX"] = biol::GetAtomTypes().HD11;
        leu_data[ "HD#"] = biol::GetAtomTypes().HD11;
        s_map[ biol::GetAATypes().LEU] = leu_data;
        // PHE
        storage::Map< std::string, biol::AtomType> phe_data;
        phe_data[ "HN"] = biol::GetAtomTypes().H;
        phe_data[ "QPB"] = biol::GetAtomTypes().HB2;
        phe_data[ "QB"] = biol::GetAtomTypes().HB2;
        phe_data[ "HD"] = biol::GetAtomTypes().HD1;
        phe_data[ "CG"] = biol::GetAtomTypes().HD1;
        phe_data[ "HE"] = biol::GetAtomTypes().HE1;
        phe_data[ "CZ"] = biol::GetAtomTypes().HE1;
        phe_data[ "HZ"] = biol::GetAtomTypes().HZ;
        phe_data[ "QD"] = biol::GetAtomTypes().HD1;
        phe_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().PHE] = phe_data;
        // PRO
        storage::Map< std::string, biol::AtomType> pro_data;
        pro_data[ "QPB"] = biol::GetAtomTypes().HB2;
        pro_data[ "QB"] = biol::GetAtomTypes().HB2;
        pro_data[ "QPG"] = biol::GetAtomTypes().HG2;
        pro_data[ "QG"] = biol::GetAtomTypes().HG2;
        pro_data[ "QPD"] = biol::GetAtomTypes().HD2;
        pro_data[ "QD"] = biol::GetAtomTypes().HD2;
        s_map[ biol::GetAATypes().PRO] = pro_data;
        // MET
        storage::Map< std::string, biol::AtomType> met_data;
        met_data[ "HN"] = biol::GetAtomTypes().H;
        met_data[ "QPB"] = biol::GetAtomTypes().HB2;
        met_data[ "QB"] = biol::GetAtomTypes().HB2;
        met_data[ "QPG"] = biol::GetAtomTypes().HG2;
        met_data[ "QG"] = biol::GetAtomTypes().HG2;
        met_data[ "ME"] = biol::GetAtomTypes().HE1;
        met_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().MET] = met_data;
        // TRP
        storage::Map< std::string, biol::AtomType> trp_data;
        trp_data[ "HN"] = biol::GetAtomTypes().H;
        trp_data[ "QPB"] = biol::GetAtomTypes().HB2;
        trp_data[ "QB"] = biol::GetAtomTypes().HB2;
        trp_data[ "HNE"] = biol::GetAtomTypes().HE1;
        trp_data[ "HE1"] = biol::GetAtomTypes().HE1;
        trp_data[ "HD"] = biol::GetAtomTypes().HD1;
        trp_data[ "HE"] = biol::GetAtomTypes().HE3;
        trp_data[ "HH"] = biol::GetAtomTypes().HH2;
        trp_data[ "N1"] = biol::GetAtomTypes().HE1;
        trp_data[ "N1H"] = biol::GetAtomTypes().HE1;
        trp_data[ "H2"] = biol::GetAtomTypes().HD1;
        trp_data[ "C2"] = biol::GetAtomTypes().HD1;
        trp_data[ "C2H"] = biol::GetAtomTypes().HD1;
        trp_data[ "H4"] = biol::GetAtomTypes().HE3;
        trp_data[ "C4"] = biol::GetAtomTypes().HE3;
        trp_data[ "C4H"] = biol::GetAtomTypes().HE3;
        trp_data[ "H5"] = biol::GetAtomTypes().HZ3;
        trp_data[ "C5"] = biol::GetAtomTypes().HZ3;
        trp_data[ "C5H"] = biol::GetAtomTypes().HZ3;
        trp_data[ "H6"] = biol::GetAtomTypes().HH2;
        trp_data[ "C6"] = biol::GetAtomTypes().HH2;
        trp_data[ "C6H"] = biol::GetAtomTypes().HH2;
        trp_data[ "H7"] = biol::GetAtomTypes().HZ2;
        trp_data[ "C7"] = biol::GetAtomTypes().HZ2;
        trp_data[ "C7H"] = biol::GetAtomTypes().HZ2;
        s_map[ biol::GetAATypes().TRP] = trp_data;
        // SER
        storage::Map< std::string, biol::AtomType> ser_data;
        ser_data[ "HN"] = biol::GetAtomTypes().H;
        ser_data[ "QPB"] = biol::GetAtomTypes().HB2;
        ser_data[ "QB"] = biol::GetAtomTypes().HB2;
        ser_data[ "OG"] = biol::GetAtomTypes().HG;
        ser_data[ "HOG"] = biol::GetAtomTypes().HG;
        s_map[ biol::GetAATypes().SER] = ser_data;
        // THR
        storage::Map< std::string, biol::AtomType> thr_data;
        thr_data[ "HN"] = biol::GetAtomTypes().H;
        thr_data[ "HOG"] = biol::GetAtomTypes().HG1;
        thr_data[ "MG"] = biol::GetAtomTypes().HG21;
        thr_data[ "QG"] = biol::GetAtomTypes().HG21;
        thr_data[ "MG2"] = biol::GetAtomTypes().HG21;
        thr_data[ "QG2"] = biol::GetAtomTypes().HG21;
        s_map[ biol::GetAATypes().THR] = thr_data;
        // ASN
        storage::Map< std::string, biol::AtomType> asn_data;
        asn_data[ "HN"] = biol::GetAtomTypes().H;
        asn_data[ "QPB"] = biol::GetAtomTypes().HB2;
        asn_data[ "QB"] = biol::GetAtomTypes().HB2;
        asn_data[ "HND1"] = biol::GetAtomTypes().HD21;
        asn_data[ "HND2"] = biol::GetAtomTypes().HD22;
        asn_data[ "ND2"] = biol::GetAtomTypes().HD21;
        asn_data[ "HND"] = biol::GetAtomTypes().HD21;
        asn_data[ "QD2"] = biol::GetAtomTypes().HD21;
        s_map[ biol::GetAATypes().ASN] = asn_data;
        // GLN
        storage::Map< std::string, biol::AtomType> gln_data;
        gln_data[ "HN"] = biol::GetAtomTypes().H;
        gln_data[ "QPB"] = biol::GetAtomTypes().HB2;
        gln_data[ "QB"] = biol::GetAtomTypes().HB2;
        gln_data[ "QPG"] = biol::GetAtomTypes().HG2;
        gln_data[ "QG"] = biol::GetAtomTypes().HG2;
        gln_data[ "HNE1"] = biol::GetAtomTypes().HE21;
        gln_data[ "HNE2"] = biol::GetAtomTypes().HE22;
        gln_data[ "HE21"] = biol::GetAtomTypes().HE21;
        gln_data[ "HE22"] = biol::GetAtomTypes().HE22;
        gln_data[ "NE2"] = biol::GetAtomTypes().HE21;
        gln_data[ "HNE"] = biol::GetAtomTypes().HE21;
        gln_data[ "QE2"] = biol::GetAtomTypes().HE21;
        s_map[ biol::GetAATypes().GLN] = gln_data;
        // TYR
        storage::Map< std::string, biol::AtomType> tyr_data;
        tyr_data[ "HN"] = biol::GetAtomTypes().H;
        tyr_data[ "QPB"] = biol::GetAtomTypes().HB2;
        tyr_data[ "QB"] = biol::GetAtomTypes().HB2;
        tyr_data[ "HD"] = biol::GetAtomTypes().HD1;
        tyr_data[ "CG"] = biol::GetAtomTypes().HD1;
        tyr_data[ "HE"] = biol::GetAtomTypes().HE1;
        tyr_data[ "CZ"] = biol::GetAtomTypes().HE1;
        tyr_data[ "HOH"] = biol::GetAtomTypes().HH;
        tyr_data[ "QD"] = biol::GetAtomTypes().HD1;
        tyr_data[ "QE"] = biol::GetAtomTypes().HE1;
        s_map[ biol::GetAATypes().TYR] = tyr_data;
        // HIS
        storage::Map< std::string, biol::AtomType> his_data;
        his_data[ "HN"] = biol::GetAtomTypes().H;
        his_data[ "QPB"] = biol::GetAtomTypes().HB2;
        his_data[ "QB"] = biol::GetAtomTypes().HB2;
        his_data[ "HND"] = biol::GetAtomTypes().HD1;
        his_data[ "HE"] = biol::GetAtomTypes().HE1;
        his_data[ "HNE"] = biol::GetAtomTypes().HE2;
        his_data[ "HD"] = biol::GetAtomTypes().HD2;
        s_map[ biol::GetAATypes().HIS] = his_data;
        // ASP
        storage::Map< std::string, biol::AtomType> asp_data;
        asp_data[ "HN"] = biol::GetAtomTypes().H;
        asp_data[ "QPB"] = biol::GetAtomTypes().HB2;
        asp_data[ "QB"] = biol::GetAtomTypes().HB2;
        s_map[ biol::GetAATypes().ASP] = asp_data;
        // GLU
        storage::Map< std::string, biol::AtomType> glu_data;
        glu_data[ "HN"] = biol::GetAtomTypes().H;
        glu_data[ "QPB"] = biol::GetAtomTypes().HB2;
        glu_data[ "QB"] = biol::GetAtomTypes().HB2;
        glu_data[ "QPG"] = biol::GetAtomTypes().HG2;
        glu_data[ "QG"] = biol::GetAtomTypes().HG2;
        s_map[ biol::GetAATypes().GLU] = glu_data;
        // LYS
        storage::Map< std::string, biol::AtomType> lys_data;
        lys_data[ "HN"] = biol::GetAtomTypes().H;
        lys_data[ "QPB"] = biol::GetAtomTypes().HB2;
        lys_data[ "QB"] = biol::GetAtomTypes().HB2;
        lys_data[ "QPG"] = biol::GetAtomTypes().HG2;
        lys_data[ "QG"] = biol::GetAtomTypes().HG2;
        lys_data[ "QPD"] = biol::GetAtomTypes().HD2;
        lys_data[ "QD"] = biol::GetAtomTypes().HD2;
        lys_data[ "QPE"] = biol::GetAtomTypes().HE2;
        lys_data[ "QE"] = biol::GetAtomTypes().HE2;
        lys_data[ "NZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "HNZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "MNZ"] = biol::GetAtomTypes().HZ1;
        lys_data[ "QNZ"] = biol::GetAtomTypes().HZ1;
        s_map[ biol::GetAATypes().LYS] = lys_data;
        // ARG
        storage::Map< std::string, biol::AtomType> arg_data;
        arg_data[ "HN"] = biol::GetAtomTypes().H;
        arg_data[ "QPB"] = biol::GetAtomTypes().HB2;
        arg_data[ "QB"] = biol::GetAtomTypes().HB2;
        arg_data[ "QPG"] = biol::GetAtomTypes().HG2;
        arg_data[ "QG"] = biol::GetAtomTypes().HG2;
        arg_data[ "QPD"] = biol::GetAtomTypes().HD2;
        arg_data[ "QD"] = biol::GetAtomTypes().HD2;
        arg_data[ "CD"] = biol::GetAtomTypes().CD;
        arg_data[ "NE"] = biol::GetAtomTypes().HE;
        arg_data[ "HNE"] = biol::GetAtomTypes().HE;
        arg_data[ "HE"] = biol::GetAtomTypes().HE;
        arg_data[ "NH1"] = biol::GetAtomTypes().HH11;
        arg_data[ "NH2"] = biol::GetAtomTypes().HH21;
        s_map[ biol::GetAATypes().ARG] = arg_data;
        // CYS
        storage::Map< std::string, biol::AtomType> cys_data;
        cys_data[ "HN"] = biol::GetAtomTypes().H;
        cys_data[ "QPB"] = biol::GetAtomTypes().HB2;
        cys_data[ "QB"] = biol::GetAtomTypes().HB2;
        cys_data[ "HSG"] = biol::GetAtomTypes().HG;
        s_map[ biol::GetAATypes().CYS] = cys_data;
      }
      return s_map;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorCoordinatesHydrogen::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Locates atoms/coordinates, including hydrogens, from a protein model\n"
        "Handles locating hydrogen atoms since they are not typically included in BCL protein models but are\n"
        "often used in experimental restraints, most notably, NMR.  Classes that use this class may elect to locate\n"
        "the hydrogens (using LocateAtomCopy), otherwise, the CB position will be returned"
      );
      parameters.AddInitializer
      (
        "chain",
        "the desired chain id",
        io::Serialization::GetAgentWithRange( &m_ChainID, 'A', 'Z')
      );
      parameters.AddInitializer
      (
        "sequence",
        "the sequence id of the desired residue",
        io::Serialization::GetAgentWithMin( &m_SeqID, 0)
      );
      parameters.AddInitializer
      (
        "atom type",
        "atom type of interest, if known",
        io::Serialization::GetAgent( &m_AtomTypeString)
      );
      parameters.AddDataMember
      (
        "actual_atom_type",
        io::Serialization::GetAgent( &m_AtomType)
      );
      return parameters;
    }

    //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool LocatorCoordinatesHydrogen::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &SERIALIZER,
      std::ostream &ERR_STREAM
    )
    {
      // make string upper case
      std::transform
      (
        m_AtomTypeString.begin(),
        m_AtomTypeString.end(),
        m_AtomTypeString.begin(),
        []( unsigned char c) { return std::toupper( c);}
      );

      // try to get the real atom type
      m_AtomType = GetAtomTypeFromString();
      return true;
    }

  } // namespace restraint

} // namespace bcl
