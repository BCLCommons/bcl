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
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesFirstSideChainAtom::s_Instance
    (
      util::Enumerated< find::LocatorCoordinatesInterface< assemble::ProteinModel> >::AddInstance
      (
        util::Enumerated< assemble::LocatorAtomCoordinatesInterface>::AddInstance
        (
          new LocatorCoordinatesFirstSideChainAtom()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom() :
      m_LocatorAA(),
      m_AtomType()
    {
    }

    //! @brief constructor from chain ID, amino acid ID
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom( const char CHAIN_ID, const int SEQ_ID, bool USE_PDB_ID) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, USE_PDB_ID),
      m_AtomType( biol::GetAtomTypes().GetUndefined())
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and residue type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param AA_TYPE AtomName which indicates the residue type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AAType &AA_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, AA_TYPE, USE_PDB_ID),
      m_AtomType( biol::GetAtomTypes().GetUndefined())
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and atom type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param ATOM_TYPE indicates the atom type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, USE_PDB_ID),
      m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief constructor from chain ID, amino acid ID, and atom type and residue type
    //! @param CHAIN_ID char which indicates the chain
    //! @param SEQ_ID int which indicates the SeqID of the amino acid
    //! @param ATOM_TYPE AtomName which indicates the atom type
    //! @param AA_TYPE AtomName which indicates the residue type
    //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
    LocatorCoordinatesFirstSideChainAtom::LocatorCoordinatesFirstSideChainAtom
    (
      const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, const biol::AAType &AA_TYPE, bool USE_PDB_ID
    ) :
      m_LocatorAA( CHAIN_ID, SEQ_ID, AA_TYPE, USE_PDB_ID),
      m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorCoordinatesFirstSideChainAtom
    LocatorCoordinatesFirstSideChainAtom *LocatorCoordinatesFirstSideChainAtom::Clone() const
    {
      return new LocatorCoordinatesFirstSideChainAtom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorCoordinatesFirstSideChainAtom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorCoordinatesFirstSideChainAtom::GetAlias() const
    {
      static std::string s_alias( "FirstSideChainAtom");
      return s_alias;
    }

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorCoordinatesFirstSideChainAtom::GetIdentification() const
    {
      return m_LocatorAA.GetIdentification();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorCoordinatesFirstSideChainAtom::ReadIdentification( std::istream &ISTREAM)
    {
      m_LocatorAA.ReadIdentification( ISTREAM);
      return ISTREAM;
    }

    //! @brief gives the chain id the locator corresponds to
    //! @return the chain id the locator corresponds to
    char LocatorCoordinatesFirstSideChainAtom::GetChainID() const
    {
      return m_LocatorAA.GetLocatorChain().GetChainID();
    }

    //! @brief gives the seq id the locator corresponds to
    //! @return the seq id the locator corresponds to
    int LocatorCoordinatesFirstSideChainAtom::GetSeqID() const
    {
      return m_LocatorAA.GetAAID();
    }

    //! @brief gives the atom type the locator corresponds to
    //! @return the atom type the locator corresponds to
    const biol::AtomType &LocatorCoordinatesFirstSideChainAtom::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief gives the aa type the locator corresponds to
    //! @return the aa type the locator corresponds to
    const biol::AAType &LocatorCoordinatesFirstSideChainAtom::GetAAType() const
    {
      return m_LocatorAA.GetAAType();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return siptr to atom that has been located from PROTEIN_MODEL
    util::SiPtr< const biol::Atom>
    LocatorCoordinatesFirstSideChainAtom::LocateAtom( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate the aa containing the atom
      const util::SiPtr< const biol::AABase> aa( m_LocatorAA.Locate( PROTEIN_MODEL));

      // true if the aa could not be located
      if( !aa.IsDefined())
      {
        // return undefined vector3d
        return util::SiPtr< const biol::Atom>();
      }

      // get the first side chain atom of the aa
      const util::SiPtr< const biol::Atom> atom( aa->GetFirstSidechainAtom());

      // return the coordinates
      return atom;
    }

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return copy of atom that has been located from PROTEIN_MODEL
    biol::Atom LocatorCoordinatesFirstSideChainAtom::LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // try to locate atom
      const util::SiPtr< const biol::Atom> atom( LocateAtom( PROTEIN_MODEL));

      // true if atom could not be found
      if( !atom.IsDefined())
      {
        // return empty atom
        return biol::Atom();
      }

      // atom is defined so return atom behind pointer
      return *atom;
    }

    //! @brief locates the desired coordinates from a protein model
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @return the coordinates of the atom denoted by the LocatorAtom
    linal::Vector3D LocatorCoordinatesFirstSideChainAtom::Locate( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate the aa containing the atom
      const util::SiPtr< const biol::AABase> aa( m_LocatorAA.Locate( PROTEIN_MODEL));

      // true if the aa could not be located
      if( !aa.IsDefined())
      {
        // return undefined vector3d
        return linal::Vector3D( util::GetUndefinedDouble());
      }

      // get the coordinates from the first side chain atom of the aa
      const linal::Vector3D coords( aa->GetFirstSidechainAtom().GetCoordinates());

      // return the coordinates
      return coords;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorCoordinatesFirstSideChainAtom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LocatorAA, ISTREAM);
      io::Serialize::Read( m_AtomType,  ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorCoordinatesFirstSideChainAtom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LocatorAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomType,  OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorCoordinatesFirstSideChainAtom::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "locates the coordinates of the first side chain atom of the indicated residue");
      parameters.AddInitializer
      (
        "locator",
        "locator for the residue of interest",
        io::Serialization::GetAgent( &m_LocatorAA)
      );
      parameters.AddOptionalInitializer
      (
        "atom_type",
        "atom type if the first side chain atom atom-type is known",
        io::Serialization::GetAgent( &m_AtomType)
      );
      return parameters;
    }

  } // namespace restraint

} // namespace bcl
