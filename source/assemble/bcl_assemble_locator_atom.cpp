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
#include "assemble/bcl_assemble_locator_atom.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorAtom::s_Instance
    (
      util::Enumerated< find::LocatorCoordinatesInterface< ProteinModel> >::AddInstance
      (
        util::Enumerated< LocatorAtomCoordinatesInterface>::AddInstance
        (
          new LocatorAtom()
        )
      )
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorAtom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorAtom::GetIdentification() const
    {
      std::stringstream write;
      write << m_LocatorAA.GetIdentification();
      io::Serialize::Write( GetAtomType().GetName(), write, 1);
      return write.str();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &LocatorAtom::ReadIdentification( std::istream &ISTREAM)
    {
      m_LocatorAA.ReadIdentification( ISTREAM);
      std::string enum_name;
      io::Serialize::Read( enum_name, ISTREAM);
      m_Atom_Type = biol::GetAtomTypes().GetEnumFromName( enum_name);
      return ISTREAM;
    }

    //! @brief gives the chain id the locator corresponds to
    //! @return the chain id the locator corresponds to
    char LocatorAtom::GetChainID() const
    {
      return m_LocatorAA.GetLocatorChain().GetChainID();
    }

    //! @brief gives the seq id the locator corresponds to
    //! @return the seq id the locator corresponds to
    int LocatorAtom::GetSeqID() const
    {
      return m_LocatorAA.GetAAID();
    }

    //! GetAtom gives the atom name
    //! @return returns "m_Atom_Type"
    const biol::AtomType &LocatorAtom::GetAtomType() const
    {
      return m_Atom_Type;
    }

    //! @brief gives the aa type the locator corresponds to
    //! @return the aa type the locator corresponds to
    const biol::AAType &LocatorAtom::GetAAType() const
    {
      return biol::GetAATypes().e_Undefined;
    }

    //! SetAtomID changes the atom name to be located
    //! @param ATOM_TYPE AtomType which indicates the new atom to be located
    void LocatorAtom::SetAtomType( const biol::AtomType &ATOM_TYPE)
    {
      // set "m_Atom_Type" to "ATOM_TYPE"
      m_Atom_Type = ATOM_TYPE;
    }

    //! returns the AA locator
    //! @return returns the const reference to the AA locator
    const LocatorAA &LocatorAtom::GetLocatorAA() const
    {
      return m_LocatorAA;
    }

    //! @brief returns the AA locator
    //! @return returns the reference to the AA locator
    void LocatorAtom::SetLocatorAA( const LocatorAA &LOCATOR_AA)
    {
      m_LocatorAA = LOCATOR_AA;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorAtom::GetAlias() const
    {
      static const std::string s_Name( "LocatorAtom");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate translates the LocatorAtom denoting an atom into the actual atom
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @param CHAIN_ID the chain id the atom is in
    //! @param SEQ_ID the seq id of the residue the atom is in
    //! @param ATOM_TYPE the atom type that should be located
    //! @return the atom denoted by the LocatorAtom
    util::SiPtr< const biol::Atom> LocatorAtom::LocateAtomFromModel
    (
      const ProteinModel &PROTEIN_MODEL, const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE
    )
    {
      util::SiPtr< const biol::AABase> amino_acid( LocatorAA( CHAIN_ID, SEQ_ID).Locate( PROTEIN_MODEL));

      if( !amino_acid.IsDefined())
      {
        BCL_MessageDbg
        (
          "Locator atom: Chain " + util::Format()( CHAIN_ID) + " and amino acid with seqid " + util::Format()( SEQ_ID)
          + " does not exist in protein model"
        );

        return util::SiPtr< const biol::Atom>();
      }

      util::SiPtr< const biol::Atom> atom( amino_acid->GetAtom( ATOM_TYPE));
      if( !atom.IsDefined())
      {
        BCL_MessageDbg
        (
          "atom " + util::Format()( ATOM_TYPE) + " in chain " + util::Format()( CHAIN_ID)
          + " and amino acid with seqid " + util::Format()( SEQ_ID)
          + "does not exist in protein model"
        );
      }

      return atom;
    }

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return siptr to atom that has been located from PROTEIN_MODEL
    util::SiPtr< const biol::Atom>
    LocatorAtom::LocateAtom( const ProteinModel &PROTEIN_MODEL) const
    {
      util::SiPtr< const biol::AABase> amino_acid( m_LocatorAA.Locate( PROTEIN_MODEL));

      if( !amino_acid.IsDefined())
      {
        BCL_MessageDbg
        (
          "Locator atom: Chain " + util::Format()( m_LocatorAA.GetLocatorChain().GetChainID()) +
          " and amino acid with seqid " + util::Format()( m_LocatorAA.GetAAID())
          + " does not exist in protein model"
        );

        return util::SiPtr< const biol::Atom>();
      }

      util::SiPtr< const biol::Atom> atom( amino_acid->GetAtom( m_Atom_Type));
      if( !atom.IsDefined())
      {
        BCL_MessageDbg
        (
          "atom " + util::Format()( m_Atom_Type) + " in chain " +
          util::Format()( m_LocatorAA.GetLocatorChain().GetChainID())
          + " and amino acid with seqid " + util::Format()( m_LocatorAA.GetAAID())
          + "does not exist in protein model"
        );
      }

      return atom;
    }

    //! @brief locates the desired coordinates from a protein model
    //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
    //! @return the coordinates of the atom denoted by the LocatorAtom
    linal::Vector3D LocatorAtom::Locate( const ProteinModel &PROTEIN_MODEL) const
    {
      // try to locate the atom
      util::SiPtr< const biol::Atom> atom( LocateAtom( PROTEIN_MODEL));

      // true if the ptr to atom is not defined
      if( !atom.IsDefined())
      {
        // return undefined vector3d
        return linal::Vector3D( util::GetUndefinedDouble());
      }

      // return the coordinates of the atom
      return atom->GetCoordinates();
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorAtom::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "This class is used for locating a specified atom from a protein model. This class uses the member "
        "LocatorAA to find the corresponding atom and then finds the corresponding atom from that amino acid with the "
        "specified atom type."
      );

      parameters.AddInitializer
      (
        "locator_aa",
        "the residue locator that should be used to get the residue with the desired atom",
        io::Serialization::GetAgent( &m_LocatorAA)
      );

      parameters.AddInitializer
      (
        "atom_type",
        "the atom type to locate",
        io::Serialization::GetAgent( &m_Atom_Type)
      );

      return parameters;
    }

    //! @brief less than operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom which will be compared against the second LocatorAtom
    //! @param RHS the second LocatorAtom which will be compared against the first LocatorAtom
    //! @return boolean true if LHS is less than RHS - false otherwise
    bool operator<( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B)
    {
      // chain inequality
      if( LOCATOR_A.GetLocatorAA().GetLocatorChain().GetChainID() < LOCATOR_B.GetLocatorAA().GetLocatorChain().GetChainID())
      {
        return true;
      }
      if( LOCATOR_A.GetLocatorAA().GetLocatorChain().GetChainID() > LOCATOR_B.GetLocatorAA().GetLocatorChain().GetChainID())
      {
        return false;
      }

      // seq id inequality
      if( LOCATOR_A.GetLocatorAA().GetAAID() < LOCATOR_B.GetLocatorAA().GetAAID())
      {
        return true;
      }
      if( LOCATOR_A.GetLocatorAA().GetAAID() > LOCATOR_B.GetLocatorAA().GetAAID())
      {
        return false;
      }

      return LOCATOR_A.GetAtomType() < LOCATOR_B.GetAtomType();
    }

    //! @brief not equal operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom
    //! @param RHS the second LocatorAtom
    //! @return bool true if LHS not equal to RHS - false otherwise
    bool operator!=( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B)
    {
      return LOCATOR_A.GetLocatorAA().GetLocatorChain().GetChainID() != LOCATOR_B.GetLocatorAA().GetLocatorChain().GetChainID() ||
             LOCATOR_A.GetLocatorAA().GetAAID()                      != LOCATOR_B.GetLocatorAA().GetAAID()                      ||
             LOCATOR_A.GetAtomType()                                 != LOCATOR_B.GetAtomType();
    }

    //! @brief equal operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom
    //! @param RHS the second LocatorAtom
    //! @return bool true if LHS equal to RHS - false otherwise
    bool operator==( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B)
    {
      return !( LOCATOR_A != LOCATOR_B);
    }

  } // namespace assemble
} // namespace bcl
