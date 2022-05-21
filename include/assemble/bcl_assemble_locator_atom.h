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

#ifndef BCL_ASSEMBLE_LOCATOR_ATOM_H_
#define BCL_ASSEMBLE_LOCATOR_ATOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_locator_aa.h"
#include "bcl_assemble_locator_atom_coordinates_interface.h"
#include "biol/bcl_biol_atom_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorAtom
    //! @brief This class is used for locating a specified atom from a protein model
    //! @details This class uses the member LocatorAA to find the corresponding atom and then finds the corresponding
    //! atom from that amino acid with the specified atom type
    //!
    //! @see @link example_assemble_locator_atom.cpp @endlink
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorAtom :
      public LocatorAtomCoordinatesInterface
    {

    private:

    //////////
    // data //
    //////////

      LocatorAA      m_LocatorAA; //!< AA Locator
      biol::AtomType m_Atom_Type; //!< atom to be located

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
      LocatorAtom() :
        m_LocatorAA( 'A', 1),
        m_Atom_Type( biol::GetAtomTypes().CA)
      {
      }

      //! @brief constructor from chain ID, amino acid ID, and atom
      //! @param CHAIN_ID char which indicates the chain
      //! @param AA_ID size_t which indicates the SeqID of the amino acid
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorAtom( const char CHAIN_ID, const int AA_ID, const bool USE_PDB_ID = false) :
        m_LocatorAA( CHAIN_ID, AA_ID, USE_PDB_ID),
        m_Atom_Type( biol::GetAtomTypes().e_Undefined)
      {
      }

      //! @brief constructor from chain ID, amino acid ID, and atom
      //! @param CHAIN_ID char which indicates the chain
      //! @param AA_ID size_t which indicates the SeqID of the amino acid
      //! @param ATOM_TYPE AtomName which indicates the atom
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorAtom( const char CHAIN_ID, const int AA_ID, const biol::AtomType &ATOM_TYPE, const bool USE_PDB_ID = false) :
        m_LocatorAA( CHAIN_ID, AA_ID, USE_PDB_ID),
        m_Atom_Type( ATOM_TYPE)
      {
      }

      //! @brief constructor from chain ID, amino acid ID, and atom
      //! @param CHAIN_ID char which indicates the chain
      //! @param AA_ID size_t which indicates the SeqID of the amino acid
      //! @param ATOM_TYPE AtomName which indicates the atom
      //! @param AA_TYPE the amino acid type; used for verification purposes
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorAtom
      (
        const char CHAIN_ID,
        const int AA_ID,
        const biol::AtomType &ATOM_TYPE,
        const biol::AAType &AA_TYPE,
        const bool USE_PDB_ID = false
      ) :
        m_LocatorAA( CHAIN_ID, AA_ID, AA_TYPE, USE_PDB_ID),
        m_Atom_Type( ATOM_TYPE)
      {
      }

      //! @brief clone constructor
      LocatorAtom *Clone() const
      {
        return new LocatorAtom( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives formatted string describing the locator
      //! @return formatted string describing the locator
      std::string GetIdentification() const;

      //! @brief reads formatted string describing the locator
      //! @return formatted string describing the locator
      std::istream &ReadIdentification( std::istream &ISTREAM);

      //! @brief gives the chain id the locator corresponds to
      //! @return the chain id the locator corresponds to
      char GetChainID() const;

      //! @brief gives the seq id the locator corresponds to
      //! @return the seq id the locator corresponds to
      int GetSeqID() const;

      //! @brief gives the atom type the locator corresponds to
      //! @return the atom type the locator corresponds to
      const biol::AtomType &GetAtomType() const;

      //! @brief gives the aa type the locator corresponds to
      //! @return the aa type the locator corresponds to
      const biol::AAType &GetAAType() const;

      //! @brief SetAtomID changes the atom name to be located
      //! @param ATOM_TYPE AtomType which indicates the new atom to be located
      void SetAtomType( const biol::AtomType &ATOM_TYPE);

      //! @brief returns the AA locator
      //! @return returns the const reference to the AA locator
      const LocatorAA &GetLocatorAA() const;

      //! @brief returns the AA locator
      //! @return returns the reference to the AA locator
      void SetLocatorAA( const LocatorAA &LOCATOR_AA);

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate translates the LocatorAtom denoting an atom into the actual atom
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
      //! @param CHAIN_ID the chain id the atom is in
      //! @param SEQ_ID the seq id of the residue the atom is in
      //! @param ATOM_TYPE the atom type that should be located
      //! @return the atom denoted by the LocatorAtom
      static util::SiPtr< const biol::Atom> LocateAtomFromModel
      (
        const ProteinModel &PROTEIN_MODEL, const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE
      );

      //! @brief locates an atom from a protein model
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return siptr to atom that has been located from PROTEIN_MODEL
      util::SiPtr< const biol::Atom>
      LocateAtom( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates the desired coordinates from a protein model
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
      //! @return the coordinates of the atom denoted by the LocatorAtom
      linal::Vector3D Locate( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class LocatorAtom

    //! @brief less than operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom which will be compared against the second LocatorAtom
    //! @param RHS the second LocatorAtom which will be compared against the first LocatorAtom
    //! @return boolean true if LHS is less than RHS - false otherwise
    bool operator<( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B);

    //! @brief not equal operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom
    //! @param RHS the second LocatorAtom
    //! @return bool true if LHS not equal to RHS - false otherwise
    bool operator!=( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B);

    //! @brief equal operator for comparing two LocatorAtom
    //! @param LHS the first LocatorAtom
    //! @param RHS the second LocatorAtom
    //! @return bool true if LHS equal to RHS - false otherwise
    bool operator==( const LocatorAtom &LOCATOR_A, const LocatorAtom &LOCATOR_B);

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_ATOM_H_
