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

#ifndef BCL_RESTRAINT_LOCATOR_COORDINATES_FIRST_SIDE_CHAIN_ATOM_H_
#define BCL_RESTRAINT_LOCATOR_COORDINATES_FIRST_SIDE_CHAIN_ATOM_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesFirstSideChainAtom
    //! @brief locates the coordinates of the first side chain atom of the indicated residue
    //! @details locates the first side chain atom of the indicated residue
    //!
    //! @see @link example_restraint_locator_coordinates_first_side_chain_atom.cpp @endlink
    //! @author alexanns
    //! @date May 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorCoordinatesFirstSideChainAtom :
      public assemble::LocatorAtomCoordinatesInterface
    {

    private:

    //////////
    // data //
    //////////

      //! locator to find the residue where the coordinates will be in
      assemble::LocatorAA m_LocatorAA;

      //! optional atom type if the first side chain atom atom-type is known
      biol::AtomType m_AtomType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorCoordinatesFirstSideChainAtom();

      //! @brief constructor from chain ID, amino acid ID
      //! @param CHAIN_ID char which indicates the chain
      //! @param SEQ_ID int which indicates the SeqID of the amino acid
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorCoordinatesFirstSideChainAtom( const char CHAIN_ID, const int SEQ_ID, const bool USE_PDB_ID = false);

      //! @brief constructor from chain ID, amino acid ID, and residue type
      //! @param CHAIN_ID char which indicates the chain
      //! @param SEQ_ID int which indicates the SeqID of the amino acid
      //! @param AA_TYPE AtomName which indicates the residue type
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorCoordinatesFirstSideChainAtom( const char CHAIN_ID, const int SEQ_ID, const biol::AAType &AA_TYPE, const bool USE_PDB_ID = false);

      //! @brief constructor from chain ID, amino acid ID, and atom type
      //! @param CHAIN_ID char which indicates the chain
      //! @param SEQ_ID int which indicates the SeqID of the amino acid
      //! @param ATOM_TYPE indicates the atom type
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorCoordinatesFirstSideChainAtom( const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, const bool USE_PDB_ID = false);

      //! @brief constructor from chain ID, amino acid ID, and atom type and residue type
      //! @param CHAIN_ID char which indicates the chain
      //! @param SEQ_ID int which indicates the SeqID of the amino acid
      //! @param ATOM_TYPE AtomName which indicates the atom type
      //! @param AA_TYPE AtomName which indicates the residue type
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorCoordinatesFirstSideChainAtom
      (
        const char CHAIN_ID, const int SEQ_ID, const biol::AtomType &ATOM_TYPE, const biol::AAType &AA_TYPE, const bool USE_PDB_ID = false
      );

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesFirstSideChainAtom
      LocatorCoordinatesFirstSideChainAtom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

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

    ////////////////
    // operations //
    ////////////////

      //! @brief locates an atom from a protein model
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return siptr to atom that has been located from PROTEIN_MODEL
      virtual util::SiPtr< const biol::Atom> LocateAtom( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates an atom from a protein model
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return copy of atom that has been located from PROTEIN_MODEL
      virtual biol::Atom LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates the desired coordinates from a protein model
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
      //! @return the coordinates of the atom denoted by the LocatorAtom
      linal::Vector3D Locate( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class LocatorCoordinatesFirstSideChainAtom

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_LOCATOR_COORDINATES_FIRST_SIDE_CHAIN_ATOM_H_
