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

#ifndef BCL_RESTRAINT_LOCATOR_COORDINATES_HYDROGEN_H_
#define BCL_RESTRAINT_LOCATOR_COORDINATES_HYDROGEN_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesHydrogen
    //! @brief Locates atoms/coordinates, including hydrogens, from a protein model
    //! @details Handles locating hydrogen atoms since they are not typically included in BCL protein models but are
    //!          often used in experimental restraints, most notably, NMR.  The LocateAtom function will build H or HA
    //!          atoms correctly.  All other protons will be moved to the CB position.
    //!
    //! @see @link example_restraint_locator_coordinates_hydrogen.cpp @endlink
    //! @author weinerbe
    //! @date May 31, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorCoordinatesHydrogen :
      public assemble::LocatorAtomCoordinatesInterface
    {

    private:

    //////////
    // data //
    //////////

      //! chain id
      char m_ChainID;

      //! sequence id
      int m_SeqID;

      //! atom type string
      std::string m_AtomTypeString;

      //! atom type (best guess from the string w/o knowing the sequence)
      biol::AtomType m_AtomType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorCoordinatesHydrogen();

      //! @brief construct from chain id, sequence id, and atom type
      //! @param CHAIN_ID chain id
      //! @param SEQ_ID sequence id
      //! @param ATOM_TYPE atom type string
      LocatorCoordinatesHydrogen( const char CHAIN_ID, const int SEQ_ID, const std::string &ATOM_TYPE);

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesHydrogen
      LocatorCoordinatesHydrogen *Clone() const;

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
      char GetChainID() const
      {
        return m_ChainID;
      }

      //! @brief gives the seq id the locator corresponds to
      //! @return the seq id the locator corresponds to
      int GetSeqID() const
      {
        return m_SeqID;
      }

      //! @brief gets the atom type string
      //! @return the atom type string
      const std::string &GetAtomTypeString() const
      {
        return m_AtomTypeString;
      }

      //! @brief gives the atom type the locator corresponds to
      //!        but not used for this class
      //! @return the atom type the locator corresponds to
      const biol::AtomType &GetAtomType() const
      {
        return m_AtomType;
      }

      //! @brief returns reference to undefined AA TYPE, this function is required by the interface but not used for
      //!        this class
      //! @return reference to undefined AA TYPE
      const biol::AAType &GetAAType() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief locates the desired coordinates from a protein model
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAtom refers to
      //! @return the coordinates of the atom denoted by the LocatorAtom
      linal::Vector3D Locate( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates an atom from a protein model, creating H or HA atoms as needed - other hydrogen atoms will have
      //!        m_AtomType as an atom type but with coords of the first side chain atom
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return copy of atom that has been located from PROTEIN_MODEL
      biol::Atom LocateAtomCopy( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief takes a string of the Atom type as well as the amino acid type and converts them to the
      //!        correct IUPAC type
      //! @param AA_TYPE the AAType of the atom from the restraint file
      //! @return the IUPAC formatted atom type
      biol::AtomType GetAtomTypeFromString( const biol::AAType &AA_TYPE = biol::GetAATypes().e_Undefined) const;

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

      //! @brief tries replacing the last character in an atom name string with 1, 2, or 3 and getting the atom type
      //! @param ATOM_NAME atom name string
      //! @param AA_TYPE the AAType of the atom from the restraint file
      //! @return atom type
      static biol::AtomType ReplacePseudoAtomChar( const std::string &ATOM_NAME, const biol::AAType &AA_TYPE);

      //! @brief gets a map of PsuedoAtoms that can be mapped back to an IUPAC type
      //! @return the map of these PsuedoAtoms
      static const storage::Map< biol::AAType, storage::Map< std::string, biol::AtomType> > &GetPseudoAtomMap();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

    }; // class LocatorCoordinatesHydrogen

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_LOCATOR_COORDINATES_HYDROGEN_H_
