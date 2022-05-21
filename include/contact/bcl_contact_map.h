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

#ifndef BCL_CONTACT_MAP_H_
#define BCL_CONTACT_MAP_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_types.h"
#include "biol/bcl_biol_aa_data.h"
#include "storage/bcl_storage_object_nd_hash_map.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Map
    //! @brief This class, derived from ObjectNDHashMap, is designed to store whether two AA's are in contact.
    //! @details  The contact information for given amino acid sequences are calculated using tertiary information and
    //! stored in a hash-map structure where it can be retrieved by passing the amino acid pair of interest.
    //!
    //! @see @link example_contact_map.cpp @endlink
    //! @author karakam, woetzen
    //! @date 15.08.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Map :
      public storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >
    {

    private:

    //////////
    // data //
    //////////

      static const std::string            s_Identifier; //!< identifier string to be used for contact maps
      util::ShPtrVector< assemble::Chain> m_Chains;     //!< ShPtrVector of chains
      int                                 m_Boundary;   //!< number of residues from both ends of sequence to avoid

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Map();

      //! @brief construct contact map from a chain
      //! @param CHAIN Chain for which sequence contacts are going to be detected
      //! @param BOUNDARY number of residues to be excluded from both sides
      Map
      (
        const util::ShPtr< assemble::Chain> &CHAIN,
        const int BOUNDARY = 8
      );

      //! @brief construct contact map from a protein model( intra-sequence contacts and inter-sequence contacts),
      //! @param PROTEIN_MODEL Protein Model for which inter and intra sequence contacts are going to be deteced
      //! @param BOUNDARY number of residues to be excluded from both sides
      Map
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const int BOUNDARY = 8
      );

      //! @brief virtual copy constructor
      Map *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns stored chains
      //! @return stored chains
      const util::ShPtrVector< assemble::Chain> &GetChains() const;

      //! @brief iterates over the stored chains and finds the one with matching CHAIN_ID
      //! @param CHAIN_ID chain id of the chain that is beings searched
      //! @return Chain  with the searched chain id
      const util::ShPtr< assemble::Chain> &GetChain( const char CHAIN_ID) const;

      //! @brief returns boundary
      //! @return boundary value
      size_t GetBoundary() const
      {
        return m_Boundary;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Returns contacts vector for the provided AA Data pointers
      //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
      //! @return the contacts vector and the contact boolean
      const storage::Pair< linal::Vector< size_t>, bool> &GetContactVector
      (
        const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief helper function to read predictions from a file
      //! @param ISTREAM input stream to be read from
      //! @return input stream from which the map was read
      std::istream &ReadMap( std::istream &ISTREAM);

      //! @brief helper function to write predictions to a file
      //! @param OSTREAM output stream to write to
      //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
      //! @return output stream which was written to
      std::ostream &WriteMap
      (
        std::ostream &OSTREAM,
        const bool WRITE_ONLY_CONTACTS
      ) const;

    private:

      //! @brief reads contacts for the provided chains
      //! @param ISTREAM input stream
      //! @param CHAIN_A first chain
      //! @param CHAIN_B second chain
      //! @return istream from which was read
      std::istream &ReadContacts
      (
        std::istream &ISTREAM,
        const assemble::Chain &CHAIN_A,
        const assemble::Chain &CHAIN_B
      );

      //! @brief Writes predictions for the provided chains
      //! @param OSTREAM output stream
      //! @param CHAIN_A first chain
      //! @param CHAIN_B second chain
      //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
      //! @return ostream which was written to
      std::ostream &WriteContacts
      (
        std::ostream &OSTREAM,
        const assemble::Chain &CHAIN_A,
        const assemble::Chain &CHAIN_B,
        const bool WRITE_ONLY_CONTACTS
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Given 2 AAs, this function returns a bool of whether these two residues
      //! @brief are in contact according to rules of the specified contact type
      //! @param AMINO_ACID_A first amino acid
      //! @param AMINO_ACID_B second amino acid
      //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with any contact type
      static bool IsInContact( const biol::AABase &AMINO_ACID_A, const biol::AABase &AMINO_ACID_B);

      //! @brief Given 2 AAs and a contact type, this function returns a bool of whether these two residues
      //! @brief are in contact according to rules of the specified contact type
      //! @param AMINO_ACID_A first amino acid
      //! @param AMINO_ACID_B second amino acid
      //! @param CONTACT_TYPE contact type
      //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with CONTACT_TYPE
      static bool IsInContact
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const Type &CONTACT_TYPE
      );

      //! @brief checks if the two SSEs are in contact in the given contact map
      //! @param CONTACTS the contact map the defines the contacts between amino acids
      //! @param FIRST first SSE
      //! @param SECOND second SSE
      //! @return true if the SSEs are in contact
      static bool IsInContact( const Map &CONTACTS, const assemble::SSE &FIRST, const assemble::SSE &SECOND);

    private:

      //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
      //! @brief each contact type for the stored chains
      void FillMap();

      //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
      //! @brief each contact type for the given chain
      //! @param CHAIN first chain
      void FillMap( const assemble::Chain &CHAIN);

      //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
      //! @brief each contact type for the given pair of chains
      //! @param CHAIN_A first chain
      //! @param CHAIN_B second chain
      void FillMap( const assemble::Chain &CHAIN_A, const assemble::Chain &CHAIN_B);

      //! @brief checks the boundary conditions of a residue to determine whether it should be included in the map
      //! @param AMINO_ACID amino acid for which the conditions are going to be checked
      //! @param CHAIN a chain
      bool CheckBoundaryCondition( const biol::AABase &AMINO_ACID, const assemble::Chain &CHAIN) const;

      //! @brief insert the provided amino acid pair and contact type and optional reverse contact type,
      //! @brief inserts the related data into hashmap and residue pair list
      //! @param AMINO_ACID_A first amino acid
      //! @param AMINO_ACID_B second amino acid
      //! @param CONTACT_TYPE contact type
      //! @param INSERT_REVERSE also insert reverse of CONTACT_TYPE
      void InsertAminoAcidPair
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const Type &CONTACT_TYPE,
        const bool INSERT_REVERSE
      );

      //! @brief this functions checks that the chains stored also have SSE information
      //! @brief so that the contact map can be generated correctly
      void CheckChainsForSSEs() const;

    }; //class Map

  } // namespace contact
} // namespace bcl

#endif //BCL_CONTACT_MAP_H_

