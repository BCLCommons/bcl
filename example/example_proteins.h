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

#ifndef EXAMPLE_PROTEINS_H_
#define EXAMPLE_PROTEINS_H_

// include forward header of this class

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @class Proteins
  //! @brief Proteins is convenience class for reading/writing proteins in example files
  //! Proteins allows convenience functions for reading and writing proteins, chains to and from pdb files. In addition,
  //! it stores a map for every protein it has read so far along with the AAClass used, so that when another example
  //! asks for the same protein, the access is very fast without the need to read the pdb again.
  //!
  //! @author karakam
  //! @date Nov 3, 2009
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class Proteins :
    public util::ObjectInterface
  {

  private:

  //////////
  // data //
  //////////

    //! static map to store the model for each protein and for each AAClass
    static
    storage::Map
    <
      storage::Pair< std::string, biol::AAClass>,
      assemble::ProteinModel
    >
    &GetProteinsMap();

  public:

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GetClassIdentifier() const;

  ////////////////
  // operations //
  ////////////////

    //! @brief gets the protein model specified by given pdb filename, chainid and AAClass and min sse sizes
    //! @param PDB_FILENAME filename of the pdb to be read
    //! @param AA_CLASS AAClass to be used in constructing protein model
    //! @param MIN_SSE_SIZES map of requested SSTypes with their corresponding minimal sizes
    //! @return requested protein model
    static assemble::ProteinModel GetModel
    (
      const std::string &PDB_FILENAME,
      const biol::AAClass &AA_CLASS = biol::GetAAClasses().e_AABackBone,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES = GetDefaultSSETypeMinSizes()
    );

    //! @brief gets the chain from the protein model specified by given pdb filename and AAClass and min sse sizes
    //! @param PDB_FILENAME filename of the pdb to be read
    //! @param CHAIN_ID chain id
    //! @param AA_CLASS AAClass to be used in constructing protein model
    //! @param MIN_SSE_SIZES map of requested SSTypes with their corresponding minimal sizes
    //! @return requested protein model
    static assemble::Chain GetChain
    (
      const std::string &PDB_FILENAME,
      const char CHAIN_ID,
      const biol::AAClass &AA_CLASS = biol::GetAAClasses().e_AABackBone,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES = GetDefaultSSETypeMinSizes()
    );

    //! @brief gets the SSEfrom the protein model specified by given pdb filename, chain ID, seqid first and last and AAClass
    //! @param PDB_FILENAME filename of the pdb to be read
    //! @param CHAIN_ID chain id
    //! @param SEQID_BEGIN seqid of the first amino acid in the requested SSE
    //! @param SEQID_END seqid of the last amino acid in the requested SSE
    //! @param AA_CLASS AAClass to be used in constructing protein model
    //! @return requested protein model
    static util::ShPtr< assemble::SSE> GetSSE
    (
      const std::string &PDB_FILENAME,
      const char CHAIN_ID,
      const int SEQID_BEGIN,
      const int SEQID_END,
      const biol::AAClass &AA_CLASS = biol::GetAAClasses().e_AABackBone
    );

    //! @brief writes the provided protein model with the given filename
    //! @param PROTEIN_MODEL ProteinModel to be written
    //! @param PDB_FILENAME filename of the pdb file to be written
    static void WriteModelToPDB
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PDB_FILENAME
    );

    //! @brief writes the provided chain with the given filename
    //! @param CHAIN Chain to be written
    //! @param PDB_FILENAME filename of the pdb file to be written
    static void WriteChainToPDB
    (
      const assemble::Chain &CHAIN,
      const std::string &PDB_FILENAME
    );

    //! @brief writes the provided SSE with the given filename
    //! @param SP_SSE SSE to be written
    //! @param PDB_FILENAME filename of the pdb file to be written
    static void WriteSSEToPDB
    (
      const util::ShPtr< assemble::SSE> &SP_SSE,
      const std::string &PDB_FILENAME
    );

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

  public:

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the map with minimal SSE sizes each set to 0
    //! @return the map with minimal SSE sizes each set to 0
    static storage::Map< biol::SSType, size_t> GetZeroSSETypeMinSizes();

    //! @brief returns the map with minimal SSE sizes each set to default values
    //! @return the map with minimal SSE sizes each set to default values
    static storage::Map< biol::SSType, size_t> GetDefaultSSETypeMinSizes();

  private:

    //! @brief retrieves the protein model specified by given pdb filename and AAClass from the set
    //! @param PDB_FILENAME filename of the pdb to be read
    //! @param AA_CLASS AAClass to be used in constructing protein model
    //! @return requested protein model
    static const assemble::ProteinModel &RetrieveModel
    (
      const std::string &PDB_FILENAME,
      const biol::AAClass &AA_CLASS
    );

  }; // class Proteins

} // namespace bcl

#endif // EXAMPLE_PROTEINS_H_
