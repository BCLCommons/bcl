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
#include "example_proteins.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{

//////////
// data //
//////////

  //! static map to store the model for each protein and for each AAClass
  storage::Map
  <
    storage::Pair< std::string, biol::AAClass>,
    assemble::ProteinModel
  > &Proteins::GetProteinsMap()
  {
    static storage::Map
    <
      storage::Pair< std::string, biol::AAClass>,
      assemble::ProteinModel
    > s_proteins;

    return s_proteins;
  }

/////////////////
// data access //
/////////////////

  //! @brief returns class name of the object behind a pointer or the current object
  //! @return the class name
  const std::string &Proteins::GetClassIdentifier() const
  {
    return GetStaticClassName( *this);
  }

////////////////
// operations //
////////////////

  //! @brief gets the protein model specified by given pdb filename and AAClass and min sse sizes
  //! @param PDB_FILENAME filename of the pdb to be read
  //! @param AA_CLASS AAClass to be used in constructing protein model
  //! @param MIN_SSE_SIZES map of requested SSTypes with their corresponding minimal sizes
  //! @return requested protein model
  assemble::ProteinModel Proteins::GetModel
  (
    const std::string &PDB_FILENAME,
    const biol::AAClass &AA_CLASS,
    const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES
  )
  {
    // retrieve the model and make a hardcopy
    util::ShPtr< assemble::ProteinModel> new_model_ptr
    (
       RetrieveModel( PDB_FILENAME, AA_CLASS).HardCopy()
    );
    assemble::ProteinModel &new_model( *new_model_ptr);

    // filter by given min_sse_sizes
    new_model.FilterByMinSSESizes( MIN_SSE_SIZES);

    // end
    return new_model;
  }

  //! @brief gets the chain from the protein model specified by given pdb filename and AAClass and min sse sizes
  //! @param PDB_FILENAME filename of the pdb to be read
  //! @param CHAIN_ID chain id
  //! @param AA_CLASS AAClass to be used in constructing protein model
  //! @param MIN_SSE_SIZES map of requested SSTypes with their corresponding minimal sizes
  //! @return requested protein model
  assemble::Chain Proteins::GetChain
  (
    const std::string &PDB_FILENAME,
    const char CHAIN_ID,
    const biol::AAClass &AA_CLASS,
    const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES
  )
  {
    // retrieve the model and make a hardcopy
    const assemble::ProteinModel &model( RetrieveModel( PDB_FILENAME, AA_CLASS));

    // look for the chain
    const util::ShPtr< assemble::Chain> sp_chain( model.GetChain( CHAIN_ID));

    // assert it is found
    BCL_Assert( sp_chain.IsDefined(), "no chain " + util::Format()( CHAIN_ID) + " found in protein " + PDB_FILENAME);

    // make a hardcopy
    util::ShPtr< assemble::Chain> new_chain( sp_chain->HardCopy());

    // filter
    new_chain->FilterByMinSSESizes( MIN_SSE_SIZES);

    // end
    return *new_chain;
  }

  //! @brief gets the SSEfrom the protein model specified by given pdb filename, chain ID, seqid first and last and AAClass
  //! @param PDB_FILENAME filename of the pdb to be read
  //! @param CHAIN_ID chain id
  //! @param SEQID_BEGIN seqid of the first amino acid in the requested SSE
  //! @param SEQID_END seqid of the last amino acid in the requested SSE
  //! @param AA_CLASS AAClass to be used in constructing protein model
  //! @param MIN_SSE_SIZES map of requested SSTypes with their corresponding minimal sizes
  //! @return requested protein model
  util::ShPtr< assemble::SSE> Proteins::GetSSE
  (
    const std::string &PDB_FILENAME,
    const char CHAIN_ID,
    const int SEQID_BEGIN,
    const int SEQID_END,
    const biol::AAClass &AA_CLASS
  )
  {
    // retrieve the model and make a hardcopy
    const assemble::ProteinModel &model( RetrieveModel( PDB_FILENAME, AA_CLASS));

    // locate the given SSE
    util::SiPtr< const assemble::SSE> sp_sse( assemble::LocatorSSE( CHAIN_ID, SEQID_BEGIN, SEQID_END).Locate( model));

    // make sure it was found
    BCL_Assert
    (
      sp_sse.IsDefined(),
      "SSE " + util::Format()( CHAIN_ID) + " " + util::Format()( SEQID_BEGIN) + "-" + util::Format()( SEQID_END) +
      " was not found in protein from file " + PDB_FILENAME
     );

    // make a hard copy
    util::ShPtr< assemble::SSE> new_sse( sp_sse->HardCopy());

    // end
    return new_sse;
  }

  //! @brief writes the provided protein model with the given filename
  //! @param PROTEIN_MODEL ProteinModel to be written
  //! @param PDB_FILENAME filename of the pdb file to be written
  void Proteins::WriteModelToPDB
  (
    const assemble::ProteinModel &PROTEIN_MODEL,
    const std::string &PDB_FILENAME
  )
  {
    // open the output stream with the filename
    io::OFStream write;
    io::File::MustOpenOFStream( write, PDB_FILENAME);

    // call the corresponding pdb::Factory function
    pdb::Factory().WriteModelToPDB( PROTEIN_MODEL, write);
    io::File::CloseClearFStream( write);
  }

  //! @brief writes the provided chain with the given filename
  //! @param CHAIN Chain to be written
  //! @param PDB_FILENAME filename of the pdb file to be written
  void Proteins::WriteChainToPDB
  (
    const assemble::Chain &CHAIN,
    const std::string &PDB_FILENAME
  )
  {
    // open the output stream with the filename
    io::OFStream write;
    io::File::MustOpenOFStream( write, PDB_FILENAME);

    // call the corresponding pdb::Factory function
    pdb::Factory::WriteChainToPDB( CHAIN, write);
    io::File::CloseClearFStream( write);
  }

  //! @brief writes the provided SSE with the given filename
  //! @param SP_SSE SSE to be written
  //! @param PDB_FILENAME filename of the pdb file to be written
  void Proteins::WriteSSEToPDB
  (
    const util::ShPtr< assemble::SSE> &SP_SSE,
    const std::string &PDB_FILENAME
  )
  {
    // construct chain and insert this SSE
    assemble::Chain chain( SP_SSE);
    chain.Insert( SP_SSE);

    // write chain
    WriteChainToPDB( chain, PDB_FILENAME);
  }

//////////////////////
// input and output //
//////////////////////

  //! @brief read from std::istream
  //! @param ISTREAM input stream
  //! @return istream which was read from
  std::istream &Proteins::Read( std::istream &ISTREAM)
  {
    // end
    return ISTREAM;
  }

  //! @brief write to std::ostream
  //! @param OSTREAM outputstream to write to
  //! @param INDENT number of indentations
  //! @return outputstream which was written to
  std::ostream &Proteins::Write( std::ostream &OSTREAM, const size_t INDENT) const
  {
    // end
    return OSTREAM;
  }

//////////////////////
// helper functions //
//////////////////////

  //! @brief returns the map with minimal SSE sizes each set to 0
  //! @return the map with minimal SSE sizes each set to 0
  storage::Map< biol::SSType, size_t> Proteins::GetZeroSSETypeMinSizes()
  {
    // construct the map
    storage::Map< biol::SSType, size_t> default_min_sse_sizes;

    // set the values
    default_min_sse_sizes[ biol::GetSSTypes().HELIX] = 0;
    default_min_sse_sizes[ biol::GetSSTypes().STRAND] = 0;
    default_min_sse_sizes[ biol::GetSSTypes().COIL] = 0;

    // end
    return default_min_sse_sizes;
  }

  //! @brief returns the map with minimal SSE sizes each set to default values
  //! @return the map with minimal SSE sizes each set to default values
  storage::Map< biol::SSType, size_t> Proteins::GetDefaultSSETypeMinSizes()
  {
    return pdb::Factory::GetCommandlineSSETypeMinSizes();
  }

  //! @brief retrieves the protein model specified by given pdb filename and AAClass from the set
  //! @param PDB_FILENAME filename of the pdb to be read
  //! @param AA_CLASS AAClass to be used in constructing protein model
  //! @return requested protein model
  const assemble::ProteinModel &Proteins::RetrieveModel
  (
    const std::string &PDB_FILENAME,
    const biol::AAClass &AA_CLASS
  )
  {
    // static min_sse_size vector
    static storage::Map< biol::SSType, size_t> s_zero_min_sse_sizes( GetZeroSSETypeMinSizes());

    // create the key
    const storage::Pair< std::string, biol::AAClass> key( PDB_FILENAME, AA_CLASS);

    // check to see if the requested protein exists
    storage::Map
    <
      storage::Pair< std::string, biol::AAClass>, assemble::ProteinModel
    >::const_iterator protein_itr( GetProteinsMap().Find( key));

    BCL_MessageVrb
    (
      "Looking for the following pdb in the map: " + PDB_FILENAME + " " + AA_CLASS.GetName()
    );

    // if not found read the protein model
    if( protein_itr == GetProteinsMap().End())
    {
      BCL_MessageVrb( "pdb not found reading it in");

      // get the protein model with the all SSEs as defined in the pdb and store it
      GetProteinsMap()[ key] =
        pdb::Factory( AA_CLASS).ProteinModelFromPDBFilename( PDB_FILENAME, s_zero_min_sse_sizes);
    }

    // end
    return GetProteinsMap()[ key];
  }

} // namespace bcl
