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

#ifndef BCL_PDB_FACTORY_H_
#define BCL_PDB_FACTORY_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_atom_types.h"
#include "command/bcl_command_flag_interface.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Factory
    //! @brief The Factory creates objects like AASequence, Model, etc, from a pdb file read through a Handler
    //! @details Factory is the creator class for ProteinModel, Chain, AASequence, SSE, amino acid. The information is
    //! used from the Handler which already has the parsed pdb file information.
    //!
    //! @see @link example_pdb_factory.cpp @endlink
    //! @author woetzen, staritrd, karakam
    //! @date 3.3.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Factory :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Type of AAClass to be used
      biol::AAClass m_AAClass;

      //! List of PDB line printers to use
      util::ShPtrList< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > > m_Printers;

      //! static bool for whether all command line defaults have been saved
      static bool s_HaveDefaults;

    public:

      //! Flag to switch between numerated AtomIDs and the original PDBAtomIDs
      //! - default will numerated the atomIDs as they are written
      static util::ShPtr< command::FlagInterface> &GetFlagWritePDBAtomID();

      //! Flag to switch between numerated ResIDs( SeqIDs) and the original PDBResIDs and PDBICodes
      //! default will write SeqIDs
      static util::ShPtr< command::FlagInterface> &GetFlagWritePDBResID();

      //! command line flag to be used to change the Defaut AAClass that is used to generate AASequences/Chains/Proteins
      static util::ShPtr< command::FlagInterface> &GetFlagAAClass();

      //! command line flag to be used to change the default minimal size of the sse types tp be included in protein models
      static util::ShPtr< command::FlagInterface> &GetFlagMinSSESize();

      //! command line flag to be used, if pdb does not have any sse definition to use the backbone conformation to determine that
      static util::ShPtr< command::FlagInterface> &GetFlagSSEsFromBackBone();

      //! command line flag to be used, if pdb SSE definitions should be reassigned by DSSP
      static util::ShPtr< command::FlagInterface> &GetFlagDSSP();

      //! command line flag to be used to allow conversion of unnatural aminoacids to their parent amino acid
      static util::ShPtr< command::FlagInterface> &GetFlagConvertToNaturalAAType();

      //! command line flag to be used to write zero coordinates for undefined amino acids
      static util::ShPtr< command::FlagInterface> &GetFlagWriteZeroCoordinatesForUndefinedAminoAcids();

      //! command line flag to be used to write hydrogen atoms if available
      static util::ShPtr< command::FlagInterface> &GetFlagWriteHydrogens();

      //! command line flag to switch between IUPAC (standard) and pdb atom name output (when flag is used)
      static util::ShPtr< command::FlagInterface> &GetFlagPDBAtomName();

      //! command line flag for specifying which biomolecule number to read in for generating multimers
      static util::ShPtr< command::FlagInterface> &GetFlagBiomolecule();

      //! function to return all flags for the pdb factory
      static util::ShPtrVector< command::FlagInterface> &GetAllFlags();

      //! @brief reset defaults on all flags; to ensure a consistent set of defaults
      //! This function should be called whenever a given example or app may depend on the normal default parameters
      //! @return true on success
      //! @note should be called during static initialization to acquire the actual defaults
      static bool ResetFlagDefaults();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Factory();

      //! @brief constructor from AAClass
      //! @param AA_CLASS AAClass of interest
      Factory( const biol::AAClass &AA_CLASS);

      //! @brief virtual copy constructor
      Factory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the map of SSTypes and their corresponding minimal size
      //! @param FLAG the flag to get ss type min sizes from
      //! @return map of ss type and associated value
      static storage::Map< biol::SSType, size_t>
      GetCommandlineSSETypeMinSizes( const command::FlagInterface &FLAG = *GetFlagMinSSESize());

      //! @brief returns the map of SSTypes and their corresponding minimal size
      //! @param MIN_HELIX_SIZE the minimum helix size
      //! @param MIN_STRAND_SIZE the minimum strand size
      //! @param MIN_COIL_SIZE the minimum coil size
      //! @return map of ss type and associated value
      static storage::Map< biol::SSType, size_t>
      GetSSETypeMinSizes( const size_t &MIN_HELIX_SIZE, const size_t &MIN_STRAND_SIZE, const size_t &MIN_COIL_SIZE);

      //! @brief get the printers
      //! @return the printers
      const util::ShPtrList
      <
        util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> >
      > &GetPrinters() const
      {
        return m_Printers;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief resets the printers
      void ResetPrinters();

      //! @brief appends printer to list
      //! @param SP_PRINTER printer to be appended
      void AppendPrinter
      (
        const util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > > &SP_PRINTER
      );

    ///////////////////////
    // operations - read //
    ///////////////////////

      //! @brief builds biol::Atom from Line
      //! @param ATOM_TYPE type of atom in line
      //! @param LINE Line which contains Atom information
      //! @return ShPtr to Atom read from Line
      static util::ShPtr< biol::Atom>
      AtomFromLine
      (
        const biol::AtomType &ATOM_TYPE,
        const Line &LINE
      );

      //! @brief builds biol::Atom from Line
      //! @param LINE Line which contains Atom information
      //! @return ShPtr to Atom read from Line
      static
      util::ShPtr< biol::Atom>
      AtomFromLine( const Line &LINE);

      //! @brief builds an AASequence from a pdb file read by the given handler
      //! @param HANDLER Handler that contains pdb file
      //! @param PDB_ID Optional PDB ID to be used in the fasta header, otherwise, it is retrieved from the handler
      //! @return ShPtrVector of AASequences that were read from the given Handler
      util::ShPtrVector< biol::AASequence>
      AASequencesFromPDB( const Handler &HANDLER, const std::string &PDB_ID = "") const;

      //! @brief cut secondary structure elements from a given sequence
      //! @param SEQEUENCE the amino acid sequence
      //! @param SSE_RESIDUE secondary structure element defined by starting and end residue
      //! @return ShPtr to SSE, undefined if error occurs (no such amino acids, wrong chain, wrong amino acid order etc.)
      static util::ShPtr< assemble::SSE> SecondaryStructureElementFromSequence
      (
        const biol::AASequence &SEQUENCE,
        const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_RESIDUE
      );

      //! @brief cut secondary structure elements from a given sequence
      //! @param SEQEUENCE the amino acid sequence
      //! @param SSE_RESIDUES secondary structure elements defined by starting and end residue
      //! @return ShPtrList of SSEs
      static util::ShPtrVector< assemble::SSE> SecondaryStructureElementsFromSequence
      (
        const biol::AASequence &SEQUENCE,
        const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES
      );

      //! @brief merge two overlapping sequences into one
      //! @param SEQUENCE_A sequence a
      //! @param SEQUENCE_B sequence b
      //! @return the new sequence - empty if sses are not overlapping/consecutive
      static biol::AASequence MergeOverlappingSequences
      (
        const biol::AASequence &SEQUENCE_A,
        const biol::AASequence &SEQUENCE_B
      );

      //! @brief remove the amino acids that are common to both sses
      //! @param SSE_A overlapping sse a
      //! @param SSE_B overlapping sse b
      //! @return two sses that do not have the amino acids that were common to the two argument sses
      static storage::VectorND< 2, util::ShPtr< assemble::SSE> > RemoveOverlappingAAs
      (
        const util::ShPtr< assemble::SSE> &SSE_A,
        const util::ShPtr< assemble::SSE> &SSE_B
      );

      //! @brief process overlapping secondary structure elements
      //! @param SECONDARY_STRUCTURE_ELEMENTS possibly overlapping secondary structure elements
      //! @param MERGE_OVERLAPPING merge overlapping sses, if they are of the same type
      //! @return storage set of non overlapping SSEs
      static storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
      ProcessOverlappingSSEs
      (
        const util::ShPtrVector< assemble::SSE> &SECONDARY_STRUCTURE_ELEMENTS,
        const bool MERGE_OVERLAPPING
      );

      //! @brief add loops to a set of secondary structure elements
      //! @param SECONDARY_STRUCTURE_ELEMENTS secondary structure elements
      //! @param SEQUENCE the full sequenece; coils will be subsequences
      //! @return the set of secondary structure elements with coils
      static storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
      AddLoops
      (
        const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &SECONDARY_STRUCTURE_ELEMENTS,
        const biol::AASequence &SEQUENCE
      );

      //! @brief builds a Chain with sequence and SSEs from given Residues and SSE Residues
      //! builds a chain of template amino acids from a chainID and a util::ShPtrVector of SSEs represented as
      //! util::ShPtrVector of PDB Residues - each secondary structure element will have a minimal length of SSE_MINSIZE
      //! @param CHAIN_ID ChainID of the chain of interest
      //! @param SSE_RESIDUES Residues belonging to SSEs
      //! @param SEQUENCE_RESIDUES Residues from the sequence
      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
      //! @param PDB_ID PdbID of the protein
      //! @return Chain with the full sequence and SSEs read from given pdb information
      util::ShPtr< assemble::Chain>
      ChainFromPDBSSElementsAndSequence
      (
        const char CHAIN_ID,
        const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES,
        const storage::List< Residue> &SEQUENCE_RESIDUES,
        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
        const std::string &PDB_ID
      ) const;

      //! @brief builds chain without sses
      //! @param CHAIN_ID chain id of given fasta
      //! @param ISTREAM fasta file stream
      util::ShPtr< assemble::Chain>
      ChainFromFastaStream( const char CHAIN_ID, std::istream &ISTREAM) const;

      //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
      //! The types and sizes of SSEs to be considered are determined by the additional parameters
      //! @param HANDLER Handler that contains pdb information
      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
      //! @return Protein Model read from given Handler
      assemble::ProteinModel
      ProteinModelFromPDB
      (
        const Handler &HANDLER,
        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE = GetCommandlineSSETypeMinSizes()
      ) const;

      //! @brief builds an ensemble from the models stored in the handler
      //! @param HANDLER Handler that contains pdb information
      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
      //! @return Protein ensemble read from given Handler
      assemble::ProteinEnsemble ProteinEnsembleFromPDB
      (
        const Handler &HANDLER,
        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE = GetCommandlineSSETypeMinSizes()
      ) const;

//      //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
//      //! The types and sizes of SSEs to be considered are determined by the additional parameters
//      //! @param HANDLER Handler that contains pdb information
//      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
//      //! @return Protein Model read from given Handler
//      assemble::Complex
//      ComplexFromPDB
//      (
//        Handler &HANDLER,
//        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE = GetCommandlineSSETypeMinSizes()
//      ) const;

      //! @brief creates and returns a protein model based on a PDB filename
      //! @param PDB_FILENAME is the pdb filename from which the protein model will be created
      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
      //! @param IGNORE_CLASH whether to ignore clashes or not
      //! @return ProteinModel which was created from "PDB_FILENAME"
      assemble::ProteinModel ProteinModelFromPDBFilename
      (
        const std::string &PDB_FILENAME,
        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE = GetCommandlineSSETypeMinSizes(),
        const bool IGNORE_CLASH = false
      ) const;

    ////////////////////////
    // operations - write //
    ////////////////////////

      //! @brief write a single Protein Model to a pdb file stream
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param OSTREAM output stream to which ProteinModel will be written to
      //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
      //! @param CLASSIFICATION string for header classification
      //! @param TIME time (date) information to be placed in the header
      //! @return std::ostream to which ProteinModel was written to
      std::ostream &WriteModelToPDB
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM,
        const bool WRITE_BODY_INFORMATION = false,
        const std::string &CLASSIFICATION = "",
        const util::Time &TIME = util::Time::GetCurrent()
      ) const;

      //! @brief write a Chain to a pdb file stream
      //! @param CHAIN Chain of interest
      //! @param OSTREAM output stream to which Chain will be written to
      //! @return std::ostream to which Chain was written to
      static std::ostream &WriteChainToPDB( const assemble::Chain &CHAIN, std::ostream &OSTREAM);

      //! @brief write SSE definition for one SSE to Line
      //! @param THIS_SSE SSE of interest
      //! @param SERIAL SSE serial id
      //! @return ShPtr to Line containing SSE definition for given SSE
      static util::ShPtr< Line> WriteSSEDefinitionToLine
      (
        const assemble::SSE &THIS_SSE,
        const size_t &SERIAL
      );

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read Factory from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write Factory to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief construct the first side chain atom from given atoms
      //! those atoms need to contain the CA, N and O
      //! @param ATOMS
      //! @param BOND_LENGTH the distance of the constructed side chain atom coordinate to the CA atom
      //! @return position of first sidechain atom for L-amino acid
      static linal::Vector3D FirstSidechainAtomCoordinateFromBackboneAtoms
      (
        const util::SiPtrVector< const biol::Atom> &ATOMS,
        const double BOND_LENGTH
      );

      //! @brief returns the current Date using the PDB format
      //! @param TIME time (date) information to be placed in the header
      //! @return string of the current Date using the PDB format
      static std::string ConvertTimeToPDBDate( const util::Time &TIME);

      //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
      //! The types and sizes of SSEs to be considered are determined by the additional parameters
      //! @param HANDLER Handler that contains pdb information
      //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
      //! @param MODEL model in the handler to build
      //! @return Protein Model read from given Handler
      assemble::ProteinModel ProcessModel
      (
        const Handler &HANDLER,
        const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
        const Model &MODEL
      ) const;

    ///////////////////////
    // operations - read //
    ///////////////////////

      //! @brief builds all Atoms for provided AtomTypes from Lines
      //! @param LINES List of Lines which contains Atom information
      //! @param AA_TYPE aa type for which atoms are desired
      //! @param AA_TYPE_IN_FILE type of AA in the actual file
      //! @return ShPtrVector of Atoms read from given Lines
      util::ShPtrVector< biol::Atom>
      AtomsFromLines
      (
        const util::ShPtrList< Line> &LINES,
        const biol::AAType &AA_TYPE,
        const biol::AAType &AA_TYPE_IN_FILE
      ) const;

      //! @brief build amino acid from pdb residue
      //! @param RESIDUE Residue that contains amino acid information
      //! @param SEQ_ID sequence id of amino acid to be read
      //! @return ShPtr to amino acid that was read from given Residue
      util::ShPtr< biol::AABase>
      AminoAcidFromResidue( const Residue &RESIDUE, const int SEQ_ID = 1) const;

      //! @brief builds AASequence from CHAINID and util::ShPtrVector of Residues
      //! @param PDB_ID PdbID of the sequence to be read
      //! @param CHAIN_ID ChainID of the sequence to be read
      //! @param RESIDUES List of Residues that contain amino acid information for the sequence
      //! @return ShPtr to AASequence that was read from Residues
      util::ShPtr< biol::AASequence>
      AASequenceFromResidues
      (
        const std::string &PDB_ID,
        const char CHAIN_ID,
        const storage::List< Residue> &RESIDUES
      ) const;

    ////////////////////////
    // operations - write //
    ////////////////////////

    public:

      //! @brief write header information to Line
      //! @param CLASSIFICATION string for header classification
      //! @param TIME time (date) information to be placed in the header
      //! @return ShPtr to Line that contains Header information
      static util::ShPtr< Line> WriteHeaderToLine( const std::string &CLASSIFICATION, const util::Time &TIME);

      //! @brief write Atom to Line of a certain Resiue and the CHAINID and the atom SERIAL
      //! @param ATOM Atom of interest
      //! @param AMINO_ACID Amino acid which atom belgongs to
      //! @param CHAIN_ID chain id of the sequence atom belongs to
      //! @param SERIAL serial id of the atom
      //! @return ShPtr to Line that contains Atom information
      static util::ShPtr< Line> WriteAtomToLine
      (
        const biol::Atom &ATOM,
        const biol::AABase &AMINO_ACID,
        const char CHAIN_ID,
        const size_t SERIAL
      );

      //! @brief write Residue to Lines of a certain CHAINID and the atom SERIAL
      //! @param AMINO_ACID Amino acid of interest
      //! @param CHAIN_ID chain id of the sequence amino acid belongs to
      //! @param SERIAL serial id of the amino acid
      //! @return ShPtr to Line that contains amino acid information
      static util::ShPtrList< Line> WriteResiduesToLines
      (
        const biol::AABase &AMINO_ACID,
        const char CHAIN_ID,
        size_t &SERIAL
      );

      //! @brief write all atoms of a sequence to a list of lines to be used eventually for writing to a pdb
      //! @param AA_SEQUENCE AASequence of interest
      //! @param SERIAL serial id
      //! @return ShPtrList of Lines that contain AASequence information
      static util::ShPtrList< Line> WriteAASequenceToLines( const biol::AASequence &AA_SEQUENCE, size_t &SERIAL);

      //! @brief write a Chain to a SharedpointerVector of Lines
      //! @param CHAIN Chain of interest
      //! @param SERIAL serial id
      //! @return ShPtrList of Lines that contain Chain information
      static util::ShPtrList< Line> WriteChainToLines( const assemble::Chain &CHAIN, size_t &SERIAL);

      //! @brief write a single Protein Model to pdb lines
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
      //! @param CLASSIFICATION string for header classification
      //! @param TIME time (date) information to be placed in the header
      //! @return lines to which ProteinModel was written to
      util::ShPtrList< Line> WriteCompleteModelToPDBLines
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const bool WRITE_BODY_INFORMATION = false,
        const std::string &CLASSIFICATION = "",
        const util::Time &TIME = util::Time::GetCurrent()
      ) const;

      //! @brief write a ProteinModel to a SharedpointerVector of Lines
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param SERIAL serial id
      //! @return ShPtrList of Lines that contain ProteinModel information
      static util::ShPtrList< Line> WriteProteinModelToLines
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        size_t &SERIAL
      );

      //! @brief write helix definiton for one SSE to Line
      //! @param THIS_SSE SSE of interest
      //! @param HELIX_SERIAL helix serial id
      //! @return ShPtr to Line containing helix SSE definition for given SSE
      static util::ShPtr< Line> WriteHelixDefinitionToLine
      (
        const assemble::SSE &THIS_SSE,
        const size_t HELIX_SERIAL
      );

      //! @brief write Helix definitions to Lines
      //! @param SSE_VECTOR SiPtrVector of SSEs of interest
      //! @param HELIX_SERIAL helix serial id
      //! @return ShPtrList of Lines containing helix SSE defintions for given vector of SSEs
      static util::ShPtrList< Line> WriteHelixDefinitionsToLines
      (
        const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
        size_t &HELIX_SERIAL
      );

      //! @brief write Strand definition for one SSE to Line
      //! @param THIS_SSE SSE of interest
      //! @param STRAND_SERIAL strand serial id
      //! @return ShPtr to Line containing strand SSE definition for given SSE
      static util::ShPtr< Line> WriteStrandDefinitionToLine
      (
        const assemble::SSE &THIS_SSE,
        const size_t STRAND_SERIAL
      );

      //! @brief write strand definitions to Lines
      //! @param SSE_VECTOR SiPtrVector of SSEs of interest
      //! @param STRAND_SERIAL strand serial id
      //! @return ShPtrList of Lines containing strand SSE definitions for given vector of SSEs
      static util::ShPtrList< Line> WriteStrandDefinitionsToLines
      (
        const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
        size_t &STRAND_SERIAL
      );

      //! @brief write assemble::SSEGeometry Information to Lines of one SSE
      //! @param THIS_SSE SSE of interest
      //! @return ShPtrList of Lines containing body information for the given SSE
      static util::ShPtrList< Line> WriteBodyInformationToLines( const assemble::SSE &THIS_SSE);

      //! @brief write assemble::SSEGeometry Information to Lines for given ProteinModel
      //! writes TransformationMatrix and assemble::SSEGeometryextension of all SSES to Lines
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return ShPtrList of Lines containing body information for the given ProteinModel
      static util::ShPtrList< Line> WriteBodyInformationToLines( const assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief write SSE definitons to Lines for a util::ShPtrVector of models
      //! @param PROTEIN_MODELS ShPtrVector of ProteinModels of interest
      //! @return ShPtrList of Lines containing SSE definition information for given ShPtrVector of ProteinModels
      static util::ShPtrList< Line> WriteSSEDefinitionsToLines
      (
        const util::ShPtrVector< assemble::ProteinModel> &PROTEIN_MODELS
      );

      //! @brief write SEQRES lines for a sequence
      //! @param SEQUENCE AASequence of interest
      //! @return ShPtrList of Lines containing SEQRES information for given AASequence
      static util::ShPtrList< Line> WriteSeqResToLines( const biol::AASequence &SEQUENCE);

      //! @brief write SEQRES information of a Chain into Lines
      //! @param CHAIN Chain of interest
      //! @return ShPtrList of Lines containing SEQRES information for given Chain
      static util::ShPtrList< Line> WriteSeqResToLines( const assemble::Chain &CHAIN);

      //! @brief write SEQRES information of a ProteinModel into Lines
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return ShPtrList of Lines containing SEQRES information for given ProteinModel
      static util::ShPtrList< Line> WriteSeqResToLines( const assemble::ProteinModel &PROTEIN_MODEL);

    }; // class Factory

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_FACTORY_H_
