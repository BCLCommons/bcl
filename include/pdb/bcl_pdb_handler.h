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

#ifndef BCL_PDB_HANDLER_H_
#define BCL_PDB_HANDLER_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_head.h"
#include "bcl_pdb_line_group_interface.h"
#include "bcl_pdb_model.h"
#include "bcl_pdb_residue_simple.h"
#include "bcl_pdb_tail.h"
#include "biol/bcl_biol_ss_types.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Handler
    //! @brief Handler reads the pdb and parses the information and stores the information in an organized way
    //! @details handler allows reading of pdbs and parsing of the information and storing the information such as
    //! SSE definitions.
    //!
    //! @see @link example_pdb_handler.cpp @endlink
    //! @author staritrd, meilerj, haenigc, woetzen, karakam
    //! @date 11/11/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Handler :
      public LineGroupInterface
    {

    private:

    //////////
    // data //
    //////////

      //! head of pdb file
      Head m_Head;

      //! all models that are contained
      storage::Vector< Model> m_Models;

      //! rest
      Tail m_Tail;

      //! stores structure of chains, residues
      storage::Map< char, storage::List< ResidueSimple> > m_SEQRES;
      //                    chainID Residues

      //! stores information of Secondary structure elements
      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > m_SSEStructure;
      // Chains  ChainID                                        SSEType  start_residue  last_residue

      //! building structure without seqres checks if two residues have a meaningful peptide bond distance - with that flag on, it doe not
      bool m_IgnoreClash;

      //! set of helix classes to consider
      storage::Set< biol::SSType> m_HelixClasses;

      //! chainids of non protein chains
      storage::Set< char> m_NonProteinChainsIDs;

    public:

      //! Flag to choose which helix classes to process
      //! default is helix class 1 (right handed alpha helix) if there are no classes given in the list
      static util::ShPtr< command::FlagInterface> &GetFlagHelixClasses();

      //! Flag to trigger the merging of overlapping SSE definitions - although they should not occur, there are still
      //! many pdb files that do not adhere to the format standard
      static util::ShPtr< command::FlagInterface> &GetFlagMergeOverlappingSSEs();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Handler( const bool IGNORE_CLASH = false);

      //! construct PDBRreader from the std::string, containing the pdb file. It reads all information and performs a
      //! SequenceCheck. NOT appropriate when there is no SEQRES information in the pdb file!!!
      //! In that case, use default constructor and call function Read().
      Handler( std::istream &ISTREAM, const bool IGNORE_CLASH = false);

      //! copy constructor
      Handler( const Handler &PDBREADER, const bool IGNORE_CLASH = false);

      //! copy constructor
      Handler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to head
      //! @return reference to head
      const Head &GetHead() const
      {
        return m_Head;
      }

      //! @brief access the sse structure
      //! @return reference on sse structure
      const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > &
      GetSSEStructure() const
      {
        return m_SSEStructure;
      }

      //! @brief linetypes within group
      //! @return set of line types
      const storage::Set< LineType> &GetTypesOfLines() const;

      //! @brief return the Lines for a given line type
      //! @param LINE_TYPE the line type of interest
      //! @return reference to the lines - can be empty if no lines are available
      util::ShPtrList< Line> GetLines( const LineType &LINE_TYPE) const;

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @return lines that are considered by criterium
      util::ShPtrList< Line> CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const;

      //returns the entire sequence as one letter code of a chain CHAINID
      std::string GetSequence( const char CHAINID) const;

      //! @brief get the helix classes to be considered as given in the command line
      //! @return set of helix classes
      static storage::Set< biol::SSType> GetHelixClassesFromCommandLine();

      //! @brief access to all models
      //! @return reference to models
      const storage::Vector< Model> &GetModels() const
      {
        return m_Models;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get all chain ids
      std::string ChainIDs() const;

      //! @brief get chain ids for protein chains
      std::string GetProteinChains() const;

      //! @brief pushback a new line into that group
      //! @param LINE ShPtr to the line
      //! @return true, if it fits into that group (line type is eligible)
      bool PushBack( const util::ShPtr< Line> &LINE);

      //! @brief reset the line group
      void Reset();

      //! @brief Appends a list of lines
      //! @param PDBLINES list of lines
      //! @return true, if all lines could be inserted
      bool AppendLines( const util::ShPtrList< Line> &PDBLINES);

      //! @brief return all ligands
      //! @return List of ligands
      util::ShPtrList< Ligand> GetLigands() const;

      //! @brief return all sites
      //! @return list of all sites
      util::ShPtrList< Site> GetSites() const;

      //! @brief check if given residue is in one of the chains (if not might be in the HETATM section)
      //! @param RESIDUE the residue to check for
      //! @return true, if that residue in found in one of the chains, false oterhwise
      bool IsInChain( const ResidueInterface &RESIDUE) const;

      //! @brief extracts the path and tag from the provided full path
      static storage::VectorND< 2, std::string> ExtractPathAndPDBTag( const std::string &FULL_PATH);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      std::ostream &WriteLines( std::ostream &OSTREAM) const;

    protected:

      //! reads Handler from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! Writes all existing lines into std::ostrem STREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Handler

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_HANDLER_H_
