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

#ifndef BCL_PDB_MODEL_H_
#define BCL_PDB_MODEL_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line.h"
#include "bcl_pdb_line_group_interface.h"
#include "bcl_pdb_residue.h"
#include "biol/bcl_biol_ss_types.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Model
    //! @brief contains all lines and functionality necessary for a pdb model (coordinate section of a pdb file)
    //! @details
    //!
    //! @see @link example_pdb_model.cpp @endlink
    //! @author woetzen
    //! @date Feb 20, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Model :
      public LineGroupInterface
    {

    private:

    //////////
    // data //
    //////////

      //! lines belonging to that model sorted by chaind
      storage::Map< char, util::ShPtrList< Line> > m_ChainLines;

      //! lines containing the hetatm lines that do not belong to a chain
      storage::Map< char, util::ShPtrList< Line> > m_HetatmLines;

      //! structured chain
      storage::Map< char, storage::List< Residue> > m_StructuredChains;

      //! this is the number of residues that must align to be sure that the seqres and atomline residues truly match
      //! 3kvnX and 3zuxA have 5 residue his tags at the start of the model with no coordinates, which leads to AA's
      //! with undefined PDB ids in the resulting model if this number is set less than 5
      enum { s_NumberRequiredAlignedResidues = 6};

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Model();

      //! @brief Clone function
      //! @return pointer to new Model
      Model *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief linetypes within group
      //! @return set of line types
      const storage::Set< LineType> &GetTypesOfLines() const;

      //! @brief get hetatm lines
      //! @return reference to HETATM lines
      const storage::Map< char, util::ShPtrList< Line> > &GetHETATMLines() const;

      //! @brief access to lines of given type
      //! @param LINE_TYPE the desire line type
      //! @return lines of given type
      util::ShPtrList< Line> GetLines( const LineType &LINE_TYPE) const;

      //! @brief count the number of lines for a given line type
      //! @param LINE_TYPE the desire line type
      //! @return the number of lines for the given line type in the model
      size_t Count( const LineType &LINE_TYPE) const;

      //! @brief access to the structure chains
      //! @return reference to map of chains with residues, containing the pdb lines
      const storage::Map< char, storage::List< Residue> > &GetStructuredChains() const
      {
        return m_StructuredChains;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @return lines that are considered by criterium
      util::ShPtrList< Line> CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const;

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @param CHAIN_ID only lines for that chain if know apriori
      //! @return lines that are considered by criterium
      util::ShPtrList< Line>
      CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const char CHAIN_ID) const;

      //! @brief pushback a new line into that group
      //! @param ShPtr to the line
      //! @return true, if it fits into that group (line type is eligible)
      bool PushBack( const util::ShPtr< Line> &LINE);

      //! @brief reset the line group
      void Reset();

      //! @brief create the chains as they were given by the atoms
      //! @return map with chainid as key and list of residues as data
      storage::Map< char, storage::List< ResidueSimple> > GetChains() const;

      //! @brief initialize the structure
      //! @param CHAIN_ID
      //! @param SEQRES sequences for the chain
      //! @param MISSING_RESIDUES residues that are missing in the ATOM section (REMARK 465)
      void InitializeStructuredChain
      (
        const char CHAIN_ID,
        const storage::List< ResidueSimple> &SEQRES,
        const storage::List< ResidueSimple> &MISSING_RESIDUES
      );

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      std::ostream &WriteLines( std::ostream &OSTREAM) const;

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

      //! @brief list that stores for each pdbid and insertion code a map of alternate location ids their residue
      //! HETATM lines are also considered, but are only stored as hardcopies within the residues as ATOM lines
      //! @param CHAINID chain for which this list is generated
      //! @return list that stores for each pdbid and insertion code a map of alternate location ids their residue
      storage::List< storage::Map< char, Residue> >
      AtomLineResiduesFromChainID( const char CHAINID) const;

      //! @brief complete first residue in a map of alternate residues
      //! @param RESIDUES map of residues where first one is to be completed
      static
      void CompleteFirstAlternateResidue( storage::Map< char, Residue> &RESIDUES);

      //! @brief add missing residues to the atom line residues
      //! @param ATOM_LINE_RESIDUES residues acquired form residues
      //! @param MISSING_RESIDUES residues that are not located in experiment REMARK 465
      //! @return list of maps (alternate location residues) merged with the missing residues
      storage::List< storage::Map< char, Residue> >
      static MergeAtomLinesAndMissingResidueLine
      (
        const storage::List< storage::Map< char, Residue> > &ATOM_LINE_RESIDUES,
        const storage::List< ResidueSimple> &MISSING_RESIDUES
      );

      //! @brief compare a map of alternate residues to a given residue if one matches by name
      //! @param ATOM_RES map of alternate residues
      //! @param RES residue to be matched
      //! @param SET_RESIDUE if match was found, set the given RES to the residue
      //! @return true, if match was found
      static
      bool CompareAtomLineResWithResByName
      (
        const storage::Map< char, Residue> &ATOM_RES,
        Residue &RESIDUE,
        const bool SET_RESIDUE
      );

      //! @brief Get the number of residues before getting to the next type
      //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
      //! @param ATOM_LINE_ITR_END total end of atom lines
      //! @param RES_ITR itr to residue that matches up by name with the atom line
      //! @param RES_ITR_END total end of residues
      //! @return min(number next residues in atom lines that have the same res type,
      //!             number next residues in seqres lines that have the same res type)
      static size_t GetLengthHomopolymerChain
      (
        const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
        const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
        const storage::List< Residue>::iterator &RES_ITR,
        const storage::List< Residue>::iterator &RES_ITR_END
      );

      //! @brief align range of residues and atom lines
      //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
      //! @param RES_ITR itr to residue that matches up by name with the atom line
      //! @param ATOM_LINE_ITR_END total end of atom lines
      //! @param RES_ITR_END total end of residues
      //! @param SET_RESIDUE_ATOM_LINES if true, atom lines are inserted into the residues
      //! @param MAX_NUMBER_RES_TO_ALIGN stop aligning after this number is exceeded
      //! @return number fo residues that have been aligned
      static size_t AlignResiduesAndAtomLines
      (
        storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
        storage::List< Residue>::iterator    &RES_ITR,
        const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
        const storage::List< Residue>::iterator &RES_ITR_END,
        const bool SET_RESIDUE_ATOM_LINES,
        const size_t MAX_NUMBER_RES_TO_ALIGN
      );

      //! @brief Inserts missing resid informations if they are in HETATM lines instead of ATOM lines
      //! @param CHAIN_ID the chain if of that chain
      //! @param CHAIN the residues for that chain
      void InsertLinesFromHETATMLines
      (
        const char CHAIN_ID,
        storage::List< Residue> &CHAIN
      ) const;

    }; // class Model

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_MODEL_H_
