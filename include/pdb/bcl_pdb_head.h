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

#ifndef BCL_PDB_HEAD_H_
#define BCL_PDB_HEAD_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line.h"
#include "bcl_pdb_line_group_interface.h"
#include "bcl_pdb_line_type_data.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Head
    //! @brief all pdb lines above the coordinate section
    //! @details store all lines above the coordinate (ATOM) section in the order of the line types.
    //!
    //! @see @link example_pdb_head.cpp @endlink
    //! @author woetzen
    //! @date Feb 21, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Head :
      public LineGroupInterface
    {

    private:

    ///////////
    // types //
    ///////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class RemarkLineLessThan
      //! @brief has operator for checking if one remark line is less than the other one
      //! @remarks example unnecessary
      //! @author weinerbe
      //! @date Feb 17, 2012
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct RemarkLineLessThan
      {
         //! @brief operator for checking if one line is less than the other one
         //! @param LINE_LHS line on left-hand side
         //! @param LINE_RHS line on right-hand side
         //! @return if one line is less than the other one
         bool operator()( const util::ShPtr< Line> &LINE_LHS, const util::ShPtr< Line> &LINE_RHS) const;

      }; // struct RemarkLineLessThan

    //////////
    // data //
    //////////

      //! Data, containing LineType and the content for each pdb line
      storage::Map< LineType, util::ShPtrList< Line> > m_Lines;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! BCL pdb identifier
      static const std::string &GetBCLPdbID();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Head();

      //! @brief Clone function
      //! @return pointer to new Head
      Head *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief linetypes within group
      //! @return set of line types
      const storage::Set< LineType> &GetTypesOfLines() const;

      //! @brief access to lines of given type
      //! @param LINE_TYPE the desire line type
      //! @return lines of given type
      util::ShPtrList< Line> GetLines( const LineType &LINE_TYPE) const;

      //! @brief count the number of lines of given TYPE used for master record
      //! @param LINE_TYPE the line type
      //! @return the number of lines of that type found
      size_t Count( const LineType &LINE_TYPE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @return lines that are considered by criterium
      util::ShPtrList< Line> CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const;

      //! @brief locate lines of given criterium and line type
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @param LINE_TYPE only consider this line type
      //! @return lines that are considered by criterium
      util::ShPtrList< Line>
      CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const LineType &LINE_TYPE) const;

      //! @brief pushback a new line into that group
      //! @param ShPtr to the line
      //! @return true, if it fits into that group (line type is eligible)
      bool PushBack( const util::ShPtr< Line> &LINE);

      //! @brief reset the line group
      void Reset();

      //! extract all transformation matrices
      storage::Vector< math::TransformationMatrix3D> GetTransformationMatrices() const;

      //! @brief the transformation matrices to generate different bio molecules by applying transformations to different chains
      //! @param CHAIN_IDS all chain ids that should be considered
      //! @return Map of biomolecule number to a vector of transformations of chainid it should be applied to, the new chain id and the transformation
      storage::Map< size_t, storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > >
      GetBioTransformationMatrices( const std::string &CHAIN_IDS) const;

      //! @brief get the membrane from the remark lines
      //! @return the membrane from the remark lines
      util::ShPtr< biol::Membrane> GetMembrane() const;

      //! @brief gets the pdb ID read in from the pdb file
      //! @return the pdb ID read in from the pdb file
      std::string GetPDBID() const;

      //! @brief full name for all het residues
      //! @return map with key res name and data fullname
      storage::Map< std::string, std::string> GetHetFullname() const;

      //! @brief formula for all het residues
      //! @return map with key res name and data formula
      storage::Map< std::string, std::string> GetHetFormula() const;

      //! @brief create the chains as they were given in the seqres with three letter code
      //! @return map with chainid as key and list of residues as data
      storage::Map< char, storage::List< ResidueSimple> > GetSEQRESProteinChains() const;

      //! @brief map of residues that could not be located in the experiments using REMARK 465
      //! @param MODEL_NUMBER missing residues for the given model number
      //! @return map of chains and List of Residues that are not in the ATOM lines
      storage::Map< char, storage::List< ResidueSimple> > GetMissingResidues( const size_t MODEL_NUMBER) const;

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

    public:

      //! @brief pdb line to sse definition form a HELIX or SHEET line, with sstype, start and end residue
      //! @return Triplet of SSType starting and end residue
      static storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
      SSEDefinitionFromPDBLine( const Line &PDB_LINE);

      //! @brief get all sse definitions as they are given in the HELIX/SHEET section
      //! @param HELIX_CLASSES set of helix types to consider
      //! @param MERGE_OVERLAPPING merge overlapping sses of the same type
      //! @return Map of chain ids and
      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
      SSEDefinitions( const storage::Set< biol::SSType> &HELIX_CLASSES, const bool MERGE_OVERLAPPING) const;

      //! @brief check sse definitions for overlap and merge
      //! @param SSE_DEFINITIONS
      //! @return Map of chain ids and sse definitions that were merged when possible
      static storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
      SSEDefinitionsMergeOverlap( const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > &SSE_DEFINITIONS);

      //! @brief compare the sse infos for overlap
      //! @param SSE_INFO_LHS left hand side sse info
      //! @param SSE_INFO_RHS right hand side sse info
      static bool SSEInformationCompare
      (
        const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_LHS,
        const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_RHS
      );

    }; // class Head

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_HEAD_H_ 
