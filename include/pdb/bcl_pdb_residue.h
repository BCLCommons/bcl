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

#ifndef BCL_PDB_RESIDUE_H_
#define BCL_PDB_RESIDUE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line.h"
#include "bcl_pdb_residue_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Residue
    //! @brief helper class for representing amino acid information for pdb reader
    //! @details This class is a helper class for the pdb reader and stores informations for one residue in the pdb
    //! if contains the residue name, chain id, pdbid, insertion code and holdes all ATOM lines for the residue
    //!
    //! @see @link example_pdb_residue.cpp @endlink
    //! @author woetzen
    //! @date 08/21/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Residue :
      public ResidueInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string            m_ResidueName; //! residue name
      char                   m_ChainID;     //! chain ID
      int                    m_PDBID;       //! pdb id
      char                   m_ICode;       //! insertion code
      util::ShPtrList< Line> m_Lines;       //! pdb lines as sharedpointervector

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Residue();

      //! construct Residue from ResidueName and ChainID
      Residue
      (
        const std::string &RESIDUENAME,
        const char &CHAINID
      );

      //! construct Residue from ResidueName, ChainID, ID, ICode
      Residue
      (
        const std::string &RESIDUENAME,
        const char &CHAINID,
        const int  &PDBID,
        const char &ICODE
      );

      //! construct Residue from ResidueName, ChainID, ID, ICode, and Lines
      Residue
      (
        const std::string &RESIDUENAME,
        const char &CHAINID,
        const int  &PDBID,
        const char &ICODE,
        const util::ShPtrList< Line> &LINES
      );

      //! @brief copy constructor
      Residue( const Residue &RESIDUE);

      //! @brief construct form Residue
      //! @param RESIDUE the residue to copy
      Residue( const ResidueInterface &RESIDUE);

      //! copy constructor
      Residue *Clone() const;

      //! destructor
      ~Residue();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to pdb ID
      //! @return pdb sequence id
      int GetPDBID() const
      {
        return m_PDBID;
      }

      //! @brief set ID
      //! @param PDB_ID new pdb id for residue
      void SetPDBID( const int PDB_ID)
      {
        m_PDBID = PDB_ID;
      }

      //! @brief return insertion code
      //! @return pdb sequence insertion code, if multiple residues have the same pdb id
      char GetICode() const
      {
        return m_ICode;
      }

      //! @brief set insertion code
      //! @param I_CODE insertion code for that residue
      void SetICode( const char I_CODE)
      {
        m_ICode = I_CODE;
      }

      //! @brief return chain id
      //! @return the chain id this residue belongs to
      char GetChainID() const
      {
        return m_ChainID;
      }

      //! @brief set chain id
      //! @param CHAIN_ID the chain id
      void SetChainID( const char CHAIN_ID)
      {
        m_ChainID = CHAIN_ID;
      }

      //! @brief return residue name
      //! @return three letter code residue name
      const std::string &GetResidueName() const
      {
        return m_ResidueName;
      }

      //! @brief set residue name
      //! @param RESIDUE_NAME new residue name for that residue
      void SetResidueName( const std::string &RESIDUE_NAME)
      {
        m_ResidueName = RESIDUE_NAME;
      }

      //! @brief return Lines
      //! @return const reference to lines
      util::ShPtrList< Line> const &GetLines() const;

      //! @brief change Lines
      //! @return changeable reference to lines
      util::ShPtrList< Line> &ChangeLines();

    ////////////////
    // operations //
    ////////////////

      //! @brief get all atom serials for that residue
      //! @return set of atom serials
      storage::Set< size_t> AtomSerials() const;

      //! @brief line criterium to locate atom/hetatom lines for that residue
      //! @param LINE_TYPE the line type of interest
      //! @return criterium to check if a line corresponds to that residue
      LineCriterium GetCriterium( const LineType &LINE_TYPE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief compare two residues if they are equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid and insertion code are equal
      bool operator ==( const ResidueInterface &RESIDUE_RHS) const;

      //! @brief compare two residues if they are not equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid or insertion code are not equal
      bool operator !=( const ResidueInterface &RESIDUE_RHS) const;

      //! @brief compare two residues if the lhs is less than the rhs
      //! @param RESIDUE_RHS rhs residue
      //! @return true if chain id or pdbid or insertion code of lhs is smaller
      bool operator <( const ResidueInterface &RESIDUE_RHS) const;

      //! @brief compare two residues if they are equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid and insertion code are equal
      bool operator ==( const Residue &RESIDUE_RHS) const;

      //! @brief compare two residues if they are not equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid or insertion code are not equal
      bool operator !=( const Residue &RESIDUE_RHS) const;

      //! @brief compare two residues if the lhs is less than the rhs
      //! @param RESIDUE_RHS rhs residue
      //! @return true if chain id or pdbid or insertion code of lhs is smaller
      bool operator <( const Residue &RESIDUE_RHS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write Residue to STREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read Residue from std::istream
      std::istream &Read( std::istream &ISTREAM);

    }; // class Residue

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_RESIDUE_H_
