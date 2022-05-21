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

#ifndef BCL_PDB_RESIDUE_INTERFACE_H_
#define BCL_PDB_RESIDUE_INTERFACE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_format.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ResidueInterface
    //! @brief helper interface for representing hetero residues within pdb files
    //! @details This class is a helper class for the pdb reader and stores informations for one residue in the pdb
    //! if contains the residue name, chain id, pdbid, insertion code and holdes all ATOM/HETATM lines for the residue
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 03/22/2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ResidueInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief access to pdb ID
      //! @return pdb sequence id
      virtual int GetPDBID() const = 0;

      //! @brief set ID
      //! @param PDB_ID new pdb id for residue
      virtual void SetPDBID( const int PDB_ID) = 0;

      //! @brief return insertion code
      //! @return pdb sequence insertion code, if multiple residues have the same pdb id
      virtual char GetICode() const = 0;

      //! @brief set insertion code
      //! @param I_CODE insertion code for that residue
      virtual void SetICode( const char I_CODE) = 0;

      //! @brief return chain id
      //! @return the chain id this residue belongs to
      virtual char GetChainID() const = 0;

      //! @brief set chain id
      //! @param CHAIN_ID the chain id
      virtual void SetChainID( const char CHAIN_ID) = 0;

      //! @brief return residue name
      //! @return three letter code residue name
      virtual const std::string &GetResidueName() const = 0;

      //! @brief set residue name
      //! @param RESIDUE_NAME new residue name for that residue
      virtual void SetResidueName( const std::string &RESIDUE_NAME) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get Identification of this amino acid
      //! @return string with identification
      std::string GetIdentification() const
      {
        return GetResidueName() + " " + GetChainID() + " " + util::Format().W( 5)( GetPDBID()) + GetICode();
      }

      //! @brief get all atom serials for that residue
      //! @return set of atom serials
      virtual storage::Set< size_t> AtomSerials() const = 0;

      //! @brief line criterium to locate atom/hetatom lines for that residue
      //! @param LINE_TYPE the line type of interest
      //! @return criterium to check if a line corresponds to that residue
      virtual LineCriterium GetCriterium( const LineType &LINE_TYPE) const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief compare two residues if they are equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid and insertion code are equal
      virtual bool operator ==( const ResidueInterface &RESIDUE_RHS) const = 0;

      //! @brief compare two residues if they are not equal
      //! @param RESIDUE_RHS rhs residue
      //! @return true if residue name, chain id, pdbid or insertion code are not equal
      virtual bool operator !=( const ResidueInterface &RESIDUE_RHS) const = 0;

      //! @brief compare two residues if the lhs is less than the rhs
      //! @param RESIDUE_RHS rhs residue
      //! @return true if chain id or pdbid or insertion code of lhs is smaller
      virtual bool operator <( const ResidueInterface &RESIDUE_RHS) const = 0;

    }; // class ResidueInterface

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_RESIDUE_INTERFACE_H_
