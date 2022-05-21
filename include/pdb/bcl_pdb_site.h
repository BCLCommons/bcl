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

#ifndef BCL_PDB_SITE_H_
#define BCL_PDB_SITE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_ligand.h"
#include "bcl_pdb_residue_simple.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Site
    //! @brief represent a SITE definition
    //! @details A site definition contains an interesting site like a ligand binding, catalytic, metal binding or other
    //!          site with a relevant function.
    //!          A ligand can also be present, but is not necesaary.
    //!
    //! @see @link example_pdb_site.cpp @endlink
    //! @author woetzen
    //! @date Jan 31, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Site :
      public util::ObjectInterface
    {
    public:

      //! @enum evidence code for site
      enum EvidenceCode
      {
        e_Software, //!< determined by software
        e_Author,   //!< determined by author
        e_Unknwon,  //!< unknown
        s_NumberEvidenceCode
      };

      //! @brief EvidenceCode as string
      //! @param EVIDENCE_CODE the EvidenceCode
      //! @return the string for the EvidenceCode
      static const std::string &GetEvidenceCodeDescriptor( const EvidenceCode &EVIDENCE_CODE);

      //! @brief PropertyTypeEnum is used for I/O of PropertyType
      typedef util::WrapperEnum< EvidenceCode, &GetEvidenceCodeDescriptor, s_NumberEvidenceCode> EvidenceCodeEnum;

    private:

    //////////
    // data //
    //////////

      //! identifier from "REMARK 800 SITE_IDENTIFIER:"
      std::string m_Identifier;

      //! evidence code from "REMARK 800 EVIDENCE_CODE:"
      EvidenceCodeEnum m_EvidenceCode;

      //! description from "REMARK 800 SITE_DESCRIPTION:"
      std::string m_Description;

      //! residues comprising the site within the chain
      storage::List< ResidueSimple> m_ChainResidues;

      //! residues comprisiing the site that are HETATM / other ligand contacts
      storage::List< ResidueSimple> m_HetResidues;

      //! ptr to the ligand
      util::ShPtr< Ligand> m_Ligand;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Site();

      //! @brief construct from identifier, evidence code and description
      Site( const std::string &NAME, const EvidenceCodeEnum &EVIDENCE_CODE, const std::string &DESCRIPTION);

      //! @brief Clone function
      //! @return pointer to new Site
      Site *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief name/identifier of the site
      //! @return the name as string
      const std::string &GetName() const;

      //! @brief get evidence code
      //! @return EvidenceCodeEnum
      const EvidenceCodeEnum &GetEvidenceCode() const;

      //! @brief site description
      //! @return description of site
      const std::string &GetDescription() const;

      //! @brief access to the residues comprising the site that are within a chain
      //! @return List of residues
      const storage::List< ResidueSimple> &GetChainResidues() const;

      //! @brief ligand that might be associated with that site
      //! @return ShPtr to Ligand; null if there is no ligand
      const util::ShPtr< Ligand> &GetLigand() const;

      //! @brief add a residue from a chain
      //! @param RESIDUE the residue to  add to the site
      void AddChainResidue( const ResidueSimple &RESIDUE);

      //! @brief add a residue from hetatm section
      //! @param RESIDUE the residue to  add to the site
      void AddHetatmResidue( const ResidueSimple &RESIDUE);

    ////////////////
    // operations //
    ////////////////

      //! @brief locators for the residues associated with the site definition
      //! @return the amino acid locators in a list
      storage::List< assemble::LocatorAA> AALocators() const;

      //! @brief find ligand for site
      //! @param LIGANDS list of ligands, that might be associated with that site
      //! @return true, if a ligand was found, false otherwise
      bool FindLigand( const util::ShPtrList< Ligand> &LIGANDS);

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

    }; // class Site

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_SITE_H_ 
