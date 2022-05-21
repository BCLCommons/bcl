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

#ifndef BCL_PDB_ENTRY_TYPES_H_
#define BCL_PDB_ENTRY_TYPES_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_pdb_entry_type_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EntryTypes
    //! @brief collecting all pdb entry types that can occur within a pdb
    //! @details This class is a collection and enumerates all possible entry types that can be found in pdb lines. The name of
    //! each EntryType begins with the line type and is followed by the entry type. The pdb defines those entry types:
    //! http://www.wwpdb.org/documentation/format23/v2.3.html
    //!
    //! @see @link example_pdb_entry_types.cpp @endlink
    //! @author woetzen
    //! @date Feb 8, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EntryTypes :
      public util::Enumerate< EntryTypeData, EntryTypes>
    {
      friend class util::Enumerate< EntryTypeData, EntryTypes>;

    public:

    //////////
    // data //
    //////////

      const EntryType HEADERClassification;
      const EntryType HEADERDate;
      const EntryType HEADERIDCode;
      const EntryType REMARK_Number;
      const EntryType REMARK_String;
      const EntryType REMARK_350_BiomoleculeIdentifier;
      const EntryType REMARK_350_BiomoleculeNumber;
      const EntryType REMARK_350_ChainsIdentifier;
      const EntryType REMARK_350_ChainsList;
      const EntryType REMARK_350_BIOMT1_Identifier;
      const EntryType REMARK_350_BIOMT1_Number;
      const EntryType REMARK_350_BIOMT1_M11;
      const EntryType REMARK_350_BIOMT1_M12;
      const EntryType REMARK_350_BIOMT1_M13;
      const EntryType REMARK_350_BIOMT1_V1;
      const EntryType REMARK_350_BIOMT2_Identifier;
      const EntryType REMARK_350_BIOMT2_Number;
      const EntryType REMARK_350_BIOMT2_M21;
      const EntryType REMARK_350_BIOMT2_M22;
      const EntryType REMARK_350_BIOMT2_M23;
      const EntryType REMARK_350_BIOMT2_V2;
      const EntryType REMARK_350_BIOMT3_Identifier;
      const EntryType REMARK_350_BIOMT3_Number;
      const EntryType REMARK_350_BIOMT3_M31;
      const EntryType REMARK_350_BIOMT3_M32;
      const EntryType REMARK_350_BIOMT3_M33;
      const EntryType REMARK_350_BIOMT3_V3;
      const EntryType REMARK_465_ModelNumber;
      const EntryType REMARK_465_ResidueName;
      const EntryType REMARK_465_ChainID;
      const EntryType REMARK_465_ResidueSequenceID;
      const EntryType REMARK_465_InsertionCode;
      const EntryType REMARK_800_SiteIdentifierString;
      const EntryType REMARK_800_SiteIdentifier;
      const EntryType REMARK_800_SiteEvidenceCodeString;
      const EntryType REMARK_800_SiteEvidenceCode;
      const EntryType REMARK_800_SiteDescriptionString;
      const EntryType REMARK_800_SiteDescription;
      const EntryType SEQRESSerial;
      const EntryType SEQRESChainID;
      const EntryType SEQRESNrOfResiduesInChain;
      const EntryType SEQRESName_1;
      const EntryType SEQRESName_2;
      const EntryType SEQRESName_3;
      const EntryType SEQRESName_4;
      const EntryType SEQRESName_5;
      const EntryType SEQRESName_6;
      const EntryType SEQRESName_7;
      const EntryType SEQRESName_8;
      const EntryType SEQRESName_9;
      const EntryType SEQRESName_10;
      const EntryType SEQRESName_11;
      const EntryType SEQRESName_12;
      const EntryType SEQRESName_13;
      const EntryType HETIdentifier;
      const EntryType HETChainID;
      const EntryType HETSequenceID;
      const EntryType HETInsertionCode;
      const EntryType HETNumberAtoms;
      const EntryType HETDescription;
      const EntryType HETNAMContinuation;
      const EntryType HETNAMIdentifier;
      const EntryType HETNAMText;
      const EntryType HETSYNContinuation;
      const EntryType HETSYNIdentifier;
      const EntryType HETSYNSynonyms;
      const EntryType FORMULComponentNumber;
      const EntryType FORMULIdentifier;
      const EntryType FORMULContinuation;
      const EntryType FORMULAsterisk;
      const EntryType FORMULChemicalFormula;
      const EntryType HELIXSerial;
      const EntryType HELIXID;
      const EntryType HELIXResidueName_Initial;
      const EntryType HELIXChainID_Initial;
      const EntryType HELIXSequenceID_Initial;
      const EntryType HELIXInsertionCode_Initial;
      const EntryType HELIXResidueName_Terminal;
      const EntryType HELIXChainID_Terminal;
      const EntryType HELIXSequenceID_Terminal;
      const EntryType HELIXInsertionCode_Terminal;
      const EntryType HELIXClass;
      const EntryType HELIXComment;
      const EntryType HELIXLength;
      const EntryType SHEETStrandID;
      const EntryType SHEETSheetID;
      const EntryType SHEETNrOfStrandsInSheet;
      const EntryType SHEETResidueName_Initial;
      const EntryType SHEETChainID_Initial;
      const EntryType SHEETSequenceID_Initial;
      const EntryType SHEETInsertionCode_Initial;
      const EntryType SHEETResidueName_Terminal;
      const EntryType SHEETChainID_Terminal;
      const EntryType SHEETSequenceID_Terminal;
      const EntryType SHEETInsertionCode_Terminal;
      const EntryType SHEETSenseToPreviousStrand;
      const EntryType SHEETReg_AtomName_CurrentStrand;
      const EntryType SHEETReg_ResidueName_CurrentStrand;
      const EntryType SHEETReg_ChainID_CurrentStrand;
      const EntryType SHEETReg_SequenceID_CurrentStrand;
      const EntryType SHEETReg_InsertionCode_CurrStrand;
      const EntryType SHEETReg_AtomName_PrevStrand;
      const EntryType SHEETReg_ResidueName_PrevStrand;
      const EntryType SHEETReg_ChainID_PrevStrand;
      const EntryType SHEETReg_SequenceID_PrevStrand;
      const EntryType SHEETReg_InsertionCode_PrevStrand;
      const EntryType SITESequenceNumber;
      const EntryType SITEName;
      const EntryType SITENumberResidues;
      const EntryType SITEResidueName1;
      const EntryType SITEChainID1;
      const EntryType SITEResidueSequenceID1;
      const EntryType SITEResidueInsertionCode1;
      const EntryType SITEResidueName2;
      const EntryType SITEChainID2;
      const EntryType SITEResidueSequenceID2;
      const EntryType SITEResidueInsertionCode2;
      const EntryType SITEResidueName3;
      const EntryType SITEChainID3;
      const EntryType SITEResidueSequenceID3;
      const EntryType SITEResidueInsertionCode3;
      const EntryType SITEResidueName4;
      const EntryType SITEChainID4;
      const EntryType SITEResidueSequenceID4;
      const EntryType SITEResidueInsertionCode4;
      const EntryType MTRIX1_Serial;
      const EntryType MTRIX1_M11;
      const EntryType MTRIX1_M12;
      const EntryType MTRIX1_M13;
      const EntryType MTRIX1_V1;
      const EntryType MTRIX1_Given;
      const EntryType MTRIX2_Serial;
      const EntryType MTRIX2_M21;
      const EntryType MTRIX2_M22;
      const EntryType MTRIX2_M23;
      const EntryType MTRIX2_V2;
      const EntryType MTRIX2_Given;
      const EntryType MTRIX3_Serial;
      const EntryType MTRIX3_M31;
      const EntryType MTRIX3_M32;
      const EntryType MTRIX3_M33;
      const EntryType MTRIX3_V3;
      const EntryType MTRIX3_Given;
      const EntryType ATOMSerial;
      const EntryType ATOMName;
      const EntryType ATOMAlternateLocationID;
      const EntryType ATOMResidueName;
      const EntryType ATOMChainID;
      const EntryType ATOMResidueSequenceID;
      const EntryType ATOMInsertionCode;
      const EntryType ATOMX;
      const EntryType ATOMY;
      const EntryType ATOMZ;
      const EntryType ATOMOccupancy;
      const EntryType ATOMTempFactor;
      const EntryType ATOMSegmentID;
      const EntryType ATOMElement;
      const EntryType ATOMCharge;
      const EntryType ANISOUSerial;
      const EntryType ANISOUName;
      const EntryType ANISOUAlternateLocationID;
      const EntryType ANISOUResidueName;
      const EntryType ANISOUChainID;
      const EntryType ANISOUResidueSequenceID;
      const EntryType ANISOUInsertionCode;
      const EntryType ANISOUU11;
      const EntryType ANISOUU22;
      const EntryType ANISOUU33;
      const EntryType ANISOUU12;
      const EntryType ANISOUU13;
      const EntryType ANISOUU23;
      const EntryType ANISOUElement;
      const EntryType ANISOUCharge;
      const EntryType TERSerial;
      const EntryType TERResidueName;
      const EntryType TERChainID;
      const EntryType TERResidueSequenceID;
      const EntryType TERInsertionCode;
      const EntryType HETATMSerial;
      const EntryType HETATMName;
      const EntryType HETATMAlternateLocationID;
      const EntryType HETATMResidueName;
      const EntryType HETATMChainID;
      const EntryType HETATMResidueSequenceID;
      const EntryType HETATMInsertionCode;
      const EntryType HETATMX;
      const EntryType HETATMY;
      const EntryType HETATMZ;
      const EntryType HETATMOccupancy;
      const EntryType HETATMTempFactor;
      const EntryType HETATMSegmentID;
      const EntryType HETATMElement;
      const EntryType HETATMCharge;
      const EntryType MODELSerial;
      const EntryType CONECTAtomSerial;
      const EntryType CONECTBondedAtomSerial1;
      const EntryType CONECTBondedAtomSerial2;
      const EntryType CONECTBondedAtomSerial3;
      const EntryType CONECTBondedAtomSerial4;
      const EntryType MASTERNumRemark;
      const EntryType MASTER0;
      const EntryType MASTERNumHet;
      const EntryType MASTERNumHelix;
      const EntryType MASTERNumSheet;
      const EntryType MASTERNumTurn;
      const EntryType MASTERNumSite;
      const EntryType MASTERNumXform;
      const EntryType MASTERNumCoord;
      const EntryType MASTERNumTer;
      const EntryType MASTERNumConect;
      const EntryType MASTERNumSeqres;

      //! @brief double precision
      static const size_t s_DoublePrecision = 3;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all entry types
      EntryTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief add a new entry type
      //! @param DESCRIPTOR name of entry type
      //! @param ENTRY_DATA data for the entry type
      //! @return the created enum for that entry type
      EntryType AddEntryType( const std::string &DESCRIPTOR, const EntryTypeData &ENTRY_DATA);

    }; // EntryTypes

    //! @brief global access to all entry types
    BCL_API const EntryTypes &GetEntryTypes();

  } // namespace pdb

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< pdb::EntryTypeData, pdb::EntryTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_PDB_ENTRY_TYPES_H_
