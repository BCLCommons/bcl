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
#include "pdb/bcl_pdb_entry_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    //! @brief default constructor
    EntryTypes::EntryTypes() :
      HEADERClassification              ( AddEntryType( "HEADERClassification"              , EntryTypeData( GetLineTypes().HEADER, 10, 40, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      HEADERDate                        ( AddEntryType( "HEADERDate"                        , EntryTypeData( GetLineTypes().HEADER, 50,  9, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      HEADERIDCode                      ( AddEntryType( "HEADERIDCode"                      , EntryTypeData( GetLineTypes().HEADER, 62,  4, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      REMARK_Number                     ( AddEntryType( "REMARK_Number"                     , EntryTypeData( GetLineTypes().REMARK,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      REMARK_String                     ( AddEntryType( "REMARK_String"                     , EntryTypeData( GetLineTypes().REMARK, 11, 69, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      REMARK_350_BiomoleculeIdentifier  ( AddEntryType( "REMARK_350_BiomoleculeIdentifier"  , EntryTypeData( GetLineTypes().REMARK, 11, 12, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_BiomoleculeNumber      ( AddEntryType( "REMARK_350_BiomoleculeNumber"      , EntryTypeData( GetLineTypes().REMARK, 24,  3, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_SizeT ))),
      REMARK_350_ChainsIdentifier       ( AddEntryType( "REMARK_350_ChainsIdentifier"       , EntryTypeData( GetLineTypes().REMARK, 34,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_ChainsList             ( AddEntryType( "REMARK_350_ChainsList"             , EntryTypeData( GetLineTypes().REMARK, 42, 38, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_BIOMT1_Identifier      ( AddEntryType( "REMARK_350_BIOMT1_Identifier"      , EntryTypeData( GetLineTypes().REMARK, 13,  6, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_BIOMT1_Number          ( AddEntryType( "REMARK_350_BIOMT1_Number"          , EntryTypeData( GetLineTypes().REMARK, 21,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      REMARK_350_BIOMT1_M11             ( AddEntryType( "REMARK_350_BIOMT1_M11"             , EntryTypeData( GetLineTypes().REMARK, 23, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT1_M12             ( AddEntryType( "REMARK_350_BIOMT1_M12"             , EntryTypeData( GetLineTypes().REMARK, 33, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT1_M13             ( AddEntryType( "REMARK_350_BIOMT1_M13"             , EntryTypeData( GetLineTypes().REMARK, 43, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT1_V1              ( AddEntryType( "REMARK_350_BIOMT1_V1"              , EntryTypeData( GetLineTypes().REMARK, 57, 11,                             5,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT2_Identifier      ( AddEntryType( "REMARK_350_BIOMT2_Identifier"      , EntryTypeData( GetLineTypes().REMARK, 13,  6, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_BIOMT2_Number          ( AddEntryType( "REMARK_350_BIOMT2_Number"          , EntryTypeData( GetLineTypes().REMARK, 21,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      REMARK_350_BIOMT2_M21             ( AddEntryType( "REMARK_350_BIOMT2_M21"             , EntryTypeData( GetLineTypes().REMARK, 23, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT2_M22             ( AddEntryType( "REMARK_350_BIOMT2_M22"             , EntryTypeData( GetLineTypes().REMARK, 33, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT2_M23             ( AddEntryType( "REMARK_350_BIOMT2_M23"             , EntryTypeData( GetLineTypes().REMARK, 43, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT2_V2              ( AddEntryType( "REMARK_350_BIOMT2_V2"              , EntryTypeData( GetLineTypes().REMARK, 57, 11,                             5,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT3_Identifier      ( AddEntryType( "REMARK_350_BIOMT3_Identifier"      , EntryTypeData( GetLineTypes().REMARK, 13,  6, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_350_BIOMT3_Number          ( AddEntryType( "REMARK_350_BIOMT3_Number"          , EntryTypeData( GetLineTypes().REMARK, 21,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      REMARK_350_BIOMT3_M31             ( AddEntryType( "REMARK_350_BIOMT3_M31"             , EntryTypeData( GetLineTypes().REMARK, 23, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT3_M32             ( AddEntryType( "REMARK_350_BIOMT3_M32"             , EntryTypeData( GetLineTypes().REMARK, 33, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT3_M33             ( AddEntryType( "REMARK_350_BIOMT3_M33"             , EntryTypeData( GetLineTypes().REMARK, 43, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      REMARK_350_BIOMT3_V3              ( AddEntryType( "REMARK_350_BIOMT3_V3"              , EntryTypeData( GetLineTypes().REMARK, 57, 11,                             5,  true, util::CPPDataTypes::e_Double))),
      REMARK_465_ModelNumber            ( AddEntryType( "REMARK_465_ModelNumber"            , EntryTypeData( GetLineTypes().REMARK, 13,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      REMARK_465_ResidueName            ( AddEntryType( "REMARK_465_ResidueName"            , EntryTypeData( GetLineTypes().REMARK, 15,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_465_ChainID                ( AddEntryType( "REMARK_465_ChainID"                , EntryTypeData( GetLineTypes().REMARK, 19,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      REMARK_465_ResidueSequenceID      ( AddEntryType( "REMARK_465_ResidueSequenceID"      , EntryTypeData( GetLineTypes().REMARK, 21,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      REMARK_465_InsertionCode          ( AddEntryType( "REMARK_465_InsertionCode"          , EntryTypeData( GetLineTypes().REMARK, 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      REMARK_800_SiteIdentifierString   ( AddEntryType( "REMARK_800_SiteIdentifierString"   , EntryTypeData( GetLineTypes().REMARK, 11, 16, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_800_SiteIdentifier         ( AddEntryType( "REMARK_800_SiteIdentifier"         , EntryTypeData( GetLineTypes().REMARK, 28,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_800_SiteEvidenceCodeString ( AddEntryType( "REMARK_800_SiteEvidenceCodeString" , EntryTypeData( GetLineTypes().REMARK, 11, 14, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_800_SiteEvidenceCode       ( AddEntryType( "REMARK_800_SiteEvidenceCode"       , EntryTypeData( GetLineTypes().REMARK, 26,  8, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_800_SiteDescriptionString  ( AddEntryType( "REMARK_800_SiteDescriptionString"  , EntryTypeData( GetLineTypes().REMARK, 11, 17, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      REMARK_800_SiteDescription        ( AddEntryType( "REMARK_800_SiteDescription"        , EntryTypeData( GetLineTypes().REMARK, 29, 52, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESSerial                      ( AddEntryType( "SEQRESSerial"                      , EntryTypeData( GetLineTypes().SEQRES,  8,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SEQRESChainID                     ( AddEntryType( "SEQRESChainID"                     , EntryTypeData( GetLineTypes().SEQRES, 11,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SEQRESNrOfResiduesInChain         ( AddEntryType( "SEQRESNrOfResiduesInChain"         , EntryTypeData( GetLineTypes().SEQRES, 13,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SEQRESName_1                      ( AddEntryType( "SEQRESResName_1"                   , EntryTypeData( GetLineTypes().SEQRES, 19,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_2                      ( AddEntryType( "SEQRESResName_2"                   , EntryTypeData( GetLineTypes().SEQRES, 23,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_3                      ( AddEntryType( "SEQRESResName_3"                   , EntryTypeData( GetLineTypes().SEQRES, 27,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_4                      ( AddEntryType( "SEQRESResName_4"                   , EntryTypeData( GetLineTypes().SEQRES, 31,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_5                      ( AddEntryType( "SEQRESResName_5"                   , EntryTypeData( GetLineTypes().SEQRES, 35,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_6                      ( AddEntryType( "SEQRESResName_6"                   , EntryTypeData( GetLineTypes().SEQRES, 39,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_7                      ( AddEntryType( "SEQRESResName_7"                   , EntryTypeData( GetLineTypes().SEQRES, 43,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_8                      ( AddEntryType( "SEQRESResName_8"                   , EntryTypeData( GetLineTypes().SEQRES, 47,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_9                      ( AddEntryType( "SEQRESResName_9"                   , EntryTypeData( GetLineTypes().SEQRES, 51,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_10                     ( AddEntryType( "SEQRESResName_10"                  , EntryTypeData( GetLineTypes().SEQRES, 55,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_11                     ( AddEntryType( "SEQRESResName_11"                  , EntryTypeData( GetLineTypes().SEQRES, 59,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_12                     ( AddEntryType( "SEQRESResName_12"                  , EntryTypeData( GetLineTypes().SEQRES, 63,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SEQRESName_13                     ( AddEntryType( "SEQRESResName_13"                  , EntryTypeData( GetLineTypes().SEQRES, 67,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETIdentifier                     ( AddEntryType( "HETIdentifier"                     , EntryTypeData( GetLineTypes().HET   ,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETChainID                        ( AddEntryType( "HETChainID"                        , EntryTypeData( GetLineTypes().HET   , 12,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETSequenceID                     ( AddEntryType( "HETSequenceID"                     , EntryTypeData( GetLineTypes().HET   , 13,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      HETInsertionCode                  ( AddEntryType( "HETInsertionCode"                  , EntryTypeData( GetLineTypes().HET   , 17,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETNumberAtoms                    ( AddEntryType( "HETNumberAtoms"                    , EntryTypeData( GetLineTypes().HET   , 20,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HETDescription                    ( AddEntryType( "HETDescription"                    , EntryTypeData( GetLineTypes().HET   , 30, 40, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      HETNAMContinuation                ( AddEntryType( "HETNAMContinuation"                , EntryTypeData( GetLineTypes().HETNAM,  8,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HETNAMIdentifier                  ( AddEntryType( "HETNAMIdentifier"                  , EntryTypeData( GetLineTypes().HETNAM, 11,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETNAMText                        ( AddEntryType( "HETNAMText"                        , EntryTypeData( GetLineTypes().HETNAM, 15, 55, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      HETSYNContinuation                ( AddEntryType( "HETSYNContinuation"                , EntryTypeData( GetLineTypes().HETSYN,  8,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HETSYNIdentifier                  ( AddEntryType( "HETSYNIdentifier"                  , EntryTypeData( GetLineTypes().HETSYN, 11,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETSYNSynonyms                    ( AddEntryType( "HETSYNSynonyms"                    , EntryTypeData( GetLineTypes().HETSYN, 15, 55, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      FORMULComponentNumber             ( AddEntryType( "FORMULComponentNumber"             , EntryTypeData( GetLineTypes().FORMUL,  8,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      FORMULIdentifier                  ( AddEntryType( "FORMULIdentifier"                  , EntryTypeData( GetLineTypes().FORMUL, 12,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      FORMULContinuation                ( AddEntryType( "FORMULContinuation"                , EntryTypeData( GetLineTypes().FORMUL, 16,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      FORMULAsterisk                    ( AddEntryType( "FORMULAsterisk"                    , EntryTypeData( GetLineTypes().FORMUL, 18,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      FORMULChemicalFormula             ( AddEntryType( "FORMULChemicalFormula"             , EntryTypeData( GetLineTypes().FORMUL, 19, 51, util::GetUndefined< size_t>(), false, util::CPPDataTypes::e_String))),
      HELIXSerial                       ( AddEntryType( "HELIXSerial"                       , EntryTypeData( GetLineTypes().HELIX ,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HELIXID                           ( AddEntryType( "HELIXID"                           , EntryTypeData( GetLineTypes().HELIX , 11,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HELIXResidueName_Initial          ( AddEntryType( "HELIXResidueName_Initial"          , EntryTypeData( GetLineTypes().HELIX , 15,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HELIXChainID_Initial              ( AddEntryType( "HELIXChainID_Initial"              , EntryTypeData( GetLineTypes().HELIX , 19,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HELIXSequenceID_Initial           ( AddEntryType( "HELIXSequenceID_Initial"           , EntryTypeData( GetLineTypes().HELIX , 21,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      HELIXInsertionCode_Initial        ( AddEntryType( "HELIXInsertionCode_Initial"        , EntryTypeData( GetLineTypes().HELIX , 25,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HELIXResidueName_Terminal         ( AddEntryType( "HELIXResidueName_Terminal"         , EntryTypeData( GetLineTypes().HELIX , 27,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HELIXChainID_Terminal             ( AddEntryType( "HELIXChainID_Terminal"             , EntryTypeData( GetLineTypes().HELIX , 31,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HELIXSequenceID_Terminal          ( AddEntryType( "HELIXSequenceID_Terminal"          , EntryTypeData( GetLineTypes().HELIX , 33,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      HELIXInsertionCode_Terminal       ( AddEntryType( "HELIXInsertionCode_Terminal"       , EntryTypeData( GetLineTypes().HELIX , 37,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HELIXClass                        ( AddEntryType( "HELIXClass"                        , EntryTypeData( GetLineTypes().HELIX , 38,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HELIXComment                      ( AddEntryType( "HELIXComment"                      , EntryTypeData( GetLineTypes().HELIX , 40, 30, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HELIXLength                       ( AddEntryType( "HELIXLength"                       , EntryTypeData( GetLineTypes().HELIX , 71,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SHEETStrandID                     ( AddEntryType( "SHEETStrandID"                     , EntryTypeData( GetLineTypes().SHEET ,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SHEETSheetID                      ( AddEntryType( "SHEETSheetID"                      , EntryTypeData( GetLineTypes().SHEET , 11,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETNrOfStrandsInSheet           ( AddEntryType( "SHEETNrOfStrandsInSheet"           , EntryTypeData( GetLineTypes().SHEET , 14,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SHEETResidueName_Initial          ( AddEntryType( "SHEETResidueName_Initial"          , EntryTypeData( GetLineTypes().SHEET , 17,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETChainID_Initial              ( AddEntryType( "SHEETChainID_Initial"              , EntryTypeData( GetLineTypes().SHEET , 21,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETSequenceID_Initial           ( AddEntryType( "SHEETSequenceID_Initial"           , EntryTypeData( GetLineTypes().SHEET , 22,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SHEETInsertionCode_Initial        ( AddEntryType( "SHEETInsertionCode_Initial"        , EntryTypeData( GetLineTypes().SHEET , 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETResidueName_Terminal         ( AddEntryType( "SHEETResidueName_Terminal"         , EntryTypeData( GetLineTypes().SHEET , 28,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETChainID_Terminal             ( AddEntryType( "SHEETChainID_Terminal"             , EntryTypeData( GetLineTypes().SHEET , 32,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETSequenceID_Terminal          ( AddEntryType( "SHEETSequenceID_Terminal"          , EntryTypeData( GetLineTypes().SHEET , 33,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SHEETInsertionCode_Terminal       ( AddEntryType( "SHEETInsertionCode_Terminal"       , EntryTypeData( GetLineTypes().SHEET , 37,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETSenseToPreviousStrand        ( AddEntryType( "SHEETSenseToPreviousStrand"        , EntryTypeData( GetLineTypes().SHEET , 38,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SHEETReg_AtomName_CurrentStrand   ( AddEntryType( "SHEETReg_AtomName_CurrentStrand"   , EntryTypeData( GetLineTypes().SHEET , 41,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETReg_ResidueName_CurrentStrand( AddEntryType( "SHEETReg_ResidueName_CurrentStrand", EntryTypeData( GetLineTypes().SHEET , 45,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETReg_ChainID_CurrentStrand    ( AddEntryType( "SHEETReg_ChainID_CurrentStrand"    , EntryTypeData( GetLineTypes().SHEET , 49,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETReg_SequenceID_CurrentStrand ( AddEntryType( "SHEETReg_SequenceID_CurrentStrand" , EntryTypeData( GetLineTypes().SHEET , 50,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SHEETReg_InsertionCode_CurrStrand ( AddEntryType( "SHEETReg_InsertionCode_CurrStrand" , EntryTypeData( GetLineTypes().SHEET , 54,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETReg_AtomName_PrevStrand      ( AddEntryType( "SHEETReg_AtomName_PrevStrand"      , EntryTypeData( GetLineTypes().SHEET , 56,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETReg_ResidueName_PrevStrand   ( AddEntryType( "SHEETReg_ResidueName_PrevStrand"   , EntryTypeData( GetLineTypes().SHEET , 60,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SHEETReg_ChainID_PrevStrand       ( AddEntryType( "SHEETReg_ChainID_PrevStrand"       , EntryTypeData( GetLineTypes().SHEET , 64,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SHEETReg_SequenceID_PrevStrand    ( AddEntryType( "SHEETReg_SequenceID_PrevStrand"    , EntryTypeData( GetLineTypes().SHEET , 65,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SHEETReg_InsertionCode_PrevStrand ( AddEntryType( "SHEETReg_InsertionCode_PrevStrand" , EntryTypeData( GetLineTypes().SHEET , 69,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITESequenceNumber                ( AddEntryType( "SITESequenceNumber"                , EntryTypeData( GetLineTypes().SITE  ,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SITEName                          ( AddEntryType( "SITEName"                          , EntryTypeData( GetLineTypes().SITE  , 11,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SITENumberResidues                ( AddEntryType( "SITENumberResidues"                , EntryTypeData( GetLineTypes().SITE  , 15,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      SITEResidueName1                  ( AddEntryType( "SITEResidueName1"                  , EntryTypeData( GetLineTypes().SITE  , 18,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SITEChainID1                      ( AddEntryType( "SITEChainID1"                      , EntryTypeData( GetLineTypes().SITE  , 22,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueSequenceID1            ( AddEntryType( "SITEResidueSequenceID1"            , EntryTypeData( GetLineTypes().SITE  , 23,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SITEResidueInsertionCode1         ( AddEntryType( "SITEResidueInsertionCode1"         , EntryTypeData( GetLineTypes().SITE  , 27,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueName2                  ( AddEntryType( "SITEResidueName2"                  , EntryTypeData( GetLineTypes().SITE  , 29,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SITEChainID2                      ( AddEntryType( "SITEChainID2"                      , EntryTypeData( GetLineTypes().SITE  , 33,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueSequenceID2            ( AddEntryType( "SITEResidueSequenceID2"            , EntryTypeData( GetLineTypes().SITE  , 34,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SITEResidueInsertionCode2         ( AddEntryType( "SITEResidueInsertionCode2"         , EntryTypeData( GetLineTypes().SITE  , 38,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueName3                  ( AddEntryType( "SITEResidueName3"                  , EntryTypeData( GetLineTypes().SITE  , 40,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SITEChainID3                      ( AddEntryType( "SITEChainID3"                      , EntryTypeData( GetLineTypes().SITE  , 44,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueSequenceID3            ( AddEntryType( "SITEResidueSequenceID3"            , EntryTypeData( GetLineTypes().SITE  , 45,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SITEResidueInsertionCode3         ( AddEntryType( "SITEResidueInsertionCode3"         , EntryTypeData( GetLineTypes().SITE  , 49,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueName4                  ( AddEntryType( "SITEResidueName4"                  , EntryTypeData( GetLineTypes().SITE  , 51,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      SITEChainID4                      ( AddEntryType( "SITEChainID4"                      , EntryTypeData( GetLineTypes().SITE  , 55,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      SITEResidueSequenceID4            ( AddEntryType( "SITEResidueSequenceID4"            , EntryTypeData( GetLineTypes().SITE  , 56,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      SITEResidueInsertionCode4         ( AddEntryType( "SITEResidueInsertionCode4"         , EntryTypeData( GetLineTypes().SITE  , 60,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      MTRIX1_Serial                     ( AddEntryType( "MTRIX1_Serial"                     , EntryTypeData( GetLineTypes().MTRIX1,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MTRIX1_M11                        ( AddEntryType( "MTRIX1_M11"                        , EntryTypeData( GetLineTypes().MTRIX1, 10, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX1_M12                        ( AddEntryType( "MTRIX1_M12"                        , EntryTypeData( GetLineTypes().MTRIX1, 20, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX1_M13                        ( AddEntryType( "MTRIX1_M13"                        , EntryTypeData( GetLineTypes().MTRIX1, 30, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX1_V1                         ( AddEntryType( "MTRIX1_V1"                         , EntryTypeData( GetLineTypes().MTRIX1, 45, 10,                             5,  true, util::CPPDataTypes::e_Double))),
      MTRIX1_Given                      ( AddEntryType( "MTRIX1_Given"                      , EntryTypeData( GetLineTypes().MTRIX1, 59,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Bool  ))),
      MTRIX2_Serial                     ( AddEntryType( "MTRIX2_Serial"                     , EntryTypeData( GetLineTypes().MTRIX2,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MTRIX2_M21                        ( AddEntryType( "MTRIX2_M21"                        , EntryTypeData( GetLineTypes().MTRIX2, 10, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX2_M22                        ( AddEntryType( "MTRIX2_M22"                        , EntryTypeData( GetLineTypes().MTRIX2, 20, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX2_M23                        ( AddEntryType( "MTRIX2_M23"                        , EntryTypeData( GetLineTypes().MTRIX2, 30, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX2_V2                         ( AddEntryType( "MTRIX2_V2"                         , EntryTypeData( GetLineTypes().MTRIX2, 45, 10,                             5,  true, util::CPPDataTypes::e_Double))),
      MTRIX2_Given                      ( AddEntryType( "MTRIX2_Given"                      , EntryTypeData( GetLineTypes().MTRIX2, 59,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Bool  ))),
      MTRIX3_Serial                     ( AddEntryType( "MTRIX3_Serial"                     , EntryTypeData( GetLineTypes().MTRIX3,  7,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MTRIX3_M31                        ( AddEntryType( "MTRIX3_M31"                        , EntryTypeData( GetLineTypes().MTRIX3, 10, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX3_M32                        ( AddEntryType( "MTRIX3_M32"                        , EntryTypeData( GetLineTypes().MTRIX3, 20, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX3_M33                        ( AddEntryType( "MTRIX3_M33"                        , EntryTypeData( GetLineTypes().MTRIX3, 30, 10,                             6,  true, util::CPPDataTypes::e_Double))),
      MTRIX3_V3                         ( AddEntryType( "MTRIX3_V3"                         , EntryTypeData( GetLineTypes().MTRIX3, 45, 10,                             5,  true, util::CPPDataTypes::e_Double))),
      MTRIX3_Given                      ( AddEntryType( "MTRIX3_Given"                      , EntryTypeData( GetLineTypes().MTRIX3, 59,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Bool  ))),
      ATOMSerial                        ( AddEntryType( "ATOMSerial"                        , EntryTypeData( GetLineTypes().ATOM  ,  6,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      ATOMName                          ( AddEntryType( "ATOMName"                          , EntryTypeData( GetLineTypes().ATOM  , 12,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ATOMAlternateLocationID           ( AddEntryType( "ATOMAlternateLocationID"           , EntryTypeData( GetLineTypes().ATOM  , 16,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ATOMResidueName                   ( AddEntryType( "ATOMResidueName"                   , EntryTypeData( GetLineTypes().ATOM  , 17,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ATOMChainID                       ( AddEntryType( "ATOMChainID"                       , EntryTypeData( GetLineTypes().ATOM  , 21,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ATOMResidueSequenceID             ( AddEntryType( "ATOMResidueSequenceID"             , EntryTypeData( GetLineTypes().ATOM  , 22,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ATOMInsertionCode                 ( AddEntryType( "ATOMInsertionCode"                 , EntryTypeData( GetLineTypes().ATOM  , 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ATOMX                             ( AddEntryType( "ATOMX"                             , EntryTypeData( GetLineTypes().ATOM  , 30,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      ATOMY                             ( AddEntryType( "ATOMY"                             , EntryTypeData( GetLineTypes().ATOM  , 38,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      ATOMZ                             ( AddEntryType( "ATOMZ"                             , EntryTypeData( GetLineTypes().ATOM  , 46,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      ATOMOccupancy                     ( AddEntryType( "ATOMOccupancy"                     , EntryTypeData( GetLineTypes().ATOM  , 54,  6,                             2,  true, util::CPPDataTypes::e_Double))),
      ATOMTempFactor                    ( AddEntryType( "ATOMTempFactor"                    , EntryTypeData( GetLineTypes().ATOM  , 60,  6,                             2,  true, util::CPPDataTypes::e_Double))),
      ATOMSegmentID                     ( AddEntryType( "ATOMSegmentID"                     , EntryTypeData( GetLineTypes().ATOM  , 72,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ATOMElement                       ( AddEntryType( "ATOMElement"                       , EntryTypeData( GetLineTypes().ATOM  , 76,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ATOMCharge                        ( AddEntryType( "ATOMCharge"                        , EntryTypeData( GetLineTypes().ATOM  , 78,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ANISOUSerial                      ( AddEntryType( "ANISOUSerial"                      , EntryTypeData( GetLineTypes().ANISOU,  6,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      ANISOUName                        ( AddEntryType( "ANISOUName"                        , EntryTypeData( GetLineTypes().ANISOU, 12,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ANISOUAlternateLocationID         ( AddEntryType( "ANISOUAlternateLocationID"         , EntryTypeData( GetLineTypes().ANISOU, 16,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ANISOUResidueName                 ( AddEntryType( "ANISOUResidueName"                 , EntryTypeData( GetLineTypes().ANISOU, 17,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ANISOUChainID                     ( AddEntryType( "ANISOUChainID"                     , EntryTypeData( GetLineTypes().ANISOU, 21,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ANISOUResidueSequenceID           ( AddEntryType( "ANISOUResidueSequenceID"           , EntryTypeData( GetLineTypes().ANISOU, 22,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUInsertionCode               ( AddEntryType( "ANISOUInsertionCode"               , EntryTypeData( GetLineTypes().ANISOU, 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      ANISOUU11                         ( AddEntryType( "ANISOUU11"                         , EntryTypeData( GetLineTypes().ANISOU, 28,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUU22                         ( AddEntryType( "ANISOUU22"                         , EntryTypeData( GetLineTypes().ANISOU, 35,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUU33                         ( AddEntryType( "ANISOUU33"                         , EntryTypeData( GetLineTypes().ANISOU, 42,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUU12                         ( AddEntryType( "ANISOUU12"                         , EntryTypeData( GetLineTypes().ANISOU, 49,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUU13                         ( AddEntryType( "ANISOUU13"                         , EntryTypeData( GetLineTypes().ANISOU, 56,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUU23                         ( AddEntryType( "ANISOUU23"                         , EntryTypeData( GetLineTypes().ANISOU, 63,  7, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      ANISOUElement                     ( AddEntryType( "ANISOUElement"                     , EntryTypeData( GetLineTypes().ANISOU, 76,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      ANISOUCharge                      ( AddEntryType( "ANISOUCharge"                      , EntryTypeData( GetLineTypes().ANISOU, 78,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      TERSerial                         ( AddEntryType( "TER_Serial"                        , EntryTypeData( GetLineTypes().TER   ,  6,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      TERResidueName                    ( AddEntryType( "TER_ResidueName"                   , EntryTypeData( GetLineTypes().TER   , 17,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      TERChainID                        ( AddEntryType( "TER_ChainID"                       , EntryTypeData( GetLineTypes().TER   , 21,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      TERResidueSequenceID              ( AddEntryType( "TER_ResidueSequenceID"             , EntryTypeData( GetLineTypes().TER   , 22,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      TERInsertionCode                  ( AddEntryType( "TER_InsertionCode"                 , EntryTypeData( GetLineTypes().TER   , 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETATMSerial                      ( AddEntryType( "HETATMSerial"                      , EntryTypeData( GetLineTypes().HETATM,  6,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      HETATMName                        ( AddEntryType( "HETATMName"                        , EntryTypeData( GetLineTypes().HETATM, 12,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETATMAlternateLocationID         ( AddEntryType( "HETATMAlternateLocationID"         , EntryTypeData( GetLineTypes().HETATM, 16,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETATMResidueName                 ( AddEntryType( "HETATMResidueName"                 , EntryTypeData( GetLineTypes().HETATM, 17,  3, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETATMChainID                     ( AddEntryType( "HETATMChainID"                     , EntryTypeData( GetLineTypes().HETATM, 21,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETATMResidueSequenceID           ( AddEntryType( "HETATMResidueSequenceID"           , EntryTypeData( GetLineTypes().HETATM, 22,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      HETATMInsertionCode               ( AddEntryType( "HETATMInsertionCode"               , EntryTypeData( GetLineTypes().HETATM, 26,  1, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Char  ))),
      HETATMX                           ( AddEntryType( "HETATMX"                           , EntryTypeData( GetLineTypes().HETATM, 30,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      HETATMY                           ( AddEntryType( "HETATMY"                           , EntryTypeData( GetLineTypes().HETATM, 38,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      HETATMZ                           ( AddEntryType( "HETATMZ"                           , EntryTypeData( GetLineTypes().HETATM, 46,  8,                             3,  true, util::CPPDataTypes::e_Double))),
      HETATMOccupancy                   ( AddEntryType( "HETATMOccupancy"                   , EntryTypeData( GetLineTypes().HETATM, 54,  6,                             2,  true, util::CPPDataTypes::e_Double))),
      HETATMTempFactor                  ( AddEntryType( "HETATMTempFactor"                  , EntryTypeData( GetLineTypes().HETATM, 60,  6,                             2,  true, util::CPPDataTypes::e_Double))),
      HETATMSegmentID                   ( AddEntryType( "HETATMSegmentID"                   , EntryTypeData( GetLineTypes().HETATM, 72,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETATMElement                     ( AddEntryType( "HETATMElement"                     , EntryTypeData( GetLineTypes().HETATM, 76,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      HETATMCharge                      ( AddEntryType( "HETATMCharge"                      , EntryTypeData( GetLineTypes().HETATM, 78,  2, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_String))),
      MODELSerial                       ( AddEntryType( "MODELSerial"                       , EntryTypeData( GetLineTypes().MODEL , 10,  4, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      CONECTAtomSerial                  ( AddEntryType( "CONECTAtomSerial"                  , EntryTypeData( GetLineTypes().CONECT,  6,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      CONECTBondedAtomSerial1           ( AddEntryType( "CONECTBondedAtomSerial1"           , EntryTypeData( GetLineTypes().CONECT, 11,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      CONECTBondedAtomSerial2           ( AddEntryType( "CONECTBondedAtomSerial2"           , EntryTypeData( GetLineTypes().CONECT, 16,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      CONECTBondedAtomSerial3           ( AddEntryType( "CONECTBondedAtomSerial3"           , EntryTypeData( GetLineTypes().CONECT, 21,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      CONECTBondedAtomSerial4           ( AddEntryType( "CONECTBondedAtomSerial4"           , EntryTypeData( GetLineTypes().CONECT, 26,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_Int   ))),
      MASTERNumRemark                   ( AddEntryType( "MASTERNumRemark"                   , EntryTypeData( GetLineTypes().MASTER, 10,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTER0                           ( AddEntryType( "MASTER0"                           , EntryTypeData( GetLineTypes().MASTER, 15,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumHet                      ( AddEntryType( "MASTERNumHet"                      , EntryTypeData( GetLineTypes().MASTER, 20,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumHelix                    ( AddEntryType( "MASTERNumHelix"                    , EntryTypeData( GetLineTypes().MASTER, 25,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumSheet                    ( AddEntryType( "MASTERNumSheet"                    , EntryTypeData( GetLineTypes().MASTER, 30,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumTurn                     ( AddEntryType( "MASTERNumTurn"                     , EntryTypeData( GetLineTypes().MASTER, 35,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumSite                     ( AddEntryType( "MASTERNumSite"                     , EntryTypeData( GetLineTypes().MASTER, 40,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumXform                    ( AddEntryType( "MASTERNumXform"                    , EntryTypeData( GetLineTypes().MASTER, 45,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumCoord                    ( AddEntryType( "MASTERNumCoord"                    , EntryTypeData( GetLineTypes().MASTER, 50,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumTer                      ( AddEntryType( "MASTERNumTer"                      , EntryTypeData( GetLineTypes().MASTER, 55,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumConect                   ( AddEntryType( "MASTERNumConect"                   , EntryTypeData( GetLineTypes().MASTER, 60,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT ))),
      MASTERNumSeqres                   ( AddEntryType( "MASTERNumSeqres"                   , EntryTypeData( GetLineTypes().MASTER, 65,  5, util::GetUndefined< size_t>(),  true, util::CPPDataTypes::e_SizeT )))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EntryTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief add a new entry type
    //! @param DESCRIPTOR name of entry type
    //! @param ENTRY_DATA data for the entry type
    //! @return the created enum for that entry type
    EntryType EntryTypes::AddEntryType( const std::string &DESCRIPTOR, const EntryTypeData &ENTRY_DATA)
    {
      // increment number entries for that line type
      const EntryType new_enum( AddEnum( DESCRIPTOR, ENTRY_DATA));
      LineTypes::GetEnums().ConsiderEntryType( new_enum->GetLineType(), new_enum);
      return new_enum;
    }

    //! @brief global access to all entry types
    const EntryTypes &GetEntryTypes()
    {
      return EntryTypes::GetEnums();
    }

  } // namespace pdb

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< pdb::EntryTypeData, pdb::EntryTypes>;

  } // namespace util
} // namespace bcl
