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
#include "pdb/bcl_pdb.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief file extension for pdb file
    //! @return the file extension of a pdb file
    const std::string &GetDefaultFileExtension()
    {
      static const std::string s_file_extension( "pdb");
      return s_file_extension;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_entry_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EntryTypeData::EntryTypeData() :
      m_LineType( GetLineTypes().e_Undefined),
      m_Start( util::GetUndefined< size_t>()),
      m_Length( util::GetUndefined< size_t>()),
      m_DataType( util::CPPDataTypes::e_Unknown),
      m_IsNumerical( false)
    {
    }

    //! @brief construct from information about entry
    //! @param LINE_TYPE
    //! @param START
    //! @param LENGTH
    //! @param PRECISION
    //! @param RIGHT_ALIGNED
    //! @param DATA_TYPE
    EntryTypeData::EntryTypeData
    (
      const LineType &LINE_TYPE,
      const size_t START,
      const size_t LENGTH,
      const size_t PRECISION,
      const bool RIGHT_ALIGNED,
      const util::CPPDataTypes::Types DATA_TYPE
    ) :
      m_LineType( LINE_TYPE),
      m_Start( START),
      m_Length( LENGTH),
      m_Precision( PRECISION),
      m_RightAligned( RIGHT_ALIGNED),
      m_DataType( DATA_TYPE),
      m_IsNumerical
      (
        m_DataType == util::CPPDataTypes::e_Double ||
        m_DataType == util::CPPDataTypes::e_Float ||
        m_DataType == util::CPPDataTypes::e_Int ||
        m_DataType == util::CPPDataTypes::e_SizeT
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new EntryTypeData
    EntryTypeData *EntryTypeData::Clone() const
    {
      return new EntryTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &EntryTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief construct a format object
    util::Format EntryTypeData::GetFormat() const
    {
      // initialize format object
      util::Format format;

      // set width
      format.W( m_Length);

      // set precision for double or float types
      if( util::IsDefined( m_Precision))
      {
        format.FFP( m_Precision);
      }

      // set alignment
      m_RightAligned ? format.R() : format.L();

      // return
      return format;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EntryTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_LineType    , ISTREAM);
      io::Serialize::Read( m_Start       , ISTREAM);
      io::Serialize::Read( m_Length      , ISTREAM);
      io::Serialize::Read( m_Precision   , ISTREAM);
      io::Serialize::Read( m_RightAligned, ISTREAM);
      io::Serialize::Read( m_DataType    , ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &EntryTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_LineType    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Start       , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Length      , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_Precision   , OSTREAM,      0) << '\t';
      io::Serialize::Write( m_RightAligned, OSTREAM,      0) << '\t';
      io::Serialize::Write( m_DataType    , OSTREAM,      0);

      // return
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_dssp.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_printer_biomatrix.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! static bool for whether all command line defaults have been saved
    bool Factory::s_HaveDefaults( Factory::ResetFlagDefaults());

    //! Flag to switch between numerated AtomIDs and the original PDBAtomIDs
    //! - default will numerate the atomIDs as they are written
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagWritePDBAtomID()
    {
      static util::ShPtr< command::FlagInterface> s_flag_write_pdb_atom_id
      (
        new command::FlagStatic
        (
          "write_pdb_atom_ids",
          "switch between numerated AtomIDs and the original PDBAtomIDs - default will numerate the atom IDs as they are written"
        )
      );

      return s_flag_write_pdb_atom_id;
    }

    //! Flag to switch between numerated ResIDs( SeqIDs) and the original PDBResIDs and PDBICodes
    //! default will write SeqIDs
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagWritePDBResID()
    {
      static util::ShPtr< command::FlagInterface> s_flag_write_pdb_res_id
      (
        new command::FlagStatic
        (
          "write_pdb_res_ids",
          "switch between numerated ResidueIDs( SeqIDs) and the original PDBResIDs and PDBICodes - if flag is not set, sequence IDs (starting with 1 for first residue) will be written"
        )
      );

      return s_flag_write_pdb_res_id;
    }

    //! command line flag to be used to change the Defaut AAClass that is used to generate AASequences/Chains/Proteins
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagAAClass()
    {
      static util::ShPtr< command::FlagInterface> s_flag_aa_class
      (
        new command::FlagStatic
        (
          "aaclass",
          "choose the amino acid class that will be written",
          command::Parameter
          (
            "amino_acid_class",
            "choice of amino acid class that is different by the amino acids used",
            command::ParameterCheckEnumerate< biol::AAClasses>(),
            biol::GetAAClasses().e_AABackBone.GetName()
          )
        )
      );

      return s_flag_aa_class;
    }

    //! command line flag to be used to change the default minimal size of the sse types tp be included in protein models
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagMinSSESize()
    {
      static util::ShPtr< command::ParameterInterface> s_min_helix_size_parameter
      (
        new command::Parameter
        (
          "min_helix_size",
          "minimal size of helices to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::ParameterInterface> s_min_strand_size_parameter
      (
        new command::Parameter
        (
          "min_strand_size",
          "minimal size of strand to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::ParameterInterface> s_min_loop_size_parameter
      (
        new command::Parameter
        (
          "min_loop_size",
          "minimal size of loops to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::FlagInterface> s_min_sse_size_flag
      (
        new command::FlagStatic
        (
          "min_sse_size",
          "change the default minimal size of each secondary structure type that is to be included when protein models "
          "are constructed from a pdb file"
        )
      );

      // insert parameters if it was constructed for the first time
      if( s_min_sse_size_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> min_sse_size_flag( s_min_sse_size_flag);
        min_sse_size_flag->PushBack( s_min_helix_size_parameter);
        min_sse_size_flag->PushBack( s_min_strand_size_parameter);
        min_sse_size_flag->PushBack( s_min_loop_size_parameter);
      }

      // end
      return s_min_sse_size_flag;
    }

    //! command line flag to be used to allow conversion of unnatural aminoacids to their parent amino acid
    util::ShPtr< command::FlagInterface> &Factory::GetFlagConvertToNaturalAAType()
    {
      static util::ShPtr< command::FlagInterface> s_flag_convert_to_natural_aa_type
      (
        new command::FlagStatic
        (
          "convert_to_natural_aa_type",
          "when creating amino acids from a pdb, unnatural amino acid types are converted to their parent amino acid: "
          "e.g. MSE -> MET"
        )
      );

      return s_flag_convert_to_natural_aa_type;
    }

    //! command line flag to be used, if pdb does not have any sse definition to use the backbone conformation to determine that
    util::ShPtr< command::FlagInterface> &Factory::GetFlagSSEsFromBackBone()
    {
      static util::ShPtr< command::FlagInterface> s_sse_from_back_bone
      (
        new command::FlagStatic
        (
          "sse_from_backbone",
          "only if no pdb SSE definitions are given for certain chain within pdb, use the backbone conformation to calculate them"
        )
      );

      //end
      return s_sse_from_back_bone;
    }

    //! command line flag to be used, if pdb SSE definitions should be reassigned by DSSP
    util::ShPtr< command::FlagInterface> &Factory::GetFlagDSSP()
    {
      static util::ShPtr< command::FlagInterface> s_sse_from_dssp
      (
        new command::FlagStatic
        (
          "sse_from_dssp",
          "reassign pdb SSE definitions using dssp",
          command::Parameter
          (
            "max_hbond_energy",
            "maximal H-Bond energy to be considered an H-Bond",
            command::ParameterCheckRanged< double>( -10, 0),
            util::Format()( biol::DSSP::GetDefaultHBondMaxEnergy())
          )
        )
      );

      //end
      return s_sse_from_dssp;
    }

    //! command line flag to be used to write zero coordinates for undefined amino acids
    util::ShPtr< command::FlagInterface> &Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()
    {
      static util::ShPtr< command::FlagInterface> s_write_zero_coordinates
      (
        new command::FlagStatic
        (
          "write_zero_coordinates",
          "If amino acids are undefined, write them to the pdb with zero coordinates"
        )
      );

      //end
      return s_write_zero_coordinates;
    }

    // command line flag to be used to write hydrogen atoms if available
    util::ShPtr< command::FlagInterface> &Factory::GetFlagWriteHydrogens()
    {
      static util::ShPtr< command::FlagInterface> s_write_hdydrogens
      (
        new command::FlagStatic
        (
          "write_hydrogens",
          "output hydrogen atoms if available"
        )
      );

      //end
      return s_write_hdydrogens;
    }

    //! command line flag to switch between IUPAC (standard) and pdb atom name output (when flag is used)
    util::ShPtr< command::FlagInterface> &Factory::GetFlagPDBAtomName()
    {
      static util::ShPtr< command::FlagInterface> s_write_pdb_atom_names
      (
        new command::FlagStatic
        (
          "write_pdb_atom_names",
          "instead of IUPAC atom names, pdb atom nameing will be used (e.g. HH21 -> 1HH2 for ARGININE)"
        )
      );

      //end
      return s_write_pdb_atom_names;
    }

    //! command line flag for specifying which biomolecule number to read in for generating multimers
    util::ShPtr< command::FlagInterface> &Factory::GetFlagBiomolecule()
    {
      static util::ShPtr< command::FlagInterface> s_flag_biomolecule
      (
        new command::FlagStatic
        (
          "biomolecule",
          "use the REMARK 350  BIOMT entries to generate the biomolecule",
          command::Parameter
          (
            "BIOMOLECULE",
            "biomolecule number",
            command::ParameterCheckRanged< size_t>( 1, 99),
            "1"
          )
        )
      );

      // end
      return s_flag_biomolecule;
    }

    //! function to return all flags for the pdb factory
    util::ShPtrVector< command::FlagInterface> &Factory::GetAllFlags()
    {
      static util::ShPtrVector< command::FlagInterface> s_all_flags
      (
        util::ShPtrVector< command::FlagInterface>::Create
        (
          GetFlagWritePDBAtomID(),
          GetFlagWritePDBResID(),
          GetFlagAAClass(),
          GetFlagMinSSESize(),
          GetFlagSSEsFromBackBone(),
          GetFlagDSSP(),
          GetFlagConvertToNaturalAAType(),
          GetFlagWriteZeroCoordinatesForUndefinedAminoAcids(),
          GetFlagWriteHydrogens(),
          GetFlagPDBAtomName(),
          GetFlagBiomolecule()
        )
      );
      return s_all_flags;
    }

    //! @brief reset defaults on all flags; to ensure a consistent set of defaults
    //! This function should be called whenever a given example or app may depend on the normal default parameters
    //! @return true on success
    //! @note should be called during static initialization to acquire the actual defaults
    bool Factory::ResetFlagDefaults()
    {
      static util::ShPtrVector< command::FlagInterface> &s_flags( GetAllFlags());

      // record of the original defaults for all flags
      static storage::Vector< storage::Vector< std::string> > s_original_defaults;

      if( s_original_defaults.IsEmpty())
      {
        // fill up the vector with the original defaults
        for
        (
          util::ShPtrVector< command::FlagInterface>::const_iterator
            itr( s_flags.Begin()), itr_end( s_flags.End());
          itr != itr_end;
          ++itr
        )
        {
          s_original_defaults.PushBack( storage::Vector< std::string>());
          for
          (
            util::ShPtrVector< command::ParameterInterface>::const_iterator
              itr_param( ( *itr)->GetParameterList().Begin()), itr_param_end( ( *itr)->GetParameterList().End());
            itr_param != itr_param_end;
            ++itr_param
          )
          {
            // store the default value
            s_original_defaults.LastElement().PushBack( ( *itr_param)->GetDefaultValue());
          }
        }
      }
      else
      {
        // perform the reset on all the flags and reset defaults on parameters
        storage::Vector< storage::Vector< std::string> >::const_iterator itr_defaults( s_original_defaults.Begin());
        for
        (
          util::ShPtrVector< command::FlagInterface>::iterator
            itr( s_flags.Begin()), itr_end( s_flags.End());
          itr != itr_end;
          ++itr, ++itr_defaults
        )
        {
          // reset the flag itself
          ( *itr)->ResetFlag();

          // reset the defaults
          storage::Vector< std::string>::const_iterator itr_default( itr_defaults->Begin());
          for
          (
            util::ShPtrVector< command::ParameterInterface>::iterator
              itr_param( ( *itr)->GetParameterList().Begin()), itr_param_end( ( *itr)->GetParameterList().End());
            itr_param != itr_param_end;
            ++itr_param, ++itr_default
          )
          {
            // test whether there is a default value
            if( ( *itr_param)->GetWasDefaultGiven())
            {
              // yes, so set the default back to its original value
              ( *itr_param)->SetDefaultParameter( *itr_default);
            }
          }
        }
      }
      return true;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Factory::Factory() :
      m_AAClass( GetFlagAAClass()->GetFirstParameter()->GetValue()),
      m_Printers
      (
        1,
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      )
    {
    }

    //! @brief constructor from AAClass
    //! @param AA_CLASS AAClass of interest
    Factory::Factory( const biol::AAClass &AA_CLASS) :
      m_AAClass( AA_CLASS),
      m_Printers
      (
        1,
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      )
    {
    }

    //! @brief virtual copy constructor
    Factory *Factory::Clone() const
    {
      return new Factory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Factory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief resets the printers
    void Factory::ResetPrinters()
    {
      m_Printers.Reset();
      m_Printers.PushBack
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      );
    }

    //! @brief appends printer to list
    //! @param SP_PRINTER printer to be appended
    void Factory::AppendPrinter
    (
      const util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > > &SP_PRINTER
    )
    {
      m_Printers.Append( SP_PRINTER);
    }

  ///////////////////////
  // operations - read //
  ///////////////////////

    //! @brief builds biol::Atom from Line
    //! @param ATOM_TYPE type of atom in line
    //! @param LINE Line which contains Atom information
    //! @return ShPtr to Atom read from Line
    util::ShPtr< biol::Atom>
    Factory::AtomFromLine
    (
      const biol::AtomType &ATOM_TYPE,
      const Line &LINE
    )
    {
      // return util::ShPtr for biol::Atom, constructed from Atom Name, PDB Atom Serial and its Position
      return util::ShPtr< biol::Atom>
      (
        new biol::Atom
        (
          LINE.RetrieveCoordinates(),
          ATOM_TYPE,
          LINE.GetNumericalValue< size_t>( GetEntryTypes().ATOMSerial),
          LINE.GetNumericalValue< double>( GetEntryTypes().ATOMTempFactor)
        )
      );
    }

    //! @brief builds biol::Atom from Line
    //! @param LINE Line which contains Atom information
    //! @return ShPtr to Atom read from Line
    util::ShPtr< biol::Atom>
    Factory::AtomFromLine( const Line &LINE)
    {
      const biol::AtomType type( biol::GetAtomTypes().TypeFromPDBAtomName( util::TrimString( LINE.GetString( GetEntryTypes().ATOMName))));

      // if the AtomName in the Line is undefined
      if( !type.IsDefined())
      {
        BCL_MessageStd
        (
          " This atom name leads to undefined biol::AtomType |" +
            util::TrimString( LINE.GetString( GetEntryTypes().ATOMName)) + "|"
        );
      }

      // return util::ShPtr for biol::Atom, constructed from Atom Name, PDB Atom Serial and its Position
      return AtomFromLine( type, LINE);
    }

    //! @brief builds all Atoms for provided AtomTypes from Lines
    //! @param LINES List of Lines which contains Atom information
    //! @param AA_TYPE aa type for which atoms are desired
    //! @return ShPtrVector of Atoms read from given Lines
    util::ShPtrVector< biol::Atom>
    Factory::AtomsFromLines
    (
      const util::ShPtrList< Line> &LINES,
      const biol::AAType &AA_TYPE,
      const biol::AAType &AA_TYPE_IN_FILE
    ) const
    {
      util::ShPtrVector< biol::Atom> atoms;

      if( LINES.IsEmpty())
      {
        return atoms;
      }

      storage::Set< biol::AtomType> missing_types( AA_TYPE->GetAllowedAtomTypes());

      // find extra atom types that the aa type in the file may have (for unnatural aa types) that will not be
      // represented in the bcl
      storage::Set< biol::AtomType> extra_types;
      if( AA_TYPE != AA_TYPE_IN_FILE)
      {
        std::set_difference
        (
          AA_TYPE_IN_FILE->GetAllowedAtomTypes().Begin(),
          AA_TYPE_IN_FILE->GetAllowedAtomTypes().End(),
          missing_types.Begin(),
          missing_types.End(),
          std::inserter( extra_types.InternalData(), extra_types.Begin())
        );
      }

      // iterate over atoms in residue lines to match this type
      for
      (
        util::ShPtrList< Line>::const_iterator line_itr( LINES.Begin()), line_itr_end( LINES.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // get the atom type from the pdb atom name
        const biol::AtomType &current_atom_type
        (
          AA_TYPE->GetAtomTypeFromAtomName( ( *line_itr)->GetString( GetEntryTypes().ATOMName))
        );

        // if the types match, store the atom and update the found status
        if( !current_atom_type.IsDefined())
        {
          const biol::AtomType real_atom_type
          (
            biol::GetAtomTypes().TypeFromPDBAtomName( ( *line_itr)->GetString( GetEntryTypes().ATOMName))
          );

          // if the atom type is in the atom in the file but is not in the aa type used in the bcl (e.g. when using
          // convert to natural aa type), quietly skip the aa
          if( extra_types.Erase( real_atom_type) == 0)
          {
            BCL_MessageCrt
            (
              "found line with atom type, that is not compatible with the residue type: " + ( *line_itr)->GetString()
            );
          }

          continue;
        }

        // try to delete the type
        if( !biol::GetAtomTypes().GetTerminalExtraAtomTypes().Has( current_atom_type) && missing_types.Erase( current_atom_type) == 0)
        {
          BCL_MessageCrt
          (
            "found additional atom types for residue with line: " + ( *line_itr)->GetString()
          );
          continue;
        }

        // if type could have been removed, atom was not seen yet, insert
        atoms.PushBack( AtomFromLine( current_atom_type, **line_itr));
      }

      // iterate over missing types and check that all non-hydrogens atoms could be constructed
      for
      (
        storage::Set< biol::AtomType>::const_iterator itr( missing_types.Begin()), itr_end( missing_types.End());
        itr != itr_end;
        ++itr
      )
      {
        // hydrogen can be missed
        if( ( *itr)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        //
        BCL_MessageVrb
        (
          "at least one non hydrogen atom line is missing for residue starting with line: " +
          LINES.FirstElement()->GetString()
        );

        break;
      }

      // return all atoms
      return atoms;
    }

    //! @brief build amino acid from pdb residue
    //! @param RESIDUE Residue that contains amino acid information
    //! @param SEQ_ID sequence id of amino acid to be read
    //! @return ShPtr to amino acid that was read from given Residue
    util::ShPtr< biol::AABase>
    Factory::AminoAcidFromResidue( const Residue &RESIDUE, const int SEQ_ID) const
    {
      // aa type for current aa
      const biol::AAType aa_type_in_file( biol::GetAATypes().AATypeFromThreeLetterCode( RESIDUE.GetResidueName()));
      const biol::AAType current_aa_type
      (
        aa_type_in_file.IsDefined() && GetFlagConvertToNaturalAAType()->GetFlag()
        ? aa_type_in_file->GetParentType()
        : aa_type_in_file
      );

      // aa data for new amino acid
      util::ShPtr< biol::AAData> sp_aa_data
      (
        new biol::AAData
        (
          current_aa_type,
          SEQ_ID,
          RESIDUE.GetPDBID(),
          RESIDUE.GetICode(),
          RESIDUE.GetChainID()
        )
      );

      // construct an amino acid without any coordinates of specified AAType
      util::ShPtr< biol::AABase> amino_acid( ( *m_AAClass)->Empty( sp_aa_data));

      // set atoms for that amino acid
      //if no atom information is available return amino acid without atoms
      if( RESIDUE.GetLines().IsEmpty())
      {
        // construct an amino acid without any coordinates
        return amino_acid;
      };

      // get the atom types for this AA
      const storage::Set< biol::AtomType> atom_types( current_aa_type->GetAllowedAtomTypes());

      // get all atoms used for this
      util::ShPtrVector< biol::Atom> atoms( AtomsFromLines( RESIDUE.GetLines(), current_aa_type, aa_type_in_file));

      // for glycine, construct a HA2, which corresponds to CB in carbon side chain, for any other AA construct CB if missing
      if( !biol::Atom::FindAtom( atoms, current_aa_type->GetFirstSidechainAtomType())->GetType().IsDefined())
      {
        // first side chain atom coordinate
        const linal::Vector3D first_sc_atom_coord
        (
          FirstSidechainAtomCoordinateFromBackboneAtoms
          (
            atoms,
            biol::GetAtomTypes().CA->GetBondLength( current_aa_type->GetFirstSidechainAtomType())
          )
        );

        // if information to construct pseudo side chain atom is not sufficient
        if( !first_sc_atom_coord.IsDefined())
        {
          // give warning
          BCL_MessageCrt
          (
            "Insufficient data to construct first side chain atom for residue with leading Line:\n" +
            RESIDUE.GetLines().FirstElement()->GetString()
           );
        }
        else
        {
          // insert first side chain atom
          atoms.PushBack( util::ShPtr< biol::Atom>( new biol::Atom( first_sc_atom_coord, current_aa_type->GetFirstSidechainAtomType())));
        }
      }

      // first side chain atom is is in place - set the atoms and return the amino acid
      amino_acid->SetAtoms( atoms);

      // end
      return amino_acid;
    }

    //! @brief builds AASequence from CHAINID and util::ShPtrVector of Residues
    //! @param PDB_ID PdbID of the sequence to be read
    //! @param CHAIN_ID ChainID of the sequence to be read
    //! @param RESIDUES List of Residues that contain amino acid information for the sequence
    //! @return ShPtr to AASequence that was read from Residues
    util::ShPtr< biol::AASequence>
    Factory::AASequenceFromResidues
    (
      const std::string &PDB_ID,
      const char CHAIN_ID,
      const storage::List< Residue> &RESIDUES
    ) const
    {
      //instantiate new aasequences
      const std::string pdb_id( ">" + PDB_ID + ":" + CHAIN_ID + "|PDBID|CHAIN|SEQUENCE");
      util::ShPtr< biol::AASequence> aasequence( new biol::AASequence( m_AAClass, 0, CHAIN_ID, pdb_id));

      //give residues increasing number for later sequence-distance calculations
      int seqid( 1);

      //loop over all pdb residues and construt util::ShPtr< AMINOACID_TYPE> of amino acids
      for
      (
        storage::List< Residue>::const_iterator residue_itr( RESIDUES.Begin()),
          residue_itr_end( RESIDUES.End());
        residue_itr != residue_itr_end;
        ++residue_itr, ++seqid
      )
      {
        aasequence->PushBack( AminoAcidFromResidue( *residue_itr, seqid));
        // check that chain id matches
        BCL_Assert
        (
          aasequence->GetLastAA()->GetChainID() == CHAIN_ID,
          std::string( "amino acid has wrong chain id, should be: ") + CHAIN_ID +
          " but is " + aasequence->GetLastAA()->GetChainID()
        );
      }

      //end
      return aasequence;
    }

    //! @brief cut secondary structure elements from a given sequence
    //! @param SEQEUENCE the amino acid sequence
    //! @param SSE_RESIDUE secondary structure element defined by starting and end residue
    //! @return ShPtr to SSE, undefined if error occurs (no such amino acids, wrong chain, wrong amino acid order etc.)
    util::ShPtr< assemble::SSE> Factory::SecondaryStructureElementFromSequence
    (
      const biol::AASequence &SEQUENCE,
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_RESIDUE
    )
    {
      biol::AASequence::const_iterator aa_itr_first( SEQUENCE.Begin()), aa_itr_last( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());

      // find first amino acid of sse
      while( aa_itr_first != aa_itr_end && SSE_RESIDUE.Second() != **aa_itr_first)
      {
        ++aa_itr_first;
      }

      // find last amino acid of sse
      while( aa_itr_last != aa_itr_end && SSE_RESIDUE.Third() != **aa_itr_last)
      {
        ++aa_itr_last;
      }

      // if either could not be found, error
      if( aa_itr_first == aa_itr_end || aa_itr_last == aa_itr_end)
      {
        BCL_MessageVrb
        (
          "unable to find SSE in given SEQUENCE: " + SSE_RESIDUE.First().GetName() + " " +
          SSE_RESIDUE.Second().GetIdentification() + " " + SSE_RESIDUE.Third().GetIdentification()
        );

        return util::ShPtr< assemble::SSE>();
      }

      // make last the end
      ++aa_itr_last;

      if( std::distance( aa_itr_first, aa_itr_last) < 0)
      {
        BCL_MessageStd( "given sse residues are of no valid sse definition!");

        BCL_MessageStd( "residues: " + util::Format()( SSE_RESIDUE));
        BCL_MessageStd( "first: " + ( *aa_itr_first)->GetIdentification());
        BCL_MessageStd( "last: " + ( *aa_itr_last)->GetIdentification());
        return util::ShPtr< assemble::SSE>();
      }

      // return new sse from subsequence
      util::ShPtr< assemble::SSE> sp_sse
        (
          new assemble::SSE
          (
            biol::AASequence( util::ShPtrVector< biol::AABase>( aa_itr_first, aa_itr_last), SEQUENCE.GetChainID(), SEQUENCE.GetFastaHeader()),
            SSE_RESIDUE.First()
          )
        );

      // iterate over the residues in the sequence
      for
      (
        biol::AASequence::iterator
          aa_itr( sp_sse->Begin()), aa_itr_end( sp_sse->End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // set the data
        ( *aa_itr)->SetSSPrediction
        (
          sspred::GetMethods().e_PDB,
          sspred::PDB( SSE_RESIDUE.First(), biol::GetEnvironmentTypes().e_Solution)
        );
      }

      // end
      return sp_sse;
    }

    //! @brief cut secondary structure elements from a given sequence
    //! @param SEQEUENCE the amino acid sequence
    //! @param SSE_RESIDUES secondary structure elements defined by starting and end residue
    //! @return ShPtrList of SSEs
    util::ShPtrVector< assemble::SSE> Factory::SecondaryStructureElementsFromSequence
    (
      const biol::AASequence &SEQUENCE,
      const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES
    )
    {
      util::ShPtrVector< assemble::SSE> sses;

      // iterate through secondary structure element definitions
      for
      (
        storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >::const_iterator
           sse_itr( SSE_RESIDUES.Begin()), sse_itr_end( SSE_RESIDUES.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const util::ShPtr< assemble::SSE> sp_sse( SecondaryStructureElementFromSequence( SEQUENCE, *sse_itr));
        if( sp_sse.IsDefined())
        {
          sses.PushBack( sp_sse);
        }
      }

      // end
      return sses;
    }

    //! @brief builds a AASequence> from a pdb file read by the given handler
    //! @param HANDLER Handler that contains pdb file
    //! @param PDB_ID Optional PDB ID to be used in the fasta header, otherwise, it is retrieved from the handler
    //! @return ShPtrVector of AASequences that were read from the given Handler
    util::ShPtrVector< biol::AASequence>
    Factory::AASequencesFromPDB( const Handler &HANDLER, const std::string &PDB_ID) const
    {
      // instantiate ShPtrVector of chains
      util::ShPtrVector< biol::AASequence> sequences;

      // get the first model
      const Model &first_model( HANDLER.GetModels().FirstElement());

      //get iterator on pdb structure
      //loop over all chains and instantiate AASequences
      for
      (
        storage::Map< char, storage::List< Residue> >::const_iterator
          chain_itr( first_model.GetStructuredChains().Begin()), chain_itr_end( first_model.GetStructuredChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sequences.PushBack
        (
          AASequenceFromResidues( PDB_ID == "" ? HANDLER.GetHead().GetPDBID() : PDB_ID, chain_itr->first, chain_itr->second)
        );
      }

      //return sequneces
      return sequences;
    }

    //! @brief merge two overlapping sequences into one
    //! @param SEQUENCE_A sequence a
    //! @param SEQUENCE_B sequence b
    //! @return the new sequence - empty if sses are not overlapping/consecutive
    biol::AASequence Factory::MergeOverlappingSequences
    (
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {
      if( !biol::DoOverlap( SEQUENCE_A, SEQUENCE_B))
      {
        return biol::AASequence();
      }

      util::SiPtr< const biol::AASequence> first_seq( SEQUENCE_A);
      util::SiPtr< const biol::AASequence> second_seq( SEQUENCE_B);

      if( first_seq->GetFirstAA()->GetSeqID() <= second_seq->GetLastAA()->GetSeqID())
      {
        std::swap( first_seq, second_seq);
      }

      biol::AASequence new_seq( *first_seq);
      // search for last amino acid of first sequence in the second sequence
      const util::ShPtr< biol::AAData> &last_aa_data( new_seq.GetLastAA()->GetData());
      biol::AASequence::const_iterator aa_itr( second_seq->Begin()), aa_itr_end( second_seq->End());
      for( ; aa_itr != aa_itr_end && ( *aa_itr)->GetData() != last_aa_data; ++aa_itr);

      if( aa_itr == aa_itr_end)
      {
        return new_seq;
      }

      // append the rest of second sequence
      ++aa_itr;
      for( ; aa_itr != aa_itr_end; ++aa_itr)
      {
        new_seq.PushBack( aa_itr->HardCopy());
      }

      // end
      return new_seq;
    }

    //! @brief remove the amino acids that are common to both sses
    //! @param SSE_A overlapping sse a
    //! @param SSE_B overlapping sse b
    //! @return two sses that do not have the amino acids that were common to the two argument sses
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > Factory::RemoveOverlappingAAs
    (
      const util::ShPtr< assemble::SSE> &SSE_A,
      const util::ShPtr< assemble::SSE> &SSE_B
    )
    {
      if( !biol::DoOverlap( *SSE_A, *SSE_B))
      {
        return storage::VectorND< 2, util::ShPtr< assemble::SSE> >( SSE_A, SSE_B);
      }

      util::ShPtr< assemble::SSE> sse_first( SSE_A);
      util::ShPtr< assemble::SSE> sse_second( SSE_B);
      if( sse_second->GetFirstAA()->GetSeqID() <= sse_first->GetFirstAA()->GetSeqID())
      {
        std::swap( sse_first, sse_second);
      }

      biol::AASequence::const_iterator aa_itr1( sse_first->Begin()), aa_itr1_end( sse_first->End());
      biol::AASequence::const_iterator aa_itr2( sse_second->Begin()), aa_itr2_end( sse_second->End());

      // fill first
      util::ShPtrVector< biol::AABase> new_first;
      while( aa_itr1 != aa_itr1_end && ( *aa_itr1)->GetData() != ( *aa_itr2)->GetData())
      {
        new_first.PushBack( aa_itr1->HardCopy());
        ++aa_itr1;
      }

      // skip common amino acids
      for
      (
        ;
        aa_itr1 != aa_itr1_end && aa_itr2 != aa_itr2_end && ( *aa_itr1)->GetData() == ( *aa_itr2)->GetData();
        ++aa_itr1, ++aa_itr2
      );

      // fill second
      util::ShPtrVector< biol::AABase> new_second;
      while( aa_itr2 != aa_itr2_end)
      {
        new_second.PushBack( aa_itr2->HardCopy());
        ++aa_itr2;
      }

      // end
      return storage::VectorND< 2, util::ShPtr< assemble::SSE> >
      (
        util::ShPtr< assemble::SSE>( new_first.IsEmpty() ? NULL : new assemble::SSE( biol::AASequence( new_first, sse_first->GetChainID()), sse_first->GetType())),
        util::ShPtr< assemble::SSE>( new_second.IsEmpty() ? NULL : new assemble::SSE( biol::AASequence( new_second, sse_second->GetChainID()), sse_second->GetType()))
      );
    }

    //! @brief process overlapping secondary structure elements
    //! @param SECONDARY_STRUCTURE_ELEMENTS possibly overlapping secondary structure elements
    //! @param MERGE_OVERLAPPING merge overlapping sses, if they are of the same type
    //! @return storage set of non overlapping SSEs
    storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
    Factory::ProcessOverlappingSSEs
    (
      const util::ShPtrVector< assemble::SSE> &SECONDARY_STRUCTURE_ELEMENTS,
      const bool MERGE_OVERLAPPING
    )
    {
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> non_overlapping_sses;

      // iterate over sses
      for
      (
        util::ShPtrVector< assemble::SSE>::const_iterator
          sse_itr( SECONDARY_STRUCTURE_ELEMENTS.Begin()), sse_itr_end( SECONDARY_STRUCTURE_ELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> sp_sse_insert( *sse_itr);

        // search overlapping sse
        std::set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::iterator overlap_itr;

        // modify the sse to insert until there is no overlapping sse left
        while( sp_sse_insert.IsDefined() && ( overlap_itr = non_overlapping_sses.Find( sp_sse_insert)) != non_overlapping_sses.End())
        {
          // copy overlapping sse and remove from set
          util::ShPtr< assemble::SSE> sp_sse_overlap( *overlap_itr);

          // identical sses
          if
          (
               sp_sse_insert->GetFirstAA()->GetData() == sp_sse_overlap->GetFirstAA()->GetData()
            && sp_sse_insert->GetLastAA()->GetData() == sp_sse_overlap->GetLastAA()->GetData()
          )
          {
            sp_sse_insert = util::ShPtr< assemble::SSE>();
            break;
          }
          non_overlapping_sses.RemoveElement( overlap_itr);

          // merge overlapping
          if( MERGE_OVERLAPPING)
          {
            // check if overlapping sse has different type
            if( sp_sse_overlap->GetType() != sp_sse_insert->GetType())
            {
              // remove the sse
              BCL_MessageCrt
              (
                "found overlapping sses of different type; cannot merge them -> will not consider either\n:" +
                sp_sse_insert->GetIdentification() + "\n" + sp_sse_overlap->GetIdentification()
              );
              sp_sse_insert = util::ShPtr< assemble::SSE>();
              break;
            }
            else
            {
              sp_sse_insert = util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( MergeOverlappingSequences( *sp_sse_insert, *sp_sse_overlap), sp_sse_insert->GetType())
              );
            }
          }
          // split overlapping by removing all overlapping amino acids
          else
          {
            // remove overlapping amino acids
            const storage::VectorND< 2, util::ShPtr< assemble::SSE> > new_sses( RemoveOverlappingAAs( sp_sse_overlap, sp_sse_insert));
            if( new_sses.First().IsDefined() && non_overlapping_sses.Find( new_sses.First()) == non_overlapping_sses.End())
            {
              non_overlapping_sses.Insert( new_sses.First());
              sp_sse_insert = new_sses.Second();
            }
            else if( new_sses.Second().IsDefined() && non_overlapping_sses.Find( new_sses.Second()) == non_overlapping_sses.End())
            {
              non_overlapping_sses.Insert( new_sses.Second());
              sp_sse_insert = new_sses.First();
            }
            else
            {
              BCL_MessageCrt
              (
                "unable to merge SSEs:\n" +
                sp_sse_insert->GetIdentification() + '\n' + sp_sse_overlap->GetIdentification()
              );
              continue;
            }
          }
        }

        // insert sse
        if( sp_sse_insert.IsDefined())
        {
          non_overlapping_sses.Insert( sp_sse_insert);
        }
      }

      // end
      return non_overlapping_sses;
    }

    //! @brief add loops to a set of secondary structure elements
    //! @param SECONDARY_STRUCTURE_ELEMENTS secondary structure elements
    //! @param SEQUENCE the full sequence; coils will be subsequences
    //! @return the set of secondary structure elements with coils
    storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
    Factory::AddLoops
    (
      const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &SECONDARY_STRUCTURE_ELEMENTS,
      const biol::AASequence &SEQUENCE
    )
    {
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> all_sses;

      // check that there is a least one sse
      if( SECONDARY_STRUCTURE_ELEMENTS.IsEmpty())
      {
        // create ShPtr to SSE "new_sse"
        util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( SEQUENCE, biol::GetSSTypes().COIL));

        // iterate over amino acids and set the prediction
        for
        (
          biol::AASequence::iterator aa_itr( new_sse->Begin()), aa_itr_end( new_sse->End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // set the data
          ( *aa_itr)->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );
        }

        all_sses.Insert( new_sse);

        return all_sses;
      }

      // create ShPtrVector iterator "seq_itr" and "seq_itr_end" for iterating over the aa sequence of this chain
      biol::AASequence::const_iterator seq_itr( SEQUENCE.Begin()), seq_itr_end( SEQUENCE.End());

      // loop over all SSEs in this chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( SECONDARY_STRUCTURE_ELEMENTS.Begin()), sse_itr_end( SECONDARY_STRUCTURE_ELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // create size_t "sse_begin_seq_id" initialize with the SeqID of the first aa in the SSE denoted by "sse_itr"
        const int sse_begin_seq_id( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        biol::AASequence new_seq;

        // set the chain ID
        new_seq.SetChainID( SEQUENCE.GetChainID());

        // set the fasta header
        new_seq.SetFastaHeader( SEQUENCE.GetFastaHeader());

        // while seq ID of aa denoted by "seq_itr" has not reached "sse_begin_seq_id" and the end of the aa sequence
        // has not been reached
        while( ( *seq_itr)->GetSeqID() != sse_begin_seq_id && seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to undefined)
          util::ShPtr< biol::AABase> new_aa( ( *seq_itr)->Clone());

          // set the data
          new_aa->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );

          // add the new aa to "new_sse"
          new_seq.PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // true if "new_sse" has some aas in it
        if( new_seq.GetSize() > 0)
        {
          // insert "new_sse" into "new_sses"
          all_sses.Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, biol::GetSSTypes().COIL)));
        }

        // move "seq_itr" to the place of the next non-sse region
        const size_t sse_size( ( *sse_itr)->GetSize());
        for( size_t sse_index( 0); sse_index != sse_size && seq_itr != seq_itr_end; ++sse_index)
        {
          ++seq_itr;
        }

        // add the current sse also to the all sses
        all_sses.Insert( all_sses.End(), *sse_itr);
      }

      // we have to make sure to add the loop after the last SSE
      const size_t last_sse_end_id( ( *SECONDARY_STRUCTURE_ELEMENTS.ReverseBegin())->GetLastAA()->GetSeqID());
      const size_t last_aa_id( SEQUENCE.GetLastAA()->GetSeqID());

      // if there are any residues at the end of the sequence that are not in the last SSE
      if( last_sse_end_id < last_aa_id)
      {
        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        biol::AASequence new_seq;

        // set the chain ID
        new_seq.SetChainID( SEQUENCE.GetChainID());

        // set the fasta header
        new_seq.SetFastaHeader( SEQUENCE.GetFastaHeader());

        // iterate until you hit the end of amino acids in the sequecne
        while( seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to 0.000)
          util::ShPtr< biol::AABase> new_aa( ( *seq_itr)->Clone());

          // set the data
          new_aa->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );

          // add the new aa to "new_sse"
          new_seq.PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // push the last loop
        all_sses.Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, biol::GetSSTypes().COIL)));
      }

      // end
      return all_sses;
    }

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
    Factory::ChainFromPDBSSElementsAndSequence
    (
      const char CHAIN_ID,
      const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES,
      const storage::List< Residue> &SEQUENCE_RESIDUES,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const std::string &PDB_ID
    ) const
    {
      // create the current chain with the sequence
      util::ShPtr< biol::AASequence> sp_sequence( AASequenceFromResidues( PDB_ID, CHAIN_ID, SEQUENCE_RESIDUES));
      util::ShPtr< assemble::Chain> current_chain;

      // assign secondary structure usign dssp
      if( GetFlagDSSP()->GetFlag())
      {
        assemble::ProteinModel model;
        util::ShPtr< assemble::Chain> tmp_chain( new assemble::Chain( sp_sequence));
        tmp_chain->Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( *sp_sequence, biol::GetSSTypes().COIL)));
        model.Insert( tmp_chain);

        const biol::DSSP dssp( GetFlagDSSP()->GetFirstParameter()->GetNumericalValue< double>());
        math::MutateResult< assemble::ProteinModel> result( dssp( model));
        if( result.GetArgument().IsDefined())
        {
          current_chain = result.GetArgument()->GetChain( sp_sequence->GetChainID());
          BCL_MessageVrb( "Used DSSP to assign SSEs; found " + util::Format()( current_chain->GetData().GetSize()) + " SSEs");
        }
      }
      // if the sses are supposed to be created from the backbone conformation ignoring and no sses are given in pdb for that chain
      else if( GetFlagSSEsFromBackBone()->GetFlag())
      {
        // create chain with sses from the conformation
        current_chain = util::ShPtr< assemble::Chain>
        (
          assemble::ConstructChainWithSSEsFromConformation( sp_sequence).Clone()
        );
      }

      else
      {
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> non_overlapping_sses
        (
          ProcessOverlappingSSEs( SecondaryStructureElementsFromSequence( *sp_sequence, SSE_RESIDUES), false)
        );

        non_overlapping_sses = AddLoops( non_overlapping_sses, *sp_sequence);

        current_chain = util::ShPtr< assemble::Chain>
        (
          new assemble::Chain( sp_sequence, util::ShPtrVector< assemble::SSE>( non_overlapping_sses.Begin(), non_overlapping_sses.End()))
        );
      }

      // for each sse type remove sses that are too short are not even considered
      current_chain->FilterByMinSSESizes( SSE_TYPE_MINSIZE);

      // return Chain
      return current_chain;
    }

    //! builds chain without sses
    //! @param CHAIN_ID chain id of given fasta
    //! @param ISTREAM fasta file stream
    //! @return chain without SSEs read from the given fasta stream
    util::ShPtr< assemble::Chain>
    Factory::ChainFromFastaStream( const char CHAIN_ID, std::istream &ISTREAM) const
    {
      // create sequence from fasta
      util::ShPtr< biol::AASequence> sp_aa_seq
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( ISTREAM, m_AAClass, CHAIN_ID))
      );

      //! create chain
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( sp_aa_seq));

      // return
      return sp_chain;
    }

    //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
    //! The types and sizes of SSEs to be considered are determined by the additional parameters
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @return Protein Model read from given Handler
    assemble::ProteinModel
    Factory::ProteinModelFromPDB
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE
    ) const
    {
      return ProcessModel( HANDLER, SSE_TYPE_MINSIZE, HANDLER.GetModels().FirstElement());
    }

    //! @brief builds an ensemble from the models stored in the handler
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @return Protein ensemble read from given Handler
    assemble::ProteinEnsemble Factory::ProteinEnsembleFromPDB
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE
    ) const
    {
      // initialize ensemble
      assemble::ProteinEnsemble ensemble;

      // iterate over the models
      for
      (
        storage::Vector< Model>::const_iterator model_itr( HANDLER.GetModels().Begin()),
          model_itr_end( HANDLER.GetModels().End());
        model_itr != model_itr_end; ++model_itr
      )
      {
        // build the model and add it to the ensemble
        ensemble.InsertElement( util::CloneToShPtr( ProcessModel( HANDLER, SSE_TYPE_MINSIZE, *model_itr)));
      }

      // end
      return ensemble;
    }

    //! @brief creates and returns a protein model based on a PDB filename
    //! @param PDB_FILENAME is the pdb filename from which the protein model will be created
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @param IGNORE_CLASH whether to ignore clashes or not
    //! @return ProteinModel which was created from "PDB_FILENAME"
    assemble::ProteinModel Factory::ProteinModelFromPDBFilename
    (
      const std::string &PDB_FILENAME,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const bool IGNORE_CLASH
    ) const
    {
      // initialize write and read stream objects
      io::IFStream read;

      // open pdb file
      io::File::MustOpenIFStream( read, PDB_FILENAME);

      // true is used to advice handler to ignore clashes in the structure and insert residues as they are suggested by
      // the numbering without regard for the bond information
      Handler pdb( read, IGNORE_CLASH);
      io::File::CloseClearFStream( read);

      // create a protein model from the handler
      assemble::ProteinModel model( ProteinModelFromPDB( pdb, SSE_TYPE_MINSIZE));

      // add filename to the protein model data
      util::ShPtr< assemble::ProteinModelData> sp_data( model.GetProteinModelData().HardCopy());
      sp_data->Insert
      (
        assemble::ProteinModelData::e_PDBFile,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( PDB_FILENAME))
      );
      model.SetProteinModelData( sp_data);
      return model;
    }

    //! @brief returns the map of SSTypes and their corresponding minimal size
    //! @param FLAG the flag to get ss type min sizes from
    //! @return map of ss type and associated value
    storage::Map< biol::SSType, size_t>
    Factory::GetCommandlineSSETypeMinSizes( const command::FlagInterface &FLAG)
    {
      return
        GetSSETypeMinSizes
        (
          FLAG.GetParameterList()( biol::GetSSTypes().HELIX)->GetNumericalValue< size_t>(),
          FLAG.GetParameterList()( biol::GetSSTypes().STRAND)->GetNumericalValue< size_t>(),
          FLAG.GetParameterList()( biol::GetSSTypes().COIL)->GetNumericalValue< size_t>()
        );
    }

    //! @brief returns the map of SSTypes and their corresponding minimal size
    //! @param MIN_HELIX_SIZE the minimum helix size
    //! @param MIN_STRAND_SIZE the minimum strand size
    //! @param MIN_COIL_SIZE the minimum coil size
    //! @return map of ss type and associated value
    storage::Map< biol::SSType, size_t>
    Factory::GetSSETypeMinSizes
    (
      const size_t &MIN_HELIX_SIZE,
      const size_t &MIN_STRAND_SIZE,
      const size_t &MIN_COIL_SIZE
    )
    {
      storage::Map< biol::SSType, size_t> sse_type_min_sizes;
      sse_type_min_sizes[ biol::GetSSTypes().HELIX] = MIN_HELIX_SIZE;
      sse_type_min_sizes[ biol::GetSSTypes().STRAND] = MIN_STRAND_SIZE;
      sse_type_min_sizes[ biol::GetSSTypes().COIL] = MIN_COIL_SIZE;

      const storage::Set< biol::SSType> helix_classes( Handler::GetHelixClassesFromCommandLine());

      // iterate over the helix types
      for
      (
        storage::Set< biol::SSType>::const_iterator
          helix_class_itr( helix_classes.Begin()), helix_class_itr_end( helix_classes.End());
        helix_class_itr != helix_class_itr_end;
        ++helix_class_itr
      )
      {
        sse_type_min_sizes[ *helix_class_itr] = MIN_HELIX_SIZE;
      }

      // end
      return sse_type_min_sizes;
    }

  ////////////////////////
  // operations - write //
  ////////////////////////

    //! @brief write header information to Line
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return ShPtr to Line that contains Header information
    util::ShPtr< Line> Factory::WriteHeaderToLine( const std::string &CLASSIFICATION, const util::Time &TIME)
    {
      // instantiate new header line
      util::ShPtr< Line> header_line( new Line( GetLineTypes().HEADER));

      // add the entries
      header_line->Put( GetEntryTypes().HEADERClassification, CLASSIFICATION);
      header_line->Put( GetEntryTypes().HEADERDate, ConvertTimeToPDBDate( TIME));
      header_line->Put( GetEntryTypes().HEADERIDCode, Head::GetBCLPdbID());

      // return the Line
      return header_line;
    }

    //! @brief write Atom to Line of a certain Residue and the CHAINID and the atom SERIAL
    //! @param ATOM Atom of interest
    //! @param AMINO_ACID Amino acid which atom belongs to
    //! @param CHAIN_ID chain id of the sequence atom belongs to
    //! @param SERIAL serial id of the atom
    //! @return ShPtr to Line that contains Atom information
    util::ShPtr< Line> Factory::WriteAtomToLine
    (
      const biol::Atom &ATOM,
      const biol::AABase &AMINO_ACID,
      const char CHAIN_ID,
      const size_t SERIAL
    )
    {
      util::ShPtr< Line> atom_line;

      if( AMINO_ACID.GetType()->IsNaturalAminoAcid())
      {
        //instantiate new atom line
        atom_line = util::ShPtr< Line>( new Line( GetLineTypes().ATOM));

        // write defined coordinates
        if( ATOM.GetCoordinates().IsDefined())
        {
          atom_line->PutCoordinates( ATOM.GetCoordinates());

          // write default occupancy of 1.00 with the precision 2
          atom_line->Put( GetEntryTypes().ATOMOccupancy, double( 1));
        }
        // true if need to write zero coordinates for undefined atoms
        else if( GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetFlag())
        {
          atom_line->PutCoordinates( linal::Vector3D( 0.0, 0.0, 0.0));

          // write default occupancy of -1.00 with the precision 2 for missing residues, compatible to rosetta
          atom_line->Put( GetEntryTypes().ATOMOccupancy, double( -1));
        }
        else
        {
          // return empty line
          return util::ShPtr< Line>();
        }

        // write atom ID
        if( GetFlagWritePDBAtomID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMSerial, ATOM.GetPdbID());
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMSerial, SERIAL);
        }

        // use the pdb atom name instead of the IUPAC name
        if( GetFlagPDBAtomName()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMName, AMINO_ACID.GetType()->GetPDBAtomName( ATOM.GetType()));
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMName, ATOM.GetType()->AtomTypeData::GetName());
        }

        // residue name and chain id
        atom_line->Put( GetEntryTypes().ATOMResidueName, AMINO_ACID.GetType()->GetThreeLetterCode());
        atom_line->Put( GetEntryTypes().ATOMChainID, CHAIN_ID);

        // write res ID
        if( GetFlagWritePDBResID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMResidueSequenceID, AMINO_ACID.GetPdbID());
          atom_line->Put( GetEntryTypes().ATOMInsertionCode, AMINO_ACID.GetPdbICode());
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMResidueSequenceID, AMINO_ACID.GetSeqID());
        }

        // write temp factor with the precision 2
        atom_line->Put
        (
          GetEntryTypes().ATOMTempFactor,
          util::IsDefined( ATOM.GetBFactor()) ? ATOM.GetBFactor() : double( 0.0)
        );

        // element type
        atom_line->Put( GetEntryTypes().ATOMElement, std::string( ATOM.GetType()->GetElementType()->GetChemicalSymbol()));
      }
      // non natural amino acid as HETATM line
      else
      {
        //instantiate new atom line
        atom_line = util::ShPtr< Line>( new Line( GetLineTypes().HETATM));

        // write defined coordinates
        if( ATOM.GetCoordinates().IsDefined())
        {
          atom_line->PutCoordinates( ATOM.GetCoordinates());

          // write default occupancy of 1.00 with the precision 2
          atom_line->Put( GetEntryTypes().HETATMOccupancy, double( 1));
        }
        // true if need to write zero coordinates for undefined atoms
        else if( GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetFlag())
        {
          atom_line->PutCoordinates( linal::Vector3D( 0.0, 0.0, 0.0));

          // write default occupancy of -1.00 with the precision 2 for missing residues, compatible to rosetta
          atom_line->Put( GetEntryTypes().HETATMOccupancy, double( -1));
        }
        else
        {
          // return empty line
          return util::ShPtr< Line>();
        }

        // write atom ID
        if( GetFlagWritePDBAtomID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMSerial, ATOM.GetPdbID());
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMSerial, SERIAL);
        }

        // use the pdb atom name instead of the IUPAC name
        if( GetFlagPDBAtomName()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMName, AMINO_ACID.GetType()->GetPDBAtomName( ATOM.GetType()));
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMName, ATOM.GetType()->AtomTypeData::GetName());
        }

        // residue name and chain id
        atom_line->Put( GetEntryTypes().HETATMResidueName, AMINO_ACID.GetType()->GetThreeLetterCode());
        atom_line->Put( GetEntryTypes().HETATMChainID, CHAIN_ID);

        // write res ID
        if( GetFlagWritePDBResID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMResidueSequenceID, AMINO_ACID.GetPdbID());
          atom_line->Put( GetEntryTypes().HETATMInsertionCode, AMINO_ACID.GetPdbICode());
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMResidueSequenceID, AMINO_ACID.GetSeqID());
        }

        // write temp factor with the precision 2
        atom_line->Put
        (
          GetEntryTypes().HETATMTempFactor,
          util::IsDefined( ATOM.GetBFactor()) ? ATOM.GetBFactor() : double( 0.0)
        );

        // element type
        atom_line->Put( GetEntryTypes().HETATMElement, std::string( ATOM.GetType()->GetElementType()->GetChemicalSymbol()));
      }

      //return
      return atom_line;
    }

    //! @brief write Residue to Lines of a certain CHAINID and the atom SERIAL
    //! @param AMINO_ACID Amino acid of interest
    //! @param CHAIN_ID chain id of the sequence amino acid belongs to
    //! @param SERIAL serial id of the amino acid
    //! @return ShPtr to Line that contains amino acid information
    util::ShPtrList< Line> Factory::WriteResiduesToLines
    (
      const biol::AABase &AMINO_ACID,
      const char CHAIN_ID,
      size_t &SERIAL
    )
    {
      //instantiate residue lines of amino acids
      util::ShPtrList< Line> residue_lines;

      //loop over all atoms in residue
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( AMINO_ACID.GetAtoms().Begin()),
          atom_itr_end( AMINO_ACID.GetAtoms().End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // in non debug mode
        if( !util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
        {
          // by default, write all atom atoms that are not hydrogens or have defined coordinates
          // write hydrogen atoms only if the corresponding flag is set
          if
          (
               !GetFlagWriteHydrogens()->GetFlag()
            && ( *atom_itr)->GetType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen
          )
          {
            continue;
          }
        }

        // ShPtr to new line
        const util::ShPtr< Line> sp_line( WriteAtomToLine( **atom_itr, AMINO_ACID, CHAIN_ID, SERIAL));
        if( sp_line.IsDefined())
        {
          // write the actual line
          residue_lines.PushBack( sp_line);
          ++SERIAL;
        }
      }

      //return
      return residue_lines;
    }

    //! @brief write all atoms of a sequence to a list of lines to be used eventually for writing to a pdb
    //! @param AA_SEQUENCE AASequence of interest
    //! @param SERIAL serial id of the sequence
    //! @return ShPtrList of Lines that contain AASequence information
    util::ShPtrList< Line> Factory::WriteAASequenceToLines( const biol::AASequence &AA_SEQUENCE, size_t &SERIAL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> sequence_lines;

      //loop over all residues
      for
      (
        biol::AASequence::const_iterator residue_itr( AA_SEQUENCE.GetData().Begin()),
          residue_itr_end( AA_SEQUENCE.GetData().End());
        residue_itr != residue_itr_end;
        ++residue_itr
      )
      {
        sequence_lines.Append( WriteResiduesToLines( **residue_itr, AA_SEQUENCE.GetChainID(), SERIAL));
      }

      //return
      return sequence_lines;
    }

    //! @brief write a Chain to a SharedpointerVector of Lines
    //! @param CHAIN Chain of interest
    //! @param SERIAL serial id
    //! @return ShPtrList of Lines that contain Chain information
    util::ShPtrList< Line> Factory::WriteChainToLines( const assemble::Chain &CHAIN, size_t &SERIAL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> chain_lines;

      //loop over all sses in chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        chain_lines.Append( WriteAASequenceToLines( ( **sse_itr), SERIAL));
      }

      // TER line
      if( chain_lines.IsEmpty())
      {
        return chain_lines;
      }

      const Line &last_line( *chain_lines.LastElement());

      util::ShPtr< Line> ter_line( new Line( GetLineTypes().TER));
      ter_line->Put( GetEntryTypes().TERSerial, SERIAL);
      ter_line->Put( GetEntryTypes().TERChainID, CHAIN.GetChainID());
      if( last_line.GetType() == GetLineTypes().ATOM)
      {
        ter_line->Put( GetEntryTypes().TERInsertionCode, last_line.GetChar( GetEntryTypes().ATOMInsertionCode));
        ter_line->Put( GetEntryTypes().TERResidueName, last_line.GetString( GetEntryTypes().ATOMResidueName));
        ter_line->Put( GetEntryTypes().TERResidueSequenceID, last_line.GetString( GetEntryTypes().ATOMResidueSequenceID));
      }
      else
      {
        ter_line->Put( GetEntryTypes().TERInsertionCode, last_line.GetChar( GetEntryTypes().HETATMInsertionCode));
        ter_line->Put( GetEntryTypes().TERResidueName, last_line.GetString( GetEntryTypes().HETATMResidueName));
        ter_line->Put( GetEntryTypes().TERResidueSequenceID, last_line.GetString( GetEntryTypes().HETATMResidueSequenceID));
      }
      chain_lines.PushBack( ter_line);
      ++SERIAL;

      //return
      return chain_lines;
    }

    //! @brief write a single Protein Model to pdb lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return lines to which ProteinModel was written to
    util::ShPtrList< Line> Factory::WriteCompleteModelToPDBLines
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const bool WRITE_BODY_INFORMATION,
      const std::string &CLASSIFICATION,
      const util::Time &TIME
    ) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // the HEADER section
      if( CLASSIFICATION.empty())
      {
        lines.PushBack( util::ShPtr< Line>( new Line( std::string( "HEADER  BCL MODEL"))));
      }
      else
      {
        lines.PushBack( WriteHeaderToLine( CLASSIFICATION, TIME));
      }

      // the SEQRES section
      lines.Append( WriteSeqResToLines( PROTEIN_MODEL));

      size_t serial( 1);
      // the HELIX section
      lines.Append( WriteHelixDefinitionsToLines( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial));

      serial = 1;
      // the STRAND section
      lines.Append( WriteStrandDefinitionsToLines( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().STRAND), serial));

      serial = 1;
      // write the residues and atoms to pdblines
      lines.Append( WriteProteinModelToLines( PROTEIN_MODEL, serial));

      //write bodyinforamtions for each sse
      if( WRITE_BODY_INFORMATION)
      {
        lines.Append( WriteBodyInformationToLines( PROTEIN_MODEL));
      }

      // iterate over printers
      for
      (
        util::ShPtrList< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >::const_iterator
          itr( m_Printers.Begin()), itr_end( m_Printers.End());
        itr != itr_end; ++itr
      )
      {
        lines.Append( ( *itr)->operator ()( PROTEIN_MODEL));
      }

      // end
      return lines;
    }

    //! @brief write a ProteinModel to a SharedpointerVector of Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param SERIAL serial id
    //! @return ShPtrList of Lines that contain ProteinModel information
    util::ShPtrList< Line> Factory::WriteProteinModelToLines
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      size_t &SERIAL
    )
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> proteinmodel_lines;

      //loop over all chains of proteinmodel
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        proteinmodel_lines.Append( WriteChainToLines( **chain_itr, SERIAL));
      }

      //return
      return proteinmodel_lines;
    }

    //! @brief write SSE definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param SERIAL SSE serial id
    //! @return ShPtr to Line containing SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteSSEDefinitionToLine( const assemble::SSE &THIS_SSE, const size_t &SERIAL)
    {
      // if Helix type over SSType
      if( biol::GetSSTypes().GetHelixTypes().Contains( THIS_SSE.GetType()))
      {
        return WriteHelixDefinitionToLine( THIS_SSE, SERIAL);
      }
      // else if STRAND
      else if( THIS_SSE.GetType() == biol::GetSSTypes().STRAND)
      {
        return WriteStrandDefinitionToLine( THIS_SSE, SERIAL);
      }
      else
      {
        return util::ShPtr< Line>();
      }
    }

    //! @brief write helix definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param HELIX_SERIAL helix serial id
    //! @return ShPtr to Line containing helix SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteHelixDefinitionToLine
    (
      const assemble::SSE &THIS_SSE,
      const size_t HELIX_SERIAL
    )
    {
      // empty sse
      if( THIS_SSE.GetSize() == 0)
      {
        return util::ShPtr< Line>();
      }

      util::ShPtr< Line> helix_line( new Line( GetLineTypes().HELIX));

      //write all information for helix definiton line
      //helix serial
      helix_line->Put( GetEntryTypes().HELIXSerial, HELIX_SERIAL);
      //helix id
      helix_line->Put( GetEntryTypes().HELIXID, HELIX_SERIAL);

      // first residue name
      helix_line->Put
      (
        GetEntryTypes().HELIXResidueName_Initial, THIS_SSE.GetFirstAA()->GetType()->GetThreeLetterCode()
      );
      // last residue name
      helix_line->Put
      (
        GetEntryTypes().HELIXResidueName_Terminal, THIS_SSE.GetLastAA()->GetType()->GetThreeLetterCode()
      );

      // for res id use either the original or the bcl numerated seqid
      if( GetFlagWritePDBResID()->GetFlag())
      {
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Initial, THIS_SSE.GetFirstAA()->GetPdbID());
        helix_line->Put( GetEntryTypes().HELIXInsertionCode_Initial, THIS_SSE.GetFirstAA()->GetPdbICode());
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Terminal, THIS_SSE.GetLastAA()->GetPdbID());
        helix_line->Put( GetEntryTypes().HELIXInsertionCode_Terminal, THIS_SSE.GetLastAA()->GetPdbICode());
      }
      else
      {
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Initial, THIS_SSE.GetFirstAA()->GetSeqID());
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Terminal, THIS_SSE.GetLastAA()->GetSeqID());
      }

      // chain ID
      helix_line->Put( GetEntryTypes().HELIXChainID_Initial, THIS_SSE.GetChainID());
      helix_line->Put( GetEntryTypes().HELIXChainID_Terminal, THIS_SSE.GetChainID());

      // helix class - currently right handed alpha helix
      helix_line->Put< size_t>( GetEntryTypes().HELIXClass, biol::GetSSTypes().PDBHelixClassFromSSType( THIS_SSE.GetType()));

      // length of helix
      helix_line->Put( GetEntryTypes().HELIXLength, THIS_SSE.GetSize());

      // return
      return helix_line;
    }

    //! @brief write Helix definitions to Lines
    //! @param SSE_VECTOR SiPtrVector of SSEs of interest
    //! @param HELIX_SERIAL helix serial id
    //! @return ShPtrList of Lines containing helix SSE defintions for given vector of SSEs
    util::ShPtrList< Line> Factory::WriteHelixDefinitionsToLines
    (
      const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
      size_t &HELIX_SERIAL
    )
    {
      util::ShPtrList< Line> helixlines;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          helix_itr( SSE_VECTOR.Begin()), helix_itr_end( SSE_VECTOR.End());
        helix_itr != helix_itr_end; ++helix_itr)
      {
        helixlines.PushBack( WriteHelixDefinitionToLine( **helix_itr, HELIX_SERIAL++));
      }

      return helixlines;
    }

    //! @brief write Strand definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param STRAND_SERIAL strand serial id
    //! @return ShPtr to Line containing strand SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteStrandDefinitionToLine
    (
      const assemble::SSE &THIS_SSE,
      const size_t STRAND_SERIAL
    )
    {
      // empty sse
      if( THIS_SSE.GetSize() == 0)
      {
        return util::ShPtr< Line>();
      }

      util::ShPtr< Line> strand_line( new Line( GetLineTypes().SHEET));

      //write all information for strand definiton line
      //strand ID
      strand_line->Put( GetEntryTypes().SHEETStrandID, STRAND_SERIAL);
      //first residue name
      strand_line->Put
      (
        GetEntryTypes().SHEETResidueName_Initial, THIS_SSE.GetFirstAA()->GetType()->GetThreeLetterCode()
      );
      // last residue name
      strand_line->Put
      (
        GetEntryTypes().SHEETResidueName_Terminal, THIS_SSE.GetLastAA()->GetType()->GetThreeLetterCode()
      );

      // for res id use either the original or the bcl numerated seqid
      if( GetFlagWritePDBResID()->GetFlag())
      {
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Initial, THIS_SSE.GetFirstAA()->GetPdbID());
        strand_line->Put( GetEntryTypes().SHEETInsertionCode_Initial, THIS_SSE.GetFirstAA()->GetPdbICode());
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Terminal, THIS_SSE.GetLastAA()->GetPdbID());
        strand_line->Put( GetEntryTypes().SHEETInsertionCode_Terminal, THIS_SSE.GetLastAA()->GetPdbICode());
      }
      else
      {
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Initial, THIS_SSE.GetFirstAA()->GetSeqID());
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Terminal, THIS_SSE.GetLastAA()->GetSeqID());
      }

      //chain ID
      strand_line->Put( GetEntryTypes().SHEETChainID_Initial, THIS_SSE.GetChainID());
      strand_line->Put( GetEntryTypes().SHEETChainID_Terminal, THIS_SSE.GetChainID());

      // end
      return strand_line;
    }

    //! @brief write strand definitions to Lines
    //! @param SSE_VECTOR SiPtrVector of SSEs of interest
    //! @param STRAND_SERIAL strand serial id
    //! @return ShPtrList of Lines containing strand SSE definitions for given vector of SSEs
    util::ShPtrList< Line> Factory::WriteStrandDefinitionsToLines
    (
      const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
      size_t &STRAND_SERIAL
    )
    {
      util::ShPtrList< Line> strandlines;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( SSE_VECTOR.Begin()), strand_itr_end( SSE_VECTOR.End());
        strand_itr != strand_itr_end; ++strand_itr
      )
      {
        strandlines.PushBack( WriteStrandDefinitionToLine( **strand_itr, STRAND_SERIAL++));
      }

      return strandlines;
    }

    //! @brief write assemble::SSEGeometry Information to Lines of one SSE
    //! @param THIS_SSE SSE of interest
    //! @return ShPtrList of Lines containing body information for the given SSE
    util::ShPtrList< Line> Factory::WriteBodyInformationToLines( const assemble::SSE &THIS_SSE)
    {
      //instantiate lines for current body
      util::ShPtrList< Line> body_lines;

      //instantiate std::string storing the temporary line
      std::string tmp;

      //write extension
      tmp =
        std::string( "REMARK  ") + std::string( "Extension in Z direction: ") +
        util::Format()( THIS_SSE.GetExtent( coord::GetAxes().e_Z));
      body_lines.PushBack( util::ShPtr< Line>( new Line( tmp)));

      // insert atoms at the beginning and end of bodies
      util::ShPtr< Line> begin_atom( new Line( GetLineTypes().ATOM));
      util::ShPtr< Line> end_atom( new Line( GetLineTypes().ATOM));
      begin_atom->PutCoordinates( THIS_SSE.BeginOfZ());
      end_atom->PutCoordinates( THIS_SSE.EndOfZ());
      body_lines.PushBack( begin_atom);
      body_lines.PushBack( end_atom);

      //return
      return body_lines;
    }

    //! @brief write assemble::SSEGeometry Information to Lines for given ProteinModel
    //! writes TransformationMatrix and assemble::SSEGeometry extension of all SSES to Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return ShPtrList of Lines containing body information for the given ProteinModel
    util::ShPtrList< Line> Factory::WriteBodyInformationToLines( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> body_lines;

      //instantiate Lines for sheet and helix definitons
      util::ShPtrList< Line> sse_helix_definition_lines;
      util::ShPtrList< Line> sse_strand_definition_lines;

      //loop over all chains of proteinmodel
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        //loop over all sses
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          //according to type write sse line
          // the HELIX lines
          if( biol::GetSSTypes().GetHelixTypes().Contains( ( *sse_itr)->GetType()))
          {
            sse_helix_definition_lines.PushBack
            (
              WriteHelixDefinitionToLine( **sse_itr, sse_helix_definition_lines.GetSize())
            );
            sse_helix_definition_lines.Append( WriteBodyInformationToLines( **sse_itr));
          }
          // the STRAND lines
          else if( ( *sse_itr)->GetType() == biol::GetSSTypes().STRAND)
          {
            sse_strand_definition_lines.PushBack
            (
              WriteStrandDefinitionToLine( **sse_itr, sse_strand_definition_lines.GetSize())
            );
            sse_strand_definition_lines.Append( WriteBodyInformationToLines( **sse_itr));
          }
        }
      }

      util::ShPtrList< Line> sse_definition_lines( sse_helix_definition_lines);
      sse_definition_lines.Append( sse_strand_definition_lines);

      //end
      return sse_definition_lines;
    }

    //! @brief write SSE definitions to Lines for a util::ShPtrVector of models
    //! @param PROTEIN_MODELS ShPtrVector of ProteinModels of interest
    //! @return ShPtrList of Lines containing SSE definition information for given ShPtrVector of ProteinModels
    util::ShPtrList< Line> Factory::WriteSSEDefinitionsToLines
    (
      const util::ShPtrVector< assemble::ProteinModel> &PROTEIN_MODELS
    )
    {
      //instantiate Lines for sheet and helix definitons
      util::ShPtrList< Line> sse_helix_definition_lines;
      util::ShPtrList< Line> sse_strand_definition_lines;
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          model_itr( PROTEIN_MODELS.Begin()), model_itr_end( PROTEIN_MODELS.End());
        model_itr != model_itr_end; ++model_itr)
      {
        size_t serial( sse_helix_definition_lines.GetSize() + 1);
        // the HELIX lines
        sse_helix_definition_lines.Append
        (
          WriteHelixDefinitionsToLines( ( *model_itr)->GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial)
        );

        serial = 1;
        // the STRAND lines
        sse_helix_definition_lines.Append
        (
          WriteStrandDefinitionsToLines( ( *model_itr)->GetSSEs( biol::GetSSTypes().STRAND), serial)
        );
      }

      //put Helix and strand definition lines in one util::ShPtrVector of Lines
      util::ShPtrList< Line> sse_definition_lines( sse_helix_definition_lines);
      sse_definition_lines.Append( sse_strand_definition_lines);

      //return helix and strand definiton lines
      return sse_definition_lines;
    }

    //! @brief write SEQRES lines for a sequence
    //! @param SEQUENCE AASequence of interest
    //! @return ShPtrList of Lines containing SEQRES information for given AASequence
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const biol::AASequence &SEQUENCE)
    {
      const char chainid( SEQUENCE.GetChainID());

      size_t number_aas( 0);
      std::stringstream string_sequence;

      //iterate over all aas in chain an write the threelettercodes in a string and count the aas
      for
      (
        biol::AASequence::const_iterator
          aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        string_sequence << ( *aa_itr)->GetType()->GetThreeLetterCode() << ' ';
        number_aas++;
      }

      //string to read in the threelettercodes for the string_sequence
      std::string threelettercode;
      //the seqreslines that will contain the complete seqrespart for the pdb in Lines
      util::ShPtrList< Line> seqreslines;
      //iterate over all seqserials each containing 13 AAs in the end
      for( size_t seqserial( 1); seqserial <= ( number_aas / 13) + 1; ++seqserial)
      {
        //new SEQRES line
        util::ShPtr< Line> seqresline( new Line( GetLineTypes().SEQRES));
        //write seqserial to Line
        seqresline->Put( GetEntryTypes().SEQRESSerial, seqserial);
        //write chainid to seqresline
        seqresline->Put( GetEntryTypes().SEQRESChainID, chainid);
        //write number of residues in complete chain
        seqresline->Put( GetEntryTypes().SEQRESNrOfResiduesInChain, number_aas);
        //write 13 aas in each line
        for( size_t seqres_pos( 0); seqres_pos < 13 && seqres_pos < ( number_aas - ( seqserial - 1) * 13); ++seqres_pos)
        {
          //read the threelettercode
          string_sequence >> threelettercode;
          //write the current three letter code to the Line
          seqresline->Put( EntryType( GetEntryTypes().SEQRESName_1 + seqres_pos), threelettercode);
        }
        //push the line to the lines
        seqreslines.PushBack( seqresline);
      }

      return seqreslines;
    }

    //! @brief write SEQRES information of a Chain into Lines
    //! @param CHAIN Chain of interest
    //! @return ShPtrList of Lines containing SEQRES information for given Chain
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const assemble::Chain &CHAIN)
    {
      //if there is a sequence given with the model write the seqresinformation from the complete sequence
      if( CHAIN.GetSequence().IsDefined())
      {
        return WriteSeqResToLines( *CHAIN.GetSequence());
      }

      //if there are no sses in the Proteinmodel return an empty Seqreslines vector
      if( CHAIN.GetData().IsEmpty())
      {
        return util::ShPtrList< Line>();
      }

      //store the chainid
      const char chainid( CHAIN.GetChainID());

      size_t number_aas( 0);
      std::stringstream string_sequence;

      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        for
        (
          biol::AASequence::const_iterator
            aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End());
          aa_itr != aa_itr_end; ++aa_itr)
        {
          number_aas++;
          while( int( number_aas) < ( *aa_itr)->GetSeqID())
          {
            string_sequence << biol::GetAATypes().e_Undefined->GetThreeLetterCode() << ' ';
            number_aas++;
          }
          string_sequence << ( *aa_itr)->GetType()->GetThreeLetterCode() << ' ';
        }
      }

      std::string threelettercode;
      util::ShPtrList< Line> seqreslines;
      for( size_t seqserial( 1); seqserial <= ( number_aas / 13) + 1; ++seqserial)
      {
        util::ShPtr< Line> seqresline( new Line( GetLineTypes().SEQRES));
        seqresline->Put( GetEntryTypes().SEQRESSerial, seqserial);
        seqresline->Put( GetEntryTypes().SEQRESChainID, chainid);
        seqresline->Put( GetEntryTypes().SEQRESNrOfResiduesInChain, number_aas);
        for( size_t seqres_pos( 0); seqres_pos < 13 && seqres_pos < ( number_aas - ( seqserial - 1) * 13); ++seqres_pos)
        {
          string_sequence >> threelettercode;
          seqresline->Put( EntryType( GetEntryTypes().SEQRESName_1 + seqres_pos), threelettercode);
        }
        seqreslines.PushBack( seqresline);
      }

      //return
      return seqreslines;
    }

    //! @brief write SEQRES information of a ProteinModel into Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return ShPtrList of Lines containing SEQRES information for given ProteinModel
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      //accumulate seqreslines for each chain
      util::ShPtrList< Line> seqreslines;

      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr)
      {
        seqreslines.Append( WriteSeqResToLines( **chain_itr));
      }

      //return
      return seqreslines;
    }

    //! @brief write a single Protein Model to a pdb file stream
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param OSTREAM output stream to which ProteinModel will be written to
    //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return std::ostream to which ProteinModel was written to
    std::ostream &Factory::WriteModelToPDB
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM,
      const bool WRITE_BODY_INFORMATION,
      const std::string &CLASSIFICATION,
      const util::Time &TIME
    ) const
    {
      // create handler
      Handler newpdb;

      // add lines
      newpdb.AppendLines
      (
        WriteCompleteModelToPDBLines( PROTEIN_MODEL, WRITE_BODY_INFORMATION, CLASSIFICATION, TIME)
      );

      // write lines
      return newpdb.WriteLines( OSTREAM);
    }

    //! @brief write a Chain to a pdb file stream
    //! @param CHAIN Chain of interest
    //! @param OSTREAM output stream to which Chain will be written to
    //! @return std::ostream to which Chain was written to
    std::ostream &Factory::WriteChainToPDB( const assemble::Chain &CHAIN, std::ostream &OSTREAM)
    {
      Handler newpdb;
      newpdb.PushBack
      (
        util::ShPtr< Line>
        (
          new Line( std::string( "HEADER  BCL MODEL"))// + util::Time::GetCurrent().GetTimeAsDate()))
        )
      );

      // the SEQRES section
      newpdb.AppendLines( WriteSeqResToLines( CHAIN));

      size_t serial( 1);
      // the HELIX section
      newpdb.AppendLines( WriteHelixDefinitionsToLines( CHAIN.GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial));

      serial = 1;
      // the STRAND section
      newpdb.AppendLines( WriteStrandDefinitionsToLines( CHAIN.GetSSEs( biol::GetSSTypes().STRAND), serial));

      serial = 1;
      // write the residues and atoms to pdblines
      newpdb.AppendLines( WriteChainToLines( CHAIN, serial));

      return newpdb.WriteLines( OSTREAM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Factory from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Factory::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AAClass, ISTREAM);
      io::Serialize::Read( m_Printers, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write Factory to given OSTREAM
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return std::ostream to which was written
    std::ostream &Factory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write members
      io::Serialize::Write( m_AAClass, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Printers, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct the first side chain atom from given atoms
    //! those atoms need to contain the CA, N and O
    //! @param ATOMS CA, N and C atom
    //! @param BOND_LENGTH the distance of the constructed side chain atom coordinate to the CA atom
    //! @return position of first side chain atom for L-amino acid
    linal::Vector3D Factory::FirstSidechainAtomCoordinateFromBackboneAtoms
    (
      const util::SiPtrVector< const biol::Atom> &ATOMS,
      const double BOND_LENGTH
    )
    {
      // initialize pointers to atoms required for building a pseudo HA2
      // CA
      const util::SiPtr< const biol::Atom> ca_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().CA));
      if( !ca_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // N
      const util::SiPtr< const biol::Atom> n_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().N));
      if( !n_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // C
      const util::SiPtr< const biol::Atom> c_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().C));
      if( !c_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // angle for tetrahedral geometry
      static const double s_tetrahedral_angle( 109.5 / 180.0 * math::g_Pi);

      // construct atom coordinate
      return linal::CoordinatesAngle
             (
               ca_ptr->GetCoordinates(),
               n_ptr->GetCoordinates(),
               c_ptr->GetCoordinates(),
               BOND_LENGTH, s_tetrahedral_angle, s_tetrahedral_angle
             );
    }

    //! @brief returns the current Date using the PDB format
    //! @param TIME time (date) information to be placed in the header
    //! @return string of the current Date using the PDB format
    std::string Factory::ConvertTimeToPDBDate( const util::Time &TIME)
    {
      // split the time as date string by " "
      storage::Vector< std::string> split_date( util::SplitString( TIME.GetTimeAsDate(), " "));

      // make the month uppercase
      for( size_t i( 0); i != 3; ++i)
      {
        split_date( 1)[ i] = toupper( split_date( 1)[ i]);
      }

      // return the formatted date
      return split_date( 2) + "-" + split_date( 1) + "-" + split_date( 4)[ 2] + split_date( 4)[ 3];
    }

    //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
    //! The types and sizes of SSEs to be considered are determined by the additional parameters
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @param MODEL model in the handler to build
    //! @return Protein Model read from given Handler
    assemble::ProteinModel Factory::ProcessModel
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const Model &MODEL
    ) const
    {
      // model
      assemble::ProteinModel protein_model;

      // structured chain
      const storage::Map< char, storage::List< Residue> > structured_protein_chains( MODEL.GetStructuredChains());

      // secondary structure element definitions
      const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
        sse_definitions( HANDLER.GetSSEStructure());

      // get iterator on pdb ssestructure
      for
      (
        storage::Map< char, storage::List< Residue> >::const_iterator
          chain_itr( structured_protein_chains.Begin()), chain_itr_end( structured_protein_chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        util::ShPtr< assemble::Chain> sp_chain;

        const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >::const_iterator
          sse_itr( sse_definitions.Find( chain_itr->first));
        if( sse_itr == sse_definitions.End())
        {
          BCL_MessageVrb( "cannot find sse definitions for chain: " + util::Format()( chain_itr->first));
          //create each chain in pdb
          sp_chain = util::ShPtr< assemble::Chain>
          (
            ChainFromPDBSSElementsAndSequence
            (
              chain_itr->first, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >(),
              chain_itr->second,
              SSE_TYPE_MINSIZE,
              HANDLER.GetHead().GetPDBID()
            )
          );
        }
        else
        {
          //create each chain in pdb
          sp_chain = util::ShPtr< assemble::Chain>
          (
            ChainFromPDBSSElementsAndSequence
            (
              chain_itr->first, sse_itr->second,
              chain_itr->second,
              SSE_TYPE_MINSIZE,
              HANDLER.GetHead().GetPDBID()
            )
          );
        }

        // insert chain into model
        protein_model.Insert( sp_chain);
      }

      // try to get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane( HANDLER.GetHead().GetMembrane());

      // if a membrane is defined
      if( sp_membrane.IsDefined())
      {
        // add to protein model data
        util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
        sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
        protein_model.SetProteinModelData( sp_data);

        // set pdb environments
        sspred::PDB::SetEnvironmentTypes( protein_model);
      }

      // if the biomolecule flag was set
      if( GetFlagBiomolecule()->GetFlag())
      {
        // get any bio transformation matrices
        const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_matrices
        (
          HANDLER.GetHead().GetBioTransformationMatrices( HANDLER.GetProteinChains())
          [
            GetFlagBiomolecule()->GetFirstParameter()->GetNumericalValue< size_t>()
          ]
        );

        // if matrices were found
        if( !bio_matrices.IsEmpty())
        {
          // construct the protein model multiplier from the first biomolecule
          util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
          (
            new assemble::ProteinModelMultiplier( bio_matrices, protein_model)
          );

          // add to protein model data
          util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
          sp_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
          protein_model.SetProteinModelData( sp_data);
        }
        else
        {
          BCL_MessageCrt
          (
            "no biomatrices in REMARK 350 in pdb given, assuming that all atoms in pdb belong to biological molecule"
          );
        }
      }

      //return model
      return protein_model;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_handler.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_site.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! Flag to choose which helix classes to process
    //! default is helix class 1 (right handed alpha helix) if there are no classes given in the list
    util::ShPtr< command::FlagInterface> &Handler::GetFlagHelixClasses()
    {
      static util::ShPtr< command::FlagInterface> s_flag_helix_class
      (
        new command::FlagDynamic
        (
          "helix_classes",
          "list of helix classes to be read in, by default, only 1-right handed helix will be considered",
          command::Parameter
          (
            "helix_class_code",
            "helix class form the pdb format",
            command::ParameterCheckRanged< size_t>( 1, 10),
            "1"
          ),
          0, 10
        )
      );

      return s_flag_helix_class;
    }

    //! Flag to trigger the merging of overlapping SSE definitions - although they should not occur, there are still
    //! many pdb files that do not adhere to the format standard
    util::ShPtr< command::FlagInterface> &Handler::GetFlagMergeOverlappingSSEs()
    {
      static util::ShPtr< command::FlagInterface> s_flag_merge
      (
        new command::FlagStatic
        (
          "merge_overlapping_sses",
          "merge overlapping sses given in a pdb file. Although overlapping SSEs are not pdb format conform, they occur"
          " quite often. If two SSE definitions are of the same type and they overlap, they will be merged into one."
        )
      );

      return s_flag_merge;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Handler::Handler( const bool IGNORE_CLASH) :
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine())
    {
    }

    //! construct PDBRreader from the std::string, containing the pdb file. It reads all information and performs a
    //! SequenceCheck. NOT appropiate when there is no SEQRES information in the pdb file!!!
    //! In that case, use default constructor and call function Read().
    Handler::Handler( std::istream &ISTREAM, const bool IGNORE_CLASH) :
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine())
    {
      Read( ISTREAM);
    }

    //! copy constructor
    Handler::Handler( const Handler &PDBREADER, const bool IGNORE_CLASH) :
      m_Head( PDBREADER.m_Head),
      m_Models( PDBREADER.m_Models),
      m_Tail( PDBREADER.m_Tail),
      m_SEQRES( PDBREADER.m_SEQRES),
      m_SSEStructure( PDBREADER.m_SSEStructure),
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine()),
      m_NonProteinChainsIDs( PDBREADER.m_NonProteinChainsIDs)
    {
    }

    //! copy constructor
    Handler *Handler::Clone() const
    {
      return new Handler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Handler::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types( GetLineTypes().Begin(), GetLineTypes().End());
      return s_line_types;
    }

    //! @brief return the Lines for a given line type
    //! @param LINE_TYPE the line type of interest
    //! @return reference to the lines - can be empty if no lines are available
    util::ShPtrList< Line> Handler::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      // from head
      lines.Append( m_Head.GetLines( LINE_TYPE));

      // from all models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        lines.Append( itr->GetLines( LINE_TYPE));
      }

      // from tail
      lines.Append( m_Tail.GetLines( LINE_TYPE));

      // end
      return lines;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Handler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Handler::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // from head
      lines.Append( m_Head.CollectLines( CRITERIUM));

      // from all models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        lines.Append( itr->CollectLines( CRITERIUM));
      }

      // from tail
      lines.Append( m_Tail.CollectLines( CRITERIUM));

      // end
      return lines;
    }

    // returns the entire sequence as one letter code of a chain CHAINID
    std::string Handler::GetSequence( const char CHAINID) const
    {
      // initialize empty sequence
      std::string sequence( "");

      // get chain for given CHAINID
      const storage::Map< char, storage::List< ResidueSimple> >::const_iterator chain_itr( m_SEQRES.Find( CHAINID));
      if( chain_itr == m_SEQRES.End())
      {
        return sequence;
      }

      // add one-letter-code of each residue in chain to sequence
      for
      (
        storage::List< ResidueSimple>::const_iterator res_itr( chain_itr->second.Begin()), res_itr_end( chain_itr->second.End());
        res_itr != res_itr_end;
        ++res_itr
      )
      {
        sequence += biol::GetAATypes().AATypeFromThreeLetterCode( res_itr->GetResidueName())->GetOneLetterCode();
      }

      // end
      return sequence;
    }

    //! @brief get the helix classes to be considered as given in the command line
    //! @return set of helix classes
    storage::Set< biol::SSType> Handler::GetHelixClassesFromCommandLine()
    {
      // get all classes specified on the commandline
      const storage::Vector< size_t> class_list( GetFlagHelixClasses()->GetNumericalList< size_t>());

      // the set of classes
      storage::Set< biol::SSType> classes;

      for
      (
        storage::Vector< size_t>::const_iterator itr( class_list.Begin()), itr_end( class_list.End());
        itr != itr_end;
        ++itr
      )
      {
        classes.Insert( biol::GetSSTypes().SSTypeFromPDBHelixClass( *itr));
      }

      // if empty, default is right handed helix
      if( classes.IsEmpty())
      {
        classes.Insert( biol::GetSSTypes().HELIX);
      }

      // end
      return classes;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all chain ids
    std::string Handler::ChainIDs() const
    {
      storage::Set< char> keys( m_SEQRES.GetKeys());
      keys.InsertElements( m_NonProteinChainsIDs);
      return std::string( keys.Begin(), keys.End());
    }

    //! @brief get chain ids for protein chains
    std::string Handler::GetProteinChains() const
    {
      storage::Set< char> all_chainids( m_SEQRES.GetKeys());

      std::string protein_chainids;
      std::set_difference
      (
        all_chainids.Begin(), all_chainids.End(),
        m_NonProteinChainsIDs.Begin(), m_NonProteinChainsIDs.End(),
        std::inserter( protein_chainids, protein_chainids.begin())
      );

      // end
      return protein_chainids;
    }

    //! @brief pushback a new line into that group
    //! @param LINE ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Handler::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( m_Head.PushBack( LINE))
      {
        return true;
      }
      if( m_Tail.PushBack( LINE))
      {
        return true;
      }

      if( m_Models.IsEmpty())
      {
        m_Models.PushBack( Model());
      }

      if( m_Models.GetSize() == 1 && m_Models.FirstElement().PushBack( LINE))
      {
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Handler::Reset()
    {
      m_Head.Reset();
      m_Models.Reset();
      m_Tail.Reset();
      m_SEQRES.Reset();
      m_SSEStructure.Reset();
      m_NonProteinChainsIDs.Reset();
    }

    //! @brief Appends a list of lines
    //! @param PDBLINES list of lines
    //! @return true, if all lines could be inserted
    bool Handler::AppendLines( const util::ShPtrList< Line> &PDBLINES)
    {
      // iterate through argument lines
      for( util::ShPtrList< Line>::const_iterator itr( PDBLINES.Begin()), itr_end( PDBLINES.End()); itr != itr_end; ++itr)
      {
        if( !PushBack( *itr))
        {
          return false;
        }
      }

      // end
      return true;
    }

    //! @brief extracts the path and tag from the provided full path
    storage::VectorND< 2, std::string> Handler::ExtractPathAndPDBTag( const std::string &FULL_PATH)
    {
      const size_t pos_slash( FULL_PATH.rfind( PATH_SEPARATOR));
      const size_t pos_point( FULL_PATH.rfind( "."));

      std::string path, pdb_tag;

      // if a dot exists
      if( pos_point != std::string::npos)
      {
        // if no slash found then just filename without any path is given
        if( pos_slash == std::string::npos)
        {
          path = ".";
          pdb_tag = FULL_PATH.substr( 0, pos_point);
        }
        // if a slash is found that is before the point then take whatever in between as the pdb tag
        else if( pos_point > pos_slash)
        {
          path = FULL_PATH.substr( 0, pos_slash);
          pdb_tag = FULL_PATH.substr( pos_slash + 1, pos_point - pos_slash - 1);
        }
        else
        {
          BCL_Exit( "Unable to deduct the path and pdb tag from the provided full path \"" + FULL_PATH + "\"", -1);
        }
      }

      return storage::VectorND< 2, std::string>( path, pdb_tag);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Handler::WriteLines( std::ostream &OSTREAM) const
    {
      return Write( OSTREAM, 0);
    }

    //! reads Handler from std::istream
    std::istream &Handler::Read( std::istream &ISTREAM)
    {
      bool unknown_linetype( false);

      // previously stored pdb lines are cleaned
      Reset();

      util::SiPtr< LineGroupInterface> current_group( &m_Head);

      while( !ISTREAM.eof())
      {
        std::string buffer;
        std::getline( ISTREAM, buffer);
        if( buffer.empty())
        {
          continue;
        }

        // create line and get line type
        const util::ShPtr< Line> current_line( new Line( buffer));
        const LineType &line_type( current_line->GetType());

        if( line_type == GetLineTypes().e_Undefined && !unknown_linetype)
        {
          BCL_MessageVrb( "File contains at least one unknown linetype! " + buffer);
          unknown_linetype = true;
          continue;
        }

        // end of pdb file
        if( line_type == GetLineTypes().END)
        {
          break;
        }

        // store the line
        if( !current_group->PushBack( current_line))
        {
          // if the current group is tail already, than the current line is found at the wrong place
          if( current_group == util::SiPtr< LineGroupInterface>( &m_Tail))
          {
            BCL_MessageCrt( "pdb file does not appear to have proper order of records!");
            BCL_MessageVrb( "ignoring line: " + current_line->GetString());
            continue;
          }

          // is it the start of a new model
          if( line_type == GetLineTypes().MODEL)
          {
            BCL_MessageStd
            (
              "reading model nr: " + current_line->GetString( GetEntryTypes().MODELSerial)
            );
            // create a new model
            m_Models.PushBack( Model());
            current_group = &m_Models.LastElement();
            continue;
          }

          // end of model - set to head
          if( line_type == GetLineTypes().ENDMDL)
          {
            current_group = &m_Head;
            continue;
          }

          // start of a model without explicit indication
          if( Model().GetTypesOfLines().Contains( line_type))
          {
            BCL_MessageDbg( "start reading model");
            // create a new model
            m_Models.PushBack( Model());
            current_group = &m_Models.LastElement();
            current_group->PushBack( current_line);
            continue;
          }

          // is it a line that should go to the tail
          if( m_Tail.GetTypesOfLines().Contains( line_type))
          {
            current_group = &m_Tail;
            m_Tail.PushBack( current_line);
            continue;
          }

          // try to insert the line back into the head group
          if( m_Head.PushBack( current_line))
          {
            BCL_MessageStd
            (
              "Header line following atom or other model-specific lines merged with header lines: "
              + current_line->GetString()
            );
            continue;
          }

          BCL_MessageVrb
          (
            "order of records in pdb is not standard conform! ignoring current line:\n" + current_line->GetString()
          );
          current_group = &m_Head;
        } // change the group
      }

      // if no model was read, insert an empty model
      if( m_Models.IsEmpty())
      {
        m_Models.PushBack( Model());
      }

      m_SEQRES = m_Head.GetSEQRESProteinChains();
      m_SSEStructure = m_Head.SSEDefinitions( m_HelixClasses, GetFlagMergeOverlappingSSEs()->GetFlag());

      // no seqres given in pdb header
      if( m_SEQRES.IsEmpty())
      {
        BCL_MessageStd
        (
          "no SEQRES information given, try to retrieve sequence from ATOM section, which requires that every residue "
          "is given in there!"
        )
        m_SEQRES = m_Models.FirstElement().GetChains();
      }

      // initialize the models
      {
        size_t model_nr( 1);
        for( storage::Vector< Model>::iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr, ++model_nr)
        {
          const storage::Map< char, storage::List< ResidueSimple> > missing_residues( m_Head.GetMissingResidues( model_nr));
          // iterate through chains
          for
          (
            storage::Map< char, storage::List< ResidueSimple> >::const_iterator
              seqres_itr( m_SEQRES.Begin()), seqres_itr_end( m_SEQRES.End());
            seqres_itr != seqres_itr_end;
            ++seqres_itr
          )
          {

            const storage::Map< char, storage::List< ResidueSimple> >::const_iterator mis_itr( missing_residues.Find( seqres_itr->first));
            itr->InitializeStructuredChain
            (
              seqres_itr->first,
              seqres_itr->second,
              mis_itr == missing_residues.End() ? storage::List< ResidueSimple>() : mis_itr->second
            );
          }
        }
      }

      //return
      return ISTREAM;
    }

    //! Writes all existing lines into std::ostrem STREAM
    std::ostream &Handler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write header
      m_Head.WriteLines( OSTREAM);

      // write all models
      const bool write_model_separators( m_Models.GetSize() != 1);
      size_t model_count( 1);

      // iterate through models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        if( write_model_separators)
        {
          Line model_line( GetLineTypes().MODEL);
          model_line.Put( GetEntryTypes().MODELSerial, model_count);
          OSTREAM << model_line.GetString() << '\n';
        }
        itr->WriteLines( OSTREAM);
        if( write_model_separators)
        {
          Line end_line( GetLineTypes().ENDMDL);
          OSTREAM << end_line.GetString() << '\n';
        }
      }

      // tail
      m_Tail.UpdateMasterRecord( m_Head);
      if( !m_Models.IsEmpty())
      {
        m_Tail.UpdateMasterRecord( m_Models.FirstElement());
      }
      m_Tail.WriteLines( OSTREAM);

      // return
      return OSTREAM;
    }

    //! @brief return a map of hetero residues
    //! @brief return all ligands
    //! @return List of ligands
    util::ShPtrList< Ligand> Handler::GetLigands() const
    {
      util::ShPtrList< Ligand> ligands;

      // connections
      const storage::Map< size_t, storage::Set< size_t> > connections( m_Tail.GetConnections());

      // fullnames
      const storage::Map< std::string, std::string> fullnames( m_Head.GetHetFullname());

      // formula
      const storage::Map< std::string, std::string> formulas( m_Head.GetHetFormula());

      // search all HET entries
      const util::ShPtrList< Line> het_lines( GetLines( GetLineTypes().HET));

      // iterate through all HET lines
      for( util::ShPtrList< Line>::const_iterator itr( het_lines.Begin()), itr_end( het_lines.End()); itr != itr_end; ++itr)
      {
        // new residues for that HET entry
        // select HETATM lines corresponding to the HET entries
        LineCriterium criteria;
        criteria.SetMeetAllCriteria( true);
        criteria.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, ( *itr)->GetString( GetEntryTypes().HETSequenceID));
        criteria.AddCriterium( GetEntryTypes().HETATMInsertionCode    , ( *itr)->GetString( GetEntryTypes().HETInsertionCode));

        // locate the lines
        const size_t expected_nr_lines( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETNumberAtoms));
        const util::ShPtrList< Line> hetatm_lines( m_Models.FirstElement().CollectLines( criteria, ( *itr)->GetChar( GetEntryTypes().HETChainID)));

        if( expected_nr_lines != hetatm_lines.GetSize())
        {
          BCL_MessageCrt
          (
            "different number of HETATM lines found than expected: " + ( *itr)->GetString() +
            util::Format()( hetatm_lines.GetSize())
          );
        }

        // new ligand
        util::ShPtr< Ligand> sp_ligand( new Ligand( hetatm_lines));

        // connections
        sp_ligand->AddConnections( connections);

        // name
        {
          const storage::Map< std::string, std::string>::const_iterator name_itr( fullnames.Find( sp_ligand->GetResidueName()));
          if( name_itr == fullnames.End())
          {
            BCL_MessageVrb( "no fullname for ligand available: " + sp_ligand->GetResidueName());
            sp_ligand->SetFullname( "NO_HETNAM");
          }
          else
          {
            sp_ligand->SetFullname( name_itr->second);
          }
        }
        // formula
        {
          const storage::Map< std::string, std::string>::const_iterator formul_itr( formulas.Find( sp_ligand->GetResidueName()));
          if( formul_itr == formulas.End())
          {
            BCL_MessageVrb( "no formula for ligand available: " + sp_ligand->GetResidueName());
          }
          else
          {
            sp_ligand->SetFormula( formul_itr->second);
          }
        }

        ligands.PushBack( sp_ligand);
      }

      // end
      return ligands;
    }

    //! @brief return all sites
    //! @return list of all sites
    util::ShPtrList< Site> Handler::GetSites() const
    {
      static const std::string s_site_identifier_string(    "SITE_IDENTIFIER:");
      static const std::string s_site_evidence_code_string( "EVIDENCE_CODE:");
      static const std::string s_site_description_string(   "SITE_DESCRIPTION:");

      util::ShPtrList< Site> sites;

      // ligands
      const util::ShPtrList< Ligand> ligands( GetLigands());

      // locate all remark 800 lines
      LineCriterium criterium_remark;
      criterium_remark.AddCriterium( GetEntryTypes().REMARK_Number, 800);
      const util::ShPtrList< Line> remark_lines( m_Head.CollectLines( criterium_remark, GetLineTypes().REMARK));

      // iterate through all lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( remark_lines.Begin()), itr_end( remark_lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // first line needs to be a site identifier string
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifierString) != s_site_identifier_string)
        {
          continue;
        }

        // site identifier / name
        const std::string site_identifier( util::TrimString( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifier)));

        // next line
        ++itr;
        if( itr == itr_end) break;
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteEvidenceCodeString) != s_site_evidence_code_string)
        {
          BCL_MessageCrt( "missing REMARK 800 entries for site identifier: " + site_identifier);
          continue;
        }

        // evidence code
        const Site::EvidenceCodeEnum evidence_code
        (
          util::TrimString( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteEvidenceCode))
        );

        // next line
        ++itr;
        if( itr == itr_end) break;
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteDescriptionString) != s_site_description_string)
        {
          BCL_MessageCrt( "missing REMARK 800 entries for site identifier: " + site_identifier);
          continue;
        }

        std::string site_description( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteDescription));

        // site description over multiple lines
        ++itr;
        while( itr != itr_end)
        {
          const std::string current_description( ( *itr)->GetString( GetEntryTypes().REMARK_String));
          if( util::TrimString( current_description).empty())
          {
            break;
          }

          // sometimes there is no empty REMARK 800 line between two SITE definitions, so they have to be put back
          if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifierString) == s_site_identifier_string)
          {
            --itr;
            break;
          }
          site_description += ' ';
          site_description += current_description;
          ++itr;
        }

        // insert new site
        util::ShPtr< Site> sp_site( new Site( site_identifier, evidence_code, site_description));
        sp_site->FindLigand( ligands);

        // find all SITE residues
        LineCriterium criterium_site;
        criterium_site.AddCriterium( GetEntryTypes().SITEName, site_identifier);

        // locate site
        const util::ShPtrList< Line> site_lines( m_Head.CollectLines( criterium_site, GetLineTypes().SITE));
        if( site_lines.IsEmpty())
        {
          BCL_MessageCrt( "required SITE lines are not found: " + site_identifier);

          // last line of REMARK 800
          if( itr == itr_end) break;

          continue;
        }
        // number of residues
        const size_t number_residues( site_lines.FirstElement()->GetNumericalValue< size_t>( GetEntryTypes().SITENumberResidues));
        const size_t max_number_residues_per_site( 4);

        // only soft assert, since the SiteNumberResidues can only have two digits
        if( ( number_residues - 1) / max_number_residues_per_site + 1 != site_lines.GetSize())
        {
          BCL_MessageCrt
          (
            "number residues in SITE " + site_identifier + " does not match the number of lines: " +
            util::Format()( number_residues) + " do not fit in nr lines: " + util::Format()( site_lines.GetSize())
          );
        }

        size_t res_nr( 0);
        size_t previous_line_nr( 0);
        bool record_error( false);
        // iterate through all site lines
        for
        (
          util::ShPtrList< Line>::const_iterator sitel_itr( site_lines.Begin()), sitel_itr_end( site_lines.End());
          sitel_itr != sitel_itr_end;
          ++sitel_itr
        )
        {
          const size_t current_line_nr( ( *sitel_itr)->GetNumericalValue< size_t>( GetEntryTypes().SITESequenceNumber));
          if( previous_line_nr + 1 != current_line_nr)
          {
            record_error = true;
            break;
          }

          // update previours line nr
          previous_line_nr = current_line_nr;

          for
          (
            EntryTypes::const_iterator ent_res_itr( GetEntryTypes().SITEResidueName1.GetIterator()), ent_res_itr_last( GetEntryTypes().SITEResidueName4.GetIterator()),
              ent_chain_itr( GetEntryTypes().SITEChainID1.GetIterator()),
              ent_seq_itr( GetEntryTypes().SITEResidueSequenceID1.GetIterator()),
              ent_icode_itr( GetEntryTypes().SITEResidueInsertionCode1.GetIterator());
            ent_res_itr <= ent_res_itr_last && res_nr < number_residues;
            ent_res_itr += 4, ent_chain_itr += 4, ent_seq_itr += 4, ent_icode_itr += 4, ++res_nr
          )
          {
            // add residue if it is part of the chain to the chain residues
            const ResidueSimple site_residue
            (
              ( *sitel_itr)->GetString( *ent_res_itr), ( *sitel_itr)->GetChar( *ent_chain_itr),
              ( *sitel_itr)->GetNumericalValue< int>( *ent_seq_itr), ( *sitel_itr)->GetChar( *ent_icode_itr)
            );

            if( IsInChain( site_residue))
            {
              sp_site->AddChainResidue( site_residue);
            }
            else
            {
              sp_site->AddHetatmResidue( site_residue);
            }
          }
        }

        // record number not continuous, i.e. if the site name is not unique
        if( record_error)
        {
          BCL_MessageCrt
          (
            "problems parsing SITE lines for " + site_identifier + " since the sequence numbering is off"
          )
          // last line of REMARK 800
          if( itr == itr_end) break;
          continue;
        }

        // insert
        sites.PushBack( sp_site);

        // last line of REMARK 800
        if( itr == itr_end) break;
      }

      return sites;
    }

    //! @brief check if given residue is in one of the chains (if not might be in the HETATM section)
    //! @param RESIDUE the residue to check for
    //! @return true, if that residue in found in one of the chains, false oterhwise
    bool Handler::IsInChain( const ResidueInterface &RESIDUE) const
    {
      // first model
      const Model &first_model( m_Models.FirstElement());

      const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( first_model.GetHETATMLines().Find( RESIDUE.GetChainID()));

      // if it is not within HETATM lines
      if( chain_itr == first_model.GetHETATMLines().End())
      {
        return true;
      }

      // collect all lines in hetatm section for that residue - if it is empty, residue must be part of chain
      return LineCriterium::Filter( chain_itr->second, RESIDUE.GetCriterium( GetLineTypes().HETATM)).IsEmpty();
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_head.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_flag_interface.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_printer_membrane.h"
#include "pdb/bcl_pdb_residue_simple.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  ///////////
  // types //
  ///////////

    //! @brief operator for checking if one line is less than the other one
    //! @param LINE_LHS line on left-hand side
    //! @param LINE_RHS line on right-hand side
    //! @return if one line is less than the other one
    bool Head::RemarkLineLessThan::operator()
    (
      const util::ShPtr< Line> &LINE_LHS,
      const util::ShPtr< Line> &LINE_RHS
    ) const
    {
      // if both types are remark
      if
      (
        LINE_LHS.IsDefined() && LINE_RHS.IsDefined() &&
        LINE_LHS->GetType() == GetLineTypes().REMARK && LINE_RHS->GetType() == GetLineTypes().REMARK
      )
      {
        // compare remark #
        return LINE_LHS->GetString( GetEntryTypes().REMARK_Number) < LINE_RHS->GetString( GetEntryTypes().REMARK_Number);
      }

      // end
      return false;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Head::s_Instance
    (
      GetObjectInstances().AddInstance( new Head())
    );

    //! BCL pdb identifier
    const std::string &Head::GetBCLPdbID()
    {
      static const std::string s_bcl_pdb_id( "BCL ");
      return s_bcl_pdb_id;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Head::Head()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Head
    Head *Head::Clone() const
    {
      return new Head( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Head::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Head::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().HEADER.GetIterator(), GetLineTypes().MTRIX3.GetIterator() + 1
      );

      // end
      return s_line_types;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Head::GetLines( const LineType &LINE_TYPE) const
    {
      // search for that line type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      if( itr != m_Lines.End())
      {
        return itr->second;
      }

      return util::ShPtrList< Line>();
    }

    //! @brief count the number of lines of given TYPE used for master record
    //! @param LINE_TYPE the line type
    //! @return the number of lines of that type found
    size_t Head::Count( const LineType &LINE_TYPE) const
    {
      // search for that line type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      if( itr != m_Lines.End())
      {
        return itr->second.GetSize();
      }

      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Head::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // iterate through map
      for
      (
        storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End());
        itr != itr_end;
        ++itr
      )
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // end
      return lines;
    }

    //! @brief locate lines of given criterium and line type
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @param LINE_TYPE only consider this line type
    //! @return lines that are considered by criterium
    util::ShPtrList< Line>
    Head::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const LineType &LINE_TYPE) const
    {
      // find lines of that type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      // no lines of that type
      if( itr == m_Lines.End())
      {
        return util::ShPtrList< Line>();
      }

      // filter the lines according to the criterium
      return LineCriterium::Filter( itr->second, CRITERIUM);
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Head::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( GetTypesOfLines().Contains( LINE->GetType()))
      {
        m_Lines[ LINE->GetType()].PushBack( LINE);
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Head::Reset()
    {
      m_Lines.Reset();
    }

    //! extract all transformation matrices
    storage::Vector< math::TransformationMatrix3D> Head::GetTransformationMatrices() const
    {
      // get all the mtrix lines
      const util::ShPtrList< Line> &matrix_lines1( GetLines( GetLineTypes().MTRIX1));
      const util::ShPtrList< Line> &matrix_lines2( GetLines( GetLineTypes().MTRIX2));
      const util::ShPtrList< Line> &matrix_lines3( GetLines( GetLineTypes().MTRIX3));

      // holds the transformations
      storage::Vector< math::TransformationMatrix3D> transformationmatrices;

      // check that there is one rwo for matrix
      if( matrix_lines1.GetSize() != matrix_lines2.GetSize() || matrix_lines1.GetSize() != matrix_lines3.GetSize())
      {
        BCL_MessageCrt( "different number of MTRIXn lines!");
        return transformationmatrices;
      }

      // iterate over each row for all matrices
      for
      (
        util::ShPtrList< Line>::const_iterator
          matrix1_line_itr( matrix_lines1.Begin()), matrix1_line_itr_end( matrix_lines1.End()),
          matrix2_line_itr( matrix_lines2.Begin()), matrix2_line_itr_end( matrix_lines2.End()),
          matrix3_line_itr( matrix_lines3.Begin()), matrix3_line_itr_end( matrix_lines3.End());
        matrix1_line_itr != matrix1_line_itr_end &&
        matrix2_line_itr != matrix2_line_itr_end &&
        matrix3_line_itr != matrix3_line_itr_end;
        ++matrix1_line_itr, ++matrix2_line_itr, ++matrix3_line_itr
      )
      {
        // construct matrix
        linal::Matrix< double> new_matrix( 4, 4, double( 0));

        // first col
        new_matrix( 0, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M11);
        new_matrix( 1, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M12);
        new_matrix( 2, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M13);
        new_matrix( 3, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_V1);

        // second col
        new_matrix( 0, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M21);
        new_matrix( 1, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M22);
        new_matrix( 2, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M23);
        new_matrix( 3, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_V2);

        // last col
        new_matrix( 0, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M31);
        new_matrix( 1, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M32);
        new_matrix( 2, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M33);
        new_matrix( 3, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_V3);

        // insert
        transformationmatrices.PushBack( math::TransformationMatrix3D( new_matrix));
      }

      // end
      return transformationmatrices;
    }

    //! read one column of a BIOMT matrix remark 350 line into a matrix
    //! @param LINE the pdb line to read from
    //! @param BIOMT_FIRST_COL_ENTRY_TYPE any of REMARK_350_BIOMT1_M11, REMARK_350_BIOMT1_M21, REMARK_350_BIOMT1_M31
    //! @param COL_NR which col is read
    //! @param MATRIX the matrix values are added to
    void ReadREMARK_350_BIOMTMatrixCol
    (
      const Line &LINE,
      const EntryType &BIOMT_FIRST_COL_ENTRY_TYPE,
      const size_t COL_NR,
      linal::Matrix< double> &MATRIX
    )
    {
      EntryTypes::const_iterator itr( BIOMT_FIRST_COL_ENTRY_TYPE.GetIterator());
      MATRIX( 0, COL_NR) = LINE.GetNumericalValue< double>( *itr);
      MATRIX( 1, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
      MATRIX( 2, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
      MATRIX( 3, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
    }

    //! @brief the transformation matrices to generate different bio molecules by applying transformations to different chains
    //! @param CHAIN_IDS all chain ids that should be considered
    //! @return Map of biomolecule number to a vector of transformations of chainid it should be applied to, the new chain id and the transformation
    storage::Map< size_t, storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > >
    Head::GetBioTransformationMatrices( const std::string &CHAIN_IDS) const
    {
      storage::Map< size_t, storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > >
        transformation_matrices;

      // search for all REMARK  350 lines
      LineCriterium criterium;
      criterium.AddCriterium( GetEntryTypes().REMARK_Number, 350);

      // acquire matrix lines
      const util::ShPtrList< Line> matrix_lines( CollectLines( criterium, GetLineTypes().REMARK));

      const std::string biomolecule_identifier( "BIOMOLECULE:");
      const std::string chains_identifier( "CHAINS:");

      // if there are no chain ids
      if( CHAIN_IDS.empty())
      {
        return transformation_matrices;
      }
      char next_chain_id( *( --CHAIN_IDS.end()));

      // create set of transformed chains
      storage::Set< char> transformed_chains;

      util::ShPtrList< Line>::const_iterator line_itr( matrix_lines.Begin()), line_itr_end( matrix_lines.End());

      // iterate over all lines
      while( line_itr != line_itr_end)
      {
        // check for BIOMOLECULE string
        while
        (
             line_itr != line_itr_end
          && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BiomoleculeIdentifier) != biomolecule_identifier
        )
        {
          ++line_itr;
        }
        if( line_itr == line_itr_end)
        {
          BCL_MessageVrb( "could not find any more " + biomolecule_identifier);
          break;
        }

        // report that BIOMOLECULE was found
        const size_t current_biomolecule
        (
          ( *line_itr)->GetNumericalValue< size_t>( GetEntryTypes().REMARK_350_BiomoleculeNumber)
        );
        BCL_MessageStd( "found " + biomolecule_identifier + ' ' + util::Format()( current_biomolecule));

        // goto next line
        ++line_itr;

        // check for apply identifier
        while( line_itr != line_itr_end)
        {
          // break if new biomolecule appears
          if( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BiomoleculeIdentifier) == biomolecule_identifier)
          {
            break;
          }

          // check that the line contains a chains_identifier
          if( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsIdentifier) != chains_identifier)
          {
            ++line_itr;
            continue;
          }

          // extract the chains the transformation has to be applied to
          std::string chains_string( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsList));

          while( ++line_itr != line_itr_end && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsIdentifier) == chains_identifier)
          {
            chains_string += ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsList);
          }

          // sort the chain ids and erase every non alpha character
          std::sort( chains_string.begin(), chains_string.end());
          chains_string.erase
          (
            std::remove_if( chains_string.begin(), chains_string.end(), std::not1( std::ptr_fun( &isalpha))),
            chains_string.end()
          );

          // remove all chain id, that are not in the argument
          {
            std::string arg_chainids( CHAIN_IDS);
            std::sort( arg_chainids.begin(), arg_chainids.end());
            std::string filtered_chainids;
            std::set_intersection
            (
              chains_string.begin(), chains_string.end(),
              arg_chainids.begin(), arg_chainids.end(),
              std::inserter( filtered_chainids, filtered_chainids.begin())
            );

            chains_string = filtered_chainids;
          }

          BCL_MessageCrt
          (
            "Chains to apply transformations to:\n" + util::Format()( chains_string)
          );

          // should be at BIOMT1
          // read all matrices (one matrix consists of BIMT1, 2 and 3)
          while
          (
               line_itr != line_itr_end
            && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT1_Identifier) == "BIOMT1"
          )
          {
            linal::Matrix< double> new_matrix( 4, 4, double( 0));
            new_matrix( 3, 3) = 1;

            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT1_Identifier) == "BIOMT1",
              "wrong identifier; should be \"BIOMT1\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT1_M11, 0, new_matrix);

            //go to next line
            ++line_itr;
            BCL_Assert( line_itr != line_itr_end, "ran out of lines");
            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT2_Identifier) == "BIOMT2",
              "wrong identifier; should be \"BIOMT2\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT2_M21, 1, new_matrix);

            //go to next line
            ++line_itr;
            BCL_Assert( line_itr != line_itr_end, "ran out of lines");
            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT3_Identifier) == "BIOMT3",
              "wrong identifier; should be \"BIOMT3\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT3_M31, 2, new_matrix);

            // iterate over all chains this matrix will be applied to
            for
            (
              std::string::const_iterator itr( chains_string.begin()), itr_end( chains_string.end());
              itr != itr_end; ++itr
            )
            {
              // initialize next chain id
              char new_chain_id;

              // if this chain has already been transformed
              if( transformed_chains.Contains( *itr))
              {
                // assign to next chain
                do
                {
                  ++next_chain_id;
                }
                while( !isalpha( next_chain_id));

                new_chain_id = next_chain_id;
              }
              // this is the first transformation for this chain
              else
              {
                // the new chain will have the same chain id
                new_chain_id = *itr;
                transformed_chains.Insert( *itr);
              }

              // insert the transformation matrix and the corresponding chains it should be applied
              transformation_matrices[ current_biomolecule].PushBack
              (
                storage::Triplet< char, char, math::TransformationMatrix3D>
                (
                  *itr,
                  new_chain_id,
                  math::TransformationMatrix3D( new_matrix)
                )
              );
            }
            //go to next line
            ++line_itr;
          } // loop over matrices
        } // loop over "apply to chain"
      } // loop over biomolecules

      return transformation_matrices;
    }

    //! @brief get the membrane from the remark lines
    //! @return the membrane from the remark lines
    util::ShPtr< biol::Membrane> Head::GetMembrane() const
    {
      // search for all membrane lines
      LineCriterium criterium;
      criterium.AddCriterium( GetEntryTypes().REMARK_Number, size_t( PrinterMembrane::s_RemarkNumber));

      // acquire matrix lines
      const util::ShPtrList< Line> membrane_lines( CollectLines( criterium, GetLineTypes().REMARK));

      // if no membrane lines or wrong size
      if( membrane_lines.GetSize() != 4)
      {
        return util::ShPtr< biol::Membrane>();
      }

      // skip the first line
      util::ShPtrList< Line>::const_iterator line_itr( membrane_lines.Begin());
      ++line_itr;

      // second line is normal
      const storage::Vector< std::string> normal_entries( util::SplitString( ( *line_itr)->GetString()));
      BCL_Assert( normal_entries.GetSize() == 7, "There should be 7 membrane normal lines in the PDB");
      const linal::Vector3D normal
      (
        util::ConvertStringToNumericalValue< double>( normal_entries( 4)),
        util::ConvertStringToNumericalValue< double>( normal_entries( 5)),
        util::ConvertStringToNumericalValue< double>( normal_entries( 6))
      );
      ++line_itr;

      // third line is center
      const storage::Vector< std::string> center_entries( util::SplitString( ( *line_itr)->GetString()));
      BCL_Assert( center_entries.GetSize() == 7, "There should be 7 membrane center lines in the PDB");
      const linal::Vector3D center
      (
        util::ConvertStringToNumericalValue< double>( center_entries( 4)),
        util::ConvertStringToNumericalValue< double>( center_entries( 5)),
        util::ConvertStringToNumericalValue< double>( center_entries( 6))
      );
      ++line_itr;

      // last line is thickness, use command line thickness if set
      storage::Vector< double> thickness;

      if( biol::Membrane::GetFlagMembrane()->GetFirstParameter()->GetWasSetInCommandLine())
      {
        thickness = biol::Membrane::GetCommandLineMembrane().GetThicknesses();
      }
      else
      {
        const storage::Vector< std::string> thickness_entries( util::SplitString( ( *line_itr)->GetString()));
        BCL_Assert( thickness_entries.GetSize() >= 9, "There should be 9 membrane thickness lines in the PDB");
        thickness = storage::Vector< double>::Create
        (
          util::ConvertStringToNumericalValue< double>( thickness_entries( 4)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 5)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 6)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 7)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 8))
        );
      }

      // end
      return util::ShPtr< biol::Membrane>( new biol::Membrane( thickness, normal, center));
    }

    //! @brief gets the pdb ID read in from the pdb file
    //! @return the pdb ID read in from the pdb file
    std::string Head::GetPDBID() const
    {
      // find header lines
      const util::ShPtrList< Line> header_lines( GetLines( GetLineTypes().HEADER));

      //
      if( header_lines.IsEmpty())
      {
        return GetBCLPdbID();
      }

      // get pdb id from header line
      const std::string pdb_id( header_lines.FirstElement()->GetString( GetEntryTypes().HEADERIDCode));

      // check that it is not empty
      if( util::TrimString( pdb_id).empty())
      {
        return GetBCLPdbID();
      }

      // return the pdb id
      return pdb_id;
    }

    //! @brief full name for all het residues
    //! @return map with key res name and data fullname
    storage::Map< std::string, std::string> Head::GetHetFullname() const
    {
      storage::Map< std::string, std::string> het_fullname;

      // locate all lines HETNAM
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator find_itr( m_Lines.Find( GetLineTypes().HETNAM));
      if( find_itr == m_Lines.End())
      {
        return het_fullname;
      }

      // store the current line number for each HETNAM short name
      storage::Map< std::string, size_t> nam_lines;

      // iterate over lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( find_itr->second.Begin()), itr_end( find_itr->second.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string current( ( *itr)->GetString( GetEntryTypes().HETNAMIdentifier));
        const size_t line_count( ++nam_lines[ current]);
        if( line_count > 1)
        {
          const size_t current_linenumber( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETNAMContinuation));
          // if the current line number is not defined, it would be the second one without a line number, which is implicitly '1'
          // or if both line numbers are the same
          if( !util::IsDefined( current_linenumber) || current_linenumber != line_count)
          {
            BCL_MessageCrt
            (
              "found HETNAM line that does not fit in continuation count: " + ( *itr)->GetString()
            )
            continue;
          }
        }
        het_fullname[ current].append( util::TrimString( ( *itr)->GetString( GetEntryTypes().HETNAMText)));
      }

      // end
      return het_fullname;
    }

    //! @brief formula for all het residues
    //! @return map with key res name and data formula
    storage::Map< std::string, std::string> Head::GetHetFormula() const
    {
      storage::Map< std::string, std::string> het_formula;

      // locate all lines FORMUL
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator find_itr( m_Lines.Find( GetLineTypes().FORMUL));
      if( find_itr == m_Lines.End())
      {
        return het_formula;
      }

      // iterate over lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( find_itr->second.Begin()), itr_end( find_itr->second.End());
        itr != itr_end;
        ++itr
      )
      {
        het_formula[ ( *itr)->GetString( GetEntryTypes().FORMULIdentifier)].append( ( *itr)->GetString( GetEntryTypes().FORMULChemicalFormula));
      }

      // end
      return het_formula;
    }

    //! @brief create the chains as they were given in the seqres
    //! @return map with chainid as key and list of residues as data
    storage::Map< char, storage::List< ResidueSimple> > Head::GetSEQRESProteinChains() const
    {
      storage::Map< char, storage::List< ResidueSimple> > chains;

      // retrieve seqres lines for all chains
      const util::ShPtrList< Line> seqreslines( GetLines( GetLineTypes().SEQRES));

      // store number of residues for that chain
      storage::Map< char, size_t> number_res_in_chain;

      // iterate over pdb seqres lines
      for
      (
        util::ShPtrList< Line>::const_iterator seq_line_itr( seqreslines.Begin()), seq_line_itr_end( seqreslines.End());
        seq_line_itr != seq_line_itr_end;
        ++seq_line_itr
      )
      {
        const char current_chain_id( ( *seq_line_itr)->GetChar( GetEntryTypes().SEQRESChainID));
        storage::Map< char, size_t>::const_iterator num_res_itr( number_res_in_chain.Find( current_chain_id));
        // check if there is already a number of residues stored for that chain
        if( num_res_itr == number_res_in_chain.End())
        {
          // store number of residues for that chain
          number_res_in_chain[ current_chain_id] =
            ( *seq_line_itr)->GetNumericalValue< size_t>( GetEntryTypes().SEQRESNrOfResiduesInChain);
        }
        else
        {
          // assert that all seqres lines for that chain have the same number of residues
          BCL_Assert
          (
            num_res_itr->second ==
              ( *seq_line_itr)->GetNumericalValue< size_t>( GetEntryTypes().SEQRESNrOfResiduesInChain),
            "number of residues in seqreslines for that chain are not all equal: \'" + util::Format()( current_chain_id)
            + "\'"
          );
        }

        // iterate over all 13 residues in a seqres line
        for
        (
          EntryTypes::const_iterator
            res_itr( GetEntryTypes().SEQRESName_1.GetIterator()),
            res_itr_end( ++GetEntryTypes().SEQRESName_13.GetIterator());
          res_itr != res_itr_end;
          ++res_itr
        )
        {
          const std::string res_name( util::TrimString( ( *seq_line_itr)->GetString( *res_itr)));

          // is this residue entry given in this SEQRES line
          if( res_name.size() < 3)
          {
            break;
          }

          chains[ current_chain_id].PushBack( ResidueSimple( res_name, current_chain_id));
        }
      }

      // check that number of res in seqresline SequenceNrOfResiduesInChain does match the actual number of residues
      for
      (
        storage::Map< char, size_t>::const_iterator
          itr( number_res_in_chain.Begin()), itr_end( number_res_in_chain.End());
        itr != itr_end;
        ++itr
      )
      {
        // find the chain with residues for that chain id
        const storage::Map< char, storage::List< ResidueSimple> >::const_iterator chain_itr( chains.Find( itr->first));

        // check if there have been residues for that chain at all
        if( chain_itr == chains.End())
        {
          BCL_MessageCrt
          (
            "Although SEQRES lines for this chain were found: \'" + util::Format()( itr->first) + "\' no residues " +
            "have been recognized. Probably it is a DNA or RNA chain which is not handled yet"
          );
        }

        // check that number of residues for that chain match
        else if( chain_itr->second.GetSize() != itr->second)
        {
          BCL_MessageCrt
          (
            "number of pretended residues in SEQRES lines for that chain \'" + util::Format()( itr->first) + "\' " +
            "does not match with actual number of residues: " + util::Format()( itr->second) + " != " +
            util::Format()( chains[ itr->first].GetSize())
          );
        }
      }

      // end
      return chains;
    }

    //! @brief map of residues that could not be located in the experiments using REMARK 465
    //! @param MODEL_NUMBER missing residues for the given model number
    //! @return map of chains and ShPtrList of Residues that are not in the ATOM lines
    storage::Map< char, storage::List< ResidueSimple> > Head::GetMissingResidues( const size_t MODEL_NUMBER) const
    {
      // residues
      storage::Map< char, storage::List< ResidueSimple> > residues;

      // search for all REMARK  465 lines
      LineCriterium criteria;
      criteria.SetMeetAllCriteria( true);
      criteria.AddCriterium( GetEntryTypes().REMARK_Number, 465);

      // acquire missing residue lines
      const util::ShPtrList< Line> remark_465_lines( CollectLines( criteria, GetLineTypes().REMARK));
      util::ShPtrList< Line> remark_465_lines_no_header;
      {
        // header identifier
        const std::string header_identifier( "M RES C SSSEQI");
        util::ShPtrList< Line>::const_iterator line_itr( remark_465_lines.Begin()), line_itr_end( remark_465_lines.End());
        // check for header string
        while
        (
             line_itr != line_itr_end
          && util::TrimString( ( *line_itr)->GetString( GetEntryTypes().REMARK_String)) != header_identifier
        )
        {
          ++line_itr;
        }
        if( line_itr != line_itr_end)
        {
          ++line_itr;
          remark_465_lines_no_header = util::ShPtrList< Line>( line_itr, line_itr_end);
        }
      }

      criteria.AddCriterium( GetEntryTypes().REMARK_465_ModelNumber, MODEL_NUMBER);
      util::ShPtrList< Line> missing_res_lines( LineCriterium::Filter( remark_465_lines_no_header, criteria));

      if( missing_res_lines.IsEmpty() && MODEL_NUMBER == 1)
      {
        criteria.Reset();
        criteria.AddCriterium( GetEntryTypes().REMARK_465_ModelNumber, "");
        missing_res_lines = LineCriterium::Filter( remark_465_lines_no_header, criteria);
      }

      if( missing_res_lines.IsEmpty())
      {
        return residues;
      }

      // as long as their are remaining lines
      for
      (
        util::ShPtrList< Line>::const_iterator
          line_itr( missing_res_lines.Begin()), line_itr_end( missing_res_lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // get residue name
        const std::string residue_name( ( *line_itr)->GetString( GetEntryTypes().REMARK_465_ResidueName));

        // check that it is valid residue name
        if( !biol::GetAATypes().AATypeFromThreeLetterCode( residue_name).IsDefined())
        {
          BCL_MessageCrt( "found undefined residues in missing residue line: " + ( *line_itr)->GetString());
          // end
          return residues;
        }

        const char chain_id( ( *line_itr)->GetChar( GetEntryTypes().REMARK_465_ChainID));
        // current residue
        const ResidueSimple current_residue
        (
          residue_name,
          chain_id,
          ( *line_itr)->GetNumericalValue< int>( GetEntryTypes().REMARK_465_ResidueSequenceID),
          ( *line_itr)->GetChar( GetEntryTypes().REMARK_465_InsertionCode)
        );

        // insert
        residues[ chain_id].PushBack( current_residue);
      }

      // end
      return residues;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Head::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through map
      for
      (
        storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // if this is a remark section
        if( itr->first == GetLineTypes().REMARK)
        {
          // copy the lines
          util::ShPtrList< Line> remark_lines( itr->second);

          // sort by remark #
          remark_lines.Sort( RemarkLineLessThan());

          // iterate through all lines
          for
          (
            util::ShPtrList< Line>::const_iterator line_itr( remark_lines.Begin()), line_itr_end( remark_lines.End());
            line_itr != line_itr_end;
            ++line_itr
          )
          {
            OSTREAM << ( *line_itr)->GetString() << '\n';
          }
        }
        else
        {
          // iterate through all lines
          for
          (
            util::ShPtrList< Line>::const_iterator line_itr( itr->second.Begin()), line_itr_end( itr->second.End());
            line_itr != line_itr_end;
            ++line_itr
          )
          {
            OSTREAM << ( *line_itr)->GetString() << '\n';
          }
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Head::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Lines, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Head::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Lines, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief pdb line to sse definition form a HELIX or SHEET line, with sstype, start and end residue
    //! @return Triplet of SSType starting and end residue
    storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
    Head::SSEDefinitionFromPDBLine( const Line &PDB_LINE)
    {
      // helix definition
      if( PDB_LINE.GetType() == GetLineTypes().HELIX)
      {
        return
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          (
            biol::GetSSTypes().SSTypeFromPDBHelixClass( PDB_LINE.GetNumericalValue< size_t>( GetEntryTypes().HELIXClass)),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().HELIXResidueName_Initial),
              PDB_LINE.GetChar( GetEntryTypes().HELIXChainID_Initial),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().HELIXSequenceID_Initial),
              PDB_LINE.GetChar( GetEntryTypes().HELIXInsertionCode_Initial)
            ),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().HELIXResidueName_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().HELIXChainID_Terminal),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().HELIXSequenceID_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().HELIXInsertionCode_Terminal)
            )
          );
      }
      // sheet definition
      else if( PDB_LINE.GetType() == GetLineTypes().SHEET)
      {
        return
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          (
            biol::GetSSTypes().STRAND,
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().SHEETResidueName_Initial),
              PDB_LINE.GetChar( GetEntryTypes().SHEETChainID_Initial),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().SHEETSequenceID_Initial),
              PDB_LINE.GetChar( GetEntryTypes().SHEETInsertionCode_Initial)
            ),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().SHEETResidueName_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().SHEETChainID_Terminal),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().SHEETSequenceID_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().SHEETInsertionCode_Terminal)
            )
          );
      }

      // cannot handle other sse definitions yet
      BCL_MessageCrt( "cannot generate sse definition for that pdbline\n" + PDB_LINE.GetString());

      // end
      return storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>();
    }

    //! @brief get all sse definitions as they are given in the HELIX/SHEET section
    //! @param HELIX_CLASSES set of helix types to consider
    //! @param MERGE_OVERLAPPING merge overlapping sses of the same type
    //! @return Map of chain ids and
    storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
    Head::SSEDefinitions( const storage::Set< biol::SSType> &HELIX_CLASSES, const bool MERGE_OVERLAPPING) const
    {
      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
        set_of_sse_definitions;

      // extract HELIX and SHEET lines form pdb
      util::ShPtrList< Line> helix_strand_lines;
      helix_strand_lines.Append( GetLines( GetLineTypes().HELIX));
      helix_strand_lines.Append( GetLines( GetLineTypes().SHEET));

      // insert helix and sheet information from pdb into ordered set
      for
      (
        util::ShPtrList< Line>::const_iterator itr( helix_strand_lines.Begin()), itr_end( helix_strand_lines.End());
        itr != itr_end; ++itr
      )
      {
        // create the sse definition
        const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          sse_definition( SSEDefinitionFromPDBLine( **itr));

        const char chain_id( sse_definition.Second().GetChainID());
        // check if starting and terminating residue are from the same chain
        if( chain_id != sse_definition.Third().GetChainID())
        {
          BCL_MessageCrt( "start and end residue for that line are from different chains\n" + ( *itr)->GetString());
          continue; // next sse definition
        };

        // check that residues are given in correct order
        if( sse_definition.Third() < sse_definition.Second())
        {
          BCL_MessageCrt
          (
            "SSE definition in this line might be given in wrong order, so skipping!\n" + ( *itr)->GetString()
          );
          continue;
        }

        // check that helix types are considered
        if( ( *itr)->GetType() == GetLineTypes().HELIX && HELIX_CLASSES.Find( sse_definition.First()) == HELIX_CLASSES.End())
        {
          BCL_MessageStd( "HELIX line is ignored, since helix class is not considered:\n" + ( *itr)->GetString());
          continue;
        }

        set_of_sse_definitions[ chain_id].PushBack( sse_definition);
      }

      if( MERGE_OVERLAPPING)
      {
        return SSEDefinitionsMergeOverlap( set_of_sse_definitions);
      }

      return set_of_sse_definitions;
    }

    //! @brief check and merge sse definitions for overlap
    //! @param SSE_DEFINITIONS
    //! @return Map of chain ids and sse definitions that were merged when possible
    storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
    Head::SSEDefinitionsMergeOverlap( const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > &SSE_DEFINITIONS)
    {
      typedef std::pointer_to_binary_function< const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &, const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &, bool> SSEInfoCompareType;

      // initiate set of sse information
      std::map< char, std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> >
        set_of_sse_definitions;

      // insert helix and sheet information from pdb into ordered set
      for
      (
        storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >::const_iterator chain_itr( SSE_DEFINITIONS.Begin()), chain_itr_end( SSE_DEFINITIONS.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> &this_set
        (
          set_of_sse_definitions.insert
          (
            std::make_pair
            (
              chain_itr->first,
              std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType>( std::ptr_fun( &SSEInformationCompare))
            )
          ).first->second
        );

        for
        (
          storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          // create the sse definition
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> sse_definition( *itr);

          // check that residues are given in correct order
          if( sse_definition.Third() < sse_definition.Second())
          {
            BCL_MessageCrt
            (
              "SSE definition might be given in wrong order!\n" + util::Format()( *itr)
            );
          }

          std::pair< std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType>::iterator, bool>
            insert_itr_success;

          // as long as it is impossible to store the sse definition, there might be an overlapping definition
          // one overlapping residue is permitted by the comparison function, more than one is not
          while( !( insert_itr_success = this_set.insert( sse_definition)).second)
          {
            BCL_MessageCrt
            (
              "similar sse definition was already found in pdb. Tried to insert:\n" + util::Format()( sse_definition) +
              "\n which is similar to already existing:\n" + util::Format()( *insert_itr_success.first) +
              "\n they will be merged into one"
            );

            // check that overlapping definitions are of the same sse type
            if( insert_itr_success.first->First() != sse_definition.First())
            {
              BCL_MessageCrt( "overlapping sse definitions have different sstypes");
              break;
            };

            // if sse definition comes before the one already in the set
            if( SSEInformationCompare( sse_definition, *insert_itr_success.first))
            {
              // the sse definition becomes the combination of both
              sse_definition =
                storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
                (
                  sse_definition.First(),
                  sse_definition.Second(),
                  insert_itr_success.first->Third()
                );
            }
            else
            {
              sse_definition =
                storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
                (
                  sse_definition.First(),
                  insert_itr_success.first->Second(),
                  sse_definition.Third()
                );
            }

            // remove the overlapping definition already in the set
            this_set.erase( insert_itr_success.first);
          }
        }
      }

      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > result;

      // convert temp storage to final result
      for
      (
        std::map< char, std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> >::const_iterator
          itr( set_of_sse_definitions.begin()), itr_end( set_of_sse_definitions.end());
        itr != itr_end;
        ++itr
      )
      {
        result[ itr->first] = storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >( itr->second.begin(), itr->second.end());
      }

      return result;
    }

    //! @brief compare the sse infos for overlap
    //! @param SSE_INFO_LHS left hand side sse info
    //! @param SSE_INFO_RHS right hand side sse info
    bool Head::SSEInformationCompare
    (
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_LHS,
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_RHS
    )
    {
      // in the case that two sses share one and the same residue
      if( SSE_INFO_LHS.Third() == SSE_INFO_RHS.Second())
      {
        // compare the initial residue of SSE_INFO_LHS with Terminating residue if SSE_INFO_RHS
        return SSE_INFO_LHS.Second() < SSE_INFO_RHS.Third();
      }

      // return if LHS comes before RHS
      return SSE_INFO_LHS.Third() < SSE_INFO_RHS.Second();
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_ligand.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Ligand::s_Instance
    (
      GetObjectInstances().AddInstance( new Ligand())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Ligand::Ligand()
    {
    }

    //! @brief constructor from lines
    Ligand::Ligand( const util::ShPtrList< Line> &LINES) :
      m_Lines( LINES)
    {
      if( LINES.IsEmpty())
      {
        return;
      }

      const Line &first_line( *m_Lines.FirstElement());
      BCL_Assert( first_line.GetType() == GetLineTypes().HETATM, "ligands can only be constructed from HETATM lines");

      m_ResidueName = first_line.GetString(               GetEntryTypes().HETATMResidueName);
      m_ChainID     = first_line.GetChar(                 GetEntryTypes().HETATMChainID);
      m_PDBID       = first_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID);
      m_ICode       = first_line.GetChar(                 GetEntryTypes().HETATMInsertionCode);
    }

    //! @brief Clone function
    //! @return pointer to new Ligand
    Ligand *Ligand::Clone() const
    {
      return new Ligand( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Ligand::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return Lines
    util::ShPtrList< Line> const &Ligand::GetLines() const
    {
      return m_Lines;
    }

    //! change Lines
    util::ShPtrList< Line> &Ligand::ChangeLines()
    {
      return m_Lines;
    }

    //! @brief access to the fullname
    //! @return the full ligand name
    const std::string &Ligand::GetFullname() const
    {
      return m_Fullname;
    }

    //! @brief access to the formula
    //! @return the ligand formula
    const std::string &Ligand::GetFormula() const
    {
      return m_Formula;
    }

    //! @brief set full name
    //! @param FULL_NAME as reported for that residue name in the HETNAME line
    void Ligand::SetFullname( const std::string &FULL_NAME)
    {
      m_Fullname = FULL_NAME;
    }

    //! @brief set the formula
    //! @brief FORMULA as reported in the FORMUL lines for the residue name
    void Ligand::SetFormula( const std::string &FORMULA)
    {
      m_Formula = FORMULA;
    }

    //! @brief add connections for an atom serial
    //! @param CONNECTIONS map of atom serials with sets of all atoms that are connected; all connections are present
    //!        twice and inverted, as they come from the pdb
    void Ligand::AddConnections( const storage::Map< size_t, storage::Set< size_t> > &CONNECTIONS)
    {
      // get the atom serial associated with that residues
      const storage::Set< size_t> serials( AtomSerials());

      // iterate over the connections
      for
      (
        storage::Map< size_t, storage::Set< size_t> >::const_iterator
          center_itr( CONNECTIONS.Begin()), center_itr_end( CONNECTIONS.End());
        center_itr != center_itr_end;
        ++center_itr
      )
      {
        if( !serials.Contains( center_itr->first))
        {
          continue;
        }

        // collect all internal connections
        {
          storage::Set< size_t> connected;
          std::set_intersection
          (
            center_itr->second.Begin(), center_itr->second.End(),
            serials.Begin(), serials.End(),
            std::inserter( connected.InternalData(), connected.End())
          );

          if( !connected.IsEmpty())
          {
            m_InternalConnections[ center_itr->first] = connected;
          }
        }

        // collect external connections
        {
          storage::Set< size_t> connected;
          std::set_difference
          (
            center_itr->second.Begin(), center_itr->second.End(),
            serials.Begin(), serials.End(),
            std::inserter( connected.InternalData(), connected.End())
          );

          if( !connected.IsEmpty())
          {
            m_ExternalConnections[ center_itr->first] = connected;
          }
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> Ligand::AtomSerials() const
    {
      storage::Set< size_t> serials;

      // iterate over all atom lines
      for( util::ShPtrList< Line>::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End()); itr != itr_end; ++itr)
      {
        serials.Insert( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETATMSerial));
      }

      // end
      return serials;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Ligand::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ResidueName        , ISTREAM);
      io::Serialize::Read( m_ChainID            , ISTREAM);
      io::Serialize::Read( m_PDBID              , ISTREAM);
      io::Serialize::Read( m_ICode              , ISTREAM);
      io::Serialize::Read( m_Lines              , ISTREAM);
      io::Serialize::Read( m_Fullname           , ISTREAM);
      io::Serialize::Read( m_Formula            , ISTREAM);
      io::Serialize::Read( m_InternalConnections, ISTREAM);
      io::Serialize::Read( m_ExternalConnections, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Ligand::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ResidueName        , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID            , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID              , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode              , OSTREAM) << '\n';
      io::Serialize::Write( m_Lines              , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Fullname           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Formula            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InternalConnections, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExternalConnections, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace pdb
  
} // namespace bcl
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
#include "pdb/bcl_pdb_line.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Line::s_Instance
    (
      GetObjectInstances().AddInstance( new Line())
    );

    //! @brief record format
    const util::Format Line::s_RecordFormat
    (
      util::Format().W( 6).Fill( ' ').L()
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Line::Line() :
      m_LineType( GetLineTypes().e_Undefined),
      m_String( 128, ' ')
    {
    }

    //! construct Line from LineType and it will write the description
    Line::Line( const LineType &LINE_TYPE) :
      m_LineType( LINE_TYPE),
      m_String( s_LineLength, ' ')
    {
      m_String.replace( s_RecordStart, s_RecordLength, s_RecordFormat( LINE_TYPE.GetName()), s_RecordStart, s_RecordLength);
    }

    //! construct Line from std::string (complete line from pdb) and determines linetype
    Line::Line( const std::string &STRING) :
      m_LineType( GetLineTypes().LineTypeFromPDBLine( STRING)),
      m_String( STRING)
    {
      //ensure that line is long enough
      m_String.resize( s_LineLength, ' ');
    }

    //! copy constructor
    Line *Line::Clone() const
    {
      return new Line( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Line::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns coordinates for ATOM or HETATM line
    //! @return Vector3D with x, Y and Z coordinate
    linal::Vector3D Line::RetrieveCoordinates() const
    {
      //function only works for coordinates in ATOM or HETATM line
      if( m_LineType == GetLineTypes().ATOM)
      {
        //return
        return linal::Vector3D
               (
                 GetNumericalValue< double>( GetEntryTypes().ATOMX),
                 GetNumericalValue< double>( GetEntryTypes().ATOMY),
                 GetNumericalValue< double>( GetEntryTypes().ATOMZ)
               );
      }
      else if( m_LineType == GetLineTypes().HETATM)
      {
        //return
        return linal::Vector3D
               (
                 GetNumericalValue< double>( GetEntryTypes().HETATMX),
                 GetNumericalValue< double>( GetEntryTypes().HETATMY),
                 GetNumericalValue< double>( GetEntryTypes().HETATMZ)
               );
      }
      BCL_MessageCrt( "Function Position called for a non-ATOM line!");
      return linal::Vector3D( util::GetUndefined< double>());
    }

    //! @brief sets x, y and z coordinates of an ATOM or HETATM line
    //! @param COORDINATES the coordinates of that atom
    void Line::PutCoordinates( const linal::Vector3D &COORDINATES)
    {
      if( m_LineType == GetLineTypes().ATOM)
      {
        Put( GetEntryTypes().ATOMX, COORDINATES.X());
        Put( GetEntryTypes().ATOMY, COORDINATES.Y());
        Put( GetEntryTypes().ATOMZ, COORDINATES.Z());
      }
      else if( m_LineType == GetLineTypes().HETATM)
      {
        Put( GetEntryTypes().HETATMX, COORDINATES.X());
        Put( GetEntryTypes().HETATMY, COORDINATES.Y());
        Put( GetEntryTypes().HETATMZ, COORDINATES.Z());
      }
      else
      {
        BCL_MessageCrt( "cannot write Coordinates in non ATOM line");
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief equal operator
    Line &Line::operator =( const Line &LINE)
    {
      // copy member
      m_LineType = LINE.m_LineType;
      m_String   = LINE.m_String;

      // end
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns bool whether Line matches certain criterium
    //! @param CRITERIUM pair of entry type and the string that should be present
    //! @return true if the line type matches the entry type and if the entry is equal to the string
    bool Line::MatchesCriteria( const storage::Pair< EntryType, std::string> &CRITERIUM) const
    {
      //return if LineType matches Entries and string comparison between criteriastring and trimmed string within ENTRIES
      return m_LineType                    == CRITERIUM.First()->GetLineType() &&
             GetString( CRITERIUM.First()) == CRITERIUM.Second();
    }

    //! get certain entry as char
    char Line::GetChar( const EntryType &ENTRY) const
    {
      BCL_Assert
      (
        m_LineType == ENTRY->GetLineType(),
        "EntryType: " + ENTRY.GetName() + " does not match LineType " + m_LineType.GetName() +
        " for pdb line:\n" + m_String
      );
      BCL_Assert( ENTRY->GetDataType() == util::CPPDataTypes::e_Char, "cannot get non-char entry");

      return m_String[ ENTRY->GetStart()];
    }

    //! get certain entry as string
    std::string Line::GetString( const EntryType &ENTRY) const
    {
      BCL_Assert
      (
        m_LineType == ENTRY->GetLineType(),
        "EntryType: " + ENTRY.GetName() + " does not match LineType " + m_LineType.GetName() +
        " for pdb line:\n" + m_String
      );

      return m_String.substr( ENTRY->GetStart(), ENTRY->GetLength());
    }

    //! @brief clear the line except for the line type and line type record
    void Line::Clear()
    {
      // replace all positions behind record locator with ' '
      std::fill( m_String.begin() + LineTypes::s_TypeRecordStart + LineTypes::s_TypeRecordLength, m_String.end(), ' ');
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Line into std::ostream
    std::ostream &Line::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LineType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_String, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Line from std::istream
    std::istream &Line::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_LineType, ISTREAM);
      io::Serialize::Read( m_String, ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief copy a HETATM to ATOM line
    //! @param HETATM_LINE HETATM line
    //! @return ATOM line with entries from given HETATM line
    Line Line::CopyHetatmToAtomLine( const Line &HETATM_LINE)
    {
      // check that line is HETATM line
      if( HETATM_LINE.GetType() != GetLineTypes().HETATM)
      {
        BCL_MessageCrt( "passed non HETATM line but: " + HETATM_LINE.GetString());

        return Line();
      }

      // copy line, set type and change record string
      Line atom_line( HETATM_LINE);
      atom_line.m_LineType = GetLineTypes().ATOM;
      atom_line.m_String.replace( s_RecordStart, s_RecordLength, s_RecordFormat( GetLineTypes().ATOM.GetName()), s_RecordStart, s_RecordLength);

      // end
      return atom_line;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_line_criterium.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_line.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LineCriterium::s_Instance
    (
      GetObjectInstances().AddInstance( new LineCriterium())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LineCriterium::LineCriterium() :
      m_MeetAllCriteria( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LineCriterium
    LineCriterium *LineCriterium::Clone() const
    {
      return new LineCriterium( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LineCriterium::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief just add a line type
    //! @param LINE_TYPE the line needs to have a specific line type
    void LineCriterium::AddCriterium( const LineType &LINE_TYPE)
    {
      m_Criteria[ LINE_TYPE];
    }

    //! @brief set meet all critiera
    //! is set to true, all critiera need be met to be true, if false, only one has to be true
    //! @param MEET_ALL_CRITIERA
    void LineCriterium::SetMeetAllCriteria( const bool MEET_ALL_CRITERIA)
    {
      m_MeetAllCriteria = MEET_ALL_CRITERIA;
    }

    //! @brief reset the criteria
    void LineCriterium::Reset()
    {
      m_Criteria.Reset();
    }

    //! @brief filter all lines that match the criterium from a given set of lines
    //! @param CRITERIUM the criterium that each copied line has to fullfill
    //! @return the list of lines that match the given criterium
    util::ShPtrList< Line>
    LineCriterium::Filter( const util::ShPtrList< Line> &LINES, const util::FunctionInterface< Line, bool> &CRITERIUM)
    {
      util::ShPtrList< Line> result;

      // iterate through all lines
      for( util::ShPtrList< Line>::const_iterator itr( LINES.Begin()), itr_end( LINES.End()); itr != itr_end; ++itr)
      {
        if( CRITERIUM( **itr))
        {
          result.PushBack( *itr);
        }
      }

      return result;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that checks of line meets all criteria
    //! @param LINE the line to check
    //! @return true, if the criteria are met (for intersect, if all are met, for union if one is met)
    bool LineCriterium::operator()( const Line &LINE) const
    {
      // check if there are criteria for that line
      const storage::Map< LineType, storage::Vector< storage::Pair< EntryType, std::string> > >::const_iterator
        line_type_itr( m_Criteria.Find( LINE.GetType()));

      if( line_type_itr == m_Criteria.End())
      {
        return false;
      }

      if( line_type_itr->second.IsEmpty())
      {
        return true;
      }

      if( m_MeetAllCriteria)
      {
        // iterate through criteria
        for
        (
          storage::Vector< storage::Pair< EntryType, std::string> >::const_iterator
            crit_itr( line_type_itr->second.Begin()), crit_itr_end( line_type_itr->second.End());
          crit_itr != crit_itr_end;
          ++crit_itr
        )
        {
          if( !LINE.MatchesCriteria( *crit_itr))
          {
            return false;
          }
        }

        return true;
      }

      // meet one criteria
      // iterate through criteria
      for
      (
        storage::Vector< storage::Pair< EntryType, std::string> >::const_iterator
          crit_itr( line_type_itr->second.Begin()), crit_itr_end( line_type_itr->second.End());
        crit_itr != crit_itr_end;
        ++crit_itr
      )
      {
        if( LINE.MatchesCriteria( *crit_itr))
        {
          return true;
        }
      }

      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LineCriterium::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Criteria       , ISTREAM);
      io::Serialize::Read( m_MeetAllCriteria, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LineCriterium::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Criteria       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MeetAllCriteria, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_line_type_data.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LineTypeData::LineTypeData() :
      m_MultipleTimes( false),
      m_MultipleLines( false),
      m_First( util::GetUndefined< size_t>()),
      m_Last( util::GetUndefined< size_t>())
    {
    }

    //! @brief construct from record type
    //! @param MULTIPLE_TIMES occurs multiple times
    //! @param MULTIPLE_LINES over multiple lines
    LineTypeData::LineTypeData( const bool MULTIPLE_TIMES, const bool MULTIPLE_LINES) :
      m_MultipleTimes( MULTIPLE_TIMES),
      m_MultipleLines( MULTIPLE_LINES),
      m_First( util::GetUndefined< size_t>()),
      m_Last( util::GetUndefined< size_t>())
    {
    }

    //! @brief Clone function
    //! @return pointer to new LineTypeData
    LineTypeData *LineTypeData::Clone() const
    {
      return new LineTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LineTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief the first entry type for that linetype
    //! @return entry type
    EntryType LineTypeData::GetFirstEntryType() const
    {
      return EntryType( m_First);
    }

    //! @brief the last entry type for that linetype
    //! @return entry type
    EntryType LineTypeData::GetLastEntryType() const
    {
      return EntryType( m_Last);
    }

    //! @brief number fo entry types for tat linetype
    //! @return number of entries
    size_t LineTypeData::GetNumberOfEntryTypes() const
    {
      return m_Last - m_First + 1;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief consider an additional EntryType
    //! @param ENTRY_TYPE
    void LineTypeData::ConsiderNewEntryType( const EntryType &ENTRY_TYPE)
    {
      // update last
      m_Last = ENTRY_TYPE.GetIndex();

      // check if first is still undefined
      if( !util::IsDefined( m_First))
      {
        m_First = m_Last;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LineTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_First, ISTREAM);
      io::Serialize::Read( m_Last, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &LineTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_First, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Last, OSTREAM, 0);

      // return
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_line_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //! construct all line types
    LineTypes::LineTypes() :
      HEADER( AddEnum( "HEADER", LineTypeData( false, false))),
      OBSLTE( AddEnum( "OBSLTE", LineTypeData( false,  true))),
      TITLE ( AddEnum( "TITLE" , LineTypeData( false,  true))),
      SPLIT ( AddEnum( "SPLIT" , LineTypeData( false,  true))),
      CAVEAT( AddEnum( "CAVEAT", LineTypeData( false,  true))),
      COMPND( AddEnum( "COMPND", LineTypeData( false,  true))),
      SOURCE( AddEnum( "SOURCE", LineTypeData( false,  true))),
      KEYWDS( AddEnum( "KEYWDS", LineTypeData( false,  true))),
      EXPDTA( AddEnum( "EXPDTA", LineTypeData( false,  true))),
      NUMMDL( AddEnum( "NUMMDL", LineTypeData( false, false))),
      MDLTYP( AddEnum( "MDLTYP", LineTypeData( false,  true))),
      AUTHOR( AddEnum( "AUTHOR", LineTypeData( false,  true))),
      REVDAT( AddEnum( "REVDAT", LineTypeData(  true, false))),
      SPRSDE( AddEnum( "SPRSDE", LineTypeData( false,  true))),
      JRNL  ( AddEnum( "JRNL"  , LineTypeData(  true,  true))), // special format
      REMARK( AddEnum( "REMARK", LineTypeData(  true,  true))), // special format
      DBREF ( AddEnum( "DBREF" , LineTypeData(  true, false))),
      DBREF1( AddEnum( "DBREF1", LineTypeData(  true, false))),
      DBREF2( AddEnum( "DBREF2", LineTypeData(  true, false))),
      SEQADV( AddEnum( "SEQADV", LineTypeData(  true, false))),
      SEQRES( AddEnum( "SEQRES", LineTypeData(  true,  true))),
      MODRES( AddEnum( "MODRES", LineTypeData(  true, false))),
      HET   ( AddEnum( "HET"   , LineTypeData(  true, false))),
      HETNAM( AddEnum( "HETNAM", LineTypeData(  true,  true))),
      HETSYN( AddEnum( "HETSYN", LineTypeData(  true,  true))),
      FORMUL( AddEnum( "FORMUL", LineTypeData(  true,  true))),
      HELIX ( AddEnum( "HELIX" , LineTypeData(  true, false))),
      SHEET ( AddEnum( "SHEET" , LineTypeData(  true, false))),
      SSBOND( AddEnum( "SSBOND", LineTypeData(  true, false))),
      LINK  ( AddEnum( "LINK"  , LineTypeData(  true, false))),
      CISPEP( AddEnum( "CISPEP", LineTypeData(  true, false))),
      SITE  ( AddEnum( "SITE"  , LineTypeData(  true,  true))),
      CRYST1( AddEnum( "CRYST1", LineTypeData( false, false))),
      CRYST2( AddEnum( "CRYST2", LineTypeData( false, false))),
      CRYST3( AddEnum( "CRYST3", LineTypeData( false, false))),
      ORIGX1( AddEnum( "ORIGX1", LineTypeData( false, false))),
      ORIGX2( AddEnum( "ORIGX2", LineTypeData( false, false))),
      ORIGX3( AddEnum( "ORIGX3", LineTypeData( false, false))),
      SCALE1( AddEnum( "SCALE1", LineTypeData( false, false))),
      SCALE2( AddEnum( "SCALE2", LineTypeData( false, false))),
      SCALE3( AddEnum( "SCALE3", LineTypeData( false, false))),
      MTRIX1( AddEnum( "MTRIX1", LineTypeData(  true, false))),
      MTRIX2( AddEnum( "MTRIX2", LineTypeData(  true, false))),
      MTRIX3( AddEnum( "MTRIX3", LineTypeData(  true, false))),
      MODEL ( AddEnum( "MODEL" , LineTypeData( false,  true))), // for grouping
      ATOM  ( AddEnum( "ATOM"  , LineTypeData(  true, false))),
      ANISOU( AddEnum( "ANISOU", LineTypeData(  true, false))),
      TER   ( AddEnum( "TER"   , LineTypeData( false,  true))), // for grouping
      HETATM( AddEnum( "HETATM", LineTypeData(  true, false))),
      ENDMDL( AddEnum( "ENDMDL", LineTypeData( false,  true))), // for grouping
      CONECT( AddEnum( "CONECT", LineTypeData(  true, false))),
      MASTER( AddEnum( "MASTER", LineTypeData( false, false))),
      END   ( AddEnum( "END"   , LineTypeData( false, false))),
      m_GroupingTypes( MODEL, TER, ENDMDL),
      m_OutOfOrderTypes( storage::Set< LineType>::Create( MODEL, TER, ENDMDL, ATOM, ANISOU, HETATM))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LineTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set of grouping line types
    //! types that are used for grouping
    //! @return set of types that are used for grouping
    const storage::Set< LineType> &LineTypes::GetGroupingTypes() const
    {
      return m_GroupingTypes;
    }

    //! @brief set of line types that can appear out of order
    //! @return set of types that can appear out of order
    const storage::Set< LineType> &LineTypes::GetOutOfOrderTypes() const
    {
      return m_OutOfOrderTypes;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief retrieve LineType from the pdb line string
    //! @param PDB_LINE string of the line in the pdb file
    //! @return LineType for that string
    LineType LineTypes::LineTypeFromPDBLine( const std::string &PDB_LINE) const
    {
      return LineType( util::TrimString( PDB_LINE.substr( s_TypeRecordStart, s_TypeRecordLength)));
    }

    //! @brief consider entry type for given line type
    //! @param LINE_TYPE line type for which the entry type is considered
    //! @param ENTRY_TYPE the entry type
    void LineTypes::ConsiderEntryType( const LineType &LINE_TYPE, const EntryType &ENTRY_TYPE)
    {
      LineType( LINE_TYPE)->ConsiderNewEntryType( ENTRY_TYPE);
    }

    const LineTypes &GetLineTypes()
    {
      return LineTypes::GetEnums();
    }

  } // namespace pdb

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< pdb::LineTypeData, pdb::LineTypes>;

  } // namespace util
} // namespace bcl
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
#include "pdb/bcl_pdb_model.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_residue_simple.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Model::s_Instance
    (
      GetObjectInstances().AddInstance( new Model())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Model::Model() :
      m_ChainLines()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Model
    Model *Model::Clone() const
    {
      return new Model( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Model::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Model::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().ATOM.GetIterator(), GetLineTypes().HETATM.GetIterator() + 1
      );

      return s_line_types;
    }

    //! @brief get hetatm lines
    //! @return reference to HETATM lines
    const storage::Map< char, util::ShPtrList< Line> > &Model::GetHETATMLines() const
    {
      return m_HetatmLines;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Model::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            lines.PushBack( *itr);
          }
        }
      }

      // iterate through all hetatms
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            lines.PushBack( *itr);
          }
        }
      }

      // end
      return lines;
    }

    //! @brief count the number of lines for a given line type
    //! @param LINE_TYPE the desire line type
    //! @return the number of lines for the given line type in the model
    size_t Model::Count( const LineType &LINE_TYPE) const
    {
      size_t count( 0);

      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            ++count;
          }
        }
      }

      // iterate through all hetatm
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            ++count;
          }
        }
      }

      // end
      return count;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Model::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      // collect lines for each chain
      util::ShPtrList< Line> lines;

      // iterate through chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator itr( m_ChainLines.Begin()), itr_end( m_ChainLines.End()); itr != itr_end; ++itr)
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // iterate through hetatms
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator itr( m_HetatmLines.Begin()), itr_end( m_HetatmLines.End()); itr != itr_end; ++itr)
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // end
      return lines;
    }

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @param CHAIN_ID only lines for that chain if know apriori
    //! @return lines that are considered by criterium
    util::ShPtrList< Line>
    Model::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const char CHAIN_ID) const
    {
      util::ShPtrList< Line> collected_lines;

      // from chain
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Find( CHAIN_ID));

        // no such chain
        if( chain_itr != m_ChainLines.End())
        {
          collected_lines = LineCriterium::Filter( chain_itr->second, CRITERIUM);
        }
      }

      // from hetatm
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Find( CHAIN_ID));

        // no such chain
        if( chain_itr != m_HetatmLines.End())
        {
          collected_lines.Append( LineCriterium::Filter( chain_itr->second, CRITERIUM));
        }
      }

      return collected_lines;
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Model::PushBack( const util::ShPtr< Line> &LINE)
    {
      char chain_id( '\0');

      // atom line
      if( LINE->GetType() == GetLineTypes().ATOM)
      {
        chain_id = LINE->GetChar( GetEntryTypes().ATOMChainID);
      }
      // hetatm line
      else if( LINE->GetType() == GetLineTypes().HETATM)
      {
        chain_id = LINE->GetChar( GetEntryTypes().HETATMChainID);
      }
      // anisou line
      else if( LINE->GetType() == GetLineTypes().ANISOU)
      {
        chain_id = LINE->GetChar( GetEntryTypes().ANISOUChainID);
      }
      // ter line
      else if( LINE->GetType() == GetLineTypes().TER)
      {
        chain_id = LINE->GetChar( GetEntryTypes().TERChainID);
        if( chain_id == ' ')
        {
          if( !m_ChainLines.IsEmpty())
          {
            storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.End());
            --chain_itr;
            if( !chain_itr->second.IsEmpty())
            {
              const Line &last_chain_line( *chain_itr->second.LastElement());
              if( last_chain_line.GetType() == GetLineTypes().ATOM)
              {
                chain_id = last_chain_line.GetChar( GetEntryTypes().ATOMChainID);
              }
              else if( last_chain_line.GetType() == GetLineTypes().ATOM)
              {
                chain_id = last_chain_line.GetChar( GetEntryTypes().HETATMChainID);
              }
              else
              {
                BCL_MessageCrt( "TER line does not make sense here!");
                return false;
              }
            }
            else
            {
              BCL_MessageCrt( "TER line does not make sense here!");
              return false;
            }
          }
          else
          {
            BCL_MessageCrt( "TER line does not make sense here!");
            return false;
          }
          BCL_MessageVrb
          (
            "TER line might not be complete: " + LINE->GetString() + "\nassuming it is terminating the highest chainid "
                "read so far: " + chain_id
          );
        }
      }
      // did not find any of the considered line types
      else
      {
        return false;
      }

      // insert the line for the correct chain id
      // check if the TER for the chain in chain lines was given already
      storage::Map< char, util::ShPtrList< Line> >::iterator chain_line_itr( m_ChainLines.Find( chain_id));
      if( chain_line_itr == m_ChainLines.End())
      {
        m_ChainLines[ chain_id].PushBack( LINE);
      }
      else if( chain_line_itr->second.IsEmpty() || chain_line_itr->second.LastElement()->GetType() != GetLineTypes().TER)
      {
        chain_line_itr->second.PushBack( LINE);
      } // TER was already given for that chain id; insert into the hetatm lines
      else
      {
        // only insert hetatm lines
        if( LINE->GetType() != GetLineTypes().HETATM)
        {
          return false;
        }
        m_HetatmLines[ chain_id].PushBack( LINE);
      }

      return true;
    }

    //! @brief reset the line group
    void Model::Reset()
    {
      m_ChainLines.Reset();
      m_HetatmLines.Reset();
      m_StructuredChains.Reset();
    }

    //! @brief create the chains as they were given by the atoms
    //! @return map with chainid as key and list of residues as data
    storage::Map< char, storage::List< ResidueSimple> > Model::GetChains() const
    {
      storage::Map< char, storage::List< ResidueSimple> > chains;

      // iterate over chains
      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
        )
        {
          // current line
          const Line &current_line( **itr);

          // normally, that should be an atom line
          if( current_line.GetType() == GetLineTypes().ATOM)
          {
            // current residue from first atom line
            const ResidueSimple current_residue
            (
              current_line.GetString( GetEntryTypes().ATOMResidueName),
              current_line.GetChar( GetEntryTypes().ATOMChainID),
              current_line.GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID),
              current_line.GetChar( GetEntryTypes().ATOMInsertionCode)
            );
            chains[ current_residue.GetChainID()].PushBack( current_residue);

            // skip all other atom lines belonging to that residue
            while( ++itr != itr_end)
            {
              const Line &next_line( **itr);
              if
              (
                   next_line.GetType() != GetLineTypes().ATOM
                || current_residue.GetResidueName() != next_line.GetString( GetEntryTypes().ATOMResidueName)
                || current_residue.GetChainID()     != next_line.GetChar( GetEntryTypes().ATOMChainID)
                || current_residue.GetPDBID()       != next_line.GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID)
                || current_residue.GetICode()       != next_line.GetChar( GetEntryTypes().ATOMInsertionCode)
              )
              {
                break;
              }
            } // advance iterator to next residue
          } // if ATOM line
          // sometime residues are not standard
          else if( current_line.GetType() == GetLineTypes().HETATM)
          {
            // current residue from first HETATM line
            const ResidueSimple current_residue
            (
              current_line.GetString( GetEntryTypes().HETATMResidueName),
              current_line.GetChar( GetEntryTypes().HETATMChainID),
              current_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID),
              current_line.GetChar( GetEntryTypes().HETATMInsertionCode)
            );
            chains[ current_residue.GetChainID()].PushBack( current_residue);

            // skip all other hetatm lines
            while( ++itr != itr_end)
            {
              const Line &next_line( **itr);
              if
              (
                   next_line.GetType() != GetLineTypes().HETATM
                || current_residue.GetResidueName() != next_line.GetString( GetEntryTypes().HETATMResidueName)
                || current_residue.GetChainID()     != next_line.GetChar( GetEntryTypes().HETATMChainID)
                || current_residue.GetPDBID()       != next_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID)
                || current_residue.GetICode()       != next_line.GetChar( GetEntryTypes().HETATMInsertionCode)
              )
              {
                break;
              }
            } // advance iterator to next residue
          } // if HETATM line
          else // other line type
          {
            ++itr;
          }
        } // iterate over pdb lines
      } // iterate through chains

      // end
      return chains;
    }

    namespace
    {
      //! @brief get the longest contiguous stretch of AAs from a given position
      //! @param ITR, ITR_END iterator positions
      size_t GetNextLongestContiguousStretch
      (
        storage::List< storage::Map< char, Residue> >::const_iterator ITR,
        const storage::List< storage::Map< char, Residue> >::const_iterator &ITR_END
      )
      {
        if( ITR == ITR_END)
        {
          return 0;
        }
        int sz( 1);
        const int initial_pdb_id( ITR->Begin()->second.GetPDBID());
        for( ++ITR; ITR != ITR_END && ITR->Begin()->second.GetPDBID() == initial_pdb_id + sz; ++ITR, ++sz)
        {
        }
        return sz;
      }
    }

    //! @brief initialize the structure
    //! @param CHAIN_ID
    //! @param SEQRES sequences for the chain
    //! @param MISSING_RESIDUES residues that are missing in the ATOM section (REMARK 465)
    void Model::InitializeStructuredChain
    (
      const char CHAIN_ID,
      const storage::List< ResidueSimple> &SEQRES,
      const storage::List< ResidueSimple> &MISSING_RESIDUES
    )
    {
      // initialize the chain from the seqres
      storage::List< Residue> &structure( m_StructuredChains[ CHAIN_ID]);
      for( storage::List< ResidueSimple>::const_iterator res_itr( SEQRES.Begin()), res_itr_end( SEQRES.End()); res_itr != res_itr_end; ++res_itr)
      {
        structure.PushBack( Residue( *res_itr));
      }

      // locate all atom lines for the current chain
      // all atom lines are collected in residues, which have a Map of alternate locations for one and the same residue
      storage::List< storage::Map< char, Residue> >
        atom_line_residues( AtomLineResiduesFromChainID( CHAIN_ID));

      // merge the atom line residues with the missing residues
      const storage::List< storage::Map< char, Residue> >
        merged_atom_line_residues( MergeAtomLinesAndMissingResidueLine( atom_line_residues, MISSING_RESIDUES));

      // iterator on list of residues from atom lines
      storage::List< storage::Map< char, Residue> >::const_iterator
        atom_itr( merged_atom_line_residues.Begin()), atom_itr_end( merged_atom_line_residues.End());

      // iterator for residues in that chain from seqres
      storage::List< Residue>::iterator res_itr( structure.Begin()), res_itr_end( structure.End());

      // align atom line and residue gaps
      while( atom_itr != atom_itr_end && res_itr != res_itr_end)
      {
        // make copies of "res_itr" and "atom_itr"
        // they point to residues and atom lines, respectively, with matching residue types
        storage::List< storage::Map< char, Residue> >::const_iterator atom_itr2( atom_itr);
        storage::List< Residue>::iterator res_itr2( res_itr);

        // create size_t "number_aligned_residues" this will hold the number of aligned residues
        size_t number_aligned_residues( 0);

        // determine the maximum # of residues that could possibly be aligned by AlignResiduesAndAtomLines
        size_t dist_res_end( std::distance( res_itr2, res_itr_end));
        size_t dist_atom_end( GetNextLongestContiguousStretch( atom_itr2, atom_itr_end));
        size_t max_residues_aligned
        (
          std::min
          (
            size_t( s_NumberRequiredAlignedResidues),
            std::min( dist_res_end, dist_atom_end)
          )
        );
        storage::List< Residue>::iterator last_aligned_residue( res_itr);

        // match up first residue from Structure with first by res name corresponding atom line
        // and check if at least max_residues_aligned can be aligned if name matches
        while
        (
             !CompareAtomLineResWithResByName( *atom_itr, *res_itr, false)
          || ( number_aligned_residues = AlignResiduesAndAtomLines( atom_itr2, res_itr2, atom_itr_end, res_itr_end, false, max_residues_aligned)) < max_residues_aligned
        )
        {
          ++res_itr;
          if( res_itr == res_itr_end)
          {
            // outer loop will also break, since res_itr points to end
            break;
          }
          else
          {
            BCL_MessageDbg
            (
              "could align " + util::Format()( number_aligned_residues) + " starting from residue\n" +
              util::Format()( *res_itr)
            );
          }
          atom_itr2 = atom_itr;
          res_itr2 = res_itr;
          dist_res_end = std::distance( res_itr2, res_itr_end);
          dist_atom_end = GetNextLongestContiguousStretch( atom_itr2, atom_itr_end);
          max_residues_aligned =
            std::min
            (
              size_t( s_NumberRequiredAlignedResidues),
              std::min( dist_res_end, dist_atom_end)
            );
        }

        // align the residues that could be aligned so far by assigning atom line to residues
        if( !AlignResiduesAndAtomLines( atom_itr, res_itr, atom_itr_end, res_itr_end, true, util::GetUndefined< size_t>()))
        {
          BCL_MessageStd
          (
            "could not find any residue in m_Structure that matches with given atom line residue:\n" +
            util::Format()( *atom_itr)
          );
          ++atom_itr;
          res_itr = last_aligned_residue;
        }
      }

      // add the pdb id to residues in the beginning and end, if they do not have some yet
      res_itr = structure.Begin();
      size_t first_res_with_defined_pdb_id( 0);
      // find the first res_itr with a defined pdb id
      while( res_itr != res_itr_end && !util::IsDefined( res_itr->GetPDBID()))
      {
        ++res_itr;
        ++first_res_with_defined_pdb_id;
      }

      // all residues had undefined pdbid
      if( res_itr == res_itr_end)
      {
        BCL_MessageCrt
        (
          "there are no alignable atom lines for that chain in the pdb: " + util::Format()( CHAIN_ID)
        );
        return;
      }

      // residues in front had undefined pdb ids
      if( res_itr != structure.Begin())
      {
        // get the current pdbid
        int current_res_id( res_itr->GetPDBID());
        --res_itr;
        --current_res_id;

        //while( res_itr != res_itr_begin)
        for
        (
          size_t current( 0); current != first_res_with_defined_pdb_id; ++current
        )
        {
          // for pdbs there is no pdb residue id of 0
          if( current_res_id == 0)
          {
            --current_res_id;
          }
          res_itr->SetPDBID( current_res_id);

          // this has to be checked since VS does not allow --Begin on forward iterators
          // and this can't be asserted in the for loop sentinel, since we still have to fix the first residue
          if( res_itr == structure.Begin())
          {
            break;
          }

          BCL_MessageDbg( "did set the resid of that residue\n" + util::Format()( *res_itr));
          --res_itr;
          --current_res_id;
        }
      }

      // iterator to rev begin of chain
      storage::List< Residue>::reverse_iterator
        res_itr_rev( structure.ReverseBegin()),
        res_itr_rev_end( structure.ReverseEnd());

      size_t number_skipped_to_first_defined_pdb_id( 0);
      // find the first res_itr_rev with a defined pdb id
      while( res_itr_rev != res_itr_rev_end && !util::IsDefined( res_itr_rev->GetPDBID()))
      {
        ++res_itr_rev;
        ++number_skipped_to_first_defined_pdb_id;
      }

      // all residues had undefined pdbid
      if( res_itr_rev == res_itr_rev_end)
      {
        BCL_MessageCrt
        (
          "there are no alignable atom lines for that chain in the pdb: " + util::Format()( CHAIN_ID)
        );
        return;
      }

      // residues in tail had undefined pdb ids
      if( res_itr_rev != structure.ReverseBegin())
      {
        // get the current pdbid
        int current_res_id( res_itr_rev->GetPDBID());
        --res_itr_rev;
        ++current_res_id;

        //while( res_itr_rev != res_itr_rev_begin)
        for( size_t current( 0); current < number_skipped_to_first_defined_pdb_id; ++current)
        {
          // for pdbs there is no pdb residue id of 0
          if( current_res_id == 0)
          {
            ++current_res_id;
          }

          res_itr_rev->SetPDBID( current_res_id);

          BCL_MessageDbg( "did set the resid of that residue\n" + util::Format()( *res_itr_rev));

          // this has to be checked since VS does not allow --Begin on forward iterators
          // and this can't be asserted in the for loop sentinel, since we still have to fix the first residue
          if( res_itr_rev == structure.ReverseBegin())
          {
            break;
          }

          --res_itr_rev;
          ++current_res_id;
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Model::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          OSTREAM << ( *itr)->GetString() << '\n';
        }
      }

      // iterate through all hetatm lines
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          OSTREAM << ( *itr)->GetString() << '\n';
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Model::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ChainLines , ISTREAM);
      io::Serialize::Read( m_HetatmLines, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Model::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainLines , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HetatmLines, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief list that stores for each pdbid and insertion code a map of alternate location ids their residue
    //! HETATM lines are also considered, but are only stored as hard copies within the residues as ATOM lines
    //! @param CHAINID chain for which this list is generated
    //! @return list that stores for each pdbid and insertion code a map of alternate location ids their residue
    storage::List< storage::Map< char, Residue> >
    Model::AtomLineResiduesFromChainID( const char CHAINID) const
    {
      util::ShPtrList< Line> chain_atom_lines;
      // from chain
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Find( CHAINID));

        // no such chain
        if( chain_itr != m_ChainLines.End())
        {
          // search for all atom and het atom lines for that chain
          LineCriterium atomline_chain_criteria;
          atomline_chain_criteria.SetMeetAllCriteria( false);
          atomline_chain_criteria.AddCriterium( GetLineTypes().ATOM);
          atomline_chain_criteria.AddCriterium( GetLineTypes().HETATM);

          chain_atom_lines = LineCriterium::Filter( chain_itr->second, atomline_chain_criteria);
        }
      }

      // list that stores for each pdbid and insertion code a map of alternate location ids their residue
      storage::List< storage::Map< char, Residue> > atom_line_residues;

      // iterator for atom lines
      util::ShPtrList< Line>::const_iterator
        atom_itr( chain_atom_lines.Begin()), atom_itr_end( chain_atom_lines.End());

      // list that stores for each pdbid and insertion code a map of alternate location ids their residue
      storage::List< storage::Map< char, Residue> >::iterator alr_res_itr( atom_line_residues.End());

      // create string "current_res_id" initialize with sequence id of atom line currently denoted by "ATOM_LINE_ITR"
      int previous_res_id( util::GetUndefined< int>());

      // create string "current_insert_id" initialize with insertion code of atom line denoted by "ATOM_LINE_ITR"
      char previous_insert_icode( util::GetUndefined< char>());

      while( atom_itr != atom_itr_end)
      {
        // shptr to current atom line
        util::ShPtr< Line> current_line( *atom_itr);
        if( current_line->GetType() == GetLineTypes().HETATM)
        {
          current_line = util::ShPtr< Line>( Line::CopyHetatmToAtomLine( *current_line).Clone());
        }

        // create string "current_res_id" initialize with sequence id of atom line currently denoted by "ATOM_LINE_ITR"
        const int current_res_id( current_line->GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID));

        // create string "current_insert_id" initialize with insertion code of atom line denoted by "ATOM_LINE_ITR"
        const char current_insert_icode( current_line->GetChar( GetEntryTypes().ATOMInsertionCode));

        // hit new residue in atom lines
        if( current_res_id != previous_res_id || previous_insert_icode != current_insert_icode)
        {
          // check that first insertion id has all side chain atom lines
          if( alr_res_itr != atom_line_residues.End())
          {
            CompleteFirstAlternateResidue( *alr_res_itr);
          }

          // insert new empty residue map
          atom_line_residues.PushBack( storage::Map< char, Residue>());

          // set iterator to the new residue map
          alr_res_itr = atom_line_residues.Last();
        }

        // alternate location id
        const char current_alternate_location_id( current_line->GetChar( GetEntryTypes().ATOMAlternateLocationID));

        // point to current residue for that alternate_location_id
        storage::Map< char, Residue>::iterator current_residue_itr( alr_res_itr->Find( current_alternate_location_id));

        // if there is no residue with this current_alternate_location_id
        if( current_residue_itr == alr_res_itr->End())
        {
          // make a new one
          alr_res_itr->operator[]( current_alternate_location_id) =
              Residue
              (
                current_line->GetString( GetEntryTypes().ATOMResidueName),
                CHAINID,
                current_res_id,
                current_insert_icode
              );
          current_residue_itr = alr_res_itr->Find( current_alternate_location_id);
        }

        // assert that residue name is identical for that atom line and that residue
        BCL_Assert
        (
          current_residue_itr->second.GetResidueName() == current_line->GetString( GetEntryTypes().ATOMResidueName),
          "two atom line with same alternate location id, name pdbid and insertion code:\n" +
          ( *atom_itr)->GetString() + "\n" + util::Format()( current_residue_itr->second)
        );

        // insert current line to the residue
        current_residue_itr->second.ChangeLines().PushBack( current_line);

        // remember pdb id and i code
        previous_res_id = current_res_id;
        previous_insert_icode = current_insert_icode;

        // go to next line
        ++atom_itr;
      }

      // end
      return atom_line_residues;
    }

    //! @brief complete first residue in a map of alternate residues
    //! @param RESIDUES map of residues where first one is to be completed
    void
    Model::CompleteFirstAlternateResidue( storage::Map< char, Residue> &RESIDUES)
    {
      // nothing to do if there is only one alternate residue
      if( RESIDUES.GetSize() <= 1)
      {
        return;
      }

      // append atoms lines from the second alternative location
      if( RESIDUES.Begin()->first == ' ')
      {
        Residue &correct_residues( RESIDUES.Begin()->second);
        correct_residues.ChangeLines().Append( ( ++RESIDUES.Begin())->second.GetLines());
      }
    }

    //! @brief add missing residues to the atom line residues
    //! @param ATOM_LINE_RESIDUES residues acquired form residues
    //! @param MISSING_RESIDUES residues that are not located in experiment REMARK 465
    //! @return list of maps (alternate location residues) merged with the missing residues
    storage::List< storage::Map< char, Residue> >
    Model::MergeAtomLinesAndMissingResidueLine
    (
      const storage::List< storage::Map< char, Residue> > &ATOM_LINE_RESIDUES,
      const storage::List< ResidueSimple> &MISSING_RESIDUES
    )
    {
      // assuming that the missing residues and the atom line residues have consistent numbering, the missing residues
      // are just added in between
      storage::List< storage::Map< char, Residue> > merged_residues;

      // iterators for the atom line residues
      storage::List< storage::Map< char, Residue> >::const_iterator
        atom_itr( ATOM_LINE_RESIDUES.Begin()), atom_itr_end( ATOM_LINE_RESIDUES.End());

      // iterator on the missing residues
      storage::List< ResidueSimple>::const_iterator
        mis_itr( MISSING_RESIDUES.Begin()), mis_itr_end( MISSING_RESIDUES.End());

      // try to inserted sorted
      while( atom_itr != atom_itr_end && mis_itr != mis_itr_end)
      {
        bool nothing_happened( true);
        // insert atom line residues till the first missing residue
        for( ; atom_itr != atom_itr_end && atom_itr->Begin()->second < *mis_itr; ++atom_itr, nothing_happened &= false)
        {
          merged_residues.PushBack( *atom_itr);
        }

        // end of atom line residues
        if( atom_itr == atom_itr_end)
        {
          break;
        }

        // insert missing residues till the first atom line
        for( ; mis_itr != mis_itr_end && *mis_itr < atom_itr->Begin()->second; ++mis_itr, nothing_happened &= false)
        {
          merged_residues.InsertElement
          (
            storage::Map< char, Residue>::Create
            (
              std::pair< char, Residue>( ' ', Residue( *mis_itr))
            )
          );
        }

        if( nothing_happened)
        {
          BCL_MessageCrt( "unable to align missing residues to ATOM section residues! Ignoring Missing residues!");
          return ATOM_LINE_RESIDUES;
        }
      }

      // fill up to the end with atom line residues
      for( ; atom_itr != atom_itr_end; ++atom_itr)
      {
        merged_residues.InsertElement( *atom_itr);
      }

      // fill up to the end with missing residues
      for( ; mis_itr != mis_itr_end; ++mis_itr)
      {
        merged_residues.InsertElement
        (
          storage::Map< char, Residue>::Create
          (
            std::pair< char, Residue>( ' ', Residue( *mis_itr))
          )
        );
      }

      // check that merged list is as long as sum of both lists
      BCL_Assert
      (
        merged_residues.GetSize() == ATOM_LINE_RESIDUES.GetSize() + MISSING_RESIDUES.GetSize(),
        "was not able to merge missing residues with residues from atom lines"
      )

      // end
      return merged_residues;
    }

    //! @brief compare a map of alternate residues to a given residue if one matches by name
    //! @param ATOM_RES map of alternate residues
    //! @param RES residue to be matched
    //! @param SET_RESIDUE if match was found, set the given RES to the residue
    //! @return true, if match was found
    bool Model::CompareAtomLineResWithResByName
    (
      const storage::Map< char, Residue> &ATOM_RES,
      Residue &RESIDUE,
      const bool SET_RESIDUE
    )
    {
      // check for match to every alternate residue
      for
      (
        storage::Map< char, Residue>::const_iterator
          atom_res_itr( ATOM_RES.Begin()), atom_res_itr_end( ATOM_RES.End());
        atom_res_itr != atom_res_itr_end;
        ++atom_res_itr
      )
      {
        // there was one residue of the ones with the alternate location id, that did match
        if( biol::GetAATypes().HaveSameParent( atom_res_itr->second.GetResidueName(), RESIDUE.GetResidueName()))
        {
          if( SET_RESIDUE)
          {
            RESIDUE = atom_res_itr->second;
          }
          return true;
        }
      }

      // no match found
      return false;
    }

    //! @brief align range of residues and atom lines
    //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
    //! @param RES_ITR itr to residue that matches up by name with the atom line
    //! @param ATOM_LINE_ITR_END total end of atom lines
    //! @param RES_ITR_END total end of residues
    //! @param SET_RESIDUE_ATOM_LINES if true, atom lines are inserted into the residues
    //! @param MAX_NUMBER_RES_TO_ALIGN stop aligning after this number is exceeded
    //! @return number fo residues that have been aligned
    size_t Model::AlignResiduesAndAtomLines
    (
      storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
      storage::List< Residue>::iterator    &RES_ITR,
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
      const storage::List< Residue>::iterator    &RES_ITR_END,
      const bool SET_RESIDUE_ATOM_LINES,
      const size_t MAX_NUMBER_RES_TO_ALIGN
    )
    {
      // create size_t "number_residue_matches" to hold the number of sequential residues that can be aligned
      size_t number_residue_matches( 0);

      // require that all following residues with the same type as the current type be aligned
      size_t max_number_res_to_align( MAX_NUMBER_RES_TO_ALIGN);
      if( util::IsDefined( MAX_NUMBER_RES_TO_ALIGN))
      {
        max_number_res_to_align
          = std::min
            (
              std::max
              (
                MAX_NUMBER_RES_TO_ALIGN,
                GetLengthHomopolymerChain( ATOM_LINE_ITR, ATOM_LINE_ITR_END, RES_ITR, RES_ITR_END)
              ),
              size_t( std::min( std::distance( RES_ITR, RES_ITR_END), std::distance( ATOM_LINE_ITR, ATOM_LINE_ITR_END)))
            );
      }

      // see how many sequential residues match up
      while
      (
        // try to align only MAX_NUMBER_RES_TO_ALIGN
          number_residue_matches < max_number_res_to_align
        // check that "ATOM_LINE_ITR" has not reached "ATOM_LINE_ITR_END"
        && ATOM_LINE_ITR != ATOM_LINE_ITR_END
        // check that "RES_ITR" has not reached "RES_ITR_END"
        && RES_ITR != RES_ITR_END
        // check that residues match up and evtl set the res to the matched one (SET_RESIDUE_ATOM_LINES flag)
        && CompareAtomLineResWithResByName( *ATOM_LINE_ITR, *RES_ITR, SET_RESIDUE_ATOM_LINES)
      )
      {
        const Residue &current_res( ATOM_LINE_ITR->Begin()->second);
        // Retrieve current res id for sequence distance check
        const int current_res_id( current_res.GetPDBID());

        // got to next atom line residue
        ++ATOM_LINE_ITR;

        // no residues in atom lines left, nothing to match up anymore
        if( ATOM_LINE_ITR == ATOM_LINE_ITR_END)
        {
          // increase the number of residue matches
          ++number_residue_matches;

          // move "RES_ITR" to next residue
          ++RES_ITR;

          break;
        }

        const Residue &next_res( ATOM_LINE_ITR->Begin()->second);

        // different between current and next atom line res id
        const int res_id_diff( next_res.GetPDBID() - current_res_id);

        // if the res id difference is unexpected, there could be a deletion, not mentioned in the REMARK 465
        // since the atom lines also contain the residues with deletions, those residues would not be empty, if they actually
        // come from the ATOM lines
        if( res_id_diff > 1 && !current_res.GetLines().IsEmpty() && !next_res.GetLines().IsEmpty())
        {
          storage::List< Residue>::iterator peek_res_itr( RES_ITR);
          BCL_MessageDbg
          (
            "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
            "two following atom line residues do have difference in res id larger than 1: " +
            util::Format()( res_id_diff) + ". Will try to skip that many residues! till atom line residue line:\n" +
            util::Format()( *ATOM_LINE_ITR)
          );

          // try to advance iterator by the number of residues that are defined by the given difference
          // peek to the next residue and confirm that it would agree at least by its residue name
          storage::AdvanceIterator( peek_res_itr, RES_ITR_END, res_id_diff);

          // if peek hits end of residues, assume that difference does not mean anything
          if( peek_res_itr == RES_ITR_END)
          {
            BCL_MessageDbg
            (
              "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
              "could not skip the residues, since the resulting residue does not exist"
            );

            // increase the number of residue matches
            ++number_residue_matches;

            // move "RES_ITR" to next residue
            ++RES_ITR;

            // next matching pair
            continue;
          }

          // if the peek residue is not at the end
          storage::List< storage::Map< char, Residue> >::const_iterator atom_itr2( ATOM_LINE_ITR);
          storage::List< Residue>::iterator res_itr2( peek_res_itr);

          const size_t number_remaining_residues
          (
            std::min
            (
              std::distance( atom_itr2, ATOM_LINE_ITR_END),
              std::distance( res_itr2, RES_ITR_END)
            )
          );
          const size_t min_aligned_peek_residues
          (
            std::min
            (
              number_remaining_residues,
              size_t( s_NumberRequiredAlignedResidues)
            )
          );

          const size_t number_aligned_peek_residues
          (
            AlignResiduesAndAtomLines
            (
              atom_itr2, res_itr2, ATOM_LINE_ITR_END, RES_ITR_END, false, min_aligned_peek_residues
            )
          );
          if( number_aligned_peek_residues < min_aligned_peek_residues)
          {
            BCL_MessageDbg
            (
              "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
              "could not skip the residues, since the resulting pair of residue did not match: " +
              util::Format()( *peek_res_itr) + "\n!=\n" + util::Format()( *ATOM_LINE_ITR)
            );

            // increase the number of residue matches
            ++number_residue_matches;

            // move "RES_ITR" to next residue
            ++RES_ITR;

            // next iteration
            continue;
          }

          // if aligned more than 3 peek residues
          // set the pdbid of the residues that do not have an atom line
          if( SET_RESIDUE_ATOM_LINES)
          {
            // go to the start of the "gapped" residues
            ++RES_ITR;
            for
            (
              int number_skipped_residues( 1);
              RES_ITR != peek_res_itr && number_skipped_residues < res_id_diff;
              ++RES_ITR, ++number_skipped_residues
            )
            {
              RES_ITR->SetPDBID( current_res_id + number_skipped_residues);
            }
          }

          BCL_MessageDbg
          (
            "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
            "successful skipping the lines, since the residues did match: " +
            util::Format()( *peek_res_itr) + "\n!=\n" + util::Format()( *ATOM_LINE_ITR) +
            "\nand " + util::Format()( number_aligned_peek_residues) + " peek residues could be aligned"
          );

          // set RES_ITR to peek res_itr position
          RES_ITR = peek_res_itr;

          // increase the number of matched residues by the difference in residue ids
          number_residue_matches += res_id_diff - 1;

          // next iteration
          continue;
        }
        // increase the number of residue matches
        ++number_residue_matches;

        // move "RES_ITR" to next residue
        ++RES_ITR;
      }
      return number_residue_matches;
    }

    //! @brief Get the number of residues before getting to the next type
    //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
    //! @param ATOM_LINE_ITR_END total end of atom lines
    //! @param RES_ITR itr to residue that matches up by name with the atom line
    //! @param RES_ITR_END total end of residues
    //! @return min(number next residues in atom lines that have the same res type,
    //!             number next residues in seqres lines that have the same res type)
    size_t Model::GetLengthHomopolymerChain
    (
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
      const storage::List< Residue>::iterator &RES_ITR,
      const storage::List< Residue>::iterator &RES_ITR_END
    )
    {
      if( RES_ITR == RES_ITR_END)
      {
        return 0;
      }
      // determine the # of AAs till the next AA type in the atom lines
      const Residue &first_residue( *RES_ITR);
      const char chainid( first_residue.GetChainID());
      size_t number_seqresidues_till_next_type( 0);
      for( storage::List< Residue>::const_iterator res_itr( RES_ITR); res_itr != RES_ITR_END; ++res_itr)
      {
        if( res_itr->GetResidueName() == first_residue.GetResidueName() && res_itr->GetChainID() == chainid)
        {
          ++number_seqresidues_till_next_type;
        }
        else
        {
          break;
        }
      }

      // determine the # of AAs till the next AA type in the seqres lines
      storage::List< storage::Map< char, Residue> >::const_iterator atom_line_itr( ATOM_LINE_ITR);
      const storage::Map< char, Residue> &first_seqres( *atom_line_itr);
      size_t number_residues_till_next_type( 1);
      for
      (
        ++atom_line_itr;
        atom_line_itr != ATOM_LINE_ITR_END;
        ++atom_line_itr, ++number_residues_till_next_type
      )
      {
        bool does_match( false);
        for
        (
          storage::Map< char, Residue>::const_iterator
            first_res_itr( first_seqres.Begin()), first_res_itr_end( first_seqres.End());
          first_res_itr != first_res_itr_end && !does_match;
          ++first_res_itr
        )
        {
          // check for match to every alternate residue
          for
          (
            storage::Map< char, Residue>::const_iterator
              atom_res_itr( atom_line_itr->Begin()), atom_res_itr_end( atom_line_itr->End());
            atom_res_itr != atom_res_itr_end;
            ++atom_res_itr
          )
          {
            // there was one residue of the ones with the alternate location id, that did match
            if
            (
              biol::GetAATypes().HaveSameParent
              (
                atom_res_itr->second.GetResidueName(),
                first_res_itr->second.GetResidueName()
              )
            )
            {
              does_match = true;
              break;
            }
          }
        }
        if( !does_match)
        {
          break;
        }
      }
      return std::max( number_residues_till_next_type, number_seqresidues_till_next_type);
    }

    //! @brief Inserts missing resid informations if they are in HETATM lines instead of ATOM lines
    //! @param CHAIN_ID the chain if of that chain
    //! @param CHAIN the residues for that chain
    void Model::InsertLinesFromHETATMLines
    (
      const char CHAIN_ID,
      storage::List< Residue> &CHAIN
    ) const
    {
      // instantiate an iterator pointing to the first residue of current chain (called res_itr,
      // because the chain structure has been build with the atom line information
      for
      (
        storage::List< Residue>::iterator res_itr( CHAIN.Begin()), res_itr_end( CHAIN.End());
        res_itr != res_itr_end;
        ++res_itr
      )
      {
        //checks for emptiness of current residue - if it is not empty -> go to next residue
        if( !res_itr->GetLines().IsEmpty())
        {
          continue;
        }

        // construct criteria vector to search for current residue to find HETATM lines
        LineCriterium residue_criteria;
        residue_criteria.SetMeetAllCriteria( true);
        residue_criteria.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, res_itr->GetPDBID());
        residue_criteria.AddCriterium( GetEntryTypes().HETATMInsertionCode, res_itr->GetICode());

        //set new type of each line to be ATOM line
        util::ShPtrList< Line> hetatm_reslines( LineCriterium::Filter( m_ChainLines.Find( CHAIN_ID)->second, residue_criteria));
        for
        (
          util::ShPtrList< Line>::iterator
            hetatom_itr( hetatm_reslines.Begin()), hetatom_itr_end( hetatm_reslines.End());
          hetatom_itr != hetatom_itr_end; ++hetatom_itr
        )
        {
          if
          (
            biol::GetAATypes().HaveSameParent
            (
              res_itr->GetResidueName(),
              ( *hetatom_itr)->GetString( GetEntryTypes().HETATMResidueName)
            )
          )
          {
            util::ShPtr< Line> atom_line( Line::CopyHetatmToAtomLine( **hetatom_itr).Clone());
            res_itr->ChangeLines().PushBack( atom_line);
          }
        }
      }
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_biomatrix.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "pdb/bcl_pdb_line.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterBiomatrix::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterBiomatrix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterBiomatrix::PrinterBiomatrix()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterBiomatrix
    PrinterBiomatrix *PrinterBiomatrix::Clone() const
    {
      return new PrinterBiomatrix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterBiomatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterBiomatrix::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines to return
      util::ShPtrList< Line> bio_molecule_lines;

      // check that there is a multiplier
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );
      if( !sp_multiplier.IsDefined())
      {
        // return empty list
        return bio_molecule_lines;
      }

      // create a map of ShPtr transformation matrices and chain ids
      storage::Map< char, storage::Vector< math::TransformationMatrix3D> > transformations;

      // iterate through the chain multipliers
      for
      (
        storage::Set< util::ShPtr< assemble::ChainMultiplier>, assemble::ChainMultiplierLessThan>::const_iterator
          multiplier_itr( sp_multiplier->GetChainMultipliers().Begin()),
          multiplier_itr_end( sp_multiplier->GetChainMultipliers().End());
        multiplier_itr != multiplier_itr_end; ++multiplier_itr
      )
      {
        // add the transformation to the map
        transformations[ ( *multiplier_itr)->GetInitialChainID()].PushBack
        (
          *( *multiplier_itr)->GetTransformationMatrix()
        );
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      bio_molecule_lines.PushBack( empty_remark_line);

      // create biomolecule # line
      util::ShPtr< Line> biomolecule_nr_line( empty_remark_line->Clone());
      biomolecule_nr_line->Put( GetEntryTypes().REMARK_350_BiomoleculeIdentifier, "BIOMOLECULE: ");
      biomolecule_nr_line->Put( GetEntryTypes().REMARK_350_BiomoleculeNumber, "1");
      bio_molecule_lines.PushBack( biomolecule_nr_line);

      // iterate through the map of transformations
      for
      (
        storage::Map< char, storage::Vector< math::TransformationMatrix3D> >::const_iterator
          transformation_itr( transformations.Begin()), transformation_itr_end( transformations.End());
        transformation_itr != transformation_itr_end; ++transformation_itr
      )
      {
        // build the chain information line
        const std::string chain_text( "APPLY THE FOLLOWING TO CHAINS: ");
        util::ShPtr< Line> chain_line( empty_remark_line->Clone());
        chain_line->Put( GetEntryTypes().REMARK_String, chain_text + transformation_itr->first);
        bio_molecule_lines.PushBack( chain_line);

        // initialize transformation counter
        size_t transformation_ctr( 1);

        // write the transformation matrix to lines and append them
        for
        (
          storage::Vector< math::TransformationMatrix3D>::const_iterator
            matrix_itr( transformation_itr->second.Begin()), matrix_itr_end( transformation_itr->second.End());
          matrix_itr != matrix_itr_end; ++matrix_itr, ++transformation_ctr
        )
        {
          bio_molecule_lines.Append( WriteTransformationToLines( *matrix_itr, transformation_ctr));
        }
      }

      // end
      return bio_molecule_lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterBiomatrix::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterBiomatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write transformation matrix to Biomolecule lines
    //! @param TRANSFORMATION transformation matrix to be written
    //! @param COUNTER counter for current transformation
    //! @return formatted lines containing transformation information
     util::ShPtrList< Line> PrinterBiomatrix::WriteTransformationToLines
    (
      const math::TransformationMatrix3D &TRANSFORMATION,
      const size_t COUNTER
    )
    {
       // initialize lines to return
       util::ShPtrList< Line> bio_molecule_lines;

       // create an empty remark line
       util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
       empty_remark_line->Put( GetEntryTypes().REMARK_Number, "350");

       // get an iterator on the entry type
       EntryTypes::const_iterator type_itr( GetEntryTypes().REMARK_350_BIOMT1_Identifier.GetIterator());

       // iterate through the number of lines to add
       for( size_t i( 0); i != 3; ++i)
       {
         // create a new line with the correct identification
         util::ShPtr< Line> matrix_line_a( empty_remark_line->Clone());
         matrix_line_a->Put( *type_itr, "BIOMT" + util::Format()( i + 1));
         ++type_itr;
         matrix_line_a->Put( *type_itr, COUNTER);
         ++type_itr;

         // iterate through the rows of the matrix
         for( size_t j( 0); j != 4; ++j, ++type_itr)
         {
           // add the value to the line
           matrix_line_a->Put( *type_itr, TRANSFORMATION.GetMatrix()( j, i));
         }

         bio_molecule_lines.PushBack( matrix_line_a);
       }

       // end
       return bio_molecule_lines;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_body_assignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "pdb/bcl_pdb_line.h"
#include "restraint/bcl_restraint_assignment.h"
#include "restraint/bcl_restraint_body.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterBodyAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterBodyAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterBodyAssignment::PrinterBodyAssignment() :
      m_Restraints()
    {
    }

    //! @brief construct from restraint information
    //! @param RESTRAINTS holds the restraints (density map)
    PrinterBodyAssignment::PrinterBodyAssignment
    (
      const util::ShPtr
      <
        util::ShPtrVector< restraint::Body>
      > &RESTRAINTS
    ) :
      m_Restraints( RESTRAINTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterBodyAssignment
    PrinterBodyAssignment *PrinterBodyAssignment::Clone() const
    {
      return new PrinterBodyAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterBodyAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterBodyAssignment::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no restraints
      if( !m_Restraints.IsDefined())
      {
        return lines;
      }

      BCL_MessageCrt
      (
        "writing PrinterProteinModelBodyAssignment::WriteBodySSEAssignmentInformation!!! "
      );

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool( PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool));
      BCL_Assert( sp_pool.IsDefined(), "No pool stored for the given model");

      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just a
      //  single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list( GenerateAssignments( PROTEIN_MODEL));

      // create SiPtrVector "protein_model_sses" and initialize with the SSEs of "PROTEIN_MODEL"
      const util::SiPtrVector< const assemble::SSE> protein_model_sses( PROTEIN_MODEL.GetSSEs());

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      lines.PushBack( empty_remark_line);

      // iterate through all of the SSEs of the protein model in order to see which SSE corresponds with which body in
      // the assignment
      // Need to do this because the information about the SSE from the protein model is lost when it is put into the
      // assignment because the assignment just keeps the body information.
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( protein_model_sses.Begin()), sse_itr_end( protein_model_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // create linal::Vector3D "sse_center" and initialize with the origin of the SSE body
        const linal::Vector3D sse_center( ( *sse_itr)->GetCenter());

        // iterate through the vector of assignments (for one density map this for loop will only execute once)
        for
        (
          storage::Vector< restraint::SSEAssignment>::const_iterator
            assignment_itr( density_sse_assignment_list.Begin()), assignment_itr_end( density_sse_assignment_list.End());
          assignment_itr != assignment_itr_end;
          ++assignment_itr
        )
        {
          // iterate through the individual bodies of the assignment denoted by "assignment_itr"
          // iterate through the bodies of the sses of the protein model (stored in Group Collection of the Assignment)
          for
          (
            restraint::GroupCollection< size_t, assemble::SSE>::const_iterator
              assignment_sses_iter( assignment_itr->GetGroupCollection().Begin()),
              assignment_sses_iter_end( assignment_itr->GetGroupCollection().End());
            assignment_sses_iter != assignment_sses_iter_end;
            ++assignment_sses_iter
          )
          {
            // true if the center of the body behind "assignment_sses_iter"
            // (this is NOT the body of the density rod but the center of a body coming from an SSE in the protein model)
            // is equal to "sse_center"
            if( math::EqualWithinTolerance( ( *assignment_sses_iter->second.Begin())->GetCenter(), sse_center))
            {
              // get the index of the corresponding restraint body (i.e. density rod)
              const size_t density_rod_index( assignment_sses_iter->first);

              // compare the orientation of the density rod and the helix in the protein model
              // initialize bool that holds relative orientation of density rod and the helix in the protein model
              bool relative_orientation_rod_helix( true);

              // get vector along helical axis of sse from protein model
              const linal::Vector3D sse_axis( ( *sse_itr)->GetAxis( coord::GetAxes().e_Z));

              // get vector along helical axis of density rod
              const linal::Vector3D density_rod_axis
              (
                assignment_itr->GetRestraint()->operator()( density_rod_index)->GetAxis( coord::GetAxes().e_Z)
              );

              // check whether both vectors are parallel or antiparallel and store result in relative_orientation_rod_helix
              const double scalar_product( linal::ScalarProduct( sse_axis, density_rod_axis));
              BCL_MessageDbg( "the scalar product is: " + util::Format()( scalar_product));
              // check sign of scalar product (positive if parallel, negative if antiparallel)
              // parallel
              if( math::EqualWithinTolerance( scalar_product, 1.0))
              {
                relative_orientation_rod_helix = true;
              }
              // antiparallel
              else if( math::EqualWithinTolerance( scalar_product, -1.0))
              {
                relative_orientation_rod_helix = false;
              }
              // they should always be parallel or antiparallel
              else
              {
                BCL_MessageCrt( "magnitude of scalar product should always be 1, but it is: " + util::Format()( scalar_product));
                relative_orientation_rod_helix = false;
              }

              // find best match from pool for this particular sse
              // (this is necessary as sses can be changed by moves from their original pool state)
              const util::SiPtr< const assemble::SSE> best_match_from_pool( sp_pool->FindBestMatchFromPool( **sse_itr).First());

              std::string assignment_string
              (
                "map: " + util::Format()( assignment_itr - density_sse_assignment_list.Begin()) +
                " rod: " + util::Format()( density_rod_index)
              );

              if( best_match_from_pool.IsDefined())
              {
                // still iterate over pool to get the pool index of best_match_from_pool sse
                // create size_t "pool_sse_index_counter" and initialize with zero
                // this variable will be used to give an index of an SSE in the pool
                size_t pool_sse_index_counter( 0);

                // now iterate through the pool to see which element in the pool the best_match_from_pool sse corresponds
                // to; the element in the pool will be denoted by its index
                for
                (
                  storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::const_iterator
                    pool_itr( sp_pool->Begin()), pool_itr_end( sp_pool->End());
                  pool_itr != pool_itr_end;
                  ++pool_itr, ++pool_sse_index_counter
                )
                {
                  if( *pool_itr == best_match_from_pool)
                  {
                    break;
                  }
                }

                // generate the assignment string
                assignment_string +=
                  " index: " + util::Format()( pool_sse_index_counter) +
                  " begin: " + util::Format()( best_match_from_pool->GetFirstAA()->GetSeqID()) +
                  " end: " + util::Format()( best_match_from_pool->GetLastAA()->GetSeqID()) +
                  " orientation: " + util::Format()( relative_orientation_rod_helix);
              }
              else
              {
                assignment_string += " no assignment found";
              }

              // add to pdb lines
              util::ShPtr< Line> new_line( empty_remark_line->Clone());
              new_line->Put( GetEntryTypes().REMARK_String, assignment_string);
              lines.PushBack( new_line);
            }
          }
        }
      }

      // only print out body sse agreement information if message level is sufficiently low
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        lines.Append( WriteBodySSEAgreementInformation( PROTEIN_MODEL));
      }

      // end
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterBodyAssignment::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterBodyAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write relative rotation and translation information between every SSE in the protein model and the body
    //! it is assigned to
    //! @param PROTEIN_MODEL protein model for which the translation and rotation agreement information is determined
    //! @return pdb lines which are written to
    util::ShPtrList< Line> PrinterBodyAssignment::WriteBodySSEAgreementInformation
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no restraints
      if( !m_Restraints.IsDefined())
      {
        return lines;
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));

      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just a
      //  single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list( GenerateAssignments( PROTEIN_MODEL));

      // iterate over all the assignments that are contained in the density_sse_assignment_list (in case there is more
      // than 1 density map)
      for
      (
        storage::Vector< restraint::SSEAssignment>::const_iterator
          assignment_itr( density_sse_assignment_list.Begin()), assignment_itr_end( density_sse_assignment_list.End());
        assignment_itr != assignment_itr_end; ++assignment_itr
      )
      {
        // iterate through the assignment denoted by assignment iterator
        for
        (
          restraint::GroupCollection< size_t, assemble::SSE>::const_iterator
            group_collection_itr( assignment_itr->GetGroupCollection().Begin()),
            group_collection_itr_end( assignment_itr->GetGroupCollection().End());
          group_collection_itr != group_collection_itr_end; ++group_collection_itr
        )
        {
          // get body of SSE in protein model
          const assemble::SSE &model_body( *group_collection_itr->second.FirstElement());

          // get body of restraint (density rod)
          const assemble::SSEGeometryInterface &native_body
          (
            *assignment_itr->GetRestraint()->operator()( group_collection_itr->first)
          );

          const math::TransformationMatrix3D difference
          (
            math::Inverse( model_body.GetOrientation())( native_body.GetOrientation())
          );

          // calculate the effective rotation angle between these two bodies
          const double rotation_angle( difference.GetRotation().EffectiveRotationAngle());
          const double translation
          (
            ( model_body.GetOrientation().GetTranslation() - native_body.GetOrientation().GetTranslation()).Norm()
          );

          // generate the assignment string
          const std::string assignment_string
          (
            "map: " + util::Format()( assignment_itr - density_sse_assignment_list.Begin()) +
            " rod: " + util::Format()( group_collection_itr->first) +
            " angle: " + util::Format()( math::Angle::Degree( rotation_angle)) +
            " translation: " + util::Format()( translation)
          );

          // add to pdb lines
          util::ShPtr< Line> new_line( empty_remark_line->Clone());
          new_line->Put( GetEntryTypes().REMARK_String, assignment_string);
          lines.PushBack( new_line);
        }
      }

      // end
      return lines;
    }

    //! @brief generate assignments from a protein model
    //! @param PROTEIN_MODEL protein model of interest
    //! @return assignments
    storage::Vector< restraint::SSEAssignment>
    PrinterBodyAssignment::GenerateAssignments( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just
      // a single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list;

      // create SiPtrVector "protein_model_sses" and initialize with the SSEs of "PROTEIN_MODEL"
      const util::SiPtrVector< const assemble::SSE> protein_model_sses( PROTEIN_MODEL.GetSSEs());

      // iterate through "m_restraints" to build assignment objects
      // An assignment contains the information about the density map and which SSEs from the protein model go with which
      // density rods.
      // In the usual case, this loop will only execute a single time (i.e. only one density map is given)
      for
      (
        util::ShPtrVector< restraint::Body>::const_iterator itr( m_Restraints->Begin()), itr_end( m_Restraints->End());
        itr != itr_end;
        ++itr
      )
      {
        density_sse_assignment_list.PushBack( ( *itr)->GenerateAssignment( protein_model_sses));
      }

      return density_sse_assignment_list;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "pdb/bcl_pdb_line.h"
#include "storage/bcl_storage_table.hpp"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterLoopClosure::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterLoopClosure())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterLoopClosure::PrinterLoopClosure() :
      m_ClosureThreshold( util::GetUndefined< double>())
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param CLOSURE_THRESHOLD the threshold for considering the loop closed
    PrinterLoopClosure::PrinterLoopClosure( const double CLOSURE_THRESHOLD) :
      m_ClosureThreshold( CLOSURE_THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterLoopClosure
    PrinterLoopClosure *PrinterLoopClosure::Clone() const
    {
      return new PrinterLoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterLoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterLoopClosure::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no threshold was set
      if( !util::IsDefined( m_ClosureThreshold))
      {
        return lines;
      }

      // get the loop domain locator from model data
      util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_loop_domains
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_LoopDomainLocators)
      );

      // make sure it's defined
      if( sp_loop_domains.IsDefined())
      {
        util::Format format_num, format_short, format_med, format_long;
        format_num.W( 6).FFP( 3);
        format_short.W( 1);
        format_med.W( 4);
        format_long.W( 15);

        // create table which will hold all the information
        storage::Table< std::string> loop_info
        (
          storage::TableHeader
          (
            storage::Vector< std::string>::Create
            (
              format_short( "T"),
              format_long( "domain"),
              format_num( "length"),
              format_num( "omega"),
              format_num( "closed"),
              format_num( "rms"),
              format_med( "rms<" + util::Format()( m_ClosureThreshold))
            )
          )
        );

        // iterate through the loop domains to fill table with their loop closure status
        for
        (
          util::ShPtrList< fold::LocatorLoopDomain>::const_iterator
            loop_itr( sp_loop_domains->Begin()), loop_itr_end( sp_loop_domains->End());
          loop_itr != loop_itr_end;
          ++loop_itr
        )
        {
          const fold::LocatorLoopDomain &current_locator( **loop_itr);

          storage::Vector< std::string> row;
          assemble::LocatorAA left_aa_loc;
          assemble::LocatorAA right_aa_loc;
          if( current_locator.IsNToCSequenceDirection())
          {
            // amino acids in peptide bond
            left_aa_loc = current_locator.GetLoopSegments().LastElement().GetLocatorSSE().EndAALocator();
            right_aa_loc = assemble::LocatorAA( left_aa_loc.GetLocatorChain().GetChainID(), left_aa_loc.GetAAID() + 1, left_aa_loc.GetUsePDBID());

            row.PushBack( format_short( "C"));
          }
          else
          {
            // amino acids in peptide bond
            right_aa_loc = current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().StartAALocator();
            left_aa_loc = assemble::LocatorAA( right_aa_loc.GetLocatorChain().GetChainID(), right_aa_loc.GetAAID() - 1, right_aa_loc.GetUsePDBID());

            row.PushBack( format_short( "N"));
          }

          const std::string flex_ident
          (
            util::Format()( current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().GetSSEID().First()) +
            "-" +
            util::Format()( current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().GetSSEID().Second()) +
            "," +
            util::Format()( current_locator.GetLoopSegments().LastElement().GetLocatorSSE().GetSSEID().First()) +
            "-" +
            util::Format()( current_locator.GetLoopSegments().LastElement().GetLocatorSSE().GetSSEID().Second())
          );
          row.PushBack( format_long( flex_ident));

          // pointer to left and right amino acid
          const util::SiPtr< const biol::AABase> aa_left( left_aa_loc.Locate( PROTEIN_MODEL));
          const util::SiPtr< const biol::AABase> aa_right( right_aa_loc.Locate( PROTEIN_MODEL));

          if( !aa_left.IsDefined() || !aa_right.IsDefined())
          {
            continue;
          }

          // row name using the amino acids
          const std::string row_name
          (
            util::Format()( aa_left->GetChainID()) + " " +
            util::Format()( aa_left->GetSeqID()) + " " + aa_left->GetType()->GetOneLetterCode() + " - " +
            util::Format()( aa_right->GetSeqID()) + " " + aa_right->GetType()->GetOneLetterCode()
          );

          // peptide bond parameters
          const storage::VectorND< 2, double> peptide_bond_length_angle( biol::AABase::PeptideBondLengthAndAngle( *aa_left, *aa_right));
          row.PushBack( format_num( peptide_bond_length_angle.First()));
          row.PushBack( format_num( peptide_bond_length_angle.Second()));
          row.PushBack( format_num( biol::AABase::AreAminoAcidsPeptideBonded( *aa_left, *aa_right, true)));

          // rms of pseudo residue
          const double rms( fold::LocatorLoopDomain::CalculateRMS( PROTEIN_MODEL, current_locator));
          row.PushBack( format_num( rms));
          row.PushBack( format_med( rms < m_ClosureThreshold));

          // add the row to the table of information
          loop_info.InsertRow( row_name, row, true);
        }

        // create an empty remark line
        util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
        empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
        lines.PushBack( empty_remark_line);

        // write to stream
        std::stringstream ss;
        loop_info.WriteFormatted( ss, format_short, "Loops");
        storage::Vector< std::string> string_lines( util::StringLineListFromIStream( ss));

        // iterate through the lines
        for
        (
          storage::Vector< std::string>::const_iterator line_itr( string_lines.Begin()),
            line_itr_end( string_lines.End());
          line_itr != line_itr_end; ++line_itr
        )
        {
          // get a new line
          util::ShPtr< Line> new_line( empty_remark_line->Clone());

          // add the entry
          new_line->Put( GetEntryTypes().REMARK_String, *line_itr);
          lines.PushBack( new_line);
        }
      }
      else
      {
        BCL_MessageVrb( "Loop domain locators are not stored with the given model");
      }

      // convert table to lines
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterLoopClosure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ClosureThreshold, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterLoopClosure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ClosureThreshold, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "pdb/bcl_pdb_line.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterMembrane::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterMembrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterMembrane::PrinterMembrane()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterMembrane
    PrinterMembrane *PrinterMembrane::Clone() const
    {
      return new PrinterMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterMembrane::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if no membrane is found
      if( !sp_membrane.IsDefined())// && sp_membrane->IsDefined())
      {
        return lines;
      }

      // if membrane is undefined
      if( !sp_membrane->IsDefined())
      {
        return lines;
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      lines.PushBack( empty_remark_line);

      // create format object
      util::Format number_format;
      number_format.FFP( 3);

      // write the normal
      util::ShPtr< Line> normal_line( empty_remark_line->Clone());
      normal_line->Put
      (
        GetEntryTypes().REMARK_String,
        "Membrane normal: " + number_format( sp_membrane->GetNormal().X()) + "\t" +
        number_format( sp_membrane->GetNormal().Y()) + "\t" + number_format( sp_membrane->GetNormal().Z())
      );
      lines.PushBack( normal_line);

      // write the center
      util::ShPtr< Line> center_line( empty_remark_line->Clone());
      center_line->Put
      (
        GetEntryTypes().REMARK_String,
        "Membrane center: " + number_format( sp_membrane->GetCenter().X()) + "\t" +
        number_format( sp_membrane->GetCenter().Y()) + "\t" + number_format( sp_membrane->GetCenter().Z())
      );
      lines.PushBack( center_line);

      // write the thickness
      util::ShPtr< Line> thickness_line( empty_remark_line->Clone());
      const linal::Vector< double> thicknesses( sp_membrane->GetThicknesses());
      std::string thickness_string( "Membrane thickness: ");
      for
      (
        linal::Vector< double>::const_iterator itr( thicknesses.Begin()), itr_end( thicknesses.End());
        itr != itr_end; ++itr
      )
      {
        thickness_string += number_format( *itr) + "\t";
      }
      thickness_line->Put( GetEntryTypes().REMARK_String, thickness_string);
      lines.PushBack( thickness_line);

      // end
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterMembrane::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_quality_docking.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_printer_score.h"

namespace bcl
{
  namespace pdb
  {
    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityDocking())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityDocking::PrinterQualityDocking() :
      m_QualityMeasures(),
      m_ChainIDs()
    {
    }

    //! @brief construct from given quality measures and chain IDs
    //! @param QUALITY_MEASURES quality measures to use
    //! @param CHAIN_IDS chain IDs of the ligand
    PrinterQualityDocking::PrinterQualityDocking
    (
      const storage::Set< quality::Measure> &QUALITY_MEASURES,
      const std::string &CHAIN_IDS
    ) :
      m_QualityMeasures( QUALITY_MEASURES),
      m_ChainIDs( CHAIN_IDS)
    {
    }

    //! @brief Clone function
    //! @return pointer to a copy of this PrinterQualityDocking object
    PrinterQualityDocking *PrinterQualityDocking::Clone() const
    {
      return new PrinterQualityDocking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &PrinterQualityDocking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get quality measures
    //! @return quality measures
    const storage::Set< quality::Measure> &PrinterQualityDocking::GetQualityMeasures() const
    {
      return m_QualityMeasures;
    }

    //! @brief set quality measures
    void PrinterQualityDocking::SetQualityMeasures
    (
      const storage::Set< quality::Measure> &QUALITY_MEASURES
    )
    {
      m_QualityMeasures = QUALITY_MEASURES;
    }

    //! @brief get chain IDs of the ligand
    //! @return chain IDs of the ligand
    const std::string &PrinterQualityDocking::GetChainIDs() const
    {
      return m_ChainIDs;
    }

    //! @brief set chain IDs to retrieve for the ligand
    void PrinterQualityDocking::SetChainIDs( const std::string &CHAIN_IDS)
    {
      m_ChainIDs = CHAIN_IDS;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief write lines containing model quality details
    //! @param MODEL the protein model of which to write quality measure lines
    //! @return lines containing model quality details
    util::ShPtrList< Line> PrinterQualityDocking::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // create a ProteinModel object from the docked ligand chains
      assemble::ProteinModel docked_ligand_model( MODEL.GetChains( m_ChainIDs));

      // create a ProteinModel object from the native ligand chains
      util::ShPtr< assemble::ProteinModel> sp_native_model
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_NativeModel)
      );

      // if native model is not defined, return nan
      if( !sp_native_model.IsDefined())
      {
        // get table headers
        storage::Table< double> undefined_table
        (
          math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
        );
        for
        (
          storage::Set< quality::Measure>::const_iterator
            measure_itr( m_QualityMeasures.Begin()), measure_itr_end( m_QualityMeasures.End());
          measure_itr != measure_itr_end;
          ++measure_itr
        )
        {
          undefined_table.InsertRow
          (
            measure_itr->GetName(),
            storage::Vector< double>::Create( 1.0, util::GetUndefinedDouble(), util::GetUndefinedDouble())
          );
        }
        return PrinterScore::WriteTableToLines( undefined_table);
      }

      // get native ligand model
      const assemble::ProteinModel native_ligand_model( sp_native_model->GetChains( m_ChainIDs));

      // add the native ligand model into docked ligand model
      util::ShPtr< assemble::ProteinModelData> sp_data( new assemble::ProteinModelData());
      sp_data->Insert
      (
        assemble::ProteinModelData::e_NativeModel,
        util::ShPtr< assemble::ProteinModel>( native_ligand_model.Clone())
      );
      sp_data->Insert
      (
        assemble::ProteinModelData::e_NativeFilteredModel,
        util::ShPtr< assemble::ProteinModel>( native_ligand_model.Clone())
      );
      docked_ligand_model.SetProteinModelData( sp_data);

      // construct table for qualities
      const assemble::QualityBatch qualities( m_QualityMeasures, biol::GetAtomTypes().CA);
      const storage::Table< double> table( qualities.ConstructTable( docked_ligand_model));

      // write table to PDB lines
      return PrinterScore::WriteTableToLines( table, false);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterQualityDocking::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityDocking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_quality_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "pdb/bcl_pdb_line.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterQualityMembrane::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityMembrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityMembrane::PrinterQualityMembrane() :
      m_Qualities(),
      m_Environments(),
      m_NativeModel()
    {
    }

    //! @brief construct from members
    //! @param QUALITIES quality measures to calculate
    //! @param ENVIRONMENTS environment types
    //! @param NATIVE native model
    PrinterQualityMembrane::PrinterQualityMembrane
    (
      const storage::Set< quality::Measure> &QUALITIES,
      const storage::Set< biol::EnvironmentType> &ENVIRONMENTS,
      const util::ShPtr< assemble::ProteinModel> &NATIVE
    ) :
      m_Qualities( QUALITIES),
      m_Environments( ENVIRONMENTS),
      m_NativeModel( NATIVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterQualityMembrane
    PrinterQualityMembrane *PrinterQualityMembrane::Clone() const
    {
      return new PrinterQualityMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterQualityMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterQualityMembrane::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if no membrae data, return empty list
      if( !sp_membrane.IsDefined() || !m_NativeModel.IsDefined())
      {
        return util::ShPtrList< Line>();
      }

      // initialize undefined CA atom
      const biol::Atom undefined_ca( linal::Vector3D( util::GetUndefined< double>()), biol::GetAtomTypes().CA);

      // hardcopy the native
      util::ShPtr< assemble::ProteinModel> native_copy( m_NativeModel->HardCopy());

      // iterate over the SSEs
      const util::SiPtrVector< const assemble::SSE> model_sses( native_copy->GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( model_sses.Begin()),
          sse_itr_end( model_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // hardcopy the sse
        util::ShPtr< assemble::SSE> sp_sse( ( *sse_itr)->HardCopy());

        // iterate over the SSE
        for
        (
          util::ShPtrVector< biol::AABase>::iterator aa_itr( sp_sse->Begin()), aa_itr_end( sp_sse->End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // get the CA
          biol::Atom ca( ( *aa_itr)->GetCA());

          // if it is not in the right environment
          if( !m_Environments.Contains( sp_membrane->DetermineEnvironmentType( ca.GetCoordinates())))
          {
            // set the CA to undefined
            ( *aa_itr)->SetAtom( undefined_ca);
          }
        }

        // replace the SSE in the model
        native_copy->Replace( sp_sse);
      }

      // construct the table
      const storage::Table< double> table
      (
        assemble::QualityBatch( m_Qualities, biol::GetAtomTypes().CA, "_mem").ConstructTable
        (
          PROTEIN_MODEL,
          *native_copy
        )
      );

      // return pdb lines
      return PrinterScore::WriteTableToLines( table, false);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterQualityMembrane::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_quality_multimer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_line.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterQualityMultimer::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityMultimer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityMultimer::PrinterQualityMultimer() :
      m_Qualities(),
      m_NativeMultimer()
    {
    }

    //! @brief construct from members
    //! @param QUALITIES quality measures to calculate
    //! @param NATIVE native multimer model
    PrinterQualityMultimer::PrinterQualityMultimer
    (
      const storage::Set< quality::Measure> &QUALITIES,
      const util::ShPtr< assemble::ProteinModel> &NATIVE
    ) :
      m_Qualities( QUALITIES),
      m_NativeMultimer( NATIVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterQualityMultimer
    PrinterQualityMultimer *PrinterQualityMultimer::Clone() const
    {
      return new PrinterQualityMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterQualityMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterQualityMultimer::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // cast a pointer to the multiplier data if any
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // if no multiplier data, return empty list
      if( !sp_multiplier.IsDefined() || !m_NativeMultimer.IsDefined())
      {
        return util::ShPtrList< Line>();
      }

      // model to hold the model multimer
      assemble::ProteinModel model_multimer( sp_multiplier->operator()( PROTEIN_MODEL));

      // set model as best multimer based on lowest RMSD
      model_multimer = assemble::PrinterProteinModelMultimer::CalculateBestMultimer
      (
        model_multimer,
        *m_NativeMultimer,
        quality::GetMeasures().e_RMSD,
        sp_multiplier
      );

      // construct the table
      const storage::Table< double> table
      (
        assemble::QualityBatch( m_Qualities, biol::GetAtomTypes().CA, "_mult").ConstructTable
        (
          model_multimer,
          *m_NativeMultimer
        )
      );

      // return pdb lines
      return PrinterScore::WriteTableToLines( table, false);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterQualityMultimer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityMultimer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_printer_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_line.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterScore::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterScore())
    );

    //! @brief get format for table column format
    //! @return format for table column format
    const util::Format &PrinterScore::GetTableColumnFormat()
    {
      // initialize static const format
      static const util::Format s_format( util::Format().W( 13).FFP( 4).R());

      // end
      return s_format;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterScore::PrinterScore() :
      m_ScoringFunction(),
      m_Qualities()
    {
    }

    //! @brief construct with member variables
    //! @param SP_SCORE optional scoring function
    //! @param QUALITY_MEASURES set of quality measures to calculate
    PrinterScore::PrinterScore
    (
      const util::ShPtr< score::ProteinModelScoreSum> &SP_SCORE,
      const storage::Set< quality::Measure> &QUALITY_MEASURES
    ) :
      m_ScoringFunction( SP_SCORE),
      m_Qualities( QUALITY_MEASURES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterScore
    PrinterScore *PrinterScore::Clone() const
    {
      return new PrinterScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief writes a storage table to PDB lines
    //! @param TABLE table to be written
    //! @param WRITE_EMPTY_REMARK_LINE bool whether to precede table with empty remark line
    //! @return PDB lines
    util::ShPtrList< Line> PrinterScore::WriteTableToLines
    (
      const storage::Table< double> &TABLE,
      const bool WRITE_EMPTY_REMARK_LINE
    )
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if the table is not empty
      if( !TABLE.IsEmpty())
      {
        // create an empty remark line
        util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
        empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
        if( WRITE_EMPTY_REMARK_LINE)
        {
          lines.PushBack( empty_remark_line);
        }

        // write to stream
        std::stringstream ss;
        TABLE.WriteFormatted( ss, GetTableColumnFormat(), "Scores");
        storage::Vector< std::string> string_lines( util::StringLineListFromIStream( ss));

        // iterate through the lines
        for
        (
          storage::Vector< std::string>::const_iterator line_itr( string_lines.Begin()),
            line_itr_end( string_lines.End());
          line_itr != line_itr_end; ++line_itr
        )
        {
          // get a new line
          util::ShPtr< Line> new_line( empty_remark_line->Clone());

          // add the entry
          new_line->Put( GetEntryTypes().REMARK_String, *line_itr);
          lines.PushBack( new_line);
        }
      }

      // end
      return lines;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterScore::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize table
      storage::Table< double> table
      (
        math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
      );

      // write score to file
      if( m_ScoringFunction.IsDefined())
      {
        table.Append( m_ScoringFunction->CreateValueTableVertical( PROTEIN_MODEL));
      }

      // add quality measures to table
      const assemble::QualityBatch quality_batch( m_Qualities, biol::GetAtomTypes().CA);
      if( !m_Qualities.IsEmpty()) // if qualities were given, add all the quality measures
      {
        table.Append( quality_batch.ConstructTable( PROTEIN_MODEL));
      }
      else // if no qualities were given, just add the completeness estimate
      {
        table.Append( quality_batch.ContructTableWithCompletenessOnly( PROTEIN_MODEL));
      }

      // return created lines
      return WriteTableToLines( table);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterScore::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_residue.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_line_criterium.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Residue::s_Instance
    (
      GetObjectInstances().AddInstance( new Residue())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Residue::Residue() :
      m_ResidueName(),
      m_ChainID( ' '),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' '),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName and ChainID
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' '),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName, ChainID, ID, ICode
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID,
      const int  &PDBID,
      const char &ICODE
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( PDBID),
      m_ICode( ICODE),
      m_Lines()
    {
    }

    //! construct Residue from ResidueName, ChainID, ID, ICOde, and Lines
    Residue::Residue
    (
      const std::string &RESIDUENAME,
      const char &CHAINID,
      const int  &PDBID,
      const char &ICODE,
      const util::ShPtrList< Line> &LINES
    ) :
      m_ResidueName( RESIDUENAME),
      m_ChainID( CHAINID),
      m_PDBID( PDBID),
      m_ICode( ICODE),
      m_Lines( LINES)
    {
    }

    //! @brief copy constructor
    Residue::Residue( const Residue &RESIDUE) :
      m_ResidueName( RESIDUE.m_ResidueName),
      m_ChainID( RESIDUE.m_ChainID),
      m_PDBID( RESIDUE.m_PDBID),
      m_ICode( RESIDUE.m_ICode),
      m_Lines( RESIDUE.m_Lines)
    {
    }

    //! @brief construct form Residue
    //! @param RESIDUE the residue to copy
    Residue::Residue( const ResidueInterface &RESIDUE) :
      m_ResidueName( RESIDUE.GetResidueName()),
      m_ChainID( RESIDUE.GetChainID()),
      m_PDBID( RESIDUE.GetPDBID()),
      m_ICode( RESIDUE.GetICode()),
      m_Lines()
    {
    }

    //! copy constructor
    Residue *Residue::Clone() const
    {
      return new Residue( *this);
    }

    //! destructor
    Residue::~Residue()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Residue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return Lines
    util::ShPtrList< Line> const &Residue::GetLines() const
    {
      return m_Lines;
    }

    //! change Lines
    util::ShPtrList< Line> &Residue::ChangeLines()
    {
      return m_Lines;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> Residue::AtomSerials() const
    {
      storage::Set< size_t> serials;

      // iterate over all atom lines
      for( util::ShPtrList< Line>::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End()); itr != itr_end; ++itr)
      {
        serials.Insert( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().ATOMSerial));
      }

      // end
      return serials;
    }

    //! @brief line criterium to locate atom/hetatom lines for that residue
    //! @param LINE_TYPE the line type of interest
    //! @return criterium to check if a line corresponds to that residue
    LineCriterium Residue::GetCriterium( const LineType &LINE_TYPE) const
    {
      LineCriterium criterium;

      if( LINE_TYPE == GetLineTypes().ATOM)
      {
        criterium.AddCriterium( GetEntryTypes().ATOMName             , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().ATOMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().ATOMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().ATOMInsertionCode    , m_ICode);
      }
      else if( LINE_TYPE == GetLineTypes().HETATM)
      {
        criterium.AddCriterium( GetEntryTypes().HETATMName             , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().HETATMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().HETATMInsertionCode    , m_ICode);
      }

      // end
      return criterium;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Residue to STREAM
    std::ostream &Residue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ResidueName, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID    , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID      , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode      , OSTREAM) << '\n';
      io::Serialize::Write( m_Lines      , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Residue from std::istream
    std::istream &Residue::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_ResidueName, ISTREAM);
      io::Serialize::Read( m_ChainID    , ISTREAM);
      io::Serialize::Read( m_PDBID      , ISTREAM);
      io::Serialize::Read( m_ICode      , ISTREAM);
      io::Serialize::Read( m_Lines      , ISTREAM);

      // return
      return ISTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool Residue::operator ==( const ResidueInterface &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.GetChainID()
             && m_PDBID   == RESIDUE_RHS.GetPDBID()
             && m_ICode   == RESIDUE_RHS.GetICode();
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool Residue::operator !=( const ResidueInterface &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool Residue::operator <( const ResidueInterface &RESIDUE_RHS) const
    {
      // compar chainid
      const char rhs_chain_id( RESIDUE_RHS.GetChainID());
      if( m_ChainID < rhs_chain_id)
      {
        return true;
      }

      if( m_ChainID > rhs_chain_id)
      {
        return false;
      }

      // compare pdb id
      const int rhs_pdb_id( RESIDUE_RHS.GetPDBID());
      if( m_PDBID < rhs_pdb_id)
      {
        return true;
      }

      if( m_PDBID > rhs_pdb_id)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.GetICode())
      {
        return true;
      }

      return false;
    }

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_LHS lhs residue
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool Residue::operator ==( const Residue &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.m_ChainID
             && m_PDBID   == RESIDUE_RHS.m_PDBID
             && m_ICode   == RESIDUE_RHS.m_ICode;
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool Residue::operator !=( const Residue &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool Residue::operator <( const Residue &RESIDUE_RHS) const
    {
      // compare chain id
      if( m_ChainID < RESIDUE_RHS.m_ChainID)
      {
        return true;
      }

      if( m_ChainID > RESIDUE_RHS.m_ChainID)
      {
        return false;
      }

      // compare pdb id
      if( m_PDBID < RESIDUE_RHS.m_PDBID)
      {
        return true;
      }

      if( m_PDBID > RESIDUE_RHS.m_PDBID)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.m_ICode)
      {
        return true;
      }

      return false;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_residue_simple.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "pdb/bcl_pdb_entry_types.h"
#include "pdb/bcl_pdb_line_criterium.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ResidueSimple::s_Instance
    (
      GetObjectInstances().AddInstance( new ResidueSimple())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ResidueSimple::ResidueSimple() :
      m_ResidueName(),
      m_ChainID( ' '),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' ')
    {
    }

    //! @brief construct ResidueSimple from ResidueName and ChainID
    //! @param RESIDUE_NAME name
    //! @param CHAIN_ID chain id
    ResidueSimple::ResidueSimple
    (
      const std::string &RESIDUE_NAME,
      const char CHAIN_ID
    ) :
      m_ResidueName( RESIDUE_NAME),
      m_ChainID( CHAIN_ID),
      m_PDBID( util::GetUndefined< int>()),
      m_ICode( ' ')
    {
    }

    //! @brief construct ResidueSimple from ResidueName, ChainID, ID, ICode
    //! @param RESIDUE_NAME name
    //! @param CHAIN_ID chain id
    //! @param PDB_ID pdb sequence id
    //! @parmam I_CODE pdb residue insertion code
    ResidueSimple::ResidueSimple
    (
      const std::string &RESIDUE_NAME,
      const char CHAIN_ID,
      const int PDB_ID,
      const char I_CODE
    ) :
      m_ResidueName( RESIDUE_NAME),
      m_ChainID( CHAIN_ID),
      m_PDBID( PDB_ID),
      m_ICode( I_CODE)
    {
    }

    //! @brief copy constructor
    ResidueSimple::ResidueSimple( const ResidueSimple &RESIDUE) :
      m_ResidueName( RESIDUE.m_ResidueName),
      m_ChainID( RESIDUE.m_ChainID),
      m_PDBID( RESIDUE.m_PDBID),
      m_ICode( RESIDUE.m_ICode)
    {
    }

    //! copy constructor
    ResidueSimple *ResidueSimple::Clone() const
    {
      return new ResidueSimple( *this);
    }

    //! destructor
    ResidueSimple::~ResidueSimple()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ResidueSimple::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to pdb ID
    //! @return pdb sequence id
    int ResidueSimple::GetPDBID() const
    {
      return m_PDBID;
    }

    //! @brief set ID
    //! @param PDB_ID new pdb id for residue
    void ResidueSimple::SetPDBID( const int PDB_ID)
    {
      m_PDBID = PDB_ID;
    }

    //! @brief return insertion code
    //! @return pdb sequence insertion code, if multiple residues have the same pdb id
    char ResidueSimple::GetICode() const
    {
      return m_ICode;
    }

    //! @brief set insertion code
    //! @param I_CODE insertion code for that residue
    void ResidueSimple::SetICode( const char I_CODE)
    {
      m_ICode = I_CODE;
    }

    //! @brief return chain id
    //! @return the chain id this residue belongs to
    char ResidueSimple::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief set chain id
    //! @param CHAIN_ID the chain id
    void ResidueSimple::SetChainID( const char CHAIN_ID)
    {
      m_ChainID = CHAIN_ID;
    }

    //! @brief return residue name
    //! @return three letter code residue name
    const std::string &ResidueSimple::GetResidueName() const
    {
      return m_ResidueName;
    }

    //! @brief set residue name
    //! @param RESIDUE_NAME new residue name for that residue
    void ResidueSimple::SetResidueName( const std::string &RESIDUE_NAME)
    {
      m_ResidueName = RESIDUE_NAME;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> ResidueSimple::AtomSerials() const
    {
      return storage::Set< size_t>();
    }

    //! @brief line criterium to locate atom/hetatom lines for that residue
    //! @param LINE_TYPE the line type of interest
    //! @return criterium to check if a line corresponds to that residue
    LineCriterium ResidueSimple::GetCriterium( const LineType &LINE_TYPE) const
    {
      LineCriterium criterium;

      if( LINE_TYPE == GetLineTypes().ATOM)
      {
        criterium.AddCriterium( GetEntryTypes().ATOMResidueName      , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().ATOMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().ATOMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().ATOMInsertionCode    , m_ICode);
      }
      else if( LINE_TYPE == GetLineTypes().HETATM)
      {
        criterium.AddCriterium( GetEntryTypes().HETATMResidueName      , m_ResidueName);
        criterium.AddCriterium( GetEntryTypes().HETATMChainID          , m_ChainID);
        criterium.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, m_PDBID);
        criterium.AddCriterium( GetEntryTypes().HETATMInsertionCode    , m_ICode);
      }

      // end
      return criterium;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write ResidueSimple to STREAM
    std::ostream &ResidueSimple::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ResidueName, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID    , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID      , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode      , OSTREAM) << '\n';

      // end
      return OSTREAM;
    }

    //! read ResidueSimple from std::istream
    std::istream &ResidueSimple::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_ResidueName, ISTREAM);
      io::Serialize::Read( m_ChainID    , ISTREAM);
      io::Serialize::Read( m_PDBID      , ISTREAM);
      io::Serialize::Read( m_ICode      , ISTREAM);

      // return
      return ISTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool ResidueSimple::operator ==( const ResidueInterface &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.GetChainID()
             && m_PDBID   == RESIDUE_RHS.GetPDBID()
             && m_ICode   == RESIDUE_RHS.GetICode();
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool ResidueSimple::operator !=( const ResidueInterface &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool ResidueSimple::operator <( const ResidueInterface &RESIDUE_RHS) const
    {
      // compare chainid
      const char rhs_chain_id( RESIDUE_RHS.GetChainID());
      if( m_ChainID < rhs_chain_id)
      {
        return true;
      }

      if( m_ChainID > rhs_chain_id)
      {
        return false;
      }

      // compare pdb id
      const int rhs_pdb_id( RESIDUE_RHS.GetPDBID());
      if( m_PDBID < rhs_pdb_id)
      {
        return true;
      }

      if( m_PDBID > rhs_pdb_id)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.GetICode())
      {
        return true;
      }

      return false;
    }

    //! @brief compare two residues if they are equal
    //! @param RESIDUE_LHS lhs residue
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid and insertion code are equal
    bool ResidueSimple::operator ==( const ResidueSimple &RESIDUE_RHS) const
    {
      return    m_ChainID == RESIDUE_RHS.m_ChainID
             && m_PDBID   == RESIDUE_RHS.m_PDBID
             && m_ICode   == RESIDUE_RHS.m_ICode;
    }

    //! @brief compare two residues if they are not equal
    //! @param RESIDUE_RHS rhs residue
    //! @return true if residue name, chain id, pdbid or insertion code are not equal
    bool ResidueSimple::operator !=( const ResidueSimple &RESIDUE_RHS) const
    {
      return !operator==( RESIDUE_RHS);
    }

    //! @brief compare two residues if the lhs is less than the rhs
    //! @param RESIDUE_RHS rhs residue
    //! @return true if chain id or pdbid or insertion code of lhs is smaller
    bool ResidueSimple::operator <( const ResidueSimple &RESIDUE_RHS) const
    {
      // compare chain id
      if( m_ChainID < RESIDUE_RHS.m_ChainID)
      {
        return true;
      }

      if( m_ChainID > RESIDUE_RHS.m_ChainID)
      {
        return false;
      }

      // compare pdb id
      if( m_PDBID < RESIDUE_RHS.m_PDBID)
      {
        return true;
      }

      if( m_PDBID > RESIDUE_RHS.m_PDBID)
      {
        return false;
      }

      // compare insertion code
      if( m_ICode < RESIDUE_RHS.m_ICode)
      {
        return true;
      }

      return false;
    }

    //! @brief compare with biol amino acid if they are equal
    //! @param AMINO_ACID right hand side amino acid to compare with
    //! @return true, if chain id, pdb sequence id and insertion code match
    bool ResidueSimple::operator ==( const biol::AABase &AMINO_ACID) const
    {
      return    m_ChainID == AMINO_ACID.GetChainID()
             && m_PDBID   == AMINO_ACID.GetPdbID()
             && m_ICode   == AMINO_ACID.GetPdbICode();
    }

    //! @brief compare with biol amino acid if they are not equal
    //! @param AMINO_ACID right hand side amino acid to compare with
    //! @return true, if either chain id, pdb sequence id and insertion code do not match
    bool ResidueSimple::operator !=( const biol::AABase &AMINO_ACID) const
    {
      return !operator==( AMINO_ACID);
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_site.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

    //! @brief EvidenceCode as string
    //! @param EVIDENCE_CODE the EvidenceCode
    //! @return the string for the EvidenceCode
    const std::string &Site::GetEvidenceCodeDescriptor( const Site::EvidenceCode &EVIDENCE_CODE)
    {
      static const std::string s_descriptors[] =
      {
          "SOFTWARE",
          "AUTHOR",
          "UNKNOWN",
          GetStaticClassName< EvidenceCodeEnum>()
      };

      return s_descriptors[ EVIDENCE_CODE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Site::s_Instance
    (
      GetObjectInstances().AddInstance( new Site())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Site::Site()
    {
    }

    //! @brief construct from identifier, evidence code and description
    Site::Site( const std::string &NAME, const EvidenceCodeEnum &EVIDENCE_CODE, const std::string &DESCRIPTION) :
      m_Identifier( NAME),
      m_EvidenceCode( EVIDENCE_CODE),
      m_Description( DESCRIPTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Site
    Site *Site::Clone() const
    {
      return new Site( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Site::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief name/identifier of the site
    //! @return the name as string
    const std::string &Site::GetName() const
    {
      return m_Identifier;
    }

    //! @brief get evidence code
    //! @return enum EvicenceCodeEnum
    const Site::EvidenceCodeEnum &Site::GetEvidenceCode() const
    {
      return m_EvidenceCode;
    }

    //! @brief site description
    //! @return description of site
    const std::string &Site::GetDescription() const
    {
      return m_Description;
    }

    //! @brief access to the residues comprising the site that are within a chain
    //! @return ShPtrList to residues
    const storage::List< ResidueSimple> &Site::GetChainResidues() const
    {
      return m_ChainResidues;
    }

    //! @brief ligand that might be associated with that site
    //! @return ShPtr to Ligand; null if there is no ligand
    const util::ShPtr< Ligand> &Site::GetLigand() const
    {
      return m_Ligand;
    }

    //! @brief add a residue from a chain
    //! @param RESIDUE the residue to  add to the site
    void Site::AddChainResidue( const ResidueSimple &RESIDUE)
    {
      m_ChainResidues.PushBack( RESIDUE);
    }

    //! @brief add a residue from hetatm section
    //! @param RESIDUE the residue to  add to the site
    void Site::AddHetatmResidue( const ResidueSimple &RESIDUE)
    {
      m_HetResidues.PushBack( RESIDUE);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locators for the residues associated with the site definition
    //! @return the amino acid locators in a list
    storage::List< assemble::LocatorAA> Site::AALocators() const
    {
      storage::List< assemble::LocatorAA> aa_locators;

      // iterate through all residues
      for
      (
        storage::List< ResidueSimple>::const_iterator itr( m_ChainResidues.Begin()), itr_end( m_ChainResidues.End());
        itr != itr_end;
        ++itr
      )
      {
        aa_locators.PushBack( assemble::LocatorAA( itr->GetChainID(), itr->GetPDBID(), true));
      }

      // end
      return aa_locators;
    }

    //! @brief find ligand for site
    //! @param LIGANDS list of ligands, that might be associated with that site
    //! @return true, if a ligand was found, false otherwise
    bool Site::FindLigand( const util::ShPtrList< Ligand> &LIGANDS)
    {
      // check if this is a binding site
      static const std::string s_binding_site_string( "BINDING SITE FOR RESIDUE");

      // check if site is binding site for ligand
      if( m_Description.find( s_binding_site_string) == std::string::npos)
      {
        return false;
      }

      // ligand information
      const size_t pos_ligand_info( s_binding_site_string.length() + 1);
      const size_t ligand_info_length( 12);

      // name, chain and seq id
      std::stringstream ligand_info( util::TrimString( m_Description.substr( pos_ligand_info, ligand_info_length)));

      // ligand name
      std::string lig_name;
      ligand_info >> lig_name;
      lig_name = util::Format().R().W( 3).Fill( ' ')( lig_name);

      // chainid
      char chain_id( ' ');
      while( ligand_info.good() && chain_id == ' ')
      {
        ligand_info.get( chain_id);
      }

      // still expecting seq id
      if( !ligand_info.good())
      {
        return false;
      }

      std::string pdb_id_string;
      ligand_info >> pdb_id_string;

      if( util::TrimString( pdb_id_string).empty())
      {
        return false;
      }

      // insertion code could be the last character
      char pdb_icode( pdb_id_string[ pdb_id_string.length() - 1]);

      // if the insertion code is a digit, it is actually part of the seqid, and no insertion code is given
      if( std::isdigit( pdb_icode))
      {
        pdb_icode = ' ';
      }
      else
      {
        // found insertion code, the seq id string is one shorter
        pdb_id_string.erase( pdb_id_string.length() - 1, 1);
      }

      // the remaining string should be numerical
      if( !util::IsNumerical( pdb_id_string))
      {
        return false;
      }

      // the actual sequence id can be extracted
      const int seq_id( util::ConvertStringToNumericalValue< int>( pdb_id_string));

      BCL_MessageVrb
      (
        "site binds ligand: " + lig_name + ' ' + chain_id + ' ' + util::Format()( seq_id) + ' ' + pdb_icode
      );

      // iterate over ligands to find the correct one
      for
      (
        util::ShPtrList< Ligand>::const_iterator itr( LIGANDS.Begin()), itr_end( LIGANDS.End()); itr != itr_end; ++itr
      )
      {
        if( ( *itr)->GetChainID()     != chain_id)  continue;
        if( ( *itr)->GetResidueName() != lig_name)  continue;
        if( ( *itr)->GetPDBID()       != seq_id)    continue;
        if( ( *itr)->GetICode()       != pdb_icode) continue;

        // ligand found
        m_Ligand = *itr;
        return true;
      }

      // nothing found
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Site::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Identifier   , ISTREAM);
      io::Serialize::Read( m_EvidenceCode , ISTREAM);
      io::Serialize::Read( m_Description  , ISTREAM);
      io::Serialize::Read( m_ChainResidues, ISTREAM);
      io::Serialize::Read( m_HetResidues  , ISTREAM);
      io::Serialize::Read( m_Ligand       , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Site::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Identifier   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EvidenceCode , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainResidues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HetResidues  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ligand       , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
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
#include "pdb/bcl_pdb_tail.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_head.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Tail::s_Instance
    (
      GetObjectInstances().AddInstance( new Tail())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Tail::Tail() :
      m_MasterRecord( new Line( GetLineTypes().MASTER)),
      m_End( new Line( GetLineTypes().END))
    {
      InitializeMasterRecord();
    }

    //! @brief Clone function
    //! @return pointer to new Tail
    Tail *Tail::Clone() const
    {
      return new Tail( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Tail::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Tail::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().CONECT.GetIterator(),
        GetLineTypes().End()
      );

      // end
      return s_line_types;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Tail::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      if( LINE_TYPE == m_MasterRecord->GetType())
      {
        lines.PushBack( m_MasterRecord);
      }
      else if( LINE_TYPE == m_End->GetType())
      {
        lines.PushBack( m_End);
      }
      else if( LINE_TYPE == GetLineTypes().CONECT)
      {
        lines.Append( m_ConectLines);
      }

      // end
      return lines;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Tail::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // check conect lines
      lines.Append( LineCriterium::Filter( m_ConectLines, CRITERIUM));

      // master
      if( CRITERIUM( *m_MasterRecord))
      {
        lines.PushBack( m_MasterRecord);
      }

      // end
      if( CRITERIUM( *m_End))
      {
        lines.PushBack( m_End);
      }

      // end
      return lines;
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Tail::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( LINE->GetType() == GetLineTypes().CONECT)
      {
        m_ConectLines.PushBack( LINE);
        return true;
      }
      else if( LINE->GetType() == m_MasterRecord->GetType())
      {
        m_MasterRecord = LINE;
        return true;
      }
      else if( LINE->GetType() == m_End->GetType())
      {
        m_End = LINE;
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Tail::Reset()
    {
      m_ConectLines.Reset();
      InitializeMasterRecord();
    }

    //! @brief connections defined through atom serials in the CONNECT lines
    //! @return Map of atom serial center atoms and a set of serials for all connected atoms
    storage::Map< size_t, storage::Set< size_t> > Tail::GetConnections() const
    {
      const util::ShPtrList< Line> &lines( GetLines( GetLineTypes().CONECT));

      storage::Map< size_t, storage::Set< size_t> > connections;

      // iterate over the connection lines
      for( util::ShPtrList< Line>::const_iterator itr( lines.Begin()), itr_end( lines.End()); itr != itr_end; ++itr)
      {
        const size_t center( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().CONECTAtomSerial));
        storage::Set< size_t> &connected( connections[ center]);

        // iterate through all connected serials
        for
        (
          EntryTypes::const_iterator
            entry_itr( GetEntryTypes().CONECTBondedAtomSerial1.GetIterator()),
            entry_itr_end( GetEntryTypes().CONECTBondedAtomSerial4.GetIterator() + 1);
          entry_itr != entry_itr_end;
          ++entry_itr
        )
        {
          const size_t serial( ( *itr)->GetNumericalValue< size_t>( *entry_itr));
          if( !util::IsDefined( serial))
          {
            continue;
          }
          connected.Insert( connected.End(), serial);
        }
      }

      // end
      return connections;
    }

    //! @brief update the MASTER record for header information
    //! @param HEAD a pdb head section will all lines before the coordinate section
    void Tail::UpdateMasterRecord( const Head &HEAD) const
    {
      // Number of REMARK records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumRemark, HEAD.Count( GetLineTypes().REMARK));

      // Number of HET records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumHet, HEAD.Count( GetLineTypes().HET));

      // Number of HELIX records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumHelix, HEAD.Count( GetLineTypes().HELIX));

      // Number of SHEET records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSheet, HEAD.Count( GetLineTypes().SHEET));

      // Number of SITE records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSite, HEAD.Count( GetLineTypes().SITE));

      // Number of coordinate transformation records  (ORIGX+SCALE+MTRIX)
      {
        size_t number_trans( 0);
        for
        (
          LineTypes::const_iterator
            itr( GetLineTypes().ORIGX1.GetIterator()), itr_end( GetLineTypes().MTRIX3.GetIterator() + 1);
          itr != itr_end;
          ++itr
        )
        {
          number_trans += HEAD.Count( *itr);
        }
        m_MasterRecord->Put( GetEntryTypes().MASTERNumXform, number_trans);
      }

      // Number of CONECT records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumConect, m_ConectLines.GetSize());

      // Number of SEQRES records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSeqres, HEAD.Count( GetLineTypes().SEQRES));
    }

    //! @brief update the MASTER record for coordinate information
    //! @param FIRST_MODEL only use the first model
    void Tail::UpdateMasterRecord( const Model &FIRST_MODEL) const
    {
      // Number of atomic coordinate records (ATOM+HETATM)
      m_MasterRecord->Put
      (
        GetEntryTypes().MASTERNumCoord,
        FIRST_MODEL.Count( GetLineTypes().ATOM) + FIRST_MODEL.Count( GetLineTypes().HETATM)
      );

      // Number of TER records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumTer, FIRST_MODEL.Count( GetLineTypes().TER));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Tail::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through all lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( m_ConectLines.Begin()), itr_end( m_ConectLines.End());
        itr != itr_end;
        ++itr
      )
      {
        OSTREAM << ( *itr)->GetString() << '\n';
      }

      // master
      OSTREAM << m_MasterRecord->GetString() << '\n';

      // end
      OSTREAM << m_End->GetString() << '\n';

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Tail::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConectLines , ISTREAM);
      io::Serialize::Read( m_MasterRecord, ISTREAM);
      io::Serialize::Read( m_End         , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Tail::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConectLines , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MasterRecord, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_End         , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize the master record with zeros
    void Tail::InitializeMasterRecord()
    {
      // iterate through all entries
      for
      (
        EntryTypes::const_iterator
          itr( GetEntryTypes().MASTERNumRemark.GetIterator()),
          itr_end( GetEntryTypes().MASTERNumSeqres.GetIterator() + 1);
        itr != itr_end;
        ++itr
      )
      {
        m_MasterRecord->Put( *itr, 0);
      }
    }

  } // namespace pdb
} // namespace bcl
