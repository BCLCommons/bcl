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
#include "sdf/bcl_sdf_mdl_entry_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief default constructor
    MdlEntryTypes::MdlEntryTypes() :
      Header_NumberAtomLines(        AddEnum( "Header_NumberAtomLines",        MdlEntryTypeData( e_HeaderLine,       0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberBondLines(        AddEnum( "Header_NumberBondLines",        MdlEntryTypeData( e_HeaderLine,       3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberAtomLists(        AddEnum( "Header_NumberAtomLists",        MdlEntryTypeData( e_HeaderLine,       6,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete1(              AddEnum( "Header_Obsolete1",              MdlEntryTypeData( e_HeaderLine,       9,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_ChiralFlag(             AddEnum( "Header_ChiralFlag",             MdlEntryTypeData( e_HeaderLine,      12,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_NumberSTextEntries(     AddEnum( "Header_NumberSTextEntries",     MdlEntryTypeData( e_HeaderLine,      15,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete2(              AddEnum( "Header_Obsolete2",              MdlEntryTypeData( e_HeaderLine,      18,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete3(              AddEnum( "Header_Obsolete3",              MdlEntryTypeData( e_HeaderLine,      21,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete4(              AddEnum( "Header_Obsolete4",              MdlEntryTypeData( e_HeaderLine,      24,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete5(              AddEnum( "Header_Obsolete5",              MdlEntryTypeData( e_HeaderLine,      27,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_Obsolete6(              AddEnum( "Header_Obsolete6",              MdlEntryTypeData( e_HeaderLine,      30,   3,    "999", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Header_CtabVersion(            AddEnum( "Header_CtabVersion",            MdlEntryTypeData( e_HeaderLine,      33,   6,  "V2000", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXNHeader_NumberReactantLines( AddEnum( "RXNHeader_NumberReactantLines", MdlEntryTypeData( e_RXNHeaderLine,    0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      RXNHeader_NumberProductLines(  AddEnum( "RXNHeader_NumberProductLines",  MdlEntryTypeData( e_RXNHeaderLine,    3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Atom_CoordinateX(              AddEnum( "Atom_CoordinateX",              MdlEntryTypeData( e_AtomLine,         0,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_CoordinateY(              AddEnum( "Atom_CoordinateY",              MdlEntryTypeData( e_AtomLine,        10,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_CoordinateZ(              AddEnum( "Atom_CoordinateZ",              MdlEntryTypeData( e_AtomLine,        20,  10, "0.0000", util::Format().R().FFP( 4), util::CPPDataTypes::e_Double))),
      Atom_Symbol(                   AddEnum( "Atom_Symbol",                   MdlEntryTypeData( e_AtomLine,        31,   3,      "X", util::Format().L(),         util::CPPDataTypes::e_String))),
      Atom_MassDifference(           AddEnum( "Atom_MassDifference",           MdlEntryTypeData( e_AtomLine,        34,   2,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_Charge(                   AddEnum( "Atom_Charge",                   MdlEntryTypeData( e_AtomLine,        36,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_StereoParity(             AddEnum( "Atom_StereoParity",             MdlEntryTypeData( e_AtomLine,        39,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_HydrogenCount(            AddEnum( "Atom_HydrogenCount",            MdlEntryTypeData( e_AtomLine,        42,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_StereoCareBox(            AddEnum( "Atom_StereoCareBox",            MdlEntryTypeData( e_AtomLine,        45,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_Valence(                  AddEnum( "Atom_Valence",                  MdlEntryTypeData( e_AtomLine,        48,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_HODesignator(             AddEnum( "Atom_HODesignator",             MdlEntryTypeData( e_AtomLine,        51,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_NotUsed1(                 AddEnum( "Atom_NotUsed1",                 MdlEntryTypeData( e_AtomLine,        54,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_NotUsed2(                 AddEnum( "Atom_NotUsed2",                 MdlEntryTypeData( e_AtomLine,        57,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_AtomMappingNumber(        AddEnum( "Atom_AtomMappingNumber",        MdlEntryTypeData( e_AtomLine,        60,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_InversionFlag(            AddEnum( "Atom_InversionFlag",            MdlEntryTypeData( e_AtomLine,        63,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Atom_ExactChangeFlag(          AddEnum( "Atom_ExactChangeFlag",          MdlEntryTypeData( e_AtomLine,        66,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_FirstAtomIndex(           AddEnum( "Bond_FirstAtomIndex",           MdlEntryTypeData( e_BondLine,         0,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Bond_SecondAtomIndex(          AddEnum( "Bond_SecondAtomIndex",          MdlEntryTypeData( e_BondLine,         3,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_SizeT))),
      Bond_Type(                     AddEnum( "Bond_Type",                     MdlEntryTypeData( e_BondLine,         6,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_Stereo(                   AddEnum( "Bond_Stereo",                   MdlEntryTypeData( e_BondLine,         9,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_NotUsed(                  AddEnum( "Bond_NotUsed",                  MdlEntryTypeData( e_BondLine,        12,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_Topology(                 AddEnum( "Bond_Topology",                 MdlEntryTypeData( e_BondLine,        15,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      Bond_ReactingCenterStatus(     AddEnum( "Bond_ReactingCenterStatus",     MdlEntryTypeData( e_BondLine,        18,   3,      "0", util::Format().R(),         util::CPPDataTypes::e_Int))),
      RXN_RXNStartLine(              AddEnum( "RXN_RXNStartLine",              MdlEntryTypeData( e_RXNStartLine,     0,   4,   "$RXN", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXN_MolStartLine(              AddEnum( "RXN_MolStartLine",              MdlEntryTypeData( e_RXNMolStartLine,  0,   4,   "$MOL", util::Format().R(),         util::CPPDataTypes::e_String))),
      RXN_RXNTerminationLine(        AddEnum( "RXN_TerminationLine",           MdlEntryTypeData( e_RXNTerminationLine, 0, 0,       "", util::Format().R(),         util::CPPDataTypes::e_String))),
      Terminator(                    AddEnum( "Terminator",                    MdlEntryTypeData( e_TerminationLine,  0,   4,   "$$$$", util::Format().R(),         util::CPPDataTypes::e_String)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MdlEntryTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief global access to all entry types
    const MdlEntryTypes &GetMdlEntryTypes()
    {
      return MdlEntryTypes::GetEnums();
    }

  } // namespace sdf

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< sdf::MdlEntryTypeData, sdf::MdlEntryTypes>;

  } // namespace util
} // namespace bcl
