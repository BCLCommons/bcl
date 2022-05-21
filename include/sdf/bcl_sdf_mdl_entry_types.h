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

#ifndef BCL_SDF_MDL_ENTRY_TYPES_H_
#define BCL_SDF_MDL_ENTRY_TYPES_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_entry_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MdlEntryTypes
    //! @brief TODO: document
    //!
    //! @see @link example_sdf_mdl_entry_types.cpp @endlink
    //! @author butkiem1
    //! @date
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MdlEntryTypes :
      public util::Enumerate< MdlEntryTypeData, MdlEntryTypes>
    {
      friend class util::Enumerate< MdlEntryTypeData, MdlEntryTypes>;

    public:

    //////////
    // data //
    //////////

      const MdlEntryType Header_NumberAtomLines;
      const MdlEntryType Header_NumberBondLines;
      const MdlEntryType Header_NumberAtomLists;
      const MdlEntryType Header_Obsolete1;
      const MdlEntryType Header_ChiralFlag;
      const MdlEntryType Header_NumberSTextEntries;
      const MdlEntryType Header_Obsolete2;
      const MdlEntryType Header_Obsolete3;
      const MdlEntryType Header_Obsolete4;
      const MdlEntryType Header_Obsolete5;
      const MdlEntryType Header_Obsolete6;
      const MdlEntryType Header_CtabVersion;
      const MdlEntryType RXNHeader_NumberReactantLines;
      const MdlEntryType RXNHeader_NumberProductLines;
      const MdlEntryType Atom_CoordinateX;
      const MdlEntryType Atom_CoordinateY;
      const MdlEntryType Atom_CoordinateZ;
      const MdlEntryType Atom_Symbol;
      const MdlEntryType Atom_MassDifference;
      const MdlEntryType Atom_Charge;
      const MdlEntryType Atom_StereoParity;
      const MdlEntryType Atom_HydrogenCount;
      const MdlEntryType Atom_StereoCareBox;
      const MdlEntryType Atom_Valence;
      const MdlEntryType Atom_HODesignator;
      const MdlEntryType Atom_NotUsed1;
      const MdlEntryType Atom_NotUsed2;
      const MdlEntryType Atom_AtomMappingNumber;
      const MdlEntryType Atom_InversionFlag;
      const MdlEntryType Atom_ExactChangeFlag;
      const MdlEntryType Bond_FirstAtomIndex;
      const MdlEntryType Bond_SecondAtomIndex;
      const MdlEntryType Bond_Type;
      const MdlEntryType Bond_Stereo;
      const MdlEntryType Bond_NotUsed;
      const MdlEntryType Bond_Topology;
      const MdlEntryType Bond_ReactingCenterStatus;
      const MdlEntryType RXN_RXNStartLine;
      const MdlEntryType RXN_MolStartLine;
      const MdlEntryType RXN_RXNTerminationLine;
      const MdlEntryType Terminator;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all single entry type items
      MdlEntryTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class MdlEntryTypes

    //! @brief global access to all entry types
    BCL_API
    const MdlEntryTypes &GetMdlEntryTypes();

  } // namespace sdf

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< sdf::MdlEntryTypeData, sdf::MdlEntryTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_SDF_MDL_ENTRY_TYPES_H_
