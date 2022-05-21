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

#ifndef BCL_PDB_LINE_TYPES_H_
#define BCL_PDB_LINE_TYPES_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line_type_data.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineTypes
    //! @brief enumeration of LineTypeData for pdb file line types
    //! @details In pdb files, first word on each line specifies the type of that line such as ATOM, SEQRES, HELIX
    //! etc. This class provides the enumeration of all the line types as defined by PDB
    //!
    //! @see @link example_pdb_line_types.cpp @endlink
    //! @author woetzen
    //! @date Feb 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineTypes :
      public util::Enumerate< LineTypeData, LineTypes>
    {
      friend class util::Enumerate< LineTypeData, LineTypes>;

    public:

    //////////
    // data //
    //////////

      LineType HEADER;
      LineType OBSLTE;
      LineType TITLE;
      LineType SPLIT;
      LineType CAVEAT;
      LineType COMPND;
      LineType SOURCE;
      LineType KEYWDS;
      LineType EXPDTA;
      LineType NUMMDL;
      LineType MDLTYP;
      LineType AUTHOR;
      LineType REVDAT;
      LineType SPRSDE;
      LineType JRNL;
      LineType REMARK;
      LineType DBREF;
      LineType DBREF1;
      LineType DBREF2;
      LineType SEQADV;
      LineType SEQRES;
      LineType MODRES;
      LineType HET;
      LineType HETNAM;
      LineType HETSYN;
      LineType FORMUL;
      LineType HELIX;
      LineType SHEET;
      LineType SSBOND;
      LineType LINK;
      LineType CISPEP;
      LineType SITE;
      LineType CRYST1;
      LineType CRYST2;
      LineType CRYST3;
      LineType ORIGX1;
      LineType ORIGX2;
      LineType ORIGX3;
      LineType SCALE1;
      LineType SCALE2;
      LineType SCALE3;
      LineType MTRIX1;
      LineType MTRIX2;
      LineType MTRIX3;
      LineType MODEL;
      LineType ATOM;
      LineType ANISOU;
      LineType TER;
      LineType HETATM;
      LineType ENDMDL;
      LineType CONECT;
      LineType MASTER;
      LineType END;

      //! types for grouping
      storage::Set< LineType> m_GroupingTypes;

      //! types that can be out of order
      storage::Set< LineType> m_OutOfOrderTypes;

      //! starting position of line type record
      static const size_t s_TypeRecordStart = 0;

      //! length of line type record
      static const size_t s_TypeRecordLength = 6;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all LineTypes
      LineTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set of grouping line types
      //! types that are used for grouping
      //! @return set of types that are used for grouping
      const storage::Set< LineType> &GetGroupingTypes() const;

      //! @brief set of line types that can appear out of order
      //! @return set of types that can appear out of order
      const storage::Set< LineType> &GetOutOfOrderTypes() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add a new linetype to the enumeration of line types
      //! @param RECORD_NAME name of record e.g "ATOM"
      //! @param LINE_TYPE_DATA LineTypeDataObject
      LineType AddLineType
      (
        const std::string &RECORD_NAME,
        const LineTypeData &LINE_TYPE_DATA
      );

      //! @brief retrieve LineType from the pdb line string
      //! @param PDB_LINE string of the line in the pdb file
      //! @return LineType for that string
      LineType LineTypeFromPDBLine( const std::string &PDB_LINE) const;

      //! @brief consider entry type for given line type
      //! @param LINE_TYPE line type for which the entry type is considered
      //! @param ENTRY_TYPE the entry type
      void ConsiderEntryType( const LineType &LINE_TYPE, const EntryType &ENTRY_TYPE);

    }; // LineTypes

    BCL_API const LineTypes &GetLineTypes();

  } // namespace pdb

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< pdb::LineTypeData, pdb::LineTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_PDB_LINE_TYPES_H_
