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
