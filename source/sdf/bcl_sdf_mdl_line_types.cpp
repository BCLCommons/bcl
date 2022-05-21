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
#include "sdf/bcl_sdf_mdl_line_types.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief LineType as string
    //! @param LINE_TYPE the line type for which a name is desired
    //! @return the line type as string
    const std::string &GetLineTypeName( const MdlLineType &LINE_TYPE)
    {
      static const std::string s_names[] =
      {
        "Header",
        "Atom",
        "Bond",
        "Terminal",
        "RXNHeader",
        "RXNStart",
        "RXNMolStart",
        "RXNTermination",
        GetStaticClassName< MdlLineType>()
      };
      return s_names[ LINE_TYPE];
    }

    //! @brief create a vector containing the default lines for each line type
    //! @return a vector containing the default lines for each line type
    storage::Vector< std::string> GetDefaultLines()
    {
      // length of each line
      storage::Vector< size_t> line_length( s_NumberLineTypes, size_t( 0));

      // iterate over all entry types
      for
      (
        MdlEntryTypes::const_iterator itr( GetMdlEntryTypes().Begin()),
        itr_end( GetMdlEntryTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        const MdlLineType type( ( *itr)->GetMdlLineType());

        // get the position this entry ends at
        const size_t end_pos( ( *itr)->GetStart() + ( *itr)->GetLength());

        // set line length to the max of this end_pos and the existing end_pos
        line_length( size_t( type)) = std::max( line_length( size_t( type)), end_pos);
      }

      storage::Vector< std::string> lines( s_NumberLineTypes);

      // initialize all line types with their length
      for( size_t line_type_number( 0); line_type_number < s_NumberLineTypes; ++line_type_number)
      {
        // create a default string full of spaces
        lines( line_type_number) = std::string( line_length( line_type_number), ' ');
      }

      // iterate over all entry types; set the default value in the corresponding string
      for
      (
        MdlEntryTypes::const_iterator itr( GetMdlEntryTypes().Begin()), itr_end( GetMdlEntryTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        const MdlEntryType &type( *itr);
        type->Set( lines( type->GetMdlLineType()), type->GetDefault());
      }
      return lines;
    }

    //! @brief gather all entry types for given line type
    //! @param LINE_TYPE line type the entry types should be collected for
    //! @return set of entry types for that line type
    const std::string &GetDefaultLine( const MdlLineType &LINE_TYPE)
    {
      static const storage::Vector< std::string> default_lines( GetDefaultLines());
      if( LINE_TYPE == s_NumberLineTypes)
      {
        static const std::string undefined_line;
        return undefined_line;
      }
      return default_lines( size_t( LINE_TYPE));
    }

  } // namespace sdf
} // namespace bcl
