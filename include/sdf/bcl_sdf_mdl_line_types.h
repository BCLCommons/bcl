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

#ifndef BCL_SDF_MDL_LINE_TYPES_H_
#define BCL_SDF_MDL_LINE_TYPES_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically
#include <string>

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_sdf_mdl_line_types.h
    //! @brief Enum of different mdl line types
    //!
    //! @see @link example_sdf_mdl_entry_types.cpp @endlink
    //! @author butkiem1, mendenjl
    //! @date May 7, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    enum MdlLineType
    {
      e_HeaderLine,         //!< Header line, should have V2000 in it
      e_AtomLine,           //!< Atom line
      e_BondLine,           //!< Bond line
      e_TerminationLine,    //!< $$$$    Terminates the molecule
      e_RXNHeaderLine,      //!< Header line for RXN files
      e_RXNStartLine,       //!< $RXN    Begins a reaction block
      e_RXNMolStartLine,    //!< $MOL    Terminates a molecule of a reaction
      e_RXNTerminationLine, //!< empty    Terminates a molecule of a reaction
      s_NumberLineTypes     //!< Unknown line type
    };

    //! @brief LineType as string
    //! @param LINE_TYPE the line type for which a name is desired
    //! @return the line type as string
    BCL_API const std::string &GetLineTypeName( const MdlLineType &LINE_TYPE);

    //! Wrapper for MdlLineType enum
    typedef util::WrapperEnum< MdlLineType, &GetLineTypeName, s_NumberLineTypes> MdlLineTypeEnum;

    //! @brief gather all entry types for given line type
    //! @param LINE_TYPE line type the entry types should be collected for
    //! @return set of entry types for that line type
    BCL_API const std::string &GetDefaultLine( const MdlLineType &LINE_TYPE);

  } // namespace sdf

} // namespace bcl

#endif // BCL_SDF_MDL_LINE_TYPES_H_
