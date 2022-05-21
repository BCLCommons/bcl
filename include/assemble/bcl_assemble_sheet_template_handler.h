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

#ifndef BCL_ASSEMBLE_SHEET_TEMPLATE_HANDLER_H_
#define BCL_ASSEMBLE_SHEET_TEMPLATE_HANDLER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_compare.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SheetTemplateHandler
    //! @brief This is class provides storage for geometrical templates for beta-sheets
    //! @details This class is the handler for storing and accessing the beta-sheet templates. If provided with a set
    //! set of pdbs to exclude, then access with GetFilteredTemplates() excludes beta-sheet templates from these given
    //! pdbs. It allows picking a beta-sheet template that has a given number of strands, in addition to being able to
    //! pick templates not only the number of strands but also on the lengths of strands. The sheet-templates are pre-
    //! processed and are stored in a simplified representation in the compressed file sheet_templates.input.bz2.
    //!
    //! @see @link example_assemble_sheet_template_handler.cpp @endlink
    //! @author karakam, fischea
    //! @date Oct 10, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SheetTemplateHandler :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! static map to store the templates
      static storage::Map< size_t, util::ShPtrVector< FoldTemplate> > s_TemplateMap;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SheetTemplateHandler
      SheetTemplateHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return command line flag for using sheet templates
      //! @return command line flag for using sheet templates
      static util::ShPtr< command::FlagInterface> &GetFlagSheetTemplates();

    ////////////////
    // operations //
    ////////////////

      //! @brief picks a random template of the appropriate size
      //! @param NR_STRANDS number of strands to be in the template
      //! @return a random template of the appropriate size
      static const FoldTemplate &GetRandomTemplate( const size_t NR_STRANDS);

      //! @brief picks a random template of the appropriate size and geometry lengths
      //! @param SSES sses used to chose subtemplate based on length
      //! @param SSE_GEOMETRY_COMPARE comparison method
      //! @return a random template
      static const FoldTemplate &GetRandomTemplate
      (
        const util::SiPtrVector< const SSE> &SSES,
        const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE =
          SSEGeometryWithinSizeTolerance()
      );

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief singleton function to return sheet templates previously generated and stored in a file
      //! @return sheet templates previously generated and stored in a file
      static storage::Map< size_t, util::ShPtrVector< FoldTemplate> > &GetTemplates();

      //! @brief singleton function to return sheet templates previously generated and stored in a file
      //! @param FILENAMES set of filenames to be read in
      //! @return sheet templates previously generated and stored in a file
      static storage::Map< size_t, util::ShPtrVector< FoldTemplate> > ReadTemplates
      (
        const storage::Set< std::string> &FILENAMES
      );

      //! @brief default input filename
      //! @param TEMPLATE_TYPE type of templates to use
      //! @return default input filename
      static const std::string &GetInputFilename( const std::string &TEMPLATE_TYPE);

      //! @brief return set of excluded pdb ids
      //! @return set of excluded pdb ids
      static storage::Set< std::string> &GetExcludedPdbs();

      //! @brief initialize the sheet template topology
      //! @param SHEET_TEMPLATE sheet template for which topology will be initialized
      static void InitializeTopology( FoldTemplate &SHEET_TEMPLATE);

    }; // class SheetTemplateHandler

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SHEET_TEMPLATE_HANDLER_H_ 
