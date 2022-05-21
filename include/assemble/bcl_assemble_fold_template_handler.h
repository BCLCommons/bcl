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

#ifndef BCL_ASSEMBLE_FOLD_TEMPLATE_HANDLER_H_
#define BCL_ASSEMBLE_FOLD_TEMPLATE_HANDLER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_compare.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FoldTemplateHandler
    //! @brief This class is designed to store and create fold templates
    //! @details This class is the handler for creating fold templates for a given set of PDBs. By default, it uses
    //! the fold_templates.input.bz2 compressed file which stores a simplified representation for fold templates for
    //! each protein in the dataset. The map of fold templates organized by the number of helices and strands is
    //! initialized only once via a singleton access function
    //!
    //! @see @link example_assemble_fold_template_handler.cpp @endlink
    //! @author weinerbe
    //! @date Feb 9, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FoldTemplateHandler :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! static map to store the templates
      static storage::Map< storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate> > s_TemplateMap;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FoldTemplateHandler();

      //! @brief Clone function
      //! @return pointer to new FoldTemplateHandler
      FoldTemplateHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return command line flag for using fold templates
      //! @return command line flag for using fold templates
      static util::ShPtr< command::FlagInterface> &GetFlagFoldTemplates();

    ////////////////
    // operations //
    ////////////////

      //! @brief picks a random template of the appropriate size
      //! @param HELICES number of helices to be in the template
      //! @param STRANDS number of strands to be in the template
      //! @return a random template of the appropriate size
      static const FoldTemplate &GetRandomTemplate( const size_t HELICES, const size_t STRANDS);

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

      //! @brief returns a random fold template (that has geometries close in size to the passed sses)
      //!        generated from a larger template
      //! @param SSES sses used to chose subtemplate based on length
      //! @param SSE_GEOMETRY_COMPARE comparison method
      //! @return a random fold template generated from a larger template
      static FoldTemplate GetRandomSubTemplate
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

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief singleton function to return fold templates previously generated and stored in a file
      //! @return fold templates previously generated and stored in a file
      static const storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      > &GetTemplates();

      //! @brief singleton function to return fold templates previously generated and stored in a file
      //! @param FILENAMES set of filenames to be read in
      //! @return fold templates previously generated and stored in a file
      static storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      > ReadTemplates( const storage::Set< std::string> &FILENAMES);

      //! @brief gets all fold templates from the requested composition up to a constant size
      //! @param HELICES helical SSEs to be used to find a suitable template
      //! @param STRANDS strand SSEs to be used to find a suitable template
      //! @param SSE_GEOMETRY_COMPARE comparison method
      //! @return all fold templates from the requested composition up to a constant size
      static util::ShPtrVector< FoldTemplate> GetLargeTemplates
      (
        const util::SiPtrList< const SSE> &HELICES,
        const util::SiPtrList< const SSE> &STRANDS,
        const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
      );

      //! @brief default input filename
      //! @param TEMPLATE_TYPE type of templates to use
      //! @return default input filename
      static const std::string &GetInputFilename( const std::string &TEMPLATE_TYPE);

      //! @brief return set of excluded pdb ids
      //! @return set of excluded pdb ids
      static storage::Set< std::string> &GetExcludedPdbs();

    }; // class FoldTemplateHandler

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_FOLD_TEMPLATE_HANDLER_H_
