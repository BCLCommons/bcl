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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_PICKERS_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_PICKERS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packer_interface.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingPickers
    //! @brief class enumerates different ways to calculate best packing between sub-geometries for a given pair of geometries
    //! @details This enumerator, enumerates different SSEGeometryPackerInterface derived classes with a return type
    //! of a single SSEGeometryPacking. For a given pair of GEOMETRY_A (composed of M sub-geometries) and
    //! GEOMETRY_B (composed of N sub-geometries), the return type is a single best SSEGeometryPacking
    //! Each enumerator is repeated with a NoCache option, which makes sure even if the caching flag is provided
    //! that instance won't use caching.
    //!
    //! @see @link example_assemble_sse_geometry_packing_pickers.cpp @endlink
    //! @author karakam
    //! @date Jun 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingPickers :
      public util::Enumerate< util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >, SSEGeometryPackingPickers>
    {

      friend class util::Enumerate< util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >, SSEGeometryPackingPickers>;

    public:

    //////////
    // data //
    //////////

      const SSEGeometryPackingPicker e_BestInteractionWeight;         //!< best interaction
      const SSEGeometryPackingPicker e_BestInteractionWeight_NoCache; //!< best interaction no cache
      const SSEGeometryPackingPicker e_BestStrandPairing;             //!< best strand
      const SSEGeometryPackingPicker e_BestStrandPairing_NoCache;     //!< best strand no cache

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief private constructor for constructing all SSEGeometryPackingPickers
      SSEGeometryPackingPickers();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class SSEGeometryPackingPickers

    //! @brief access function to a static instance of SSEGeometryPackingPickers
    //! @return reference to a static instance of SSEGeometryPackingPickers
    BCL_API
    const SSEGeometryPackingPickers &GetSSEGeometryPackingPickers();

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< assemble::SSEGeometryPackerInterface< assemble::SSEGeometryPacking> >, assemble::SSEGeometryPackingPickers>;

  } // namespace util
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_PICKERS_H_ 
