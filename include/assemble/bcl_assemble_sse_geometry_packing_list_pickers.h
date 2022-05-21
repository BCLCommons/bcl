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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_LIST_PICKERS_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_LIST_PICKERS_H_

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
    //! @class SSEGeometryPackingListPickers
    //! @brief class enumerates different ways to calculate best packing lists between sub-geometries for a given pair of geometries
    //! @details This enumerator, enumerates different SSEGeometryPackerInterface derived classes with a return type
    //! of list of SSEGeometryPackings. For a given pair of GEOMETRY_A (composed of M sub-geometries) and
    //! GEOMETRY_B (composed of N sub-geometries), the return type is a list with M packings, where each packing
    //! was used after comparison of N packings.
    //! Each enumerator is repeated with a NoCache option, which makes sure even if the caching flag is provided
    //! that instance won't use caching.
    //!
    //! @see @link example_assemble_sse_geometry_packing_list_pickers.cpp @endlink
    //! @author karakam
    //! @date Jun 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingListPickers :
      public util::Enumerate
      <
        util::ShPtr< SSEGeometryPackerInterface< storage::List< SSEGeometryPacking> > >,
        SSEGeometryPackingListPickers
      >
    {

      friend class
      util::Enumerate
      <
        util::ShPtr< SSEGeometryPackerInterface< storage::List< SSEGeometryPacking> > >,
        SSEGeometryPackingListPickers
      >;

    public:

    //////////
    // data //
    //////////

      const SSEGeometryPackingListPicker e_BestInteractionWeight;  //!< best interaction
      const SSEGeometryPackingListPicker e_BestDistance;           //!< best distance
      const SSEGeometryPackingListPicker e_SSEClash;               //!< SSE clash
      const SSEGeometryPackingListPicker e_SSEClash_NoCache;       //!< SSE clash no cache

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief private constructor for constructing all SSEGeometryPackingListPickers
      SSEGeometryPackingListPickers();

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

    }; // class SSEGeometryFragmentsPacker

    //! @brief access static instance of SSEGeometryPackingListPickers enums
    //! @return static instance of SSEGeometryPackingListPickers enums
    BCL_API
    const SSEGeometryPackingListPickers &GetSSEGeometryPackingListPickers();

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< assemble::SSEGeometryPackerInterface< storage::List< assemble::SSEGeometryPacking> > >, assemble::SSEGeometryPackingListPickers>;

  } // namespace util
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_LIST_PICKERS_H_
