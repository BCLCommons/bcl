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
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_sheet.h"
#include "assemble/bcl_assemble_sse_geometry_packer_best_fragment_pair.h"
#include "assemble/bcl_assemble_sse_geometry_packing_compare.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief private constructor for constructing all SSEGeometryPackingPickers
    SSEGeometryPackingPickers::SSEGeometryPackingPickers() :
        e_BestInteractionWeight
        (
          AddEnum
          (
            "BestInteractionWeight",
            SSEGeometryPackerInterface< SSEGeometryPacking>::WrapPacker
            (
              util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >
              (
                new SSEGeometryPackerBestFragmentPair
                (
                  GetSSEGeometryPackers().e_Default,
                  SSEGeometryPackingCompareInteractionWeight()
                )
              ), true, false
            )
          )
        ),
        e_BestInteractionWeight_NoCache
        (
          AddEnum
          (
            "BestInteractionWeight_NoCache",
            SSEGeometryPackerInterface< SSEGeometryPacking>::WrapPacker
            (
              util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >
              (
                new SSEGeometryPackerBestFragmentPair
                (
                  GetSSEGeometryPackers().e_Default_NoCache,
                  SSEGeometryPackingCompareInteractionWeight()
                )
              ), false, false
            )
          )
        ),
        e_BestStrandPairing
        (
          AddEnum
          (
            "BestStrandPairing",
            SSEGeometryPackerInterface< SSEGeometryPacking>::WrapPacker
            (
              util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >
              (
                new SSEGeometryPackerBestFragmentPair
                (
                  GetSSEGeometryPackers().e_Default,
                  *CollectorTopologySheet::GetPackingCriteria(),
                  *CollectorTopologySheet::GetPackingComparison()
                )
              ),
              true, false
            )
          )
        ),
        e_BestStrandPairing_NoCache
        (
          AddEnum
          (
            "BestStrandPairing_NoCache",
            SSEGeometryPackerInterface< SSEGeometryPacking>::WrapPacker
            (
              util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >
              (
                new SSEGeometryPackerBestFragmentPair
                (
                  GetSSEGeometryPackers().e_Default_NoCache,
                  *CollectorTopologySheet::GetPackingCriteria(),
                  *CollectorTopologySheet::GetPackingComparison()
                )
              ), false, false
            )
          )
        )
      {
      }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPackingPickers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access function to a static instance of SSEGeometryPackingPickers
    //! @return reference to a static instance of SSEGeometryPackingPickers
    const SSEGeometryPackingPickers &GetSSEGeometryPackingPickers()
    {
      return SSEGeometryPackingPickers::GetEnums();
    }

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< assemble::SSEGeometryPackerInterface< assemble::SSEGeometryPacking> >, assemble::SSEGeometryPackingPickers>;

  } // namespace util
} // namespace bcl
