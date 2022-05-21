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

#ifndef BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_INTERFACE_H_
#define BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorTopologyInterface
    //! @brief Topology collector interface
    //! @details Interface class that determines topologies of SSEs.  In addition to collecting all unconnected
    //!          topologies, derived classes can also calculate the combined topology from the given SSEs.
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Oct 24, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorTopologyInterface :
      public find::CollectorInterface< util::ShPtrVector< Topology>, util::SiPtrVector< const SSEGeometryInterface> >
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief Clone function.
      //! @return pointer to a new CollectorTopologyInterface
      virtual CollectorTopologyInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @calculate the topology, which may contain isolated sub-topologies
      //! @param GEOMETRY_VECTOR Vector of geometries of interest
      //! @return calculated topology
      virtual Topology CalculateTopology
      (
        const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
      ) const = 0;

    }; // class CollectorTopologyInterface

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_INTERFACE_H_
