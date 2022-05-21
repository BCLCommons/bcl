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

#ifndef BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_SHEET_H_
#define BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_SHEET_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_collector_topology_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorTopologySheet
    //! @brief Collector class for identifying beta-sheet topology for a vector of geometries
    //! @details This class uses the given list of SSEGeometryInterface derived classes to identify beta-sheet
    //! topologies, where the order, connectivity information and other relavant information is stored.
    //!
    //! @see @link example_assemble_collector_topology_sheet.cpp @endlink
    //! @author karakam
    //! @date Sep 10, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorTopologySheet :
      public CollectorTopologyInterface
    {

    public:

    //////////
    // data //
    //////////

      //! maximum distance to identify as strand strand pairing
      static const double s_MaximumDistance;

      //! minimum relative position weight needed to be identified as strand strand pairing
      static const double s_MinimumStrandStrandPairingWeight;

      //! @brief return packing criteria function
      //! @return packing criteria function
      static
      const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
      &GetPackingCriteria();

      //! @brief return packing comparison function
      //! @return packing comparison function
      static
      const util::ShPtr< math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> >
      &GetPackingComparison();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorTopologySheet();

      //! @brief Clone function
      //! @return pointer to new CollectorTopologySheet
      CollectorTopologySheet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect returns all beta-sheet topologies in the given geometry vector
      //! @param GEOMETRY_VECTOR Vector of geometries of interest
      //! @return returns ShPtrVector of beta-sheets found in the given geometry vector
      util::ShPtrVector< Topology> Collect
      (
        const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
      ) const;

      //! @calculate the topology, which may contain isolated sub-topologies
      //! @param GEOMETRY_VECTOR Vector of geometries of interest
      //! @return calculated topology
      Topology CalculateTopology
      (
        const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorTopologySheet

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_TOPOLOGY_SHEET_H_ 
