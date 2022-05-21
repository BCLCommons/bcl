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
#include "assemble/bcl_assemble_collector_topology_sheet.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packing_compare.h"
#include "assemble/bcl_assemble_sse_geometry_packing_criteria.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "assemble/bcl_assemble_topology.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! maximum distance to identify as strand strand pairing
    const double CollectorTopologySheet::s_MaximumDistance( 7.0);

    //! minimum relative position weight needed to be identified as strand strand pairing
    const double CollectorTopologySheet::s_MinimumStrandStrandPairingWeight( 0.5);

    //! @brief return packing criteria function
    //! @return packing criteria function
    const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &
    CollectorTopologySheet::GetPackingCriteria()
    {
      // construct the vector of criteria vector
      static const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > s_packing_criteria
      (
        new SSEGeometryPackingCriteriaCombine
        (
          util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >::Create
          (
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaContactType( contact::GetTypes().STRAND_STRAND)
            ),
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaDistance( s_MaximumDistance, math::Comparisons< double>::GetEnums().e_LessEqual)
            ),
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaStrandWeight
              (
                s_MinimumStrandStrandPairingWeight, math::Comparisons< double>::GetEnums().e_GreaterEqual
              )
            )
          )
        )
      );

      // end
      return s_packing_criteria;
    }

    //! @brief return packing comparison function
    //! @return packing comparison function
    const util::ShPtr< math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> > &
    CollectorTopologySheet::GetPackingComparison()
    {
      // initialize the static instance of the packer
      static const util::ShPtr< math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> >
      s_packing_comparison
      (
        new SSEGeometryPackingCompareInteractionWeight()
      );

      // end
      return s_packing_comparison;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorTopologySheet::s_Instance
    (
      util::Enumerated< CollectorTopologyInterface>::AddInstance( new CollectorTopologySheet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorTopologySheet::CollectorTopologySheet()
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorSheetTopology
    CollectorTopologySheet *CollectorTopologySheet::Clone() const
    {
      return new CollectorTopologySheet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorTopologySheet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorTopologySheet::GetAlias() const
    {
      static const std::string s_alias( "CollectorTopologySheet");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorTopologySheet::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collect sheet topologies");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns all beta-sheet topologies in the given geometry vector
    //! @param GEOMETRY_VECTOR Vector of geometries of interest
    //! @return returns ShPtrVector of beta-sheets found in the given geometry vector
    util::ShPtrVector< Topology> CollectorTopologySheet::Collect
    (
      const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
    ) const
    {
      // construct sheet vector
      util::ShPtrVector< Topology> topology_vector;

      // get full topology
      const Topology full_topology( CalculateTopology( GEOMETRY_VECTOR));

      // return empty vector if no strands
      if( full_topology.GetElements().IsEmpty())
      {
        return topology_vector;
      }

      // now call the topology function
      topology_vector = full_topology.GetSubTopologies();

      // for each topology
      for
      (
        util::ShPtrVector< Topology>::iterator
          topology_itr( topology_vector.Begin()), topology_itr_end( topology_vector.End());
        topology_itr != topology_itr_end; ++topology_itr
      )
      {
        // create reference on this topology
        Topology &this_topology( **topology_itr);

        // try to order the geometries and store the return
        const bool is_valid_geometry( this_topology.OrderElements());

        // if not successful then skip this topology
        if( !is_valid_geometry)
        {
          continue;
        }

        // set the type
        this_topology.SetType
        (
          this_topology.GetElements().GetSize() > 2 &&
          this_topology.GetGraph().AreVerticesConnected
          (
            this_topology.GetElements().FirstElement(), this_topology.GetElements().LastElement()
          ) ?
            Topology::e_BetaBarrel : Topology::e_Sheet
        );

        // update the geometry
        this_topology.SetOrientationFromType();
      }

      // end
      return topology_vector;
    }

    //! @calculate the topology, which may contain isolated sub-topologies
    //! @param GEOMETRY_VECTOR Vector of geometries of interest
    //! @return calculated topology
    Topology CollectorTopologySheet::CalculateTopology
    (
      const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
    ) const
    {
      // initialize strand vector
      util::SiPtrVector< const SSEGeometryInterface> strand_vector;

      // iterate over given geometries
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          geometry_itr( GEOMETRY_VECTOR.Begin()), geometry_itr_end( GEOMETRY_VECTOR.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // if of type strand
        if( ( *geometry_itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          // add it to the strand vector
          strand_vector.PushBack( *geometry_itr);
        }
      }

      // construct the topology and return
      return *Topology::BuildTopologyFromGeometries
      (
        strand_vector, GetSSEGeometryPackingPickers().e_BestStrandPairing, *GetPackingCriteria()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
