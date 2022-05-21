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
#include "assemble/bcl_assemble_collector_topology_combined.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_sheet.h"
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

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorTopologyCombined::s_Instance
    (
      util::Enumerated< CollectorTopologyInterface>::AddInstance( new CollectorTopologyCombined())
    );

    //! @brief return packing criteria function
    //! @return packing criteria function
    const util::ShPtr
    <
      math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    > &CollectorTopologyCombined::GetPackingCriteria()
    {
      // initialize static criteria function
      static const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > s_default_criteria
      (
        new SSEGeometryPackingCriteriaCombine
        (
          util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >::Create
          (
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaDistancePerType()
            ),
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaInteractionWeight
              (
                0.5, math::Comparisons< double>::GetEnums().e_GreaterEqual
              )
            )
          )
        )
      );

      // end
      return s_default_criteria;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param INCLUDE_UNPACKED_GEOMETRIES whether to include unpacked geometries in the topologies
    CollectorTopologyCombined::CollectorTopologyCombined( const bool INCLUDE_UNPACKED_GEOMETRIES) :
      m_IncludeUnpackedGeometries( INCLUDE_UNPACKED_GEOMETRIES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorCombinedTopology
    CollectorTopologyCombined *CollectorTopologyCombined::Clone() const
    {
      return new CollectorTopologyCombined( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorTopologyCombined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorTopologyCombined::GetAlias() const
    {
      static const std::string s_alias( "CollectorTopologyCombined");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorTopologyCombined::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Determines complete topology of a protein");
      serializer.AddInitializer
      (
        "include unpacked",
        "include unpacked geometries in the topologies",
        io::Serialization::GetAgent( &m_IncludeUnpackedGeometries),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns all topologies in the given geometry vector
    //! @param GEOMETRY_VECTOR Vector of geometries of interest
    //! @return returns ShPtrVector of topologies found in the given geometry vector
    util::ShPtrVector< Topology> CollectorTopologyCombined::Collect
    (
      const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
    ) const
    {
      return CalculateTopology( GEOMETRY_VECTOR).GetSubTopologies();
    }

    //! @calculate the topology, which may contain isolated sub-topologies
    //! @param GEOMETRY_VECTOR Vector of geometries of interest
    //! @return calculated topology
    Topology CollectorTopologyCombined::CalculateTopology
    (
      const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
    ) const
    {
      // get the helix and strand geometries
      util::SiPtrVector< const SSEGeometryInterface> helices;
      util::SiPtrVector< const SSEGeometryInterface> strands;
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator itr( GEOMETRY_VECTOR.Begin()),
          itr_end( GEOMETRY_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        if( ( *itr)->GetType() == biol::GetSSTypes().HELIX)
        {
          helices.PushBack( *itr);
        }
        else if( ( *itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          strands.PushBack( *itr);
        }
      }

      // build a graph using sheet packing pickers
      Topology::GraphType graph_sheet
      (
        Topology::BuildTopologyGraphFromGeometries
        (
          strands,
          GetSSEGeometryPackingPickers().e_BestStrandPairing,
          *CollectorTopologySheet::GetPackingCriteria()
        )
      );

      // if removing unpacked geometries
      if( !m_IncludeUnpackedGeometries && strands.GetSize() > 1)
      {
        util::SiPtrVector< const SSEGeometryInterface> filtered_strands;

        // iterate over the strands
        for
        (
          util::SiPtrVector< const SSEGeometryInterface>::const_iterator strand_itr( strands.Begin()),
            strand_itr_end( strands.End());
          strand_itr != strand_itr_end; ++strand_itr
        )
        {
          // find the corresponding vertex
          const util::ShPtr< Topology::GraphType::VertexType> vertex
          (
            graph_sheet.FindVertex( *strand_itr)
          );

          // if the vertex has no edges
          if( vertex->GetDegree() == 0)
          {
            // delete it from the graph
            graph_sheet.DeleteVertex( vertex);
          }
          // if it has edges
          else
          {
            // add to filtered vector
            filtered_strands.PushBack( *strand_itr);
          }
        }

        strands = filtered_strands;
      }

      // combine the geometries
      util::SiPtrVector< const SSEGeometryInterface> combined_geometries( strands);
      combined_geometries.Append( helices);

      // build a graph from Fold template setting
      Topology::GraphType graph_fold
      (
        Topology::BuildTopologyGraphFromGeometries
        (
          combined_geometries,
          GetSSEGeometryPackingPickers().e_BestInteractionWeight,
          *GetPackingCriteria()
        )
      );

      // if there are no edges then return graph_fold
      if( graph_sheet.GetNumberEdges() == 0)
      {
        return Topology( graph_fold);
      }

      // now we need to overwrite the edges in the graph_fold for which there is an edge in the sheet graph
      // for this, we need to iterate over all the edges in the sheet
      for
      (
        Topology::GraphType::VertexContainerType::const_iterator
          vertex_itr( graph_sheet.GetVertices().Begin()), vertex_itr_end( graph_sheet.GetVertices().End());
        vertex_itr != vertex_itr_end; ++vertex_itr
      )
      {
        // create ref on vertex
        const Topology::GraphType::VertexType &vertex_a( **vertex_itr);

        // iterate over the edges for this vertex
        for
        (
          Topology::GraphType::EdgeContainerType::const_iterator
            edge_itr( vertex_a.GetEdges().Begin()), edge_itr_end( vertex_a.GetEdges().End());
          edge_itr != edge_itr_end; ++edge_itr
        )
        {
          // create ref on the vertex_b
          const Topology::GraphType::VertexType &vertex_b( *edge_itr->GetTarget());

          // now delete the edge in the graph_fold, if the edge does not exists it won't do anything
          graph_fold.DeleteEdge( vertex_a.GetData(), vertex_b.GetData());

          // now insert the sheet packing information with a new edge
          graph_fold.AddEdge( vertex_a.GetData(), vertex_b.GetData(), edge_itr->GetData());
        }
      }

      // if removing unpacked geometries
      if( !m_IncludeUnpackedGeometries)
      {
        graph_fold.DeleteUnconnectedVertices();
      }

      // end
      return Topology( graph_fold);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
