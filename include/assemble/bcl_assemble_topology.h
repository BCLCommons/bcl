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

#ifndef BCL_ASSEMBLE_TOPOLOGY_H_
#define BCL_ASSEMBLE_TOPOLOGY_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packing.h"
#include "coord/bcl_coord_orientation_interface.h"
#include "graph/bcl_graph_graph_with_data.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Topology
    //! @brief class provides a graph class for representing topologies of ProteinModels and similar classes
    //! @details Topology class provides the functionality to represent the connectivity thus the topology of SSEs in a
    //! Protein Model or SSEGeometries in a FoldTemplate or similar classes. It can be used with any class derived
    //! from SSEGeometryInterface.
    //!
    //! @see @link example_assemble_topology.cpp @endlink
    //! @author karakam, weinerbe
    //! @date Jul 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Topology :
      public coord::OrientationInterface
    {

    public:

    ///////////
    // enums //
    ///////////

      //! enumerator for topology type
      enum Type
      {
        e_Sheet,
        e_BetaBarrel,
        e_HelixBundle,
        s_NumberTypes
      };

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

      //! Typedef for graph type
      typedef graph::GraphWithData< util::SiPtr< const SSEGeometryInterface>, SSEGeometryPacking> GraphType;

    private:

    //////////
    // data //
    //////////

      //! geometry type
      TypeEnum m_Type;

      //! graph that stores the SSEs and corresponding connections
      GraphType m_Graph;

      //! storage for all the element SSEGeometries
      util::SiPtrVector< const SSEGeometryInterface> m_Elements;

      //! transformation matrix that defines the geometry of the topology
      math::TransformationMatrix3D m_Geometry;

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
      Topology();

      //! @brief constructor from vector of SSE Geometries
      //! @param GEOMETRY_VECTOR SiPtrVector of Geometries
      Topology( const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR);

      //! @brief constructor from an orientation
      //! @param ORIENTATION math::TransformationMatrix3D that defines the orientation
      Topology( const math::TransformationMatrix3D &ORIENTATION);

      //! @brief constructor from a graph
      //! @param GRAPH graph to be used
      Topology( const GraphType &GRAPH);

      //! @brief Clone function
      //! @return pointer to new Topology
      Topology *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get topology type
      //! @return topology type
      Type GetType() const
      {
        return m_Type;
      }

      //! @brief set type to given TOPOLOGY_TYPE
      //! @param TOPOLOGY_TYPE Topology::Type to be used
      void SetType( const Type TOPOLOGY_TYPE)
      {
        m_Type = TOPOLOGY_TYPE;
      }

      //! @brief set the graph
      //! @param GRAPH Graph to be copied
      void SetGraph( const GraphType &GRAPH);

      //! @brief return the graph that represents the connectivity of the topology
      //! @return graph that represents the connectivity of the topology
      const GraphType &GetGraph() const
      {
        return m_Graph;
      }

      //! @brief get all the represented geometries
      //! @return vector of all geometries
      const util::SiPtrVector< const SSEGeometryInterface> &GetElements() const
      {
        return m_Elements;
      }

      //! @brief set get all the represented geometries
      //! @param ELEMENTS vector of all geometries
      void SetElements( const util::SiPtrVector< const SSEGeometryInterface> &ELEMENTS)
      {
        m_Elements = ELEMENTS;
      }

      //! @brief get packing object for two SSEGeometries
      //! @param SP_GEOMETRY_A SiPtr to first SSEGeometry of interest
      //! @param SP_GEOMETRY_B SiPtr to second SSEGeometry of interest
      //! @return packing between SSEGeometry behind SP_GEOMETRY_A and SP_GEOMETRY_B
      const SSEGeometryPacking &GetPackingForSSEGeometryPair
      (
        const util::SiPtr< const SSEGeometryInterface> &SP_GEOMETRY_A,
        const util::SiPtr< const SSEGeometryInterface> &SP_GEOMETRY_B
      ) const;

      //! @brief returns the geometric center of the Sheet
      //! @return the geometric center of the Sheet
      linal::Vector3D GetCenter() const;

      //! @brief return the orientation of the Sheet
      //! @return orientation
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! @brief return the orientation and Position as TransformationMatrix3D
      //! @return TransformationMatrix3D that defines orientation and position
      const math::TransformationMatrix3D GetOrientation() const;

      //! @brief set the orientation of this topology
      //! @param ORIENTATION Orientation to be used
      void SetOrientation( const math::TransformationMatrix3D &ORIENTATION)
      {
        m_Geometry = ORIENTATION;
      }

      //! @brief initialize the orientation using type information
      void SetOrientationFromType();

      //! @brief get the identification of this topology with detailed connectivity map information
      //! @return the identification of this topology  with detailed connectivity map information
      std::string GetIdentification() const;

      //! @brief get Identification of this topology in a ordered fashion
      //! @return string with identification
      std::string GetOrderedIdentification() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief finds unconnected sub-topologies within this topology and returns them in a vector
      //! @return vector of sub-topologies
      util::ShPtrVector< Topology> GetSubTopologies() const;

      //! @brief creates a subtopology composed of the given SSEGeometryInterfaces
      //! @param ELEMENTS vector of elements that should be added to the subtopology
      //! @return ShPtr to a new SubTopology
      util::ShPtr< Topology> GetSubTopology( const util::SiPtrVector< const SSEGeometryInterface> &ELEMENTS) const;

      //! @brief orders the elements in the topology
      //! @return whether ordering succeeded
      bool OrderElements();

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

    public:

      //! @brief write graphviz script
      //! @param OSTREAM the stream the script is written to
      //! @return the stream the script was written to
      std::ostream &WriteGraphVizScript
      (
        std::ostream &OSTREAM,
        const math::SumFunction< SSEGeometryPacking, double> &COLOR_EDGE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief builds a topology graph from the given vector of SSEGeometryInterface derived object
      //! @param SSE_GEOMETRY_VECTOR SiPtrVector of SSEGeometryInterface derived objects
      //! @param PACKER SSEGeometryPackingPicker to be used to calculate the edges
      //! @param PACKING_CRITERIA criteria to decide whether the SSEPacking can be considered a connected edge
      //! @return topology graph constructed
      static Topology::GraphType BuildTopologyGraphFromGeometries
      (
        const util::SiPtrVector< const SSEGeometryInterface> &SSE_GEOMETRY_VECTOR,
        const SSEGeometryPackingPicker &PACKER,
        const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &PACKING_CRITERIA
      );

      //! @brief builds a topology from the given vector of SSEGeometryInterface derived object, packer and criteria
      //! @param SSE_GEOMETRY_VECTOR SiPtrVector of SSEGeometryInterface derived objects
      //! @param PACKER SSEGeometryPackingPicker to be used to calculate the edges
      //! @param PACKING_CRITERIA criteria to decide whether the SSEPacking can be considered a connected edge
      static util::ShPtr< Topology> BuildTopologyFromGeometries
      (
        const util::SiPtrVector< const SSEGeometryInterface> &SSE_GEOMETRY_VECTOR,
        const SSEGeometryPackingPicker &PACKER,
        const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &PACKING_CRITERIA
      );

      //! @brief builds a topology from the given protein model based on its amino acid contacts
      //! @param PROTEIN the protein model to use
      //! @return the topology
      static util::ShPtr< Topology> BuildTopologyFromContacts( const ProteinModel &PROTEIN);

      //! @brief collects the SiPtrVector of SSEGeometryInterfaces from the Graph vertex data
      //! @param GRAPH reference to the graph from which the elements vector will be collected
      //! @return the SiPtrVector of SSEGeometryInterfaces from the Graph vertex data
      static util::SiPtrVector< const SSEGeometryInterface> GetElementsVectorFromGraph
      (
        const GraphType &GRAPH
      );

      //! @brief static function that calculates the center of a vector SSEGeometries
      //! @param GEOMETRY_VECTOR SiPtrVector of SSE Geometries
      //! @return linal::Vector3D that represents the center of the given geometries
      linal::Vector3D CalculateCenterOfGeometries
      (
        const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
      );

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief expands the given subgraph starting from origin_vertex and updates the remaining vertices
      //! @param SP_ORIGIN_VERTEX SiPtr to the origin vertex from the complete graph
      //! @param SUB_GRAPH non-const reference to subgraph to be expanded
      //! @param REMAINING_VERTICES non-const reference to container of remaining vertices
      void ExpandSubGraph
      (
        const util::SiPtr< const GraphType::VertexType> &SP_ORIGIN_VERTEX,
        GraphType &SUB_GRAPH,
        GraphType::VertexContainerType &REMAINING_VERTICES
      ) const;

      //! @brief finds and removes the vertex from a vector of vertices
      //! @param VERTEX vertex to be searched
      //! @param VERTEX_VECTOR Vertex vector to be searched
      static void RemoveVertexFromVector
      (
        const GraphType::VertexType &VERTEX,
        GraphType::VertexContainerType &VERTEX_VECTOR
      );

    }; // class Topology

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_TOPOLOGY_H_
