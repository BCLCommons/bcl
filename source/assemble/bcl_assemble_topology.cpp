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
#include "assemble/bcl_assemble_topology.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "contact/bcl_contact_map.h"
#include "math/bcl_math_running_min_max.h"
#include "math/bcl_math_sum_function.h"
#include "util/bcl_util_color_gradient.h"
#include "util/bcl_util_colors.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  ///////////
  // enums //
  ///////////

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &Topology::GetTypeName( const Type &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "beta-sheet",
        "beta-barrel",
        "helix-bundle",
        "undefined"
      };
      return s_descriptors[ size_t( TYPE)];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Topology::s_Instance
    (
      GetObjectInstances().AddInstance( new Topology())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Topology::Topology() :
      m_Type( s_NumberTypes),
      m_Graph(),
      m_Elements(),
      m_Geometry()
    {
    }

    //! @brief constructor from vector of SSE Geometries
    //! @param GEOMETRY_VECTOR SiPtrVector of Geometries
    Topology::Topology( const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR) :
      m_Type( s_NumberTypes),
      m_Graph(),
      m_Elements( GEOMETRY_VECTOR),
      m_Geometry( CalculateCenterOfGeometries( GEOMETRY_VECTOR))
    {
    }

    //! @brief constructor from an orientation
    //! @param ORIENTATION math::TransformationMatrix3D that defines the orientation
    Topology::Topology( const math::TransformationMatrix3D &ORIENTATION) :
      m_Type( s_NumberTypes),
      m_Graph(),
      m_Elements(),
      m_Geometry( ORIENTATION)
    {
    }

    //! @brief private constructor from data members
    //! @param GRAPH graph to be used
    Topology::Topology
    (
      const GraphType &GRAPH
    ) :
      m_Type( s_NumberTypes),
      m_Graph( GRAPH),
      m_Elements( GetElementsVectorFromGraph( GRAPH)),
      m_Geometry()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Topology
    Topology *Topology::Clone() const
    {
      return new Topology( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Topology::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the graph
    //! @param GRAPH Graph to be copied
    void Topology::SetGraph( const GraphType &GRAPH)
    {
      m_Graph = GRAPH;
      m_Elements = GetElementsVectorFromGraph( m_Graph);
    }

    //! @brief get packing object for two SSEGeometries
    //! @param SP_GEOMETRY_A SiPtr to first SSEGeometry of interest
    //! @param SP_GEOMETRY_B SiPtr to second SSEGeometry of interest
    //! @return packing between SSEGeometry behind SP_GEOMETRY_A and SP_GEOMETRY_B
    const SSEGeometryPacking &Topology::GetPackingForSSEGeometryPair
    (
      const util::SiPtr< const SSEGeometryInterface> &SP_GEOMETRY_A,
      const util::SiPtr< const SSEGeometryInterface> &SP_GEOMETRY_B
    ) const
    {
      // initialize a static undefined SSEPacking
      static const SSEGeometryPacking s_undefined_packing;

      // if any of the pointers are not defined
      if( !SP_GEOMETRY_A.IsDefined() || !SP_GEOMETRY_B.IsDefined())
      {
        return s_undefined_packing;
      }

      // find the corresponding vertices
      const util::ShPtr< GraphType::VertexType> vertex_a( m_Graph.FindVertex( SP_GEOMETRY_A));
      const util::ShPtr< GraphType::VertexType> vertex_b( m_Graph.FindVertex( SP_GEOMETRY_B));

      // if any of the vertices can't be found
      if( !vertex_a.IsDefined() || !vertex_b.IsDefined())
      {
        return s_undefined_packing;
      }

      // now find the corresponding edge
      const GraphType::EdgeType &edge( vertex_a->FindEdge( *vertex_b));

      // if the edge can't be found
      if( !edge.IsDefined())
      {
        return s_undefined_packing;
      }

      // return the data
      return edge.GetData();
    }

    //! @brief returns the geometric center of the Sheet
    //! @return the geometric center of the Sheet
    linal::Vector3D Topology::GetCenter() const
    {
      return m_Geometry.GetOrigin();
    }

    //! @brief return the orientation of the Sheet
    //! @return orientation
    linal::Vector3D Topology::GetAxis( const coord::Axis &AXIS) const
    {
      return m_Geometry.GetAxis( AXIS);
    }

    //! @brief return the orientation and Position as TransformationMatrix3D
    //! @return TransformationMatrix3D that defines orientation and position
    const math::TransformationMatrix3D Topology::GetOrientation() const
    {
      return m_Geometry;
    }

    //! @brief initialize the orientation using type information
    void Topology::SetOrientationFromType()
    {
      // if type is sheet or beta barrel
      if( m_Type == e_Sheet || m_Type == e_BetaBarrel)
      {
        // store number of strands
        const size_t nr_strands( m_Elements.GetSize());

        // decide on the center strand index
        const size_t center_strand_index( ( nr_strands - 1) / 2);

        // initialize transformation to the center strand
        math::TransformationMatrix3D transformation( m_Elements( center_strand_index)->GetOrientation());

        // if there are even number of strands
        // then calculate the translation that will move the origin to the middle point of the two center strands
        if( nr_strands % 2 == 0)
        {
          // get the translation from selected middle strand to the other one
          linal::Vector3D translation
          (
            m_Elements( center_strand_index + 1)->GetCenter() - m_Elements( center_strand_index)->GetCenter()
          );

          // update the geometry with half of this translation
          transformation( translation / 2);
        }

        // make sure the y axis points from left to right according to order vector
        // first calculate the axis from left to right
        linal::Vector3D left_to_right_vector
        (
           m_Elements.LastElement()->GetCenter() - m_Elements.FirstElement()->GetCenter()
        );

        // store the y axis
        linal::Vector3D y_axis( transformation.GetAxis( coord::GetAxes().e_Y));

        // calculate the angle between the y axis and the left_to_right_vector
        const double angle( linal::ProjAngle( left_to_right_vector, y_axis));

        // if they are pointing in opposite directions, thus the angle is larger than 90 degrees
        if( angle >= math::g_Pi / 2)
        {
          // set the y axis to the opposite
          transformation.SetAxis( coord::GetAxes().e_Y, -y_axis);
        }

        // set the geometry
        m_Geometry = transformation;
      }
      else
      {
        BCL_MessageStd( "Topology can set orientation only for sheet or beta-barrels at this point");
      }

    }

    //! @brief get the identification of this topology with detailed connectivity map information
    //! @return the identification of this topology  with detailed connectivity map information
    std::string Topology::GetIdentification() const
    {
      // initialize identification to show number of strands and the topology type
      std::string identification( "number_elements: ");
      identification += util::Format()( m_Elements.GetSize()) +
                        " type: " + util::Format()( m_Type) + "\n";
      identification += "Connectivity Map\n";

      // iterate over strands in the given order
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          element_itr( m_Elements.Begin()), element_itr_end( m_Elements.End());
        element_itr != element_itr_end; ++element_itr
      )
      {
        // add the identification of this strand
        identification += ( *element_itr)->GetIdentification() + "\n";

        // get the corresponding vertex
        util::ShPtr< GraphType::VertexType> vertex( m_Graph.FindVertex( *element_itr));

        // iterate over the edges this vertex is connected to
        for
        (
          GraphType::EdgeContainerType::const_iterator
            edge_itr( vertex->GetEdges().Begin()), edge_itr_end( vertex->GetEdges().End());
          edge_itr != edge_itr_end; ++edge_itr
        )
        {
          // add the additional information
          identification +=
            "\t"               + edge_itr->GetTarget()->GetData()->GetIdentification() +
            "\n\t" + edge_itr->GetData().GetIdentification() + "\n";
        }
      }
      identification += "\n";

      // end
      return identification;
    }

    //! @brief get Identification of this topology in a ordered fashion
    //! @return string with identification
    std::string Topology::GetOrderedIdentification() const
    {
      // initialize identification to show number of strands and the topology type
      std::string identification;
      identification +=
        "number elements " + util::Format()( m_Elements.GetSize()) +
         "\ttype: " + util::Format()( m_Type) + "\n";

      // initialize element ctr
      size_t element_ctr( 0);

      // iterate over strands in the given order
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          element_itr( m_Elements.Begin()), element_itr_end( m_Elements.End());
        element_itr != element_itr_end; ++element_itr
      )
      {
        // increment strand ctr
        ++element_ctr;

        // add identification of the strand
        identification += ( *element_itr)->GetIdentification();

        // if not the first SSE
        if( element_itr != m_Elements.Begin())
        {
          // add the orientation with respect to the previous strand
          identification += "\t";

          // get the packing
          const SSEGeometryPacking &this_packing( GetPackingForSSEGeometryPair( *( element_itr - 1), *element_itr));
          identification += this_packing.GetIdentification();
        }

        identification += "\n";
      }
      // end
      return identification;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief finds unconnected sub-topologies within this topology and returns them in a vector
    //! @return vector of sub-topologies
    util::ShPtrVector< Topology> Topology::GetSubTopologies() const
    {
      // instantiate the vector of sub-topologies
      util::ShPtrVector< Topology> sub_topologies_vector;

      // make a copy of the vertex data
      GraphType::VertexContainerType remaining_vertices( m_Graph.GetVertices());

      BCL_MessageDbg( "Creating sub-topologies");

      // while there are vertices remaining that are not part of any subgraph
      while( !remaining_vertices.IsEmpty())
      {
        // construct a vector of vertices that are reachable from the first element, remove that element from remaining
        util::ShPtr< GraphType> sp_subgraph( new GraphType( m_Graph.IsDirected()));

        // get the first vertex in the remaining vertices
        util::ShPtr< GraphType::VertexType> sp_first_vertex( remaining_vertices.FirstElement());

        // add this node and remove from the remaining vertices
        sp_subgraph->AddVertex( sp_first_vertex->GetData());
        remaining_vertices.Remove( remaining_vertices.Begin());

        // expand the subgraph as far as it can go
        ExpandSubGraph( sp_first_vertex, *sp_subgraph, remaining_vertices);

        // create the new topology
        util::ShPtr< Topology> sp_sub_topology( new Topology( *sp_subgraph));

        // output subgraph
        BCL_MessageDbg( "Subgraph found " + sp_sub_topology->GetIdentification());

        // add it to the vector
        sub_topologies_vector.PushBack( sp_sub_topology);
      }

      // end
      return sub_topologies_vector;
    }

    //! @brief creates a subtopology composed of the given SSEGeometryInterfaces
    //! @param ELEMENTS vector of elements that should be added to the subtopology
    //! @return ShPtr to a new SubTopology
    util::ShPtr< Topology> Topology::GetSubTopology
    (
      const util::SiPtrVector< const SSEGeometryInterface> &ELEMENTS
    ) const
    {
      // make sure all given elements are part of the m_Elements
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator itr( ELEMENTS.Begin()), itr_end( ELEMENTS.End());
        itr != itr_end; ++itr
      )
      {
        // if this is not found
        if( std::find( m_Elements.Begin(), m_Elements.End(), *itr) == m_Elements.End())
        {
          BCL_MessageStd
          (
            "Given element is not found in the topology" + ( *itr)->GetIdentification()
          );
        }
      }

      // create new topology
      util::ShPtr< Topology> sp_topology( new Topology( ELEMENTS));

      // set type and orientation
      sp_topology->SetType( m_Type);
      sp_topology->SetOrientationFromType();

      // construct a new subgraph
      GraphType subgraph;

      // iterate over the elements
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator itr( ELEMENTS.Begin()), itr_end( ELEMENTS.End());
        itr != itr_end; ++itr
      )
      {
        // add this vertex
        subgraph.AddVertex( *itr);
      }

      // iterate over the elements
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator itr_a( ELEMENTS.Begin()), itr_end( ELEMENTS.End());
        itr_a != itr_end; ++itr_a
      )
      {
        // iterate over all the other elements
        for
        (
          util::SiPtrVector< const SSEGeometryInterface>::const_iterator itr_b( itr_a + 1); itr_b != itr_end; ++itr_b
        )
        {
          // get the edge
          GraphType::EdgeType this_edge( m_Graph.FindEdge( *m_Graph.FindVertex( *itr_a), *m_Graph.FindVertex( *itr_b)));

          // if the edge is defined
          if( this_edge.IsDefined())
          {
            // then add it to the subgraph
            subgraph.AddEdge( *itr_a, *itr_b, this_edge.GetData());
            SSEGeometryPacking reverse_packing( this_edge.GetData());
            reverse_packing.Reverse();
            subgraph.AddEdge( *itr_b, *itr_a, reverse_packing);
          }
        }
      }

      // set the graph and return
      sp_topology->SetGraph( subgraph);
      return sp_topology;
    }

    //! @brief orders the elements in the topology
    //! @return whether ordering succeeded
    bool Topology::OrderElements()
    {
      // if there is only one element
      if( m_Elements.GetSize() == 1)
      {
        // then it's already ordered
        return true;
      }

      // first determine if this is a barrel or not
      // find all vertices that have one edge
      util::SiPtrVector< const GraphType::VertexType> terminal_vertices
      (
        m_Graph.GetVerticesWithNEdges( 1)
      );

      // also check the number of center_vertices with connectivity of 2
      const size_t nr_center_vertices( m_Graph.GetNumberVerticesWithNEdges( 2));

      // if the sum of nr_center_vertices and terminal_vertices is not equal to number of SSEs then we have a problem
      if( nr_center_vertices + terminal_vertices.GetSize() != m_Elements.GetSize())
      {
        BCL_MessageVrb
        (
          "This topology can't be ordered, #edge vertices: " +
          util::Format()( terminal_vertices.GetSize()) +
          " + #center_vertices: " + util::Format()( nr_center_vertices) +
          " != #total vertices: " + util::Format()( m_Elements.GetSize())
        )
        return false;
      }

      // store the first vertex
      const util::SiPtr< const GraphType::VertexType> first_vertex( *m_Graph.GetVertices().Begin());

      // initialize two vertex pointer, one for previous one for this one
      util::SiPtrVector< const GraphType::VertexType> seed_vertices;
      util::SiPtrVector< const GraphType::VertexType> next_seed_vertices( first_vertex);

      // initialize order vector
      util::SiPtrVector< const SSEGeometryInterface> ordered_elements;
      ordered_elements.PushBack( first_vertex->GetData());

      // while not all SSEGeometriess are collected in the order vector
      while( ordered_elements.GetSize() != m_Elements.GetSize())
      {
        // update seed vertices
        seed_vertices = next_seed_vertices;
        // reset the next seed vertices
        next_seed_vertices = util::SiPtrVector< GraphType::VertexType>();

        // iterate over seed vertices
        for
        (
          util::SiPtrVector< const GraphType::VertexType>::const_iterator
            seed_itr( seed_vertices.Begin()), seed_itr_end( seed_vertices.End());
          seed_itr != seed_itr_end; ++seed_itr
        )
        {
          BCL_MessageDbg( "seed is: " + ( *seed_itr)->GetData()->GetIdentification());

          // iterate over all connect edges
          for
          (
            GraphType::EdgeContainerType::const_iterator
              edge_itr( ( *seed_itr)->GetEdges().Begin()), edge_itr_end( ( *seed_itr)->GetEdges().End());
            edge_itr != edge_itr_end; ++edge_itr
          )
          {
            // store the target
            util::SiPtr< const GraphType::VertexType> sp_target_vertex( edge_itr->GetTarget());
            BCL_MessageDbg
            (
              "looking at edge to: " + sp_target_vertex->GetData()->GetIdentification()
            );

            // see if this vertex is already in the data
            util::SiPtrVector< const SSEGeometryInterface>::const_iterator search_itr
            (
              std::find( ordered_elements.Begin(), ordered_elements.End(), sp_target_vertex->GetData())
            );

            // if it's already in the ordered elements
            if( search_itr != ordered_elements.End())
            {
              // then skip
              continue;
            }

            // put this into next targets
            next_seed_vertices.PushBack( sp_target_vertex);

            // find the index of this_vertex in the order_vector
            const size_t index_this_vertex( ordered_elements.GetIndex( ( *seed_itr)->GetData()));

            // if it's first element
            if( index_this_vertex == 0)
            {
              // insert at the beginning
              ordered_elements.InsertElement( 0, sp_target_vertex->GetData());
            }
            // else if it's the last element
            else if( index_this_vertex == ordered_elements.GetSize() - 1)
            {
              // insert at the end
              ordered_elements.PushBack( sp_target_vertex->GetData());
            }
            // otherwise the order vector is wrong
            else
            {
              // there is a problem with the order vector warn the user and die
              BCL_MessageStd( "The ordered geometries is not correct for this topology!");
              return false;
            }
          } // edge_itr
        } // seed_itr
      }

      // if this point is reached the elements were ordered successfully
      // so update the internal elements and return
      m_Elements = ordered_elements;
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Topology::Read( std::istream &ISTREAM)
    {
      // read members

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Topology::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // end
      return OSTREAM;
    }

    //! @brief write graphviz script
    //! @param OSTREAM the stream the script is written to
    //! @return the stream the script was written to
    std::ostream &Topology::WriteGraphVizScript
    (
      std::ostream &OSTREAM,
      const math::SumFunction< SSEGeometryPacking, double> &COLOR_EDGE
    ) const
    {
      OSTREAM << "Graph G {\n";
      OSTREAM << "overlap=scale\n"; // prevent nodes from overlapping
      OSTREAM << "splines=true\n"; // prevent edges from overlapping with nodes

      const util::Format color_format( util::Format().W( 5).FFP( 3));
      const util::Format seqid_format( util::Format().W( 3).R());
      const util::Format pos_format( util::Format().FFP( 1));

      // write vertices
      std::string nodes_helices( "{node [shape=circle,style=\"setlinewidth(5)\"]\n");
      std::string nodes_strand( "{node [shape=box,height=0.4,width=0.8,style=\"setlinewidth(5)\"]\n");

      int max_seqid( 0);
      storage::VectorND< 3, math::RunningMinMax< double> > min_max_coord;
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          ele_itr( m_Elements.Begin()), ele_itr_end( m_Elements.End());
        ele_itr != ele_itr_end;
        ++ele_itr
      )
      {
        min_max_coord( 0) += ( *ele_itr)->GetCenter()( 0);
        min_max_coord( 1) += ( *ele_itr)->GetCenter()( 1);
        min_max_coord( 2) += ( *ele_itr)->GetCenter()( 2);

        util::SiPtr< const SSE> ptr_sse( *ele_itr);

        // try to cast to sse
        if( ptr_sse.IsDefined())
        {
          max_seqid = std::max( max_seqid, ptr_sse->GetLastAA()->GetSeqID());
        }
      }

      // coordinate with largest range is outside of plane
      coord::Axis axis_a( coord::GetAxes().e_X);
      coord::Axis axis_b( coord::GetAxes().e_Y);
      {
        storage::Map< double, coord::Axis> width_axis;
        width_axis[ min_max_coord( 0).GetRange()] = coord::GetAxes().e_X;
        width_axis[ min_max_coord( 1).GetRange()] = coord::GetAxes().e_Y;
        width_axis[ min_max_coord( 2).GetRange()] = coord::GetAxes().e_Z;

        if( width_axis.GetSize() >= 2)
        {
          axis_a = width_axis.ReverseBegin()->second;
          axis_b = ( ++width_axis.ReverseBegin())->second;
          if( axis_a < axis_b)
          {
            std::swap( axis_a, axis_b);
          }
        }
      }

      const math::Range< double> seqid_range( 0.0, double( max_seqid));
      const util::ColorGradient color_gradient( seqid_range, util::GetColors().GetRainbow());
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          ele_itr( m_Elements.Begin()), ele_itr_end( m_Elements.End());
        ele_itr != ele_itr_end;
        ++ele_itr
      )
      {
        util::SiPtr< const SSE> ptr_sse( *ele_itr);
        std::string node;
        node += " \"" + ( *ele_itr)->GetIdentification() + "\" [label=\"";

        // try to cast to sse
        if( ptr_sse.IsDefined())
        {
          // label
          node += std::string( 1, ptr_sse->GetChainID()) + ' ';
          node += seqid_format( ptr_sse->GetFirstAA()->GetSeqID()) + '-';
          node += seqid_format( ptr_sse->GetLastAA()->GetSeqID()) + '\"';

          // color
          const double mid_seq_id( ( ptr_sse->GetFirstAA()->GetSeqID() + ptr_sse->GetLastAA()->GetSeqID()) / 2.0);
          const linal::Vector3D color( util::Colors::ConvertRGBToHSV( color_gradient( mid_seq_id)));
          node += ",color=\"" + color_format( color.X()) + ' ' + color_format( color.Y()) + ' ' + color_format( color.Z()) + '\"';
        }
        else
        {
          node += ( *ele_itr)->GetIdentification() + '\"';
        }
        node += ",pos=\"" + pos_format( ( *ele_itr)->GetCenter()( axis_a)) + ", " + pos_format( ( *ele_itr)->GetCenter()( axis_b)) + '\"';
        node += "]\n";

        if( ( *ele_itr)->GetType() == biol::GetSSTypes().HELIX)
        {
          nodes_helices += node;
        }
        else if( ( *ele_itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          nodes_strand += node;
        }
      }
      nodes_helices += "}\n";
      nodes_strand += "}\n";

      OSTREAM << nodes_helices << '\n';
      OSTREAM << nodes_strand << '\n';

      // write edges
      storage::Set< std::string> edges;
      const util::ColorGradient edge_color_gradient( math::Range< double>( -1.0, 1.0), storage::Vector< util::Color>::Create( util::GetColors().e_Blue, util::GetColors().e_White, util::GetColors().e_Red));
      for
      (
        GraphType::VertexContainerType::const_iterator
          vert_itr( m_Graph.GetVertices().Begin()), vert_itr_end( m_Graph.GetVertices().End());
        vert_itr != vert_itr_end;
        ++vert_itr
      )
      {
        const util::SiPtr< const SSE> ptr_sse_a( ( *vert_itr)->GetData());
        double mid_seq_id_a( util::GetUndefined< double>());
        if( ptr_sse_a.IsDefined())
        {
          mid_seq_id_a = ( ptr_sse_a->GetFirstAA()->GetSeqID() + ptr_sse_a->GetLastAA()->GetSeqID()) * 0.5;
        }

        const std::string identification_a( ( *vert_itr)->GetData()->GetIdentification());

        for
        (
          GraphType::EdgeContainerType::const_iterator
            edge_itr( ( *vert_itr)->GetEdges().Begin()), edge_itr_end( ( *vert_itr)->GetEdges().End());
          edge_itr != edge_itr_end;
          ++edge_itr
        )
        {
          const std::string identification_b( edge_itr->GetTarget()->GetData()->GetIdentification());
          // draw only one connection
          if( identification_a > edge_itr->GetTarget()->GetData()->GetIdentification())
          {
            continue;
          }

          const std::string style
          (
            edge_itr->GetData().GetContactType() == contact::GetTypes().STRAND_STRAND ?
              "bold" : "dashed"
          );
          const util::SiPtr< const SSE> ptr_sse_b( edge_itr->GetTarget()->GetData());
          double width( 1.0);
          if( ptr_sse_b.IsDefined() && util::IsDefined( mid_seq_id_a))
          {
            const double mid_seq_id_b( ( ptr_sse_b->GetFirstAA()->GetSeqID() + ptr_sse_b->GetLastAA()->GetSeqID()) * 0.5);
            math::Range< double> width_range( 1.0, 10.0);
            width = width_range( math::Absolute( mid_seq_id_b - mid_seq_id_a), seqid_range);
          }

          const double score( COLOR_EDGE.operator()( edge_itr->GetData()));
          if( math::EqualWithinTolerance( score, 0.0)) // skip interactions with neutral score
          {
            continue;
          }

          const linal::Vector3D edge_color( util::Colors::ConvertRGBToHSV( edge_color_gradient( score))); // color

          std::stringstream ss;
          ss << '\"' << identification_a << "\" -- \"" << identification_b
             << "\" ["
             << "len=" << util::Format()( edge_itr->GetData().GetDistance() / 10.0)
             << ", weight=" << util::Format()( -edge_itr->GetData().GetDistance())
             << ", style=" << style
             << ", penwidth=" << width
             << ", color=\"" + color_format( edge_color.X()) + ' ' + color_format( edge_color.Y()) + ' ' + color_format( edge_color.Z()) + '\"'
             << "];";

          edges.Insert( ss.str());
        }
      }

      // iterate over edge strings
      for( storage::Set< std::string>::const_iterator itr( edges.Begin()), itr_end( edges.End()); itr != itr_end; ++itr)
      {
        OSTREAM << *itr << '\n';
      }

      OSTREAM << "}\n";

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief builds a topology from the given vector of SSEGeometryInterface derived object, packer and criteria
    //! @param SSE_GEOMETRY_VECTOR SiPtrVector of SSEGeometryInterface derived objects
    //! @param PACKER SSEGeometryPackingPicker to be used to calculate the edges
    //! @param PACKING_CRITERIA criteria to decide whether the SSEPacking can be considered a connected edge
    //! @return topology graph constructed
    Topology::GraphType Topology::BuildTopologyGraphFromGeometries
    (
      const util::SiPtrVector< const SSEGeometryInterface> &SSE_GEOMETRY_VECTOR,
      const SSEGeometryPackingPicker &PACKER,
      const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &PACKING_CRITERIA
    )
    {
      // initialize a directed graph
      GraphType graph( true);

      BCL_MessageDbg( "adding SSEGeometries to the topology");

      // iterate over all the geometries to add the vertices
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          geometry_itr( SSE_GEOMETRY_VECTOR.Begin()), geometry_itr_end( SSE_GEOMETRY_VECTOR.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // add the vertices for these geometries
        graph.AddVertex( *geometry_itr);
      }

      // iterate over all the geometries this time to add the edges
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          geometry_itr_a( SSE_GEOMETRY_VECTOR.Begin()), geometry_itr_end( SSE_GEOMETRY_VECTOR.End());
        geometry_itr_a != geometry_itr_end; ++geometry_itr_a
      )
      {
        for
        (
          util::SiPtrVector< const SSEGeometryInterface>::const_iterator geometry_itr_b( geometry_itr_a + 1);
          geometry_itr_b != geometry_itr_end; ++geometry_itr_b
        )
        {
          // calculate the packing
          SSEGeometryPacking this_packing( ( *PACKER)->operator()( **geometry_itr_a, **geometry_itr_b));

          // if this is a valid packing
          if( PACKING_CRITERIA( this_packing))
          {
            // add the edge
            graph.AddEdge( *geometry_itr_a, *geometry_itr_b, this_packing);
            // reverse the packing
            this_packing.Reverse();
            // add the other edge
            graph.AddEdge( *geometry_itr_b, *geometry_itr_a, this_packing);
          }
        }
      }

      return graph;
    }

    //! @brief builds a topology from the given vector of SSEGeometryInterface derived object, packer and criteria
    //! @param SSE_GEOMETRY_VECTOR SiPtrVector of SSEGeometryInterface derived objects
    //! @param PACKER SSEGeometryPackingPicker to be used to calculate the edges
    //! @param PACKING_CRITERIA criteria to decide whether the SSEPacking can be considered a connected edge
    //! @return topology constructed
    util::ShPtr< Topology> Topology::BuildTopologyFromGeometries
    (
      const util::SiPtrVector< const SSEGeometryInterface> &SSE_GEOMETRY_VECTOR,
      const SSEGeometryPackingPicker &PACKER,
      const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &PACKING_CRITERIA
    )
    {
      // construct a topology and return it
      return
        util::ShPtr< Topology>
        (
          new Topology( BuildTopologyGraphFromGeometries( SSE_GEOMETRY_VECTOR, PACKER, PACKING_CRITERIA))
        );
    }

    //! @brief builds a topology from the given protein model based on its amino acid contacts
    //! @param PROTEIN the protein model to use
    //! @return the topology
    util::ShPtr< Topology> Topology::BuildTopologyFromContacts( const ProteinModel &PROTEIN)
    {
      GraphType topology_graph( true);

      // iterate over all SSEs to add the vertices
      util::SiPtrVector< const SSE> sses( PROTEIN.GetSSEs());
      for
      (
        util::SiPtrVector< const SSE>::const_iterator itr( sses.Begin()), itr_end( sses.End());
        itr != itr_end; ++itr
      )
      {
        topology_graph.AddVertex( *itr);
      }

      contact::Map contact_map( PROTEIN, 0);
      for
      (
        util::SiPtrVector< const SSE>::const_iterator itr_sse_a( sses.Begin()), itr_sse_end( sses.End());
        itr_sse_a != itr_sse_end; ++itr_sse_a
      )
      {
        for
        (
          util::SiPtrVector< const SSE>::const_iterator itr_sse_b( itr_sse_a + 1);
          itr_sse_b != itr_sse_end; ++itr_sse_b
        )
        {
          if( contact::Map::IsInContact( contact_map, **itr_sse_a, **itr_sse_b))
          {
            SSEGeometryPackingPicker packer( GetSSEGeometryPackingPickers().e_BestInteractionWeight);
            SSEGeometryPacking this_packing( ( *packer)->operator()( **itr_sse_a, **itr_sse_b));
            topology_graph.AddEdge( *itr_sse_a, *itr_sse_b, this_packing);
            this_packing.Reverse();
            topology_graph.AddEdge( *itr_sse_b, *itr_sse_a, this_packing);
          }
        }
      }

      return util::ShPtr< Topology>( new Topology( topology_graph));
    }

    //! @brief collects the SiPtrVector of SSEGeometryInterfaces from the Graph vertex data
    //! @param GRAPH reference to the graph from which the elements vector will be collected
    //! @return the SiPtrVector of SSEGeometryInterfaces from the Graph vertex data
    util::SiPtrVector< const SSEGeometryInterface> Topology::GetElementsVectorFromGraph
    (
      const Topology::GraphType &GRAPH
    )
    {
      // initialize vector
      util::SiPtrVector< const SSEGeometryInterface> elements;

      // iterate over vertices
      for
      (
        Topology::GraphType::VertexContainerType::const_iterator
          vertex_itr( GRAPH.GetVertices().Begin()), vertex_itr_end( GRAPH.GetVertices().End());
        vertex_itr != vertex_itr_end; ++vertex_itr
      )
      {
        // pushback into the elements vector
        elements.PushBack( ( *vertex_itr)->GetData());
      }

      // end
      return elements;
    }

    //! @brief static function that calculates the center of a vector SSEGeometries
    //! @param GEOMETRY_VECTOR SiPtrVector of SSE Geometries
    //! @return linal::Vector3D that represents the center of the given geometries
    linal::Vector3D Topology::CalculateCenterOfGeometries
    (
      const util::SiPtrVector< const SSEGeometryInterface> &GEOMETRY_VECTOR
    )
    {
      // initialize center
      linal::Vector3D center;

      // iterate over geometries
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          geometry_itr( GEOMETRY_VECTOR.Begin()), geometry_itr_end( GEOMETRY_VECTOR.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // sum up the center
        center += ( *geometry_itr)->GetCenter();
      }

      // normalize by the number of geometries
      center /= double( GEOMETRY_VECTOR.GetSize());

      // end
      return center;
    }

    //! @brief expands the given subgraph starting from origin_vertex and updates the remaining vertices
    //! @param SP_ORIGIN_VERTEX SiPtr to the origin vertex from the complete graph
    //! @param SUB_GRAPH non-const reference to subgraph to be expanded
    //! @param REMAINING_VERTICES non-const reference to container of remaining vertices
    void Topology::ExpandSubGraph
    (
      const util::SiPtr< const Topology::GraphType::VertexType> &SP_ORIGIN_VERTEX,
      Topology::GraphType &SUB_GRAPH,
      Topology::GraphType::VertexContainerType &REMAINING_VERTICES
    ) const
    {
      // find the corresponding vertex to SP_ORIGIN_VERTEX from this sub-graph
      util::ShPtr< GraphType::VertexType> origin_vertex_subgraph
      (
        SUB_GRAPH.FindVertex( SP_ORIGIN_VERTEX->GetData())
      );

      BCL_Assert( origin_vertex_subgraph.IsDefined(), "couldn't find SP_ORIGIN_VERTEX in the subgraph");

      // iterate over the edges of the given origin vertex
      for
      (
        GraphType::EdgeContainerType::const_iterator
          edge_itr( SP_ORIGIN_VERTEX->GetEdges().Begin()), edge_itr_end( SP_ORIGIN_VERTEX->GetEdges().End());
        edge_itr != edge_itr_end; ++edge_itr
      )
      {
        // find the corresponding vertex in this sub-graph that has the same data as
        // the vertex pointed by the edge in the topology graph
        util::ShPtr< GraphType::VertexType> target_vertex_subgraph
        (
          SUB_GRAPH.FindVertex( edge_itr->GetTarget()->GetData())
        );

        // create a bool to store whether the target vertex was found
        bool target_vertex_found( target_vertex_subgraph.IsDefined());

        // if such a vertex does not exist yet in our sub-graph
        if( !target_vertex_found)
        {
          // create a corresponding vertex
          target_vertex_subgraph = util::ShPtr< GraphType::VertexType>
          (
            new GraphType::VertexType( edge_itr->GetTarget()->GetData())
          );

          // add this vertex to the subgraph
          BCL_Assert
          (
            SUB_GRAPH.AddVertex( target_vertex_subgraph),
            "The vertex was not added to the graph"
          );

          // remove it from the remaining vertices
          RemoveVertexFromVector( *edge_itr->GetTarget(), REMAINING_VERTICES);

          // make a copy of the packing
          SSEGeometryPacking this_packing( edge_itr->GetData());

          // add the edge from origin to target
          origin_vertex_subgraph->AddEdge( target_vertex_subgraph, this_packing);

          // reverse the packing and add the edge from target to origin
          this_packing.Reverse();
          target_vertex_subgraph->AddEdge( origin_vertex_subgraph, this_packing);

          // keep expanding the subgraph
          ExpandSubGraph( edge_itr->GetTarget(), SUB_GRAPH, REMAINING_VERTICES);

        }
        // else if the vertex was found but the edges are not found
        else if( !origin_vertex_subgraph->FindEdge( *target_vertex_subgraph).IsDefined())
        {
          // make a copy of the packing
          SSEGeometryPacking this_packing( edge_itr->GetData());

          // add the edge from origin to target
          origin_vertex_subgraph->AddEdge( target_vertex_subgraph, this_packing);

          // reverse the packing and add the edge from target to origin
          this_packing.Reverse();
          target_vertex_subgraph->AddEdge( origin_vertex_subgraph, this_packing);
        }
      }
    }

    //! @brief finds and removes the vertex from a vector of vertices
    //! @param VERTEX vertex to be searched
    //! @param VERTEX_VECTOR Vertex vector to be searched
    void Topology::RemoveVertexFromVector
    (
      const Topology::GraphType::VertexType &VERTEX,
      Topology::GraphType::VertexContainerType &VERTEX_VECTOR
    )
    {
      // iterate over the vertex vector
      for
      (
        GraphType::VertexContainerType::iterator
          vertex_itr( VERTEX_VECTOR.Begin()), vertex_itr_end( VERTEX_VECTOR.End());
        vertex_itr != vertex_itr_end; ++vertex_itr
      )
      {
        // if has the same vertex data
        if( vertex_itr->GetPointer() == &VERTEX)
        {
          // remove this element and break out of loop
          VERTEX_VECTOR.RemoveElement( vertex_itr);
          return;
        }
      }
    }

  } // namespace assemble
} // namespace bcl
