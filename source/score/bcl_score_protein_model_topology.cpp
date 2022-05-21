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
#include "score/bcl_score_protein_model_topology.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_sse_pool_agreement.h"
#include "assemble/bcl_assemble_topology.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModelTopology::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelTopology())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor allowing to define the score result type, and default constructor
    ProteinModelTopology::ProteinModelTopology( score_result RESULT) :
      m_Result( RESULT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelTopology
    ProteinModelTopology *ProteinModelTopology::Clone() const
    {
      return new ProteinModelTopology( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelTopology::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate and return the percentage of recovered contacts
    //! @param TEMPLATE the known correct protein model
    //! @param MODEL the protein model to be evaluated against the TEMPLATE
    //! @return the percentage of recovered contacts
    double ProteinModelTopology::operator()
    (
      const assemble::ProteinModel &TEMPLATE,
      const assemble::ProteinModel &MODEL
    ) const
    {
      assemble::CollectorTopologyCombined collector;
      const assemble::Topology template_topology( collector.CalculateTopology( TEMPLATE.GetSSEs()));
      const assemble::Topology model_topology( collector.CalculateTopology( MODEL.GetSSEs()));
      util::SiPtrVector< const assemble::SSEGeometryInterface> template_sses( template_topology.GetElements()),
          model_sses( model_topology.GetElements());

      SSEMap mapping( CalculateSSEMapping( template_sses, model_sses));

      storage::Vector< storage::VectorND< 2, assemble::Topology::GraphType::VertexType> >
          template_edge_vertices( template_topology.GetGraph().GetEdgeVertices()),
          model_edge_vertices( model_topology.GetGraph().GetEdgeVertices());

      BCL_MessageDbg
      (
        "template_edge_vertices.size=" + util::Format()( template_edge_vertices.GetSize())
        + " model_edge_vertices.size=" + util::Format()( model_edge_vertices.GetSize())
      );

      // calculate TP and FP rates
      size_t tp( 0); // correctly identified, i.e. contacts that exist between template SSEs and  between model SSEs
      for
      (
        storage::Vector< storage::VectorND< 2, assemble::Topology::GraphType::VertexType> >::const_iterator
          itr( template_edge_vertices.Begin()), itr_end( template_edge_vertices.End());
        itr != itr_end; ++itr
      )
      {
        const storage::VectorND< 2, assemble::Topology::GraphType::VertexType> &template_edge_verteces( *itr);
        const assemble::Topology::GraphType::VertexType &template_edge_vertex_a( template_edge_verteces.First()),
            &template_edge_vertex_b( template_edge_verteces.Second());
        const util::SiPtr< const assemble::SSEGeometryInterface>
            &template_edge_vertex_a_data( template_edge_vertex_a.GetData()),
            &template_edge_vertex_b_data( template_edge_vertex_b.GetData());
        BCL_MessageDbg
        (
          "Find sse contact for template SSEs " + template_edge_vertex_a_data->GetIdentification()
          + "||" + template_edge_vertex_b_data->GetIdentification()
        );

        const util::SiPtr< const assemble::SSEGeometryInterface>
            &model_edge_vertex_a_data( mapping[ template_edge_vertex_a_data]),
            &model_edge_vertex_b_data( mapping[ template_edge_vertex_b_data]);

        if( !model_edge_vertex_a_data.IsDefined() || !model_edge_vertex_b_data.IsDefined())
        {
          BCL_MessageDbg( "No mapping for template SSEs to model SSEs");
          continue;
        }

        const util::ShPtr< assemble::Topology::GraphType::VertexType>
            &model_edge_vertex_a( model_topology.GetGraph().FindVertex( model_edge_vertex_a_data)),
            &model_edge_vertex_b( model_topology.GetGraph().FindVertex( model_edge_vertex_b_data));
        // if the model_edge_vertex_{a,b}_data are defined, then we will find the model_edge_vertex_{a,b} in the graph

        if( !model_edge_vertex_a->FindEdge( *model_edge_vertex_b).IsDefined())
        {
          BCL_MessageDbg
          (
            "No contact for model SSEs" + model_edge_vertex_a_data->GetIdentification()
            + "||" + model_edge_vertex_b_data->GetIdentification()
          );
          continue;
        }

        BCL_MessageDbg
        (
          "Found contact for model SSEs " + model_edge_vertex_a_data->GetIdentification()
          + "||" + model_edge_vertex_b_data->GetIdentification()
        );
        ++tp;
      }

      size_t fp( model_edge_vertices.GetSize() - tp); // incorrectly identified = all identified - correctly identified
      size_t fn( template_edge_vertices.GetSize() - tp); // incorrectly rejected = all correct - correctly identified
      BCL_MessageVrb
      (
        "tp=" + util::Format()( tp) + " fp=" + util::Format()( fp) + " fn=" + util::Format()( fn)
        + " correct=tp+fn=" + util::Format()( template_edge_vertices.GetSize())
        + " predicted=tp+fp=" + util::Format()( model_edge_vertices.GetSize())
      );

      double tpr( ( double)tp / template_edge_vertices.GetSize()); // true positive rate or sensitivity
      double ppv( ( double)tp / model_edge_vertices.GetSize()); // positive predictive value or precision
      double f( 2.0 * tp / ( 2 * tp + fp + fn)); // f-score
      BCL_MessageVrb( "f=" + util::Format()( f) + " tpr=" + util::Format()( tpr) + " ppv=" + util::Format()( ppv));

      switch( m_Result)
      {
        case score_tpr: return tpr;
        case score_ppv: return ppv;
        case score_f  : break;
        default       : break;
      }
      return f;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelTopology::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelTopology::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  /////////////////////
  // helper funtions //
  /////////////////////

    //! @brief calculates a mapping template_sse->model_sse based on the Q3 score
    //! @param TEMPLATE_SSES all template SSEs
    //! @param MODEL_SSES all model SSEs
    //! @return the mapping
    ProteinModelTopology::SSEMap ProteinModelTopology::CalculateSSEMapping
    (
      util::SiPtrVector< const assemble::SSEGeometryInterface> TEMPLATE_SSES,
      util::SiPtrVector< const assemble::SSEGeometryInterface> MODEL_SSES
    )
    {
      SSEMap mapping; // result
      assemble::SSEPoolAgreement ssepool_score; // create score before loop

      // try to find a matching model SSE for every template SSE; stop if one of the vectors is empty
      while( !TEMPLATE_SSES.IsEmpty() && !MODEL_SSES.IsEmpty())
      {
        // remove from the end b/c these are vectors
        util::SiPtr< const assemble::SSEGeometryInterface> current_template_sse_ptr( TEMPLATE_SSES.LastElement());
        TEMPLATE_SSES.RemoveElement( --TEMPLATE_SSES.End());
        BCL_MessageDbg( "Find match for template_sse: " + current_template_sse_ptr->GetIdentification());

        // initialize everything before the loop
        util::SiPtrVector< const assemble::SSE> template_sse_vector( current_template_sse_ptr); // create sse vector
        assemble::SSEPool template_pool( template_sse_vector, true, false); // create pool from sse vector with template sse
        double max_score( 0.0); // best Q3 score found so far
        util::SiPtrVector< const assemble::SSEGeometryInterface>::iterator max_score_model_sse_itr( MODEL_SSES.End());

        for
        (
          util::SiPtrVector< const assemble::SSEGeometryInterface>::iterator
            itr( MODEL_SSES.Begin()), itr_end( MODEL_SSES.End());
          itr != itr_end; ++itr
        )
        {
          util::SiPtrVector< const assemble::SSE> model_sse_vector( *itr);
          assemble::SSEPool model_pool( model_sse_vector, true, false);
          double score( ssepool_score.Q3Score( template_pool, model_pool));
          if( score > max_score)
          {
            max_score = score;
            max_score_model_sse_itr = itr;
          }
          BCL_MessageDbg
          (
            "Compare with model_sse: " + ( **itr).GetIdentification() + " score=" + util::Format()( score)
            + " max_score=" + util::Format()( max_score)
          );
        }

        // if we have a match, create mapping, and remove model sses from the vector (the template sse is removed above)
        if( max_score_model_sse_itr != MODEL_SSES.End() && max_score_model_sse_itr->IsDefined())
        {
          BCL_MessageDbg
          (
            "Create mapping template_sse->model_sse: " + current_template_sse_ptr->GetIdentification()
            + "->" + ( **max_score_model_sse_itr).GetIdentification()
          );
          mapping.Insert
          (
            std::pair
            <
              util::SiPtr< const assemble::SSEGeometryInterface>, util::SiPtr< const assemble::SSEGeometryInterface>
            >( current_template_sse_ptr, *max_score_model_sse_itr)
          );
          MODEL_SSES.RemoveElement( max_score_model_sse_itr);
        }
        else
        {
          BCL_MessageVrb( "No mapping for template_sse " + current_template_sse_ptr->GetIdentification());
        }
      }

      return mapping;
    }

  } // namespace score
} // namespace bcl
