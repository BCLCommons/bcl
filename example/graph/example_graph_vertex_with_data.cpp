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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "graph/bcl_graph_vertex_with_data.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_graph_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_vertex_with_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphVertexWithData :
    public ExampleInterface
  {
  public:

    ExampleGraphVertexWithData *Clone() const
    {
      return new ExampleGraphVertexWithData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      graph::GraphWithData< std::string, size_t>::VertexType graph_vertex_a;

      // test constructor from data and color
      graph::GraphWithData< std::string, size_t>::VertexType graph_vertex_b( "vertex b");

      // test copy constructor
      graph::GraphWithData< std::string, size_t>::VertexType graph_vertex_c( graph_vertex_b);

      // test Clone
      util::ShPtr< graph::GraphWithData< std::string, size_t>::VertexType> sp_graph_vertex_d( graph_vertex_c.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetData
      BCL_ExampleCheck( ( graph::GraphWithData< std::string, size_t>::VertexType( "vertex b").GetData()), "vertex b");

      // test Edge methods so first create some more vertices and edges
      graph::GraphWithData< std::string, size_t>::VertexType graph_vertex_e( "vertex e");
      util::ShPtr< graph::GraphWithData< std::string, size_t>::VertexType> sp_graph_vertex_e
      (
        new graph::GraphWithData< std::string, size_t>::VertexType( graph_vertex_e)
      );
      graph::GraphWithData< std::string, size_t>::VertexType graph_vertex_f( "vertex f");
      util::ShPtr< graph::GraphWithData< std::string, size_t>::VertexType> sp_graph_vertex_f
      (
        new graph::GraphWithData< std::string, size_t>::VertexType( graph_vertex_f)
      );
      graph_vertex_b.AddEdge( sp_graph_vertex_e, 1);
      graph_vertex_b.AddEdge( sp_graph_vertex_f, 2);

      // call GetEdges and make sure we have them
      graph::GraphWithData< std::string, size_t>::EdgeContainerType graph_b_edges
      (
        graph_vertex_b.GetEdges()
      );

      // make sure adding the edges and getting them worked
      BCL_ExampleCheck( graph_vertex_b.GetEdges().FirstElement().GetData(), 1);
      BCL_ExampleCheck( graph_vertex_b.GetEdges().FirstElement().GetTarget()->GetData(), "vertex e");
      BCL_ExampleCheck( graph_vertex_b.GetEdges().LastElement().GetData(), 2);
      BCL_ExampleCheck( graph_vertex_b.GetEdges().LastElement().GetTarget()->GetData(), "vertex f");

      // test degree, should be 2
      BCL_ExampleCheck( graph_vertex_b.GetDegree(), 2);

      // test delete edge
      graph_vertex_b.DeleteEdge( "vertex f");
      BCL_ExampleIndirectCheck( graph_vertex_b.GetDegree(), 1, "graph_vertex_b.DeleteEdge( \"vertex f\")");

      // test find edge
      BCL_ExampleCheck( graph_vertex_b.FindEdge( "vertex e").GetTarget()->GetData(), "vertex e");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphVertexWithData

  const ExampleClass::EnumType ExampleGraphVertexWithData::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphVertexWithData())
  );

} // namespace bcl
