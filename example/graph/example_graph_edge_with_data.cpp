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
#include "graph/bcl_graph_edge_with_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_edge_with_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphEdgeWithData :
    public ExampleInterface
  {
  public:

    ExampleGraphEdgeWithData *Clone() const
    {
      return new ExampleGraphEdgeWithData( *this);
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
      graph::EdgeWithData< std::string, size_t> graph_edge_with_data_a;

      // test constructor from vertex and color
      graph::VertexWithData< std::string, size_t> graph_vertex( "vertex b");
      util::ShPtr< graph::VertexWithData< std::string, size_t> > sp_graph_vertex
      (
        new graph::VertexWithData< std::string, size_t>( "vertex b")
      );
      graph::EdgeWithData< std::string, size_t> graph_edge_with_data_b( sp_graph_vertex, size_t( 2));

      // test copy constructor
      graph::EdgeWithData< std::string, size_t> graph_edge_with_data_c( graph_edge_with_data_b);

      // test Clone
      util::ShPtr< graph::EdgeWithData< std::string, size_t> > sp_graph_edge_with_data_d( graph_edge_with_data_c.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetData
      BCL_ExampleCheck
      (
        ( graph::EdgeWithData< std::string, size_t>( sp_graph_vertex, size_t( 2)).GetData()), 2
      );

      // test GetTarget
      BCL_ExampleCheck
      (
        ( graph::EdgeWithData< std::string, size_t>( sp_graph_vertex, size_t( 2)).GetTarget()->GetData()),
        sp_graph_vertex->GetData()
      );

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

  }; //end ExampleGraphEdgeWithData

  const ExampleClass::EnumType ExampleGraphEdgeWithData::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphEdgeWithData())
  );

} // namespace bcl
