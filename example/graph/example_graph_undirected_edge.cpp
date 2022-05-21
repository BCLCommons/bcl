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
#include "graph/bcl_graph_undirected_edge.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_undirected_edge.cpp
  //!
  //! @author mendenjl
  //! @date Jun 06, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphUndirectedEdge :
    public ExampleInterface
  {
  public:

    ExampleGraphUndirectedEdge *Clone() const
    {
      return new ExampleGraphUndirectedEdge( *this);
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

      // constructors from members
      graph::UndirectedEdge< size_t> two_one_three( 2, 1, 3);
      graph::UndirectedEdge< size_t> zero_one_three( 0, 4, 3);

    /////////////////
    // data access //
    /////////////////

      // check that constructors always order the indices correctly
      BCL_ExampleCheck( two_one_three.GetIndexHigh(), 2);
      BCL_ExampleCheck( two_one_three.GetIndexLow(),  1);
      BCL_ExampleCheck( two_one_three.GetEdgeData(),  3);
      BCL_ExampleCheck( zero_one_three.GetIndexHigh(), 4);
      BCL_ExampleCheck( zero_one_three.GetIndexLow(),  0);
      BCL_ExampleCheck( zero_one_three.GetEdgeData(),  3);
      BCL_ExampleCheck( two_one_three, two_one_three);
      BCL_ExampleCheck( zero_one_three, zero_one_three);
      BCL_ExampleCheck( zero_one_three == two_one_three, false);
      BCL_ExampleCheck( zero_one_three < two_one_three, true);
      BCL_ExampleCheck( two_one_three < zero_one_three, false);

      // test io
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( two_one_three, graph::UndirectedEdge< size_t>()), true);
      BCL_ExampleCheck( TestBCLObjectOutputDiffers( zero_one_three, graph::UndirectedEdge< size_t>()), true);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphUndirectedEdge

  const ExampleClass::EnumType ExampleGraphUndirectedEdge::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphUndirectedEdge())
  );

} // namespace bcl
