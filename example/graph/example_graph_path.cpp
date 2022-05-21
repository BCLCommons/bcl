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
#include "graph/bcl_graph_path.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_path.cpp
  //!
  //! @author mendenjl
  //! @date Oct 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphPath :
    public ExampleInterface
  {
  public:

    ExampleGraphPath *Clone() const
    {
      return new ExampleGraphPath( *this);
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

      // make a default-constructed path
      graph::Path default_path;

      // check that the default path is empty and non-cyclical
      BCL_ExampleCheck( default_path.GetSize(), 0);
      BCL_ExampleCheck( default_path.Contains( 0), false);

      // create an undirected path from vertex 1 to 3
      graph::Path one_three_undirected
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 1, 3),
        true
      );

      // create an undirected path from vertex 3 to 1
      // this should be identical to one_three, since the path is undirected
      graph::Path three_one_undirected
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 3, 1),
        true
      );

      // check that the paths were identical
      BCL_ExampleCheck( one_three_undirected.Identical( three_one_undirected), true);

      // check that the right vertices were set
      BCL_ExampleCheck( one_three_undirected.Contains( 0), false);
      BCL_ExampleCheck( one_three_undirected.Contains( 1), true);
      BCL_ExampleCheck( one_three_undirected.Contains( 3), true);
      BCL_ExampleCheck( one_three_undirected.IsUndirected(), true);
      BCL_ExampleCheck( one_three_undirected.IsDirected(), false);
      BCL_ExampleCheck( one_three_undirected.EndsAt( 0), false);
      BCL_ExampleCheck( one_three_undirected.EndsAt( 1), true);
      BCL_ExampleCheck( one_three_undirected.EndsAt( 3), true);

      // create an directed path from vertex 1 to 3
      graph::Path one_three_directed
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 1, 3),
        false
      );

      // create an directed path from vertex 3 to 1
      // this should be different from one_three, since the path is directed
      graph::Path three_one_directed
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 3, 1),
        false
      );

      // check that the paths were identical
      BCL_ExampleCheck( one_three_directed.Identical( three_one_directed), false);

      // check that the directedness was set
      BCL_ExampleCheck( one_three_directed.IsUndirected(), false);
      BCL_ExampleCheck( one_three_directed.IsDirected(), true);

      // try reversing three_one_directed, now it should be identical to 1-3
      three_one_directed.Reverse();
      BCL_ExampleCheck( one_three_directed.Identical( three_one_directed), true);
      three_one_directed.Reverse();

      // create an undirected path from vertex 3 to 4
      graph::Path three_four_undirected
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 3, 4),
        true
      );

      // create an undirected path from vertex 0 to 2
      graph::Path zero_two_undirected
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 0, 2),
        true
      );

      // create an undirected path from vertex 0 to 1
      graph::Path zero_one_undirected
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 0, 1),
        true
      );

      // combine 1-3 and 3-5 to make 1-3-4
      graph::Path one_three_four_undirected
      (
        one_three_undirected,
        three_four_undirected,
        graph::Path( 5, storage::Vector< size_t>( 1, 3), true) // combine them at vertex 3
      );

      // check that the paths are not identical
      BCL_ExampleCheck( one_three_four_undirected.Identical( three_four_undirected), false);

      // check that the path contains 1, 3, and 4
      BCL_ExampleCheck( one_three_four_undirected.Contains( 0), false);
      BCL_ExampleCheck( one_three_four_undirected.Contains( 1), true);
      BCL_ExampleCheck( one_three_four_undirected.Contains( 3), true);
      BCL_ExampleCheck( one_three_four_undirected.Contains( 4), true);
      BCL_ExampleCheck( one_three_four_undirected.IsUndirected(), true);
      BCL_ExampleCheck( one_three_four_undirected.IsDirected(), false);
      BCL_ExampleCheck( one_three_four_undirected.EndsAt( 0), false);
      BCL_ExampleCheck( one_three_four_undirected.EndsAt( 1), true);
      BCL_ExampleCheck( one_three_four_undirected.EndsAt( 3), false);
      BCL_ExampleCheck( one_three_four_undirected.EndsAt( 4), true);

      // check that one_three_four_undirected crosses itself and the paths that it is made up of
      // it should also meet them
      BCL_ExampleCheck( one_three_four_undirected.Crosses( one_three_undirected), true);
      BCL_ExampleCheck( one_three_four_undirected.Crosses( three_four_undirected), true);
      BCL_ExampleCheck( one_three_four_undirected.Crosses( zero_two_undirected), false);

      // test covers
      BCL_ExampleCheck( one_three_four_undirected.Covers( one_three_undirected), true);
      BCL_ExampleCheck( one_three_four_undirected.Covers( three_four_undirected), true);
      BCL_ExampleCheck( one_three_four_undirected.Covers( zero_one_undirected), false);

      // test Connects
      BCL_ExampleCheck( one_three_undirected.Connects( three_four_undirected), true);
      BCL_ExampleCheck( one_three_undirected.Connects( zero_two_undirected), false);

      // test CountVertexRepetitions on cases where it should return 0
      BCL_ExampleCheck( one_three_undirected.CountVertexRepetitions(), 0);
      BCL_ExampleCheck( one_three_four_undirected.CountVertexRepetitions(), 0);

      // test CountVertexRepetitions on some cases where the path does repeat itself
      BCL_ExampleCheck
      (
        graph::Path
        (
          5, // the graph is of size 5
          storage::Vector< size_t>::Create( 1, 2, 1),
          true
        ).CountVertexRepetitions(),
        1
      );

      BCL_ExampleCheck
      (
        graph::Path
        (
          5, // the graph is of size 5
          storage::Vector< size_t>::Create( 0, 4, 0, 4, 0),
          true
        ).CountVertexRepetitions(),
        3
      );

      // test the EquivalentTour function, which tells us whether the paths tour the same vertices and have the same start and end
      graph::Path one_two_three_four
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 1, 2, 3, 4),
        true
      );

      graph::Path one_three_two_four
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 1, 3, 2, 4),
        true
      );

      graph::Path one_four_two_three
      (
        5, // the graph is of size 5
        storage::Vector< size_t>::Create( 1, 4, 2, 3),
        true
      );

      BCL_ExampleCheck( one_two_three_four.EquivalentTour( one_three_two_four), true);
      BCL_ExampleCheck( one_two_three_four.EquivalentTour( one_four_two_three), false);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( one_two_three_four, default_path),
        true,
        "Path I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphPath

  const ExampleClass::EnumType ExampleGraphPath::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphPath())
  );

} // namespace bcl
