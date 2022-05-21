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
#include "graph/bcl_graph_ring.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_ring.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 12, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphRing :
    public ExampleInterface
  {
  public:

    ExampleGraphRing *Clone() const
    {
      return new ExampleGraphRing( *this);
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

      // make a default-constructed Ring
      graph::Ring default_Ring;

      // check that the default Ring is empty and non-cyclical
      BCL_ExampleCheck( default_Ring.GetSize(), 0);
      BCL_ExampleCheck( default_Ring.Contains( 0), false);

      // create a Ring from vertex 1 - 3 - 5
      graph::Ring one_three_five( 6, storage::Vector< size_t>::Create( 1, 3, 5));

      // create an Ring from vertex 3 - 5 - 1
      // this should be identical to one_three, since the Ring is undirected
      graph::Ring three_five_one( 6, storage::Vector< size_t>::Create( 3, 5, 1));

      // check that the Rings were identical
      BCL_ExampleCheck( one_three_five.Identical( three_five_one), true);

      // check that the right vertices were set
      BCL_ExampleCheck( one_three_five.Contains( 0), false);
      BCL_ExampleCheck( one_three_five.Contains( 1), true);
      BCL_ExampleCheck( one_three_five.Contains( 3), true);
      BCL_ExampleCheck( one_three_five.GetVertices()( 0), size_t( 1));
      BCL_ExampleCheck( one_three_five.GetVertices()( 1), size_t( 3));
      BCL_ExampleCheck( one_three_five.GetVertices()( 2), size_t( 5));
      BCL_ExampleCheck( one_three_five.FirstElement(), size_t( 1));
      BCL_ExampleCheck( one_three_five.LastElement(), size_t( 5));
      BCL_ExampleCheck( one_three_five.GetSize(), size_t( 3));

      // create an directed Ring from vertex 1 - 2 - 3 - 4
      graph::Ring one_four( 5, storage::Vector< size_t>::Create( 1, 2, 3, 4));

      // create an directed Ring from vertex 4 to 1
      // this should be different from one_three, since the Ring is directed
      graph::Ring four_one( 5, storage::Vector< size_t>::Create( 4, 3, 2, 1));

      // check that the Rings were identical
      BCL_ExampleCheck( one_four.Identical( four_one), true);

      // create an undirected Ring from vertex 1 - 4 - 6 - 5
      graph::Ring one_four_six_five( 7, storage::Vector< size_t>::Create( 1, 4, 6, 5));
      graph::Ring four_one_six_five( 7, storage::Vector< size_t>::Create( 4, 1, 6, 5));
      graph::Ring five_six_one_four( 7, storage::Vector< size_t>::Create( 5, 6, 1, 4));

      // test that the two rings are not identical
      BCL_ExampleCheck( four_one_six_five.Identical( one_four_six_five), false);
      BCL_ExampleCheck( five_six_one_four.Identical( four_one_six_five), true);

      // test GetOverlap
      BCL_ExampleCheck( one_four_six_five.GetOverlap( one_three_five).GetVertices().GetSize(), size_t( 2));
      BCL_ExampleCheck( one_four_six_five.GetOverlap( one_three_five).FirstElement(), size_t( 5));
      BCL_ExampleCheck( one_four_six_five.GetOverlap( one_three_five).LastElement(), size_t( 1));

      // test self overlap
      BCL_ExampleCheck( one_four_six_five.GetOverlap( one_four_six_five).GetSize(), size_t( 4));

      // test overlap between two rings
      graph::Path overlap( one_four_six_five.GetOverlap( one_three_five));
      graph::Path removed_path( one_four_six_five.Remove( overlap));
      BCL_ExampleCheck( removed_path.GetVertices().GetSize(), size_t( 2));
      BCL_ExampleCheck( removed_path.GetVertices()( 0), size_t( 4));
      BCL_ExampleCheck( removed_path.GetVertices()( 1), size_t( 6));

      // fused ring
      graph::Ring fused_ring( 10, storage::Vector< size_t>::Create( 0, 4, 7, 9, 8, 5, 1, 3, 6, 2));

      // create component rings
      // check all possible types of numbering schemes i.e. numbering can start at fused vertices or in the middle of a ring
      // numbering may be in clockwise or counterclockwise direction
      graph::Ring large_ring_a( 10, storage::Vector< size_t>::Create( 0, 4, 7, 9, 8, 5, 1));
      graph::Ring large_ring_b( 10, storage::Vector< size_t>::Create( 0, 1, 5, 8, 9, 7, 4));
      graph::Ring large_ring_c( 10, storage::Vector< size_t>::Create( 9, 8, 5, 1, 0, 4, 7));
      graph::Ring small_ring_a( 10, storage::Vector< size_t>::Create( 0, 1, 3, 6, 2));
      graph::Ring small_ring_b( 10, storage::Vector< size_t>::Create( 0, 2, 6, 3, 1));
      graph::Ring small_ring_c( 10, storage::Vector< size_t>::Create( 6, 3, 1, 0, 2));

      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_a, small_ring_a)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_a, small_ring_b)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_a, small_ring_c)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_b, small_ring_a)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_b, small_ring_b)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_b, small_ring_c)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_c, small_ring_a)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_c, small_ring_b)), true );
      BCL_ExampleCheck( fused_ring.Identical( graph::Ring::FuseRings( large_ring_c, small_ring_c)), true );

      // checking with a bridge ring system
      graph::Ring bridge_ring( 9, storage::Vector< size_t>::Create( 0, 4, 7, 1, 5, 2));
      graph::Ring ring_a( 9, storage::Vector< size_t>::Create( 0, 4, 7, 1, 6, 3));
      graph::Ring ring_b( 9, storage::Vector< size_t>::Create( 0, 3, 6, 1, 5, 2));
      BCL_ExampleCheck( bridge_ring.Identical( graph::Ring::FuseRings( ring_a, ring_b)), true);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphRing

  const ExampleClass::EnumType ExampleGraphRing::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphRing())
  );

} // namespace bcl
