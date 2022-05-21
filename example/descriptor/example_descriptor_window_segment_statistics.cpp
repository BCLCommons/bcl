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
#include "descriptor/bcl_descriptor_window_segment_statistics.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_window_segment_statistics.cpp
  //!
  //! @author mendenjl
  //! @date Mar 05, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWindowSegmentStatistics :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWindowSegmentStatistics *Clone() const
    {
      return new ExampleDescriptorWindowSegmentStatistics( *this);
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

      // retrieve segment statistics for the current segment and neighboring segments
      util::Implementation< descriptor::Base< char, float> >
        window_segments
        (
          "WindowSegmentStatistics("
            "radius=1,"
            "descriptor=AlphabeticNumber,"
            "condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3)),"
            "statistics("
            "  Condition,Length,DescriptorAve,DescriptorSD,SegmentsBefore,SegmentsAfter,SegmentsTillNearestEnd,"
            "  ElementsFromSegmentStart,ElementsFromSegmentEnd,ElementsTillNearestSegmentEnd,"
            "  ConditionalSegmentsBefore,ConditionalSegmentsAfter,ConditionalSegmentsTillNearestEnd)"
          ")"
        );

      BCL_ExampleAssert( window_segments.IsDefined(), true);

      // ensure that all statistics were added by default
      BCL_ExampleCheck( window_segments->GetSizeOfFeatures(), 39);

      // compute segment statistics that are relative positions within the sequence and/or segment
      util::Implementation< descriptor::Base< char, float> >
        window_segments_relative_positions
        (
          "WindowSegmentStatistics("
            "radius=0,"
            "descriptor=AlphabeticNumber,"
            "condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3)),"
            "statistics(SegmentsBefore,SegmentsAfter,SegmentsTillNearestEnd,"
            "           ElementsFromSegmentStart,ElementsFromSegmentEnd,ElementsTillNearestSegmentEnd,"
            "           ConditionalSegmentsBefore,ConditionalSegmentsAfter,ConditionalSegmentsTillNearestEnd)"
          ")"
        );

      // try these segment positions on a sequence with 3 segments; e.g.
      // abcdefgAA has segment numbers 001111122 because the condition is >= c
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_segments_relative_positions, "abcdefgAA", 2),
        //S-  S+   STrm EStr EEnd ETrm CSB  CSA  CSNear
        "0.00 2.00 0.00 0.00 1.00 0.00 0.00 1.00 0.00 ; " // a
        "0.00 2.00 0.00 1.00 0.00 0.00 0.00 1.00 0.00 ; " // b
        "1.00 1.00 1.00 0.00 4.00 0.00 0.00 0.00 0.00 ; " // c
        "1.00 1.00 1.00 1.00 3.00 1.00 0.00 0.00 0.00 ; " // d
        "1.00 1.00 1.00 2.00 2.00 2.00 0.00 0.00 0.00 ; " // e
        "1.00 1.00 1.00 3.00 1.00 1.00 0.00 0.00 0.00 ; " // f
        "1.00 1.00 1.00 4.00 0.00 0.00 0.00 0.00 0.00 ; " // g
        "2.00 0.00 0.00 0.00 1.00 0.00 1.00 0.00 0.00 ; " // A
        "2.00 0.00 0.00 1.00 0.00 0.00 1.00 0.00 0.00 ; " // A
      );

      // compute segment statistics that can be returned for multiple segments
      util::Implementation< descriptor::Base< char, float> >
        window_segments_neighbor_segment_info
        (
          "WindowSegmentStatistics("
            "radius=1,"
            "descriptor=AlphabeticNumber,"
            "condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3)),"
            "statistics(Condition, Length, DescriptorAve, DescriptorSD)"
          ")"
        );

      // try these segment statistics on a sequence with 3 segments; e.g.
      // abcdefgAA has segment numbers 001111122 because the condition is >= c
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_segments_neighbor_segment_info, "abcdefgAA", 2),
        //  C    L  Ave   SD   C-   L- Ave-  SD-
        "0.00 2.00 1.50 0.50 1.00 5.00 5.00 1.41 1.00 5.00 5.00 1.41 ; " // a
        "0.00 2.00 1.50 0.50 1.00 5.00 5.00 1.41 1.00 5.00 5.00 1.41 ; " // b
        "1.00 5.00 5.00 1.41 0.00 2.00 1.50 0.50 0.00 2.00 1.00 0.00 ; " // c
        "1.00 5.00 5.00 1.41 0.00 2.00 1.50 0.50 0.00 2.00 1.00 0.00 ; " // d
        "1.00 5.00 5.00 1.41 0.00 2.00 1.50 0.50 0.00 2.00 1.00 0.00 ; " // e
        "1.00 5.00 5.00 1.41 0.00 2.00 1.50 0.50 0.00 2.00 1.00 0.00 ; " // f
        "1.00 5.00 5.00 1.41 0.00 2.00 1.50 0.50 0.00 2.00 1.00 0.00 ; " // g
        "0.00 2.00 1.00 0.00 1.00 5.00 5.00 1.41 1.00 5.00 5.00 1.41 ; " // A
        "0.00 2.00 1.00 0.00 1.00 5.00 5.00 1.41 1.00 5.00 5.00 1.41 ; " // A
      );

      // test the edge case where there is only one segment, with all the same values
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_segments, "aaaa", 2),
        "0.00 4.00 1.00 0.00 0.00 0.00 0.00 0.00 3.00 0.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 0.00 3.00 "
        "0.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 0.00 3.00 0.00 0.00 0.00 0.00 ; "
        "0.00 4.00 1.00 0.00 0.00 0.00 0.00 1.00 2.00 1.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 1.00 2.00 "
        "1.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 1.00 2.00 1.00 0.00 0.00 0.00 ; "
        "0.00 4.00 1.00 0.00 0.00 0.00 0.00 2.00 1.00 1.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 2.00 1.00 "
        "1.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 2.00 1.00 1.00 0.00 0.00 0.00 ; "
        "0.00 4.00 1.00 0.00 0.00 0.00 0.00 3.00 0.00 0.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 3.00 0.00 "
        "0.00 0.00 0.00 0.00 0.00 4.00 1.00 0.00 0.00 0.00 0.00 3.00 0.00 0.00 0.00 0.00 0.00 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWindowSegmentStatistics

  const ExampleClass::EnumType ExampleDescriptorWindowSegmentStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWindowSegmentStatistics())
  );
} // namespace bcl
