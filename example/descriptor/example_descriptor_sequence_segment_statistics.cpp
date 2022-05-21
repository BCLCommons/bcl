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
#include "descriptor/bcl_descriptor_sequence_segment_statistics.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "descriptor/bcl_descriptor_segment_finder.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_sequence_segment_statistics.cpp
  //!
  //! @author mendenjl
  //! @date Mar 06, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSequenceSegmentStatistics :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSequenceSegmentStatistics *Clone() const
    {
      return new ExampleDescriptorSequenceSegmentStatistics( *this);
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
        sequence_segments
        (
          "SequenceSegmentStatistics("
            "descriptor=AlphabeticNumber,"
            "condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3))"
          ")"
        );

      BCL_ExampleAssert( sequence_segments.IsDefined(), true);

      // ensure that all statistics were added by default
      BCL_ExampleCheck( sequence_segments->GetSizeOfFeatures(), size_t( descriptor::SegmentFinder::s_NumberStatistics));

      // compute # segments, fraction segments, elements, length segments ave, etc.
      util::Implementation< descriptor::Base< char, float> >
        sequence_segments_sequence_stats
        (
          "SequenceSegmentStatistics("
          "  descriptor=AlphabeticNumber,"
          "  condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3)),"
          "  statistics(NumSegments,FractionSegments,FractionElements,"
          "             LengthSegmentsAve,LengthSegmentsElementwiseAve)"
          ")"
        );
      BCL_ExampleCheck( sequence_segments_sequence_stats->GetSizeOfFeatures(), 5);

      // Compute some of the basic segment statistics on a string with 3 logical segments.  Remember that the condition
      // for the segment is that the letter is >= c, so
      // abcdefgAAAAAABB has sequence numbers:
      // 001111122222222
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sequence_segments_sequence_stats, "abcdefgAAAAAABB", 2),
        //NSg  fSg fEl  LSA  LSEA
        "3.00 1.00 1.00 5.00 6.20 ; " // overall statistics
      );

      // Compute all statistics on the same string
      // CSEA calculation: ( 0 x 2 + 0 x 6 + 1 x 5) / (2 + 6 + 5) = 5 / 13
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sequence_segments, "abcdefgAAAAAA", 2),
        //NSg  fSg fEl  LSA  LSSD LSEA LSESD DSA DSSD DSEA DSESD CSA CSSD CSEA CSESD
        "3.00 1.00 1.00 4.33 1.70 5.00 1.36 2.50 1.78 2.62 1.89 0.33 0.47 0.38 0.49 ; " // overall statistics
      );

      // compute segment statistics using conditions
      util::Implementation< descriptor::Base< char, float> >
        sequence_segments_conditions
        (
          "SequenceSegmentStatistics("
            "descriptor=AlphabeticNumber,"
            "condition=GreaterEqual(lhs=AlphabeticNumber,rhs=Constant(3)),"
            "conditions((0),(1),(2))"
          ")"
        );

      // note that the condition 2 does not exist in the sequence; (since greater equal only returns 0 or 1), but
      // statistics class has no way of knowing this and so should return the correct numbers
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( sequence_segments_conditions, "abcdefgAAAAAABB", 2),
        //NSg  fSg fEl  LSA  LSSD LSEA LSESD DSA DSSD DSEA DSESD CSA CSSD CSEA CSESD
        "3.00 1.00 1.00 5.00 2.45 6.20 2.14 2.58 1.71 2.53 1.75 0.33 0.47 0.33 0.47 "   // overall
        "2.00 0.67 0.67 5.00 3.00 6.80 2.40 1.38 0.12 1.30 0.10 0.00 0.00 0.00 0.00 "   // condition = 0
        "1.00 0.33 0.33 5.00 0.00 5.00 0.00 5.00 0.00 5.00 0.00 1.00 0.00 1.00 0.00 "   // condition = 1
        "0.00 0.00 0.00 0.00 0.00 0.00 0.00 2.58 1.71 2.53 1.75 2.00 0.00 2.00 0.00 ; " // condition = 2
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorSequenceSegmentStatistics

  const ExampleClass::EnumType ExampleDescriptorSequenceSegmentStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSequenceSegmentStatistics())
  );
} // namespace bcl
