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
#include "math/bcl_math_range_set.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_range_set.cpp
  //!
  //! @author mendenjl
  //! @date Jan 31, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRangeSet :
    public ExampleInterface
  {
  public:

    ExampleMathRangeSet *Clone() const
    {
      return new ExampleMathRangeSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::RangeSet< double> range_set_default;
      BCL_ExampleCheck( math::RangeSet< double>().IsEmpty(), true);

      // construct from range with both side closed interval and min and max
      math::RangeSet< double> range_set_cm5_p5c( math::Range< double>( -5.0, 5.0));
      BCL_ExampleCheck( range_set_cm5_p5c.IsWithin( -5.0), true);
      BCL_ExampleCheck( range_set_cm5_p5c.IsWithin( 5.0), true);
      BCL_ExampleCheck( range_set_cm5_p5c.IsWithin( -5.1), false);
      BCL_ExampleCheck( range_set_cm5_p5c.IsWithin( -4.9), true);
      BCL_ExampleCheck( range_set_cm5_p5c.IsEmpty(), false);

      // construct from a vector of ranges
      storage::Vector< math::Range< double> > overlapping_ranges;
      // the ranges need not be in any particular order.  If ranges overlap, they will be combined
      overlapping_ranges.PushBack( math::Range< double>( 5.0, 10.0));
      overlapping_ranges.PushBack( math::Range< double>( -5.0, 5.0));
      overlapping_ranges.PushBack( math::Range< double>( 10.0, 11.0));
      if
      (
        BCL_ExampleCheck
        (
          math::RangeSet< double>( overlapping_ranges.Begin(), overlapping_ranges.End()).GetRanges().GetSize(),
          1
        )
      )
      {
        util::ShPtr< math::Range< double> > overlapping_range
        (
          math::RangeSet< double>( overlapping_ranges.Begin(), overlapping_ranges.End()).GetRanges().Begin()->Clone()
        );

        // ensure that the range is correct
        BCL_ExampleIndirectCheck
        (
          *overlapping_range, math::Range< double>( -5.0, 11.0),
          "math::RangeSet< double> ability to combine ranges"
        );
      }

      // only the ranges that overlap should be combined
      storage::Vector< math::Range< double> > two_distinct_ranges( overlapping_ranges);
      two_distinct_ranges.PushBack( math::Range< double>( -70.0, -20.0));
      BCL_ExampleCheck
      (
        math::RangeSet< double>( two_distinct_ranges.Begin(), two_distinct_ranges.End()).GetRanges().GetSize(),
        2
      );

      math::RangeSet< double> n70_n20_n5_p5( two_distinct_ranges.Begin(), two_distinct_ranges.End());
      // check is within
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( -5.0), true);
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( 5.0), true);
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( -5.1), false);
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( -4.9), true);
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( -14.0), false);
      BCL_ExampleCheck( n70_n20_n5_p5.IsWithin( -35.0), true);
      BCL_ExampleCheck( n70_n20_n5_p5.IsEmpty(), false);

      // make sure get members from data labels gives us back a correct string
      BCL_ExampleCheck( n70_n20_n5_p5.GetLabel().ToString(), "\"[-70,-20]+[-5,11]\"");

      // check GetMappedSubset
      math::RangeSet< size_t> from6_to_9_from19_to_80;
      from6_to_9_from19_to_80 += math::Range< size_t>( 6, 9);
      from6_to_9_from19_to_80 += math::Range< size_t>( 19, 80);
      BCL_ExampleCheck
      (
        from6_to_9_from19_to_80.GetMappedSubset( math::Range< size_t>( 5, 15)).GetLabel().ToString(),
        "\"[20,30]\""
      );
      BCL_ExampleCheck
      (
        from6_to_9_from19_to_80.GetMappedSubset( math::Range< size_t>( 2, 10)).GetLabel().ToString(),
        "\"[8,9]+[19,25]\""
      );

      math::RangeSet< size_t> powers_of_two1_2_4_8_16_32_64;
      powers_of_two1_2_4_8_16_32_64.FromString( "[1]+[2]+[4]+[8]+[16]+[32]+[64]", util::GetLogger());
      BCL_ExampleCheck
      (
        powers_of_two1_2_4_8_16_32_64.GetMappedSubset( math::Range< size_t>( 2, 4)).GetLabel().ToString(),
        "\"[4,4]+[8,8]+[16,16]\""
      );

    /////////////////
    // data access //
    /////////////////

      // check from string
      math::RangeSet< double> range_from_string;
      range_from_string.FromString( "[ -70.0, -20.0 ]+[ -5.0, 11.0 ]", util::GetLogger());

      // check that the ranges are equal
      BCL_ExampleCheck
      (
        range_from_string.AsString(),
        n70_n20_n5_p5.AsString()
      );

      BCL_ExampleCheck
      (
        range_from_string.AsString(),
        "[-70,-20]+[-5,11]"
      );

      // check from string
      math::RangeSet< double> range_from_neg_string;
      range_from_neg_string.FromString( "-[ -70.0, -20.0 ]-[ -5.0, 11.0 ]", util::GetLogger());
      BCL_ExampleIndirectCheck
      (
        range_from_neg_string.IsWithin( -25.0),
        false,
        "Constructor from -[ -70.0, -20.0 ]-[ -5.0, 11.0 ]"
      );
      BCL_ExampleIndirectCheck
      (
        range_from_neg_string.IsWithin( -75.0),
        true,
        "Constructor from -[ -70.0, -20.0 ]-[ -5.0, 11.0 ]"
      );

    ////////////////
    // operations //
    ////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathRangeSet

  const ExampleClass::EnumType ExampleMathRangeSet::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRangeSet())
  );

} // namespace bcl
