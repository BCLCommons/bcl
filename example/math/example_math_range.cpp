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
#include "math/bcl_math_range.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_range.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRange :
    public ExampleInterface
  {
  public:

    ExampleMathRange *Clone() const
    {
      return new ExampleMathRange( *this);
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
      math::Range< double> range_default;
      BCL_Example_Check
      (
        !util::IsDefined( range_default.GetMin()) && !util::IsDefined( range_default.GetMax()),
        "default constructed range should be undefined"
      );

      // construct both side closed interval from min and max
      math::Range< double> range_cm5_p5c( double( -5), double( 5));
      BCL_Example_Check
      (
        range_cm5_p5c.GetLeftCondition() == math::RangeBorders::e_LeftClosed
        && range_cm5_p5c.GetMin() == double( -5)
        && range_cm5_p5c.GetMax() == double( 5)
        && range_cm5_p5c.GetRightCondition() == math::RangeBorders::e_RightClosed,
        "range [-5,5] was not constructed properly: " + util::Format()( range_cm5_p5c)
      );

      // construct both side open interval from min and max
      math::Range< double> range_om5_p5o( math::RangeBorders::e_LeftOpen, double( -5), double( 5), math::RangeBorders::e_RightOpen);
      BCL_Example_Check
      (
        range_om5_p5o.GetLeftCondition() == math::RangeBorders::e_LeftOpen
        && range_om5_p5o.GetMin() == double( -5)
        && range_om5_p5o.GetMax() == double( 5)
        && range_om5_p5o.GetRightCondition() == math::RangeBorders::e_RightOpen,
        "range (-5,5) was not constructed properly: " + util::Format()( range_om5_p5o)
      );

      // copy constructor
      math::Range< double> range_copy( range_cm5_p5c);
      BCL_Example_Check
      (
        range_cm5_p5c.GetLeftCondition() == range_copy.GetLeftCondition()
        && range_cm5_p5c.GetMin() == range_copy.GetMin()
        && range_cm5_p5c.GetMax() == range_copy.GetMax()
        && range_cm5_p5c.GetRightCondition() == range_copy.GetRightCondition(),
        "range [-5,5] was not copied properly: " + util::Format()( range_cm5_p5c) + " != " + util::Format()( range_copy)
      );

      // clone
      math::Range< double> *ptr( range_copy.Clone());
      BCL_Example_Check
      (
        range_cm5_p5c.GetLeftCondition() == ptr->GetLeftCondition()
        && range_cm5_p5c.GetMin() == ptr->GetMin()
        && range_cm5_p5c.GetMax() == ptr->GetMax()
        && range_cm5_p5c.GetRightCondition() == ptr->GetRightCondition()
        && range_cm5_p5c.GetClassIdentifier() == ptr->GetClassIdentifier(),
        "range [-5,5] was not cloned properly: " + util::Format()( range_cm5_p5c) + " != " + util::Format()( *ptr)
      );

    /////////////////
    // data access //
    /////////////////

      BCL_Example_Check
      (
        GetStaticClassName< math::Range< double> >() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

      // left border condition
      BCL_MessageStd( "left border of range [-5,5]: " + std::string( 1, math::RangeBorders::GetConditionLeftChars()[ range_cm5_p5c.GetLeftCondition()]));
      // min
      BCL_MessageStd( "minimum in range [-5,5]: " + util::Format()( range_cm5_p5c.GetMin()));
      // max
      BCL_MessageStd( "maximum in range [-5,5]: " + util::Format()( range_cm5_p5c.GetMax()));
      // right border condition
      BCL_MessageStd( "right border of range [-5,5]: " + std::string( 1, math::RangeBorders::GetConditionRightChars()[ range_cm5_p5c.GetRightCondition()]));
      // width
      BCL_MessageStd( "width of range [-5,5]: " + util::Format()( range_cm5_p5c.GetWidth()));
      BCL_Example_Check
      (
        range_cm5_p5c.GetWidth() == double( 10),
        "range [-5,5] did not have width 10"
      )

    ////////////////
    // operations //
    ////////////////

      // rescale value from one scale to another
      math::Range< double> range_m2_m1( double( -2), double( -1));
      const double two_range_m5_p5( 2);
      const double rescaled_two_m5_p5_m2_m1_expected( -1.3);
      const double rescaled_two_m5_p5_m2_m1_computed( range_m2_m1.Rescale( two_range_m5_p5, range_cm5_p5c));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rescaled_two_m5_p5_m2_m1_expected, rescaled_two_m5_p5_m2_m1_computed),
        "rescaled value " + util::Format()( two_range_m5_p5) + " in range [-5,5] is not " +
        util::Format()( rescaled_two_m5_p5_m2_m1_expected) + " in range [-2,-1] but: " +
        util::Format()( rescaled_two_m5_p5_m2_m1_computed)
      );

      // check if value is within range
      BCL_Example_Check
      (
        !range_cm5_p5c.IsWithin( -5.1)
        &&  range_cm5_p5c.IsWithin( -5)
        &&  range_cm5_p5c.IsWithin( -3)
        &&  range_cm5_p5c.IsWithin( 0)
        &&  range_cm5_p5c.IsWithin( 5)
        && !range_cm5_p5c.IsWithin( 5.1),
        "IsWithin( VALUE) does not work properly for -5.1, -5, -3, 0, 5, 5.1 in range [-5,5]"
      );

      // check if value is within range
      BCL_Example_Check
      (
        !range_om5_p5o.IsWithin( -5.1)
        && !range_om5_p5o.IsWithin( -5)
        &&  range_om5_p5o.IsWithin( -3)
        &&  range_om5_p5o.IsWithin( 0)
        && !range_om5_p5o.IsWithin( 5)
        && !range_om5_p5o.IsWithin( 5.1),
        "IsWithin( VALUE) does not work properly for -5.1, -5, -3, 0, 5, 5.1 in range (-5,5)"
      );

      // create ranges to check if one range overlaps with another
      const math::Range< double> example_a( math::RangeBorders::e_LeftOpen,   -15.0,  2.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_b( math::RangeBorders::e_LeftOpen,   -15.0, -6.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_c( math::RangeBorders::e_LeftOpen,   -15.0, -5.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_d( math::RangeBorders::e_LeftOpen,     5.0, 10.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_e( math::RangeBorders::e_LeftOpen,    -4.0, 10.2, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_f( math::RangeBorders::e_LeftOpen,     5.0,  5.1, math::RangeBorders::e_RightOpen);
      const math::Range< double> example_g( math::RangeBorders::e_LeftClosed,  -5.0, -5.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_h( math::RangeBorders::e_LeftClosed, -15.0, -6.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_i( math::RangeBorders::e_LeftClosed, -15.0, -5.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_j( math::RangeBorders::e_LeftClosed,   5.0, 10.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_k( math::RangeBorders::e_LeftClosed,  -4.0, 10.2, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_l( math::RangeBorders::e_LeftClosed,   6.0, 10.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_m( math::RangeBorders::e_LeftClosed,  -6.0, 10.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_n( math::RangeBorders::e_LeftClosed,  -5.0, 10.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> example_o( math::RangeBorders::e_LeftClosed,  -6.0,  5.0, math::RangeBorders::e_RightClosed);

      // check ranges to see if they overlap in [-5,5] by passing example a - o
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_a),
            "DoesOverlap( VALUE) does not work properly for example a in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_b),
            "DoesOverlap( VALUE) does not work properly for example b in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_c),
            "DoesOverlap( VALUE) does not work properly for example c in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_d),
            "DoesOverlap( VALUE) does not work properly for example d in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_e),
            "DoesOverlap( VALUE) does not work properly for example e in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_f),
            "DoesOverlap( VALUE) does not work properly for example f in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_g),
            "DoesOverlap( VALUE) does not work properly for example g in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_h),
            "DoesOverlap( VALUE) does not work properly for example h in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_i),
            "DoesOverlap( VALUE) does not work properly for example i in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_j),
            "DoesOverlap( VALUE) does not work properly for example j in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_k),
            "DoesOverlap( VALUE) does not work properly for example k in range [-5,5]"
      );
      BCL_Example_Check
      (
        !range_cm5_p5c.DoesOverlap( example_l),
            "DoesOverlap( VALUE) does not work properly for example l in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_m),
            "DoesOverlap( VALUE) does not work properly for example m in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_n),
            "DoesOverlap( VALUE) does not work properly for example n in range [-5,5]"
      );
      BCL_Example_Check
      (
        range_cm5_p5c.DoesOverlap( example_o),
            "DoesOverlap( VALUE) does not work properly for example o in range [-5,5]"
      );

      // check ranges to see if they overlap in (-5,5) by passing example a - o
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_a),
            "DoesOverlap( VALUE) does not work properly for example a in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_b),
            "DoesOverlap( VALUE) does not work properly for example b in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_c),
            "DoesOverlap( VALUE) does not work properly for example c in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_d),
            "DoesOverlap( VALUE) does not work properly for example d in range (-5,5)"
      );
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_e),
            "DoesOverlap( VALUE) does not work properly for example e in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_f),
            "DoesOverlap( VALUE) does not work properly for example f in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_g),
            "DoesOverlap( VALUE) does not work properly for example g in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_h),
            "DoesOverlap( VALUE) does not work properly for example h in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_i),
            "DoesOverlap( VALUE) does not work properly for example i in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_j),
            "DoesOverlap( VALUE) does not work properly for example j in range (-5,5)"
      );
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_k),
            "DoesOverlap( VALUE) does not work properly for example k in range (-5,5)"
      );
      BCL_Example_Check
      (
        !range_om5_p5o.DoesOverlap( example_l),
            "DoesOverlap( VALUE) does not work properly for example l in range (-5,5)"
      );
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_m),
            "DoesOverlap( VALUE) does not work properly for example m in range (-5,5)"
      );
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_n),
            "DoesOverlap( VALUE) does not work properly for example n in range (-5,5)"
      );
      BCL_Example_Check
      (
        range_om5_p5o.DoesOverlap( example_o),
            "DoesOverlap( VALUE) does not work properly for example o in range (-5,5)"
      );

      // test comparison operator <
      BCL_Example_Check
      (
        example_a < example_d,
        "(-15,2) < (5,10) returned false"
      );
      BCL_Example_Check
      (
        example_c < example_g,
        "(-15,-5) < [-5,-5] returned false"
      );
      BCL_Example_Check
      (
        !( example_a < example_b),
        "(-15,-2) < (-15,-6) returned true"
      );
      BCL_Example_Check
      (
        !( example_b < example_a),
        "(-15,-6) < (-15,-2) returned true"
      );
      BCL_Example_Check
      (
        !( example_o < example_j),
        "[-6,-5] < [5,10] returned true"
      );

      // check that we can read and write out the range using GetString and FromStream
      std::istringstream range_om5_p5o_stream( range_om5_p5o.GetString());
      math::Range< double> range_om5_p5o_from_stream;
      range_om5_p5o_from_stream.FromStream( range_om5_p5o_stream, util::GetLogger());

      BCL_Example_Check
      (
        range_om5_p5o_from_stream.GetMin() == range_om5_p5o.GetMin()
        && range_om5_p5o_from_stream.GetMax() == range_om5_p5o.GetMax()
        && range_om5_p5o_from_stream.GetLeftCondition() == range_om5_p5o.GetLeftCondition()
        && range_om5_p5o_from_stream.GetRightCondition() == range_om5_p5o.GetRightCondition(),
        "range written to GetString() was not read properly: "
        + util::Format()( range_om5_p5o_from_stream)
        + " != " + util::Format()( range_om5_p5o)
      );

      // Check combine range
      BCL_Example_Check
      (
        math::Range< double>::CombineRanges( example_b, example_a) == example_a,
        "(-15,-6) combined with (-15,-2) should return (-15,-2), instead returned " +
        math::Range< double>::CombineRanges( example_b, example_a).GetString()
      );
      BCL_Example_Check
      (
        math::Range< double>::CombineRanges( example_a, range_om5_p5o)
          == math::Range< double>( math::RangeBorders::e_LeftOpen, -15.0, 5.0, math::RangeBorders::e_RightOpen),
        "(-15,-2) combined with (-5,5) should return (-15,5), instead returned " +
        math::Range< double>::CombineRanges( example_a, range_om5_p5o).GetString()
      );
      BCL_Example_Check
      (
        math::Range< double>::CombineRanges( example_a, range_cm5_p5c)
          == math::Range< double>( math::RangeBorders::e_LeftOpen, -15.0, 5.0, math::RangeBorders::e_RightClosed),
        "(-15,-2) combined with [-5,5] should return (-15,5], instead returned " +
        math::Range< double>::CombineRanges( example_a, range_cm5_p5c).GetString()
      );

    ///////////////
    // operators //
    ///////////////

      // convert a number from one range to another range
      // the operator is using the rescale function, so we just have to compare the results
      BCL_ExampleCheck
      (
         range_m2_m1.Rescale( two_range_m5_p5, range_cm5_p5c),
         range_m2_m1( two_range_m5_p5, range_cm5_p5c)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write range to file
      WriteBCLObject( range_m2_m1);
      // read range
      ReadBCLObject( range_default);

      BCL_Example_Check
      (
        range_m2_m1.GetMin() == range_default.GetMin()
        && range_m2_m1.GetMax() == range_default.GetMax(),
        "range was not read properly: " + util::Format()( range_m2_m1) + " != " + util::Format()( range_default)
      );

      // clean up
      delete ptr;

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathRange

  const ExampleClass::EnumType ExampleMathRange::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRange())
  );

} // namespace bcl
