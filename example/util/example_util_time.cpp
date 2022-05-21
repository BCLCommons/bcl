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
#include "util/bcl_util_time.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_time.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilTime :
    public ExampleInterface
  {
  public:

    ExampleUtilTime *Clone() const
    { return new ExampleUtilTime( *this);}

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
    //////////////////
    // constructors //
    //////////////////

      // default constructor - initializes everything to 0
      util::Time time_zero;
      BCL_MessageStd( "default initialized time to 0: " + util::Format()( time_zero));

      // check that seconds and micro seconds are initialized to 0
      BCL_ExampleIndirectCheck( time_zero.GetSeconds(), size_t( 0), "zero time constructor");
      BCL_ExampleIndirectCheck( time_zero.GetMicroSeconds(), size_t( 0), "zero time constructor");

      // initialize to curren time
      util::Time time_current( util::Time::GetCurrent());
      BCL_MessageStd
      (
        "time initialized to current system time: " + time_current.GetTimeAsDate()
      );

      // initialize to 5 seconds and 400 microseconds
      util::Time time_5_400( 5, 400);
      BCL_MessageStd
      (
        "time initialized from seconds and micro seconds: " + time_current.GetTimeAsHourMinuteSecondMilliSeconds()
      );

      // check that if there are more than 1000000 microseconds, that they are converted to seconds
      BCL_ExampleCheck( util::Time( 4, 1000400), time_5_400);

      // check that time in milliseconds is correct
      BCL_ExampleCheck( time_5_400.GetTotalMilliseconds(), 5000);

      // initialize from compiler macro
      util::Time time_compilation( util::Time::CreateTimeFromCompilerMacro( __DATE__, __TIME__));
      BCL_MessageStd
      (
        "time initialized from compilation time: " + time_compilation.GetTimeAsDate()
      );

      // initialize from number of hours
      util::Time time_from_hours( util::Time::CreateTimeFromHours( 2));
      BCL_ExampleCheck( time_from_hours.GetSeconds(), 7200);

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd
      (
        "this is the      seconds part of the current time: " + util::Format()( time_current.GetSeconds())
      );
      BCL_MessageStd
      (
        "this is the microseconds part of the current time: " + util::Format()( time_current.GetMicroSeconds())
      );

      BCL_ExampleCheck( time_5_400.GetSeconds() == 5 && time_5_400.GetMicroSeconds() == 400, true);

    ////////////////
    // operations //
    ////////////////

      // time has three functions to retrieve the current time
      const std::string time_current_date( util::Time::GetCurrent().GetTimeAsDate());
      const std::string time_current_ddhhmmss( util::Time::GetCurrent().GetTimeAsDayHourMinuteSecond());
      const std::string time_current_hhmmss( util::Time::GetCurrent().GetTimeAsHourMinuteSecond());

      // write current time as time format
      BCL_MessageStd( "this is the current time in normal format: " + util::Format()( time_current));

      // write current time as a date format
      BCL_MessageStd( "this is the current time in readable format: " + time_current_date);

      // write current time as hh:mm:ss format
      BCL_MessageStd( "this is the current time in normal format: " + util::Format()( time_current_hhmmss));

      // write current time as dd hh:mm:ss format
      BCL_MessageStd( "this is the current time in normal format: " + util::Format()( time_current_ddhhmmss));

      // one can set a time object always to the current time
      time_zero.SetToCurrentTime();

      // delay the executable by 1 second
      util::Time::Delay( util::Time( 1, 0));
      const util::Time elapsed_time( util::Time::GetCurrent() - time_zero);
      BCL_MessageStd( "delayed program by 1000 milli seconds, actual delay: " + util::Format()( elapsed_time));

      // check that at least one second did elapse
      BCL_ExampleIndirectCheck
      (
        elapsed_time >= util::Time( 0, util::Time::s_MicroSecondsPerSecond - 10 * util::Time::s_MicroSecondsPerMiliSecond), // add tolerance due to mingw inaccuracy in delay
        true,
        "delay should be at least one second but was " + util::Format()( util::Time::GetCurrent() - time_zero)
      );

      // check zero functions
      util::Time zero_time;
      BCL_ExampleCheck( zero_time.IsZero(), true);

      zero_time.SetToCurrentTime();
      BCL_ExampleIndirectCheck( zero_time.IsZero(), false, "time should have been set to the current time, but is still zero");

      zero_time.SetZero();
      BCL_ExampleIndirectCheck( zero_time.IsZero(), true, "time should be zero, after setting to zero");

    ///////////////
    // operators //
    ///////////////

      // current time in hh:mm:ss
      const std::string time_current_hms( time_zero.GetTimeAsHourMinuteSecond());
      BCL_MessageStd( "set default initialized time to current time");
      BCL_MessageStd( "this is the current time in hour::minutes::seconds format: " + time_current_hms);
      // we can verify that by making sure that the most recent current times are close enough
      const size_t seconds_difference( time_zero.GetSeconds() - time_current.GetSeconds());
      BCL_ExampleIndirectCheck
      (
        seconds_difference < 3, true,
        "the two current times deviate by more than 2 second: " + util::Format()( seconds_difference)
      );

      // initialize two different times
      const util::Time time_5_500( 5, 500), time_6_400( 6, 400);
      util::Time time_sum;

      // adding
      time_sum = time_5_500 + time_6_400;
      BCL_MessageStd
      (
        "sum: " + util::Format()( time_5_500) + " + " + util::Format()( time_6_400) + " = " + util::Format()( time_sum)
      );
      BCL_ExampleIndirectCheck
      (
           time_sum.GetSeconds() == time_5_500.GetSeconds() + time_6_400.GetSeconds()
        && time_sum.GetMicroSeconds() == time_5_500.GetMicroSeconds() + time_6_400.GetMicroSeconds(),
        true,
        "sum is incorrect: " + time_5_500.GetTimeAsDayHourMinuteSecond() + " + "
        + time_6_400.GetTimeAsDayHourMinuteSecond() + " != " + time_sum.GetTimeAsDayHourMinuteSecond()
      );

      // subtraction
      time_sum -= time_6_400;
      BCL_ExampleIndirectCheck
      (
        time_sum, time_5_500,
        "result is incorrect after -=: should be " + time_5_500.GetTimeAsDayHourMinuteSecond() + " but is " +
        time_sum.GetTimeAsDayHourMinuteSecond()
      );

      // addition
      time_sum += time_6_400;
      BCL_ExampleIndirectCheck( time_sum, time_5_500 + time_6_400, "operator +=");

      // comparison
      // equal
      BCL_ExampleIndirectCheck( time_5_500 == time_5_500, true, "times should be equal");

      // not equal
      BCL_ExampleIndirectCheck( time_5_500 != time_6_400, true, "times should be not equal");

      // smaller
      BCL_ExampleIndirectCheck( time_5_500 < time_6_400, true, "time should be smaller");

      // smaller equal
      BCL_ExampleIndirectCheck( time_5_500 <= time_5_500 && time_5_500 <= time_6_400, true, "time should be smaller equal");

      // larger
      BCL_ExampleIndirectCheck( time_6_400 > time_5_500, true, "time should be larger");

      // larger equal
      BCL_ExampleIndirectCheck( time_6_400 >= time_6_400 && time_6_400 >= time_5_500, true, "time should be larger equal");

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( time_6_400);
      // read from file
      util::Time read_time;
      ReadBCLObject( read_time);

      BCL_ExampleIndirectCheck( time_6_400, read_time, "written time to read time");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilTime

  const ExampleClass::EnumType ExampleUtilTime::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilTime())
  );

} // namespace bcl
