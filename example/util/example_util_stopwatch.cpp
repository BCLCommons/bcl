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
#include "util/bcl_util_stopwatch.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_stopwatch.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilStopwatch :
    public ExampleInterface
  {
  public:

    ExampleUtilStopwatch *Clone() const
    {
      return new ExampleUtilStopwatch( *this);
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

    //! @return true if get static class name is what it should be
    bool CheckGetStaticClassName() const
    {
      return ( GetStaticClassName( *this) == "bcl::ExampleUtilStopwatch");
    }

    //! @return true if get class identifier is what it should be
    bool CheckGetClassIdentifier() const
    {
      return ( GetClassIdentifier() == "bcl::ExampleUtilStopwatch");
    }

    int Run() const
    {
    /////////////////
    // constructor //
    /////////////////

      //initialize a stopwatch object
      util::Stopwatch timer( "running " + GetStaticClassName( *this), util::Time(), util::Message::e_Standard);
      timer.Reset();
      BCL_MessageStd( "this is the stopped timer: " + util::Format()( timer));

      BCL_Example_Check
      (
        timer.GetProcessDuration().IsZero(),
        "timer should not have started yet"
      );

      timer.Start();
      BCL_Example_Check
      (
        timer.GetDescription() == ( "running " + GetStaticClassName( *this)),
        "timer should have this description: running "
        + GetStaticClassName( *this) + " but instead is described as " + timer.GetDescription()
      );

      timer.Stop();

      const util::Time total_duration_so_far( timer.GetTotalTime());
      util::Time::Delay( util::Time( 0, 1000));

      BCL_Example_Check
      (
        total_duration_so_far == timer.GetTotalTime(),
        "timer should have stopped, but total time progressed"
      );

      // copy the timer
      util::Stopwatch *ptr( timer.Clone());

      BCL_Example_Check
      (
        timer.GetDescription() == ptr->GetDescription()
        && timer.GetTotalTime() == ptr->GetTotalTime(),
        "timer copy is different: " + util::Format()( timer) + " != " + util::Format()( *ptr)
      );

    /////////////////
    // data access //
    /////////////////

      // ask for static class name
      BCL_MessageStd
      (
        "this is the static class name of a Stopwatch: " + GetStaticClassName< util::Stopwatch>()
      );

      // ask for class identifier
      BCL_MessageStd
      (
        "this is the class identifier of a Stopwatch: " + ptr->GetClassIdentifier()
      );
      BCL_Example_Check
      (
        GetStaticClassName< util::Stopwatch>() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

      // ask the timer for its lifetime (current time - stopwatch starting time)
      BCL_MessageStd
      (
        "this is the process duration (lifetime of timer): " + util::Format()( timer.GetProcessDuration())
      );

      // ask the timer for its description
      BCL_MessageStd
      (
        "this is the description: " + timer.GetDescription()
      );

    ////////////////
    // operations //
    ////////////////

      //to time something, you can initialize the util::Stopwatch
      //it will save the beginning time, the destructor can report the duration time
      {
        util::Stopwatch timer_del( "deleting and constructing objects");
        BCL_MessageStd( " create and delete 5 times 1'000'000 times a pointer to a double");
        for( size_t loop( 0); loop < 5; ++loop)
        {
          for( size_t i( 1); i < 1000000; ++i)
          {
            double *ptr_double = new double( 5);
            ( *ptr_double)++;

            delete ptr_double;
          }

          //report the current duration time
          BCL_MessageStd( "loop " + util::Format()( loop + 1) + ": " + util::Format()( timer_del));
        }
      }

      // to time a function that will be called many times, just call start and stop with a static stopwatch
      // at the beginning and end of each function call of interest.
      // You will be automatically notified whenever the timer is stopped if enough time accumulated on the stopwatch
      // Specifically, the second parameter of the cumulative stopwatch's constructor is how many microseconds must
      // have passed before a message is output.  This makes Stopwatch ideal for timing functions that are
      // repeatedly called
      util::Stopwatch static_name_timer( "running GetStaticClassname", util::Time( 0, 10000), util::Message::e_Standard, false);
      util::Stopwatch name_timer( "running GetClassDescription", util::Time( 0, 10000), util::Message::e_Standard, false);
      for( size_t i = 0; i < 1000000; ++i)
      {
        static_name_timer.Start();
        CheckGetStaticClassName();
        static_name_timer.Stop();
        name_timer.Start();
        CheckGetClassIdentifier();
        name_timer.Stop();
      }

      // reset the timer
      timer.Reset();
      BCL_MessageStd( "reset timer: " + util::Format()( timer));
      BCL_Example_Check
      (
        timer.GetProcessDuration().GetSeconds() == 0, "after resetting timer, the seconds should be 0"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( timer);
      // read from file
      util::Stopwatch read_timer;
      ReadBCLObject( read_timer);

      BCL_Example_Check
      (
        timer.GetDescription() == read_timer.GetDescription()
        && timer.GetLastStartTime() == read_timer.GetLastStartTime(),
        "written time is not equal to read time " + util::Format()( timer) + " != " + util::Format()( read_timer)
      );

      // delete
      delete ptr;

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilStopwatch

  const ExampleClass::EnumType ExampleUtilStopwatch::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilStopwatch())
  );

} // namespace bcl
