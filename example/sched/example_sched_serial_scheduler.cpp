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
#include "sched/bcl_sched_serial_scheduler.h"

// includes from bcl - sorted alphabetically
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_time.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sched_serial_scheduler.cpp
  //!
  //! @author nobody
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSchedSerialScheduler :
    public ExampleInterface
  {
  public:

    ExampleSchedSerialScheduler *Clone() const
    {
      return new ExampleSchedSerialScheduler( *this);
    }

    class OutputNumber
    {

    public:

    ///////////////
    // operations //
    ///////////////

      //! @brief function taking a size_t and an int.  Wait and cout a string.
      //! @param JOBNUMBER the job being ran
      //! @param SECONDS sleep for this many seconds
      void OutputAndSleep( const size_t &JOBNUMBER, size_t &SECONDS)
      {
        BCL_MessageStd( "Sequential Job number " + util::Format()( JOBNUMBER) + " running!");
        SECONDS++;
        util::Time::Delay( util::Time( 0, SECONDS * 50000));
      }

    }; // class OutputNumber

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      //initialize the serial scheduler
      sched::SerialScheduler serial_scheduler;

      //class object
      OutputNumber new_instance;

      // not really used for this example
      const size_t group_id( 1);

      size_t wait_time( 0);

    ////////////////////////////////////////////////////////////
    // sequential scheduler:  run all five jobs consecutively //
    ////////////////////////////////////////////////////////////

      BCL_MessageStd( "Will now execute 5 jobs consecutively");

      const size_t s_numbers[] = { 0, 1, 2, 3, 4};
      util::ShPtrVector< sched::JobInterface> schedule;
      // create 5 jobs to schedule
      for( size_t loop_counter = 0; loop_counter < 5; loop_counter++)
      {
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::BinaryFunctionJobWithData< const size_t, size_t, void, OutputNumber>
            (
              group_id, new_instance,
              &OutputNumber::OutputAndSleep,
              s_numbers[ loop_counter], wait_time,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
      }

      //submit all jobs to the scheduler
      for( size_t loop_counter = 0; loop_counter < 5; loop_counter++)
      {
        serial_scheduler.SubmitJob( schedule( loop_counter));
      }

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        BCL_MessageStd( "this is the output of SerialScheduler");
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType ExampleSchedSerialScheduler_Instance;

  }; //end ExampleSchedSerialScheduler

  const ExampleClass::EnumType ExampleSchedSerialScheduler::ExampleSchedSerialScheduler_Instance
  (
    GetExamples().AddEnum( ExampleSchedSerialScheduler())
  );

} // namespace bcl

