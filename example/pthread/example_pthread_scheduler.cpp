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
//#include "pthread/bcl_pthread_scheduler.cpp"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_parameter.h"
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_mutex.h"
#include "sched/bcl_sched_schedulers.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pthread_scheduler.cpp
  //!
  //! @author riddeljs, woetzen, mendenjl
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePthreadScheduler :
    public ExampleInterface
  {
  public:

    ExamplePthreadScheduler *Clone() const
    {
      return new ExamplePthreadScheduler( *this);
    }

    static const size_t s_NumberJobs = 5;
    static storage::Vector< sched::Mutex> s_Mutexes;

    class OutputNumber :
      public util::ObjectInterface
    {
    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      OutputNumber *Clone() const
      {
        return new OutputNumber( *this);
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

    ///////////////
    // operations //
    ///////////////

      //! @brief query taking a size_t and an int.  Wait and cout a string.
      //! @param JOBNUMBER the job being ran
      //! @param SECONDS sleep for this many seconds
      void OutputAndSleepConsecutive( size_t &JOBNUMBER, int &SECONDS) const
      {
        BCL_MessageStd( "Sequential Job number " + util::Format()( JOBNUMBER) + " running!");
        SECONDS++;
        util::Time::Delay( util::Time( 0, SECONDS * 100000));

        // for all jobs but the last, unlock the next mutex to allow the next job to run
        if( JOBNUMBER != s_NumberJobs - 1)
        {
          s_Mutexes( JOBNUMBER + 1).Unlock();
        }
      } // OutputAndSleep

      //! @brief function taking a size_t and an int.  Wait and cout a string.
      //! @param JOBNUMBER the job being ran
      //! @param SECONDS sleep for this many seconds
      void OutputAndSleepConcurrent( size_t &JOBNUMBER, int &SECONDS)
      {
        BCL_MessageStd( "Concurrent Job number " + util::Format()( JOBNUMBER) + " running!");
        SECONDS++;
        util::Time::Delay( util::Time( 0, SECONDS * 100000));
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM << "double to int conversion function job" << '\n';
      }

    }; // class OutputNumber

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // ensure that pthreads were actually linked, otherwise there is nothing for this example to test
      if( !sched::GetSchedulers().HaveEnumWithName( "PThread"))
      {
        BCL_MessageCrt
        (
          "Could not test " + GetClassIdentifier() + " because pthreads are not available"
        );
        return 0;
      }

      //initialize schedule to have s_NumberJobs CPUS
      // set the number of cpus to s_NumberJobs for this example
      const storage::Vector< std::string> original_scheduler
      (
        sched::GetSchedulers().GetFlagSchedulerCPUS()->GetStringList()
      );
      sched::GetSchedulers().GetFlagSchedulerCPUS()->ResetFlag();
      sched::GetSchedulers().GetFlagSchedulerCPUS()->ReadFromList
      (
        storage::Vector< std::string>::Create( "PThread", util::Format()( s_NumberJobs)),
        util::GetLogger()
      );
      sched::GetSchedulers().UpdateCurrentSchedulerFromCommandLineFlag();
      sched::SchedulerInterface &scheduler( sched::GetScheduler());

      //class object
      OutputNumber output_number;

      // not really used for this example
      const size_t group_id( 1);

      int wait_time( 0);

    ////////////////////////////////////////////////////////////
    // sequential scheduler:  run all five jobs consecutively //
    ////////////////////////////////////////////////////////////

      BCL_MessageStd
      (
        "Will now execute " + util::Format()( s_NumberJobs) + " jobs consecutively"
      );

      storage::Vector< size_t> ids( s_NumberJobs);

      util::ShPtrVector< sched::JobInterface> schedule;
      // create s_NumberJobs jobs to schedule
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        ids( loop_counter) = loop_counter;
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::BinaryFunctionJobWithData< size_t, int, void, OutputNumber>
            (
              group_id,
              output_number,
              &OutputNumber::OutputAndSleepConsecutive,
              ids( loop_counter),
              wait_time,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
      }

      // lock the last four mutexes
      for( size_t loop_counter = 1; loop_counter < s_NumberJobs; loop_counter++)
      {
        s_Mutexes( loop_counter).Lock();
      }

      // submit all jobs to the scheduler
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        scheduler.SubmitJob( schedule( loop_counter));
      }

      // join all jobs ( need to be finished)
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        scheduler.Join( schedule( loop_counter));
      }

      BCL_MessageStd
      (
        "Will now execute " + util::Format()( s_NumberJobs) + " jobs concurrently"
      );

      wait_time = 0;

      schedule.Reset();
      // create s_NumberJobs jobs to schedule
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::BinaryFunctionJobWithData< size_t, int, void, OutputNumber>
            (
              group_id,
              output_number,
              &OutputNumber::OutputAndSleepConcurrent,
              ids( loop_counter),
              wait_time,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
      }

      //submit all jobs to the scheduler
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        scheduler.SubmitJob( schedule( loop_counter));
      }

      // join all jobs ( need to be finished)
      for( size_t loop_counter = 0; loop_counter < s_NumberJobs; loop_counter++)
      {
        scheduler.Join( schedule( loop_counter));
      }

      BCL_MessageStd( "Value of variable wait_time: " + util::Format()( wait_time));

    //////////////////////
    // helper functions //
    //////////////////////

      // reset the number of cpus to 1 after this example
      sched::GetSchedulers().GetFlagSchedulerCPUS()->ResetFlag();
      sched::GetSchedulers().GetFlagSchedulerCPUS()->ReadFromList
      (
        original_scheduler,
        util::GetLogger()
      );
      sched::GetSchedulers().UpdateCurrentSchedulerFromCommandLineFlag();

      return 0;
    } // Run

    static const ExampleClass::EnumType ExamplePthreadScheduler_Instance;

  }; //end ExamplePthreadScheduler

  const size_t ExamplePthreadScheduler::s_NumberJobs;

  storage::Vector< sched::Mutex> ExamplePthreadScheduler::s_Mutexes( ExamplePthreadScheduler::s_NumberJobs);

  const ExampleClass::EnumType ExamplePthreadScheduler::ExamplePthreadScheduler_Instance
  (
    GetExamples().AddEnum( ExamplePthreadScheduler())
  );

} // namespace bcl
