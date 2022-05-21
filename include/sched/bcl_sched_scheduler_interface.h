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

#ifndef BCL_SCHED_SCHEDULER_INTERFACE_H_
#define BCL_SCHED_SCHEDULER_INTERFACE_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SchedulerInterface
    //! @brief The scheduler interface is created by a class instance to parallelize various jobs spawned in a certain order.
    //! @details Job interfaces are submitted to the scheduler individually.
    //! If a batch of jobs needs to be submitted, the user should run JobSubmit an appropriate number of times.
    //!
    //! @see @link example_sched_jobs_with_data.cpp @endlink
    //! @author riddeljs, woetzen
    //! @date 11.18.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SchedulerInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      virtual SchedulerInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief How many unused cpus are there?
      //! @return a size_t representing the number of unused CPUS
      virtual size_t GetNumberUnusedCPUS() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief reserve available CPUs for jobs. Useful for multi-level threading, because threads should not be
      //!        launched for the inner-most processes.  The reserved CPUs should not be available to other threads
      //! @param MAX maximum number of CPUs to reserve
      //! @return number of CPUs actually reserved (==0 if no CPUs are available)
      virtual size_t ReserveAvailableCPUs( const size_t &MAX) = 0;

      //! @brief release a reservation held by this thread
      //! @param NR_TO_RELEASE the number of previously reserved CPUs to release
      virtual void ReleaseReservedCPUs( const size_t &NR_TO_RELEASE) = 0;

      //! @brief Enqueue a job for the scheduler
      //! @param SP_JOB ShPtr to a job
      virtual void SubmitJob( util::ShPtr< JobInterface> &SP_JOB) const = 0;

      //! @brief Like SubmitJob, but runs the job with the current process if no additional processors are available
      //!        rather than waiting for them to become available.  Has no effect on the serial scheduler
      //! @param SP_JOB ShPtr to a job
      //! @note this functionality is necessary to avoid a race condition between requesting the # of unused processors
      //!       and running a job; since multiple threads may call GetNumberUnusedCPUS at the same time, then
      //!       submit that many jobs, only to have a deadlock when none of them completes
      //!       When such behavior is possible, and the jobs should not just be enqueued (as with SubmitJob),
      //!       call Runjob instead of SubmitJob
      virtual void RunJob( util::ShPtr< JobInterface> &SP_JOB) const = 0;

      //! @brief Join method to call when closing threads
      //! @param SP_JOB job to join
      virtual void Join( util::ShPtr< JobInterface> &SP_JOB) const = 0;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! Parameter for setting the number of cpus (or possible parallel jobs) over the command line
      static util::ShPtr< command::ParameterInterface> &GetParameterNumberCPUs();

    private:

      friend class Schedulers; //!< Must be a friend to access Initialize() and Terminate()

      //! @brief initialize the scheduler
      virtual void Initialize() = 0;

      //! @brief terminate the scheduler; waits for any currently running jobs to finish
      virtual void Terminate() = 0;

    }; // class SchedulerInterface

    //! @brief get currently used Scheduler
    BCL_API SchedulerInterface &GetScheduler();

    //! @brief return the number of cpus as given on the command line
    BCL_API size_t GetNumberCPUs();

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_SCHEDULER_INTERFACE_H_
