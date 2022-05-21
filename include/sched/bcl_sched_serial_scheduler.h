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

#ifndef BCL_SCHED_SERIAL_SCHEDULER_H_
#define BCL_SCHED_SERIAL_SCHEDULER_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_interface.h"
#include "bcl_sched_scheduler_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerialScheduler
    //! @brief for environment which presumably has just one CPU.
    //! @details Jobs are submitted to the scheduler and run on a first come first serve basis.
    //!
    //! @see @link example_sched_serial_scheduler.cpp @endlink
    //! @author riddeljs
    //! @date 1.20.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SerialScheduler :
      public SchedulerInterface
    {
    private:

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! Default constructor
      SerialScheduler()
      {
      }

      //! Copy Constructor
      SerialScheduler *Clone() const
      {
        return new SerialScheduler( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
         return GetStaticClassName( *this);
      }

      // The number of jobs is either 0 or 1, and only one job can execute!
      // Thus, most of these methods aren't implemented.

      //! @brief How many unused cpus are there?
      //! @return a size_t representing the number of unused CPUS
      size_t GetNumberUnusedCPUS() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reserve available CPUs for jobs. Useful for multi-level threading, because threads should not be
      //!        launched for the inner-most processes.  The reserved CPUs should not be available to other threads
      //! @param MAX maximum number of CPUs to reserve
      //! @return number of CPUs actually reserved (==0 if no CPUs are available)
      size_t ReserveAvailableCPUs( const size_t &MAX)
      {
        return 1;
      }

      //! @brief release a reservation held by this thread
      //! @param NR_TO_RELEASE the number of previously reserved CPUs to release
      void ReleaseReservedCPUs( const size_t &NR_TO_RELEASE)
      {
      }

      //! @brief Submit a job to the scheduler
      //! @param SP_JOB ShPtr to a job
      void SubmitJob( util::ShPtr< JobInterface> &SP_JOB) const
      {
        SP_JOB->Run();
      }

      //! @brief Like SubmitJob, but runs the job with the current process if no additional processors are available
      //!        rather than waiting for them to become available.  Has no effect on the serial scheduler
      //! @param SP_JOB ShPtr to a job
      void RunJob( util::ShPtr< JobInterface> &SP_JOB) const
      {
        SP_JOB->Run();
      }

      //! @brief Join method to call when closing threads
      //! @param SP_JOB job to join
      void Join( util::ShPtr< JobInterface> &SP_JOB) const
      {
      }

    protected:

      //! @brief Standard Read and Write
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! @brief initialize the scheduler
      void Initialize()
      {
      }

      //! @brief terminate the scheduler; waits for any currently running jobs to finish
      void Terminate()
      {
      }

    }; // class SerialScheduler

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_SERIAL_SCHEDULER_H_
