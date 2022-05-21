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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class

// includes from bcl - sorted alphabetically
#include "pthread/bcl_pthread.h"
#include "sched/bcl_sched_job_interface.h"
#include "sched/bcl_sched_schedulers.h"

// external includes - sorted alphabetically
#include <map>
#include <errno.h>
#include <pthread.h>

namespace bcl
{
  namespace pthread
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ThreadQueue
    //! @brief a singleton, class that manages a job and thread queue
    //!        Threads pull jobs from the queue; enqueueing a job launches an addition thread, if additional
    //!        CPUs were specified over the command line
    //!
    //! @author mendenjl
    //! @remarks example unnecessary
    //! @date Jun 14, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ThreadQueue
    {
    private:

    //////////
    // data //
    //////////

      //! List of running threads
      std::vector< pthread_t> m_Threads;

      //! Count of number of threads currently running (faster than calling m_RunningThreads.size())
      size_t                  m_NumberThreadsTotal;    //!< sched::GetNumberCPUs() == m_Threads.size()
      size_t                  m_NumberThreadsWaiting;  //!< # of threads that are waiting on jobs
      size_t                  m_NumberThreadsReserved; //!< # of threads reserved

      //! Jobs available for running
      std::list< util::SiPtr< sched::JobInterface> > m_AvailableJobs;

      //! Map of jobs (currently running or available) to condition variables, and bool that indicates whether the
      // condition was already signaled
      std::map
      <
        util::SiPtr< const sched::JobInterface>,
        std::pair< pthread_cond_t, bool>
      > m_JobsToConditionVariables;

      //! @brief Default mutex and cv, primarily for number_of_jobs
      // since the pthread routines must take regular pointers as inputs, they are regular pointers here.
      pthread_mutex_t m_ThreadMutex;     //!< Access to m_Threads, m_NumberThreadsTotal
      pthread_mutex_t m_QueueMutex;      //!< Access to m_AvailableJobs, m_NumberThreadsWaiting
      pthread_mutex_t m_ConditionMutex;  //!< Access to m_JobsToConditionVariables
      pthread_cond_t  m_NewJobCondition; //!< Condition that is signaled whenever new jobs are available
      bool            m_Terminate;       //!< Threads should terminate rather than waiting for another job

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ThreadQueue();

      //! @brief destructor
      ~ThreadQueue();

    public:

      static ThreadQueue s_ThreadQueue; //!< Sole instance of this class

    /////////////////
    // data access //
    /////////////////

      //! @brief Get # of threads currently available
      size_t GetNumberThreadsAvailable();

      //! @brief reserve available CPUs for jobs. Useful for multi-level threading, because threads should not be
      //!        launched for the inner-most processes.  The reserved CPUs should not be available to other threads
      //! @param MAX maximum number of CPUs to reserve
      //! @return number of CPUs actually reserved (==0 if no CPUs are available)
      size_t ReserveAvailableCPUs( const size_t &MAX);

      //! @brief release a reservation held by this thread
      //! @param NR_TO_RELEASE the number of previously reserved CPUs to release
      void ReleaseReservedCPUs( const size_t &NR_TO_RELEASE);

    ////////////////
    // operations //
    ////////////////

      //! @brief join the job; really just wait for it to complete
      //! @param JOB the job to wait for
      void Join( const util::ShPtr< sched::JobInterface> &JOB);

      //! @brief enqueue a job
      //! @param JOB the job to enqueue
      //! @param LOCKED_THREADS true if this process has already locked m_ThreadMutex
      void EnqueueJob( util::ShPtr< sched::JobInterface> &JOB, const bool &LOCKED_THREADS = false);

      //! @brief Like EnqueueJob, but utilizes the current thread to run the job if it cannot be passed off onto a newly
      //!        created thread
      //! @param JOB the job to run
      void RunJob( util::ShPtr< sched::JobInterface> &JOB);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the next job, called by theads to retrieve the next job
      util::SiPtr< sched::JobInterface> GetJobToRunOrDie();

      //! @brief signal that a job is done
      void SetJobIsDone( const util::SiPtr< sched::JobInterface> &JOB);

      //! @brief This auxiliary function can't be a non-static class member, must be void *MyFun( void *MyArg)
      //! @param EMPTY This is a void pointer, points to NULL
      static void *InitiateThread( void *EMPTY);

      //! @brief convert pthread error code to a string
      //! @param ERROR_CODE the error code
      //! @return a descriptive string for the error code
      static std::string PthreadErrorCodeToString( const int ERROR_CODE);

      friend class Scheduler; //!< Friend to access Initialize and Terminate

      //! @brief initialize the scheduler
      void Initialize();

      //! @brief terminate the scheduler; waits for any currently running jobs to finish
      void Terminate();

    }; // class ThreadQueue

    //! Instantiate the lone instance of this class
    ThreadQueue ThreadQueue::s_ThreadQueue;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Scheduler
    //! Header for a scheduler providing parallelization via pthreads.
    //! Jobs are submitted to the scheduler, placed into threads, and ran according to priorities.
    //! It is up to the user to decide how and when to schedule jobs to maximize the benefits of parallelization.
    //!
    //! @author woetzen, riddeljs
    //!
    //! @see @link example_pthread_scheduler.cpp @endlink
    //!
    //! @date 17.06.2009
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Scheduler :
      public sched::SchedulerInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone
      Scheduler *Clone() const
      {
        return new Scheduler( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief How many unused cpus are there?
      //! @return a size_t representing the number of unused CPUS
      size_t GetNumberUnusedCPUS() const
      {
        return ThreadQueue::s_ThreadQueue.GetNumberThreadsAvailable();
      }

      //! @brief reserve available CPUs for jobs. Useful for multi-level threading, because threads should not be
      //!        launched for the inner-most processes.  The reserved CPUs should not be available to other threads
      //! @param MAX maximum number of CPUs to reserve
      //! @return number of CPUs actually reserved (==0 if no CPUs are available)
      size_t ReserveAvailableCPUs( const size_t &MAX)
      {
        return ThreadQueue::s_ThreadQueue.ReserveAvailableCPUs( MAX);
      }

      //! @brief release a reservation held by this thread
      //! @param NR_TO_RELEASE the number of previously reserved CPUs to release
      void ReleaseReservedCPUs( const size_t &NR_TO_RELEASE)
      {
        ThreadQueue::s_ThreadQueue.ReleaseReservedCPUs( NR_TO_RELEASE);
      }

      //! @brief Submit a job to the scheduler
      //! @param SP_JOB ShPtr to a job
      void SubmitJob( util::ShPtr< sched::JobInterface> &SP_JOB) const
      {
        ThreadQueue::s_ThreadQueue.EnqueueJob( SP_JOB);
      }

      //! @brief Like SubmitJob, but runs the job with the current process if no additional processors are available
      //!        rather than waiting for them to become available.  Has no effect on the serial scheduler
      //! @param SP_JOB ShPtr to a job
      //! @note this functionality is necessary to avoid a race condition between requesting the # of unused processors
      //!       and running a job; since multiple threads may call GetNumberUnusedCPUS at the same time, then
      //!       submit that many jobs, only to have a deadlock when none of them completes
      //!       When such behavior is possible, and the jobs should not just be enqueued (as with SubmitJob),
      //!       call Runjob instead of SubmitJob
      void RunJob( util::ShPtr< sched::JobInterface> &SP_JOB) const
      {
        ThreadQueue::s_ThreadQueue.RunJob( SP_JOB);
      }

      //! @brief Join method to call when closing threads
      //! @param SP_JOB job to join
      void Join( util::ShPtr< sched::JobInterface> &SP_JOB) const
      {
        ThreadQueue::s_ThreadQueue.Join( SP_JOB);
      }

      //! returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
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

    private:

      //! @brief initialize the scheduler
      void Initialize()
      {
        ThreadQueue::s_ThreadQueue.Initialize();
      }

      //! @brief terminate the scheduler; waits for any currently running jobs to finish
      void Terminate()
      {
        ThreadQueue::s_ThreadQueue.Terminate();
      }

      //! instance as a Scheduler enum
      static const sched::Scheduler e_PThreadScheduler;

    }; // class Scheduler

    //! instance as a Scheduler enum
    const sched::Scheduler Scheduler::e_PThreadScheduler
    (
      sched::GetSchedulers().AddEnum( "PThread", util::ShPtr< sched::SchedulerInterface>( new Scheduler()))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ThreadQueue::ThreadQueue() :
      m_NumberThreadsTotal( 0),
      m_NumberThreadsWaiting( 0),
      m_NumberThreadsReserved( 0),
      m_Terminate( false)
    {
      BCL_Assert( pthread_mutex_init( &m_ThreadMutex, NULL) == 0, "unable to initialize thread mutex");
      BCL_Assert( pthread_mutex_init( &m_QueueMutex, NULL) == 0, "unable to initialize queue mutex");
      BCL_Assert( pthread_mutex_init( &m_ConditionMutex, NULL) == 0, "unable to initialize conditions mutex");
      BCL_Assert( pthread_cond_init( &m_NewJobCondition, NULL) == 0, "unable to initialize new job condition");
    }

    //! @brief destructor
    ThreadQueue::~ThreadQueue()
    {
      Terminate();
      BCL_Assert( pthread_mutex_destroy( &m_ThreadMutex) == 0, "unable to destroy thread mutex");
      pthread_mutex_destroy( &m_QueueMutex);
      pthread_mutex_destroy( &m_ConditionMutex);
      pthread_cond_destroy( &m_NewJobCondition);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief Get # of threads currently running (this is no guarantee that this many are available, however)
    size_t ThreadQueue::GetNumberThreadsAvailable()
    {
      // lock the thread mutex
      // This is not strictly necessary : we could return the variable directly without segfaults,
      // however, locking the mutex ensures that that the current number is accurate at the time that this function
      // is called and that m_NumberThreadsRunning is not in the process of being updated
      pthread_mutex_lock( &m_QueueMutex);
      const size_t num_waiting( std::min( m_NumberThreadsWaiting, m_NumberThreadsTotal - m_NumberThreadsReserved));
      pthread_mutex_unlock( &m_QueueMutex);
      return num_waiting;
    }

    //! @brief reserve available CPUs for jobs. Useful for multi-level threading, because threads should not be
    //!        launched for the inner-most processes.  The reserved CPUs decrease the number returned by GetNumberThreadsAvailable
    //! @param MAX maximum number of CPUs to reserve
    //! @return number of CPUs actually reserved (==0 if no CPUs are available)
    size_t ThreadQueue::ReserveAvailableCPUs( const size_t &MAX)
    {
      pthread_mutex_lock( &m_QueueMutex);
      const size_t nr_reserved( std::min( m_NumberThreadsWaiting - m_NumberThreadsReserved, MAX));
      m_NumberThreadsReserved += nr_reserved;
      pthread_mutex_unlock( &m_QueueMutex);
      return nr_reserved;
    }

    //! @brief release a reservation held by this thread
    //! @param NR_TO_RELEASE the number of previously reserved CPUs to release
    void ThreadQueue::ReleaseReservedCPUs( const size_t &NR_TO_RELEASE)
    {
      if( NR_TO_RELEASE)
      {
        pthread_mutex_lock( &m_QueueMutex);
        BCL_Assert( m_NumberThreadsReserved >= NR_TO_RELEASE, "Tried to release more CPUs than were reserved");
        m_NumberThreadsReserved -= NR_TO_RELEASE;
        pthread_mutex_unlock( &m_QueueMutex);
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief enqueue a job
    //! @param JOB the job to enqueue
    //! @param LOCKED_THREADS true if this process has already locked m_ThreadMutex
    void ThreadQueue::EnqueueJob( util::ShPtr< sched::JobInterface> &JOB, const bool &LOCKED_THREADS)
    {
      // create a new mutex and condition variable for this job
      pthread_cond_t job_condition;
      BCL_Assert( pthread_cond_init( &job_condition, NULL) == 0, "unable to initialize job condition");

      // insert the job into the corresponding condition variables
      pthread_mutex_lock( &m_ConditionMutex);
      m_JobsToConditionVariables[ JOB] = std::make_pair( job_condition, false);
      pthread_mutex_unlock( &m_ConditionMutex);

      // place the job in the queue
      if( !LOCKED_THREADS)
      {
        pthread_mutex_lock( &m_QueueMutex);
      }
      m_AvailableJobs.push_back( JOB);

      // signal a thread to wakeup and take the job
      pthread_cond_signal( &m_NewJobCondition);
      pthread_mutex_unlock( &m_QueueMutex);
    }

    //! @brief Like EnqueueJob, but utilizes the current thread to run the job if it cannot be passed off onto a newly
    //!        created thread
    //! @param JOB the job to run
    void ThreadQueue::RunJob( util::ShPtr< sched::JobInterface> &JOB)
    {
      // test whether there are any processes remaining to launch this job
      pthread_mutex_lock( &m_QueueMutex);
      if( m_NumberThreadsWaiting)
      {
        // another process can be created to run this job, do so
        // Keep the thread mutex locked so that another process does not steal this process' slot
        EnqueueJob( JOB, true);
      }
      else
      {
        pthread_mutex_unlock( &m_QueueMutex);
        // no more processes available, run the job directly
        JOB->Run();
      }
    }

    //! @brief join the job; really just wait for it to complete
    //! @param JOB the job to wait for
    void ThreadQueue::Join( const util::ShPtr< sched::JobInterface> &JOB)
    {
      if( !JOB.IsDefined())
      {
        BCL_MessageStd( "Null job given to join!");
        return;
      }

      pthread_mutex_lock( &m_ConditionMutex);

      // look for the job in the condition map
      std::map
      <
        util::SiPtr< const sched::JobInterface>,
        std::pair< pthread_cond_t, bool>
      >::iterator itr( m_JobsToConditionVariables.find( JOB));

      // case: job was never submitted or was already joined
      if( itr == m_JobsToConditionVariables.end())
      {
        pthread_mutex_unlock( &m_ConditionMutex);
        return;
      }

      // necessary to have a while loop here in case of spurious wakeups from pthread_cond_wait
      while( !itr->second.second)
      {
        // wait for the condition to signal
        pthread_cond_wait( &itr->second.first, &m_ConditionMutex);
      }

      // destroy the condition variable and mutex
      BCL_Assert( pthread_cond_destroy( &itr->second.first) == 0, "unable to destroy job condition");

      m_JobsToConditionVariables.erase( itr);
      pthread_mutex_unlock( &m_ConditionMutex);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the next job, called by theads to retrieve the next job
    util::SiPtr< sched::JobInterface> ThreadQueue::GetJobToRunOrDie()
    {
      util::SiPtr< sched::JobInterface> job;
      pthread_mutex_lock( &m_QueueMutex);

      ++m_NumberThreadsWaiting;

      // check whether the scheduler has been terminated, in which case, time to exit
      while( !m_Terminate && m_AvailableJobs.empty())
      {
        // no jobs, no reason to terminate, wait for another job
        pthread_cond_wait( &m_NewJobCondition, &m_QueueMutex);
      }

      // check whether all threads are supposed to terminate
      if( m_Terminate)
      {
        --m_NumberThreadsTotal;
        --m_NumberThreadsWaiting;
        pthread_mutex_unlock( &m_QueueMutex);
        pthread_detach( pthread_self());
        pthread_exit( NULL);
      }

      // take the last job off the stack (last-in, first-out; this prevents stalling)
      job = m_AvailableJobs.back();
      m_AvailableJobs.pop_back();
      --m_NumberThreadsWaiting;
      pthread_mutex_unlock( &m_QueueMutex);
      return job;
    }

    //! @brief signal that a job is done
    void ThreadQueue::SetJobIsDone( const util::SiPtr< sched::JobInterface> &JOB)
    {
      pthread_mutex_lock( &m_ConditionMutex);

      // look for the job in the condition map
      std::map
      <
        util::SiPtr< const sched::JobInterface>,
        std::pair< pthread_cond_t, bool>
      >::iterator itr( m_JobsToConditionVariables.find( JOB));
      BCL_Assert( itr != m_JobsToConditionVariables.end(), "Could not find job in condition map!");

      // get the pair out of the iterator
      std::pair< pthread_cond_t, bool> &condition_mutex( itr->second);

      condition_mutex.second = true;

      // signal that the job is complete
      pthread_cond_signal( &condition_mutex.first);

      // unlock the condition map mutex
      pthread_mutex_unlock( &m_ConditionMutex);
    }

    //! @brief convert pthread error code to a string
    //! @param ERROR_CODE the error code
    //! @return a descriptive string for the error code
    std::string ThreadQueue::PthreadErrorCodeToString( const int ERROR_CODE)
    {
      switch( ERROR_CODE)
      {
        case EAGAIN: return "EAGAIN";
        case EINVAL: return "EINVAL";
        case EPERM: return "EPERM";
        case ENOMEM: return "ENOMEM";
        default: return "unknown error " + util::Format()( ERROR_CODE);
      }
    }

    //! @brief This auxiliary function can't be a non-static class member, must be void *MyFun( void *MyArg)
    //! @param EMPTY This is a void pointer, points to NULL
    void *ThreadQueue::InitiateThread( void *EMPTY)
    {
      // Try to retrieve a job
      util::SiPtr< sched::JobInterface> job( s_ThreadQueue.GetJobToRunOrDie());

      // so long as there are jobs to run
      while( job.IsDefined())
      {
        // run the job
        job->Run();

        // remove the job from the active jobs
        s_ThreadQueue.SetJobIsDone( job);

        // get the next job
        job = s_ThreadQueue.GetJobToRunOrDie();
      }

      // added just to avoid VS to complain about no return value
      return ( void *) 0;
    }

    //! @brief initialize the scheduler
    void ThreadQueue::Initialize()
    {
      pthread_mutex_lock( &m_ThreadMutex);
      m_NumberThreadsTotal = sched::GetNumberCPUs();
      m_NumberThreadsWaiting = 0;
      m_Terminate = false;
      m_Threads.resize( m_NumberThreadsTotal);
      for( size_t thread_id( 0); thread_id < m_NumberThreadsTotal; ++thread_id)
      {
        // attempt to create the thread, pthread_create requires a regular pointer
        const int create_error
        (
          pthread_create
          (
            &m_Threads[ thread_id],
            NULL, // no need to create it joinable since the thread will rejoin itself
            &InitiateThread,
            NULL
          )
        );

        BCL_Assert
        (
          create_error == 0,
          "Unable to create thread! Currently active threads: " +
          util::Format()( m_NumberThreadsTotal) + " error: " + PthreadErrorCodeToString( create_error)
        );
      }
      pthread_mutex_unlock( &m_ThreadMutex);
    }

    //! @brief terminate the scheduler; waits for any currently running jobs to finish
    void ThreadQueue::Terminate()
    {
      pthread_mutex_lock( &m_ThreadMutex);
      m_Terminate = true;

      // wakeup all the waiting threads so that they can notice that they are to terminate
      while( m_NumberThreadsTotal)
      {
        pthread_mutex_lock( &m_QueueMutex);
        pthread_cond_broadcast( &m_NewJobCondition);
        pthread_mutex_unlock( &m_QueueMutex);
      }
      m_Threads.clear();

      pthread_mutex_unlock( &m_ThreadMutex);
    }

  } // namespace pthread
} // namespace bcl
