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

#ifndef BCL_SCHED_THUNK_JOB_H_
#define BCL_SCHED_THUNK_JOB_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_with_result_base.h"
#include "util/bcl_util_function_interface_nonconst.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_thunk_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ThunkJob
    //! @brief Create instances of this class to prepare thunk jobs for submission to Scheduler.
    //! @details These jobs don't take arguments or output anything, and likely just work on the class object.
    //! Because of the template arguments, operations are implemented in this header.
    //!
    //! @see @link example_sched_jobs_with_data.cpp @endlink
    //! @author riddeljs
    //! @date 12.31.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_FunctionClass, typename t_ResultType>
    class ThunkJob :
      public JobWithResultBase< t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      size_t m_GroupID; //! Priority of the job for concurrency control

      //! Pointer to instance of class
      t_FunctionClass *m_Class;

      //! @brief ShPtr to the function corresponding to this job
      util::ShPtr< util::FunctionInterfaceNonConst< t_FunctionClass, void> > m_Thunk;

      //! Status of the job
      JobInterface::JobStatus m_JobStatus;

      typedef t_ResultType ( t_FunctionClass::*PtrToMemberFunction)();
      typedef t_ResultType ( t_FunctionClass::*PtrToConstMemberFunction)() const;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ThunkJob()
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, and static function.
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a static thunk
      //! @param JOBSTATUS
      //! @param RESULT
      ThunkJob
      (
        const size_t GROUPID,
        t_FunctionClass &CLASS,
        t_ResultType ( *FUNCTION_TO_SUBMIT)(),
        const JobInterface::JobStatus &JOBSTATUS = JobInterface::e_READY,
        t_ResultType *RESULT = NULL
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_Thunk( new util::ThunkWrapper< t_FunctionClass, t_ResultType>( FUNCTION_TO_SUBMIT)),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, and query.
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param QUERY_TO_SUBMIT In this case, a thunk query
      //! @param JOBSTATUS Unique identifier for the job
      //! @param RESULT
      ThunkJob
      (
        const size_t GROUPID,
        const t_FunctionClass &CLASS,
        PtrToConstMemberFunction QUERY_TO_SUBMIT,
        const JobInterface::JobStatus &JOBSTATUS = JobInterface::e_READY,
        t_ResultType *RESULT = NULL
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( const_cast< t_FunctionClass *>( &CLASS)),
        m_Thunk( new util::ThunkWrapper< t_FunctionClass, t_ResultType>( QUERY_TO_SUBMIT)),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, and function.
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a thunk query
      //! @param JOBSTATUS Unique identifier for the job
      //! @param RESULT
      ThunkJob
      (
        const size_t GROUPID,
        t_FunctionClass &CLASS,
        PtrToMemberFunction FUNCTION_TO_SUBMIT,
        const JobInterface::JobStatus &JOBSTATUS = JobInterface::e_READY,
        t_ResultType *RESULT = NULL
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_Thunk( new util::ThunkWrapper< t_FunctionClass, t_ResultType>( FUNCTION_TO_SUBMIT)),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! Copy Constructor
      ThunkJob< t_FunctionClass, t_ResultType> *Clone() const
      {
        return new ThunkJob< t_FunctionClass, t_ResultType>( *this);
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

      //! @brief Get the GroupID for the job, the priority level
      //! @return m_GroupID
      size_t GetGroupID() const
      {
        return m_GroupID;
      }

      //! @brief Get the status of the job
      //! @return The enumerator type JobStatus
      JobInterface::JobStatus GetJobStatus() const
      {
        return m_JobStatus;
      }

      //! @brief Set the status of the job
      //! @param JOBSTATUS new status for the job
      void SetJobStatus( const JobInterface::JobStatus &JOBSTATUS)
      {
        m_JobStatus = JOBSTATUS;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

      //! @brief Run the function associated with this job
      t_ResultType RunFunction()
      {
        return m_Thunk->operator()( *m_Class);
      }

    }; // class ThunkJob

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_THUNK_JOB_H_
