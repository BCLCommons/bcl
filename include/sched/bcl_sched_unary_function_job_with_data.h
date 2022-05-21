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

#ifndef BCL_SCHED_UNARY_FUNCTION_JOB_WITH_DATA_H_
#define BCL_SCHED_UNARY_FUNCTION_JOB_WITH_DATA_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_interface.h"
#include "bcl_sched_job_with_result_base.h"
#include "util/bcl_util_function_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UnaryFunctionJobWithData
    //! @brief Create instances of this class to prepare tertiary function jobs for submission to Scheduler.
    //! Because of the template arguments, operations are implemented in this header.
    //!
    //! @see @link example_sched_unary_function_job_with_data.cpp @endlink
    //! @author riddeljs
    //! @date 12.04.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ResultType, typename t_FunctionClass>
    class UnaryFunctionJobWithData :
      public JobWithResultBase< t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      size_t m_GroupID; //! Priority of the job for concurrency control

      //! Pointer to instance of class
      t_FunctionClass *m_Class;

      //! @brief Sharepointer to the tertiary function corresponding to this job
      util::FunctionWrapper< t_ArgumentType1, t_ResultType, t_FunctionClass> m_UnaryFunction;

      //! Pointers to arguments to the tertiary function
      t_ArgumentType1 *m_Argument;

      //! Status of job
      JobInterface::JobStatus m_JobStatus;

    public:

      typedef t_ResultType ( t_FunctionClass::*PtrToMemberFunction)( t_ArgumentType1 &);
      typedef t_ResultType ( t_FunctionClass::*PtrToConstMemberFunction)( t_ArgumentType1 &) const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      UnaryFunctionJobWithData()
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Static Function, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a static function
      //! @param ARGUMENT1 First argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      UnaryFunctionJobWithData
      (
        const size_t &GROUPID,
        t_FunctionClass &CLASS,
        t_ResultType ( *FUNCTION_TO_SUBMIT)( t_ArgumentType1 &),
        t_ArgumentType1 &ARGUMENT1,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_UnaryFunction( FUNCTION_TO_SUBMIT),
        m_Argument( &ARGUMENT1),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Query, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param QUERY_TO_SUBMIT In this case, a tertiary query
      //! @param ARGUMENT1 First argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      UnaryFunctionJobWithData
      (
        const size_t &GROUPID,
        const t_FunctionClass &CLASS,
        PtrToConstMemberFunction QUERY_TO_SUBMIT,
        t_ArgumentType1 &ARGUMENT1,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( const_cast< t_FunctionClass *>( &CLASS)),
        m_UnaryFunction( QUERY_TO_SUBMIT),
        m_Argument( &ARGUMENT1),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Function, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a tertiary function
      //! @param ARGUMENT1 First argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      UnaryFunctionJobWithData
      (
        const size_t &GROUPID,
        t_FunctionClass &CLASS,
        PtrToMemberFunction FUNCTION_TO_SUBMIT,
        t_ArgumentType1 &ARGUMENT1,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_UnaryFunction( FUNCTION_TO_SUBMIT),
        m_Argument( &ARGUMENT1),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! Copy Constructor
      UnaryFunctionJobWithData *Clone() const
      {
        return new UnaryFunctionJobWithData( *this);
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

      //! @brief Run the function associated with this job and return the result, if any
      //! @return The result of the job, if any
      t_ResultType RunFunction()
      {
        return m_UnaryFunction( *m_Class, *m_Argument);
      }

    }; // class UnaryFunctionJobWithData

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_UNARY_FUNCTION_JOB_WITH_DATA_H_
