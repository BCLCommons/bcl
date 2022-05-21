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

#ifndef BCL_SCHED_TERTIARY_FUNCTION_JOB_WITH_DATA_H_
#define BCL_SCHED_TERTIARY_FUNCTION_JOB_WITH_DATA_H_

// include the namespace header
#include "bcl_sched.h"

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_with_result_base.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_tertiary_function_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TertiaryFunctionJobWithData
    //! @brief Create instances of this class to prepare tertiary function jobs for submission to Scheduler.
    //! @details Because of the template arguments, operations are implemented in this header.
    //!
    //! @see @link example_sched_jobs_with_data.cpp @endlink
    //! @author riddeljs
    //! @date 12.04.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ResultType, typename t_FunctionClass>
    class TertiaryFunctionJobWithData :
      public JobWithResultBase< t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! Priority of the job for concurrency control
      size_t m_GroupID;

      //! Pointer to instance of class
      t_FunctionClass *m_Class;

      //! @brief Sharepointer to the tertiary function corresponding to this job
      util::ShPtr< util::FourInputFunctionInterface
        <
          t_FunctionClass,
          t_ArgumentType1,
          t_ArgumentType2,
          t_ArgumentType3,
          t_ResultType
        > > m_TertiaryFunction;

      t_ArgumentType1 *m_Argument1; //!< Pointer to 1st argument of the tertiary function
      t_ArgumentType2 *m_Argument2; //!< Pointer to 2nd argument of the tertiary function
      t_ArgumentType3 *m_Argument3; //!< Pointer to 3rd argument of the tertiary function

      //! Status of job
      JobInterface::JobStatus m_JobStatus;

    public:

      typedef t_ResultType ( t_FunctionClass::*PtrToMemberFunction)( t_ArgumentType1 &, t_ArgumentType2 &, t_ArgumentType3 &);
      typedef t_ResultType ( t_FunctionClass::*PtrToConstMemberFunction)( t_ArgumentType1 &, t_ArgumentType2 &, t_ArgumentType3 &) const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>()
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Static Function, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a static function
      //! @param ARGUMENT1 First argument
      //! @param ARGUMENT2 Second argument
      //! @param ARGUMENT3 Third argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
      (
        const size_t &GROUPID,
        t_FunctionClass &CLASS,
        t_ResultType ( *FUNCTION_TO_SUBMIT)( t_ArgumentType1 &, t_ArgumentType2 &, t_ArgumentType3 &),
        t_ArgumentType1 &ARGUMENT1,
        t_ArgumentType2 &ARGUMENT2,
        t_ArgumentType3 &ARGUMENT3,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_TertiaryFunction
        (
          new util::TertiaryFunctionWrapper< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
          ( FUNCTION_TO_SUBMIT)
        ),
        m_Argument1( &ARGUMENT1),
        m_Argument2( &ARGUMENT2),
        m_Argument3( &ARGUMENT3),
        m_JobStatus( JOBSTATUS)
      {}

      //! @brief Construct the job, given the JobID, GroupID, Query, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param QUERY_TO_SUBMIT In this case, a tertiary query
      //! @param ARGUMENT1 First argument
      //! @param ARGUMENT2 Second argument
      //! @param ARGUMENT3 Third argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
      (
        const size_t &GROUPID,
        const t_FunctionClass &CLASS,
        PtrToConstMemberFunction QUERY_TO_SUBMIT,
        t_ArgumentType1 &ARGUMENT1,
        t_ArgumentType2 &ARGUMENT2,
        t_ArgumentType3 &ARGUMENT3,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( const_cast< t_FunctionClass *>( &CLASS)),
        m_TertiaryFunction
        (
          new util::TertiaryFunctionWrapper< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
          ( QUERY_TO_SUBMIT)
        ),
        m_Argument1( &ARGUMENT1),
        m_Argument2( &ARGUMENT2),
        m_Argument3( &ARGUMENT3),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Function, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param CLASS
      //! @param FUNCTION_TO_SUBMIT In this case, a tertiary function
      //! @param ARGUMENT1 First argument
      //! @param ARGUMENT2 Second argument
      //! @param ARGUMENT3 Third argument
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT
      TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
      (
        const size_t &GROUPID,
        t_FunctionClass &CLASS,
        PtrToMemberFunction FUNCTION_TO_SUBMIT,
        t_ArgumentType1 &ARGUMENT1,
        t_ArgumentType2 &ARGUMENT2,
        t_ArgumentType3 &ARGUMENT3,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_Class( &CLASS),
        m_TertiaryFunction
        (
          new util::TertiaryFunctionWrapper< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>
          ( FUNCTION_TO_SUBMIT)
        ),
        m_Argument1( &ARGUMENT1),
        m_Argument2( &ARGUMENT2),
        m_Argument3( &ARGUMENT3),
        m_JobStatus( JOBSTATUS)
      {
      }

      //! returns class name
      //! the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! Copy Constructor
      TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass> *Clone() const
      {
        return new TertiaryFunctionJobWithData< t_ArgumentType1, t_ArgumentType2, t_ArgumentType3, t_ResultType, t_FunctionClass>( *this);
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
        return m_TertiaryFunction->operator()( *m_Class, *m_Argument1, *m_Argument2, *m_Argument3);
      }

    }; // class TertiaryFunctionJobWithData

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_TERTIARY_FUNCTION_JOB_WITH_DATA_H_
