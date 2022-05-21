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

#ifndef BCL_SCHED_FUNCTION_JOB_WITH_DATA_H_
#define BCL_SCHED_FUNCTION_JOB_WITH_DATA_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_interface.h"
#include "bcl_sched_job_with_result_base.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionJobWithData
    //! @brief Create instances of this class to prepare a job for submission to Scheduler.
    //!
    //! @see @link example_sched_function_job_with_data.cpp @endlink
    //! @author woetzen
    //! @date 21.08.2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionJobWithData :
      public JobWithResultBase< t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      size_t m_GroupID; //! Priority of the job for concurrency control

      //! @brief Sharepointer to the tertiary function corresponding to this job
      util::SiPtr< const util::FunctionInterface< t_ArgumentType, t_ResultType> > m_UnaryFunction;

      //! Pointers to arguments to the tertiary function
      const t_ArgumentType *m_Argument;

      //! Status of job
      JobInterface::JobStatus m_JobStatus;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      FunctionJobWithData()
      {
      }

      //! @brief Construct the job, given the JobID, GroupID, Static Function, and arguments
      //! @param GROUPID Priority of the job for concurrency control, jobs in same group are concurrent
      //! @param SP_FUNCTION_TO_SUBMIT spointer to function
      //! @param ARGUMENT argument to pass to function
      //! @param JOBSTATUS Initial status of the job
      //! @param RESULT pointer where result is written to
      FunctionJobWithData
      (
        const size_t &GROUPID,
        const util::SiPtr< const util::FunctionInterface< t_ArgumentType, t_ResultType> > &SP_FUNCTION_TO_SUBMIT,
        const t_ArgumentType &ARGUMENT,
        const JobInterface::JobStatus &JOBSTATUS,
        t_ResultType *RESULT
      ) :
        JobWithResultBase< t_ResultType>( RESULT),
        m_GroupID( GROUPID),
        m_UnaryFunction( SP_FUNCTION_TO_SUBMIT),
        m_Argument( &ARGUMENT),
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
      FunctionJobWithData *Clone() const
      {
        return new FunctionJobWithData( *this);
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
        return m_UnaryFunction->operator()( *m_Argument);
      }

    }; // class FunctionJobWithData

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_FUNCTION_JOB_WITH_DATA_H_
