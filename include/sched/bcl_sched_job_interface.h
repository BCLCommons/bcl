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

#ifndef BCL_SCHED_JOB_INTERFACE_H_
#define BCL_SCHED_JOB_INTERFACE_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class JobInterface
    //! @brief is the bcl job interface.
    //! @details The job interface is created by a class instance to submit to the scheduler for parallelization.
    //! Acknowledgments go to woetzen for help with the design
    //!
    //! @remarks example unnecessary
    //! @author riddeljs
    //! @date 11.18.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API JobInterface :
      public util::ObjectInterface
    {

    public:

    ///////////
    // types //
    ///////////

      //! @enum JobStatus
      //! @brief Custom type to determine status of job
      enum JobStatus { e_READY, e_WAITING, e_RUNNING, e_FINISHED, e_ERROR};

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      JobInterface()
      {
      }

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      virtual JobInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief Run the underlying function for the job
      virtual void Run() = 0;

      //! @brief Get the unique identifier for the group (priority) of the job
      //! @return an integer representing the group id
      virtual size_t GetGroupID() const = 0;

      //! @brief Get the status of the job
      //! @return The enumerator type JobStatus
      virtual JobStatus GetJobStatus() const = 0;

      //! @brief Set the status of the job
      //! @param JOBSTATUS new status for the job
      virtual void SetJobStatus( const JobStatus &JOBSTATUS) = 0;

    }; // class JobInterface

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_JOB_INTERFACE_H_
