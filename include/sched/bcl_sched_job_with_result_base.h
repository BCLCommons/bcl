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

#ifndef BCL_SCHED_JOB_WITH_RESULT_BASE_H_
#define BCL_SCHED_JOB_WITH_RESULT_BASE_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_job_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class JobWithResultBase
    //! @brief abstracts the handling of jobs that do/do not return values
    //!
    //! @see @link example_sched_job_with_result_base.cpp @endlink
    //! @author mendenjl
    //! @date 06.14.2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Result>
    class JobWithResultBase :
      public JobInterface
    {

    private:

    //////////
    // data //
    //////////

      t_Result *m_Data; //!< pointer to the result of the operation

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief setup the result
      //! @param DATA a pointer to the result object, which will be assigned the result of the job
      JobWithResultBase( t_Result *DATA) :
        m_Data( DATA)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Run the underlying function for the job
      void Run()
      {
        SetJobStatus( JobInterface::e_RUNNING);
        *m_Data = RunFunction(); // run the function and do the assignment
        SetJobStatus( JobInterface::e_FINISHED);
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief run the internal function without setting up the status
      //! @return the result of running the function
      virtual t_Result RunFunction() = 0;

    }; // class JobWithResultBase

    //! specialization of JobWithResultBase for void return types (jobs that do not return results)
    template<>
    class BCL_API JobWithResultBase< void> :
      public JobInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief setup the result
      //! This makes it unnecessary to specialize for void vs. non-void results
      JobWithResultBase( void *)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Run the underlying function for the job
      void Run()
      {
        SetJobStatus( JobInterface::e_RUNNING);
        RunFunction(); // run the function, which has no return value
        SetJobStatus( JobInterface::e_FINISHED);
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief run the internal function without setting up the status
      virtual void RunFunction() = 0;

    }; // class JobWithResultBase

  } // namespace sched
} // namespace bcl

#endif //BCL_SCHED_JOB_WITH_RESULT_BASE_H_

