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

#ifndef BCL_SCHED_SCHEDULERS_H_
#define BCL_SCHED_SCHEDULERS_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sched_scheduler_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Schedulers
    //! @brief manages the system-wide list of enumerated schedulers
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date March 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Schedulers :
      public util::Enumerate< util::ShPtr< SchedulerInterface>, Schedulers>
    {
      friend class util::Enumerate< util::ShPtr< SchedulerInterface>, Schedulers>;

    public:

    //////////
    // data //
    //////////

      Scheduler e_Default; //!< default scheduler takes the job and executes it directly - serial scheduler

    ///////////
    // flags //
    ///////////

      //! Flag to switch between numerated Schedulers and apss number of possible parallel jobs
      static util::ShPtr< command::FlagInterface> &GetFlagSchedulerCPUS();

      //! @brief Initialize the logger from the command line flag
      static void UpdateCurrentSchedulerFromCommandLineFlag();

    private:

      util::SiPtr< SchedulerInterface> m_CurrentScheduler;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      Schedulers();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Initialize and return the current scheduler
      //! @return one and only reference to one of the schedulers
      SchedulerInterface &GetCurrentScheduler();

    }; // class Schedulers

    //! @brief get enumerated list of Schedulers
    BCL_API Schedulers &GetSchedulers();

  } // namespace sched

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< sched::SchedulerInterface>, sched::Schedulers>;

  } // namespace util
} // namespace bcl

#endif // BCL_SCHED_SCHEDULERS_H_
