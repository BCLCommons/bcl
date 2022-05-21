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
#include "sched/bcl_sched_schedulers.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "sched/bcl_sched_serial_scheduler.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

  //////////
  // data //
  //////////

    //! @brief enum, that adds Platform flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_SchedulerFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag
      (
        Schedulers::GetFlagSchedulerCPUS(),
        command::e_Pthread
      )
    );

    //! Flag to switch between numerated Schedulers and apss number of possible parallel jobs
    util::ShPtr< command::FlagInterface> &Schedulers::GetFlagSchedulerCPUS()
    {
      static util::ShPtr< command::FlagInterface> s_flag_scheduler_spus
      (
        new command::FlagStatic
        (
          "scheduler",
          "choice of scheduler and number of cpus",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "scheduler",
                "type of scheduler",
                command::ParameterCheckEnumerate< Schedulers>(),
                GetSchedulers().e_Default.GetName()
              )
            ),
            SchedulerInterface::GetParameterNumberCPUs()
          ),
          &Schedulers::UpdateCurrentSchedulerFromCommandLineFlag
        )
      );

      // end
      return s_flag_scheduler_spus;
    }

    //! @brief Initialize the logger from the command line flag
    void Schedulers::UpdateCurrentSchedulerFromCommandLineFlag()
    {
      // stop any current jobs from running
      if( GetSchedulers().m_CurrentScheduler.IsDefined())
      {
        GetSchedulers().m_CurrentScheduler->Terminate();
      }
      GetSchedulers().m_CurrentScheduler = **Scheduler( GetFlagSchedulerCPUS()->GetFirstParameter()->GetValue());

      // initialize the new scheduler
      if( GetSchedulers().m_CurrentScheduler.IsDefined())
      {
        GetSchedulers().m_CurrentScheduler->Initialize();
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Default constructor
    Schedulers::Schedulers() :
      util::Enumerate< util::ShPtr< SchedulerInterface>, Schedulers>( false),
      e_Default( AddEnum( "Serial", util::ShPtr< SchedulerInterface>( new SerialScheduler()))),
      m_CurrentScheduler( **e_Default)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Schedulers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Initialize and return the current scheduler
    //! @return one and only reference to one of the schedulers
    SchedulerInterface &Schedulers::GetCurrentScheduler()
    {
      // return the scheduler
      return *m_CurrentScheduler;
    }

    //! @brief get enumerated list of Schedulers
    Schedulers &GetSchedulers()
    {
      return Schedulers::GetEnums();
    }

  } // namespace sched

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< sched::SchedulerInterface>, sched::Schedulers>;

  } // namespace util
} // namespace bcl
