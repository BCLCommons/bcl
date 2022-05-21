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
#include "sched/bcl_sched_scheduler_interface.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "sched/bcl_sched_schedulers.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    //! Parameter for setting the number of cpus (or possible parallel jobs) over the command line
    util::ShPtr< command::ParameterInterface> &SchedulerInterface::GetParameterNumberCPUs()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter_number_cpu
      (
        new command::Parameter
        (
          "number_cpus",
          "number of cpus for a multi job scheduler",
          command::ParameterCheckRanged< size_t>( 1, 1000),
          "1"
        )
      );

      // end
      return s_parameter_number_cpu;
    };

    //! @brief get currently used Scheduler
    SchedulerInterface &GetScheduler()
    {
      return GetSchedulers().GetCurrentScheduler();
    }

    //! @brief return the number of cpus as given on the command line
    size_t GetNumberCPUs()
    {
      return SchedulerInterface::GetParameterNumberCPUs()->GetNumericalValue< size_t>();
    }

  } // namespace sched
} // namespace bcl
