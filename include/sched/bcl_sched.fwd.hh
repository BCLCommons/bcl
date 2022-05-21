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

#ifndef BCL_SCHED_FWD_HH_
#define BCL_SCHED_FWD_HH_

// include the dependency file for this header
#include "bcl_sched.depends.fwd.hh"

// This file contains forward declarations for the sched namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace sched
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class JobInterface;
    class Mutex;
    class SchedulerInterface;
    class Schedulers;
    class ScopeLock;
    class SerialScheduler;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_FunctionClass>
    class BinaryFunctionJobWithData;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionJobWithData;

    template< typename t_Result>
    class JobWithResultBase;

    template< typename t_ArgumentType, typename t_ResultType>
    class SumFunction;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ResultType, typename t_FunctionClass>
    class TertiaryFunctionJobWithData;

    template< typename t_FunctionClass, typename t_ResultType>
    class ThunkJob;

    template< typename t_ArgumentType1, typename t_ResultType, typename t_FunctionClass>
    class UnaryFunctionJobWithData;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< SchedulerInterface>, Schedulers> Scheduler;

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_FWD_HH_
