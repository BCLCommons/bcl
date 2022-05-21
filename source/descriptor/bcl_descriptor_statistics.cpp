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
#include "descriptor/bcl_descriptor_statistics.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! Anonymous namespace, keeps these variables local
    namespace
    {
      // create all instances of interest
      static const bool s_InstancesCreated =
        math::RunningStatistics< float>::AddInstances< Base< biol::AABase, float>, biol::AABase, Statistics>() &&
        math::RunningStatistics< float>::AddInstances
        <
          Base< chemistry::AtomConformationalInterface, float>,
          chemistry::AtomConformationalInterface,
          Statistics
        >() &&
        math::RunningStatistics< float>::AddInstances< Base< char, float>, char, Statistics>()
        && math::RunningStatistics< float>::AddInstances< Base< biol::Mutation, float>, biol::Mutation, Statistics>();
    } // namespace

  } // namespace descriptor
} // namespace bcl
