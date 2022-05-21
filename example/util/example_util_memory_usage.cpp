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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "util/bcl_util_memory_usage.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_memory_usage.cpp
  //!
  //! @author mendenjl
  //! @date Sep 26, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilMemoryUsage :
    public ExampleInterface
  {
  public:

    ExampleUtilMemoryUsage *Clone() const
    {
      return new ExampleUtilMemoryUsage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      util::MemoryUsage::WritePeakMemoryUsageInfo( util::GetLogger());
      util::MemoryUsage::WriteCurrentMemoryUsageInfo( util::GetLogger());

      // allocate and use a large array, check that memory increased
      const size_t mb_to_allocate( 64);
      const size_t bytes_per_mb( 1 << 20);

      // determine how many floats should go into the linal::Vector< float> to equal mb_to_allocate
      const size_t number_floats
      (
        ( mb_to_allocate * bytes_per_mb - sizeof( linal::Vector< float>)) / sizeof( float)
      );

      // get the memory usage before allocating the vector
      util::MemoryUsage mem_usage_before;

      // allocate the vector
      linal::Vector< float> test_vector( number_floats, 0.0);

      // get the current memory usage
      util::MemoryUsage mem_usage_after;

      // ensure that the class is internally consistent
      BCL_ExampleCheck
      (
        util::IsDefined( mem_usage_before.GetRAMInUse()),
        util::IsDefined( mem_usage_after.GetRAMInUse())
      );
      if( util::IsDefined( mem_usage_before.GetRAMInUse()))
      {
        // determine how many megabytes the ram usage increased by
        const long mb_increase( mem_usage_after.GetRAMInUse() - mem_usage_before.GetRAMInUse());

        // make sure that ram in use did not decrease
        BCL_ExampleIndirectCheck
        (
          mb_increase >= long( 0),
          true,
          "allocating a large array should not decrease ram usage"
        );

        // write out how much ram actually changed from allocating the vector
        BCL_MessageStd
        (
          "Allocating a " + util::Format()( mb_to_allocate) + "MB array increased RAM used by "
          + util::Format()( mb_increase)
          + "MB; the two numbers should be close to equal if the system is not memory bound"
        );
      }

      BCL_ExampleCheck
      (
        util::IsDefined( mem_usage_before.GetVirtualMemoryInUse()),
        util::IsDefined( mem_usage_after.GetVirtualMemoryInUse())
      );
      if( util::IsDefined( mem_usage_before.GetVirtualMemoryInUse()))
      {
        // make sure that virtual memory use changed by at least the amount allocated
        // the os is free to allocate extra memory if it so chooses, so the actual amount allocated may be slightly
        // larger than what was requested
        // This should be true even on a memory bound system; the allocated memory will just be virtual
        const long mem_allocated( mem_usage_after.GetVirtualMemoryInUse() - mem_usage_before.GetVirtualMemoryInUse());
        BCL_ExampleIndirectCheck
        (
          mem_allocated >= long( mb_to_allocate - 1),
          true,
          "allocating " + util::Format()( mb_to_allocate)
          + "MB changes virtual memory usage by at least the same amount, actual change was "
          + util::Format()( mem_allocated) + "MB"
        );
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilMemoryUsage

  const ExampleClass::EnumType ExampleUtilMemoryUsage::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilMemoryUsage())
  );

} // namespace bcl
