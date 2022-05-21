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
#include "sched/bcl_sched_mutex.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sched_mutex.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSchedMutex :
    public ExampleInterface
  {
  public:

    ExampleSchedMutex *Clone() const
    {
      return new ExampleSchedMutex( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a mutex
      sched::Mutex default_mutex;

      // ensure that new muteces are unlocked
      BCL_ExampleAssert( default_mutex.TestLock(), false);

      // create a shared pointer to a mutex
      util::ShPtr< sched::Mutex> cloned_mutex_sp_a( default_mutex.Clone());
      util::ShPtr< sched::Mutex> cloned_mutex_sp_b( cloned_mutex_sp_a);

      // ensure that cloned muteces are unlocked
      BCL_ExampleAssert( cloned_mutex_sp_a->TestLock(), false);

      // copy a mutex (this should reinitialize it)
      sched::Mutex mutex_copy( default_mutex);

      // copy a mutex via operator =
      sched::Mutex mutex_assigned;
      mutex_assigned = default_mutex;

    ////////////////
    // operations //
    ////////////////

      // ensure that locking the default mutex has no effect on the cloned muteces or the mutex copy

      // test Lock
      default_mutex.Lock();
      BCL_ExampleIndirectCheck( default_mutex.TestLock(), true, "Lock()");

      // test Unlock
      default_mutex.Unlock();
      BCL_ExampleIndirectCheck( default_mutex.TestLock(), false, "Unlock()");

      // test TryLock
      BCL_ExampleIndirectCheck( default_mutex.TryLock(), true, "TryLock on unlocked mutex");
      BCL_ExampleIndirectCheck( default_mutex.TryLock(), false, "TryLock on locked mutex");
      BCL_ExampleIndirectCheck( default_mutex.TryLock(), false, "TryLock on locked mutex");

      // ensure that the mutex is still locked
      BCL_ExampleIndirectCheck
      (
        default_mutex.TestLock(),
        true,
        "TryLock should keep mutex locked if it was already locked"
      );

      // make sure that none of the other muteces are locked, since copying a mutex should always just create a new one
      BCL_ExampleIndirectCheck( cloned_mutex_sp_a->TestLock(), false, "Clone should not copy the mutex");
      BCL_ExampleIndirectCheck( mutex_copy.TestLock(), false, "Copy constructor should not copy internal mutex");
      BCL_ExampleIndirectCheck( mutex_assigned.TestLock(), false, "Assignment operator should not copy internal mutex");

      // check that locking either shared pointer to the mutex does in fact lock the other mutex, since they should be
      // pointing to the same mutex
      cloned_mutex_sp_a->Lock();
      BCL_ExampleIndirectCheck( cloned_mutex_sp_b->TestLock(), true, "locking a mutex through a shared pointer");

      // unlock all muteces since locked muteces cannot be deleted without triggering an assert
      cloned_mutex_sp_b->Unlock();
      default_mutex.Unlock();

      return 0;
    } // Run

    static const ExampleClass::EnumType ExampleSchedMutex_Instance;

  }; //end ExampleSchedMutex

  const ExampleClass::EnumType ExampleSchedMutex::ExampleSchedMutex_Instance
  (
    GetExamples().AddEnum( ExampleSchedMutex())
  );

} // namespace bcl
