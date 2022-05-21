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
#include "sched/bcl_sched_mutex.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_assert.h"
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

// define macros to handle mutex initialization, deletion, lock, trylock, and unlock
// the logic is the same no matter which of these function calls are used, the macros are used
// to abstract away the differences between different systems and/or libraries
#ifdef BCL_WITHOUT_PTHREAD
  //! single threaded mutex; just keeps track of whether or not it is locked
  //! this allows primitive error checking; e.g. that a locked mutex is eventually unlocked
  # define MUTEX_TYPE           bool
  # define INIT_MUTEX(x)        *x = false
  # define DELETE_MUTEX_DATA(x) *x = false
  # define LOCK_MUTEX(x)        BCL_Assert( *x != true, "cannot lock a mutex that was already locked"); *x = true;
  # define TRYLOCK_MUTEX(x)     *x = true
  # define UNLOCK_MUTEX(x)      *x = false
#elif defined(_WIN32)
  # include <windows.h>
  //! the windows critical section struct - similar to pthread_mutex_t
  # define MUTEX_TYPE            CRITICAL_SECTION
  # define INIT_MUTEX(x)         InitializeCriticalSection( x)
  # define DELETE_MUTEX_DATA(x)  DeleteCriticalSection( x)
  # define LOCK_MUTEX(x)         EnterCriticalSection( x)
  # define TRYLOCK_MUTEX(x)      ( TryEnterCriticalSection( x) != 0)
  # define UNLOCK_MUTEX(x)       LeaveCriticalSection( x)
#else
  // posix muteces
  # include <pthread.h>
  //! pthread mutex
  # define MUTEX_TYPE            pthread_mutex_t
  # define INIT_MUTEX(x)         pthread_mutex_init( x, NULL)
  # define DELETE_MUTEX_DATA(x)  pthread_mutex_destroy( x)
  # define LOCK_MUTEX(x)         pthread_mutex_lock( x)
  # define TRYLOCK_MUTEX(x)      ( pthread_mutex_trylock( x) == 0)
  # define UNLOCK_MUTEX(x)       pthread_mutex_unlock( x)
#endif

namespace bcl
{
  namespace sched
  {

    //! @brief default constructor; initializes new mutex
    Mutex::Mutex() :
      m_InternalMutex( new MUTEX_TYPE),
      m_IsLocked( false)
    {
      INIT_MUTEX( static_cast< MUTEX_TYPE*>( m_InternalMutex));
    }

    //! @brief copy constructor; initializes new mutex
    //! @param MUTEX parent mutex
    Mutex::Mutex( const Mutex &MUTEX) :
      m_InternalMutex( new MUTEX_TYPE),
      m_IsLocked( false)
    {
      // ensure that the mutex isn't currently locked by using testlock
      BCL_Assert( !MUTEX.TestLock(), "cannot copy a locked mutex");
      INIT_MUTEX( static_cast< MUTEX_TYPE *>( m_InternalMutex));
    }

    //! @brief delete operator
    //! @details signals an error if any threads are currently blocked on the mtuex
    Mutex::~Mutex()
    {
      // ensure that the mutex isn't currently locked
      BCL_Assert( !m_IsLocked, "cannot delete a mutex that is still locked");
      // now destroy the mutex
      DELETE_MUTEX_DATA( static_cast< MUTEX_TYPE *>( m_InternalMutex));
      delete static_cast< MUTEX_TYPE *>( m_InternalMutex);
    }

    //! @brief Clone function
    //! @return new Pointer to a copy of the actual object behind the pointer
    Mutex *Mutex::Clone() const
    {
      return new Mutex();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Mutex::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief lock the Mutex; waits for other threads to unlock
    void Mutex::Lock()
    {
      LOCK_MUTEX( static_cast< MUTEX_TYPE *>( m_InternalMutex));
      m_IsLocked = true;
    }

    //! @brief try to lock the mutex
    //! @return bool; whether the mutex is now locked by this thread, which only happens if it is not currently locked
    //! Use trylock if a thread can perform additional work if the mutex is currently locked by another thread
    bool Mutex::TryLock()
    {
      if( !m_IsLocked) // only try to lock if the mutex is not already locked
      {
        // try to lock the mutex
        const bool lock_success( TRYLOCK_MUTEX( static_cast< MUTEX_TYPE *>( m_InternalMutex)));

        // it may happen that another thread locked the mutex first, in which case the above call failed
        // so only update m_IsLocked if we succeeded in making the lock
        if( lock_success)
        {
          m_IsLocked = true;
        }

        // return whether this thread succeeded in locking the mutex
        return lock_success;
      }

      // mutex was already locked, just return false
      return false;
    }

    //! @brief unlock the Mutex
    void Mutex::Unlock()
    {
      UNLOCK_MUTEX( static_cast< MUTEX_TYPE *>( m_InternalMutex));
      m_IsLocked = false;
    }

    //! @brief check whether the mutex is locked
    //! @return bool; whether the mutex is already locked
    bool Mutex::TestLock() const
    {
      return m_IsLocked;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator, deletes internal mutex and creates a new one
    //! @param MUTEX parent mutex
    Mutex &Mutex::operator=( const Mutex &MUTEX)
    {
      // ensure that the mutex isn't currently locked by using TestLock
      BCL_Assert( !MUTEX.TestLock(), "cannot call operator= with a locked mutex");
      Reinititialize();
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mutex::Read( std::istream &ISTREAM)
    {
      Reinititialize();
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &Mutex::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // warn the user if the mutex is still locked while writing it.  We can't lock the mutex
      // when reading it in because we don't know which thread has access to the mutex and which
      // are currently blocked by it.  In general, it would be dangerous to assume that the thread
      // that reads the mutex should be the one to lock it.
      if( m_IsLocked)
      {
        BCL_MessageVrb
        (
          "Writing a mutex that was still locked! When read, the mutex will be unlocked"
        );
      }
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function used internally to delete any existing mutex and create a new one
    void Mutex::Reinititialize()
    {
      // make sure the mutex isn't already locked
      BCL_Assert( !m_IsLocked, "cannot reinitialize a mutex (via operator= or Read) that is still locked");
      // destroy the mutex data
      DELETE_MUTEX_DATA( static_cast< MUTEX_TYPE *>( m_InternalMutex));
      // and reinitialize it
      INIT_MUTEX( static_cast< MUTEX_TYPE *>( m_InternalMutex));
    }

  } // namespace sched
} // namespace bcl

// These undefs are not strictly necessary because this is a source file
// However, we follow the practice that preprocessor defines/macros used only in one
// file should always be undef'ed at the end of that file to prevent collisions / unexpected behavior elsewhere
#undef MUTEX_TYPE
#undef INIT_MUTEX
#undef DELETE_MUTEX_DATA
#undef LOCK_MUTEX
#undef TRYLOCK_MUTEX
#undef UNLOCK_MUTEX
