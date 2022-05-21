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

#ifndef BCL_SCHED_MUTEX_H_
#define BCL_SCHED_MUTEX_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mutex
    //! @brief a mutex class; a pthread mutex if pthread libraries are linked, or a critical section on windows
    //!
    //! @see @link example_sched_mutex.cpp @endlink
    //! @author mendenjl
    //! @date Jun 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Mutex :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      //! pointer to internal mutex type
      //! on windows, a critical section
      //! on systems with pthreads, a pthread_mutex_t
      //! on non-windows systems that do not have pthreads, a bool
      //! While we could use preprocessor commands to insert the correct type in the header,
      //! this would result in windows.h being included when building on mingw
      //! windows.h has preprocessor defines and macros that conflict with other functions in the bcl (particularly in chemistry)
      void *m_InternalMutex;

      //! keep track of whether the mutex is currently locked by any thread
      //! this is necessary to allow consistent behavior across architectures, because win32 critical sections allow
      //! re-locking a critical section, whereas POSIX pthread behavior is to forbid it
      bool m_IsLocked;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor; initializes new mutex
      Mutex();

      //! @brief copy constructor; initializes new mutex
      //! @param MUTEX parent mutex; asserts that MUTEX is not locked, since there is no way to copy a locked mutex
      Mutex( const Mutex &MUTEX);

      //! @brief delete operator
      //! @details signals an error if any threads are currently blocked on the mtuex
      virtual ~Mutex();

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      Mutex *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief lock the Mutex; waits for other threads to unlock
      void Lock();

      //! @brief try to lock the mutex
      //! @return bool; whether the mutex is now locked by this thread, which only happens if it is not currently locked
      //! Use trylock if a thread can perform additional work if the mutex is currently locked by another thread
      bool TryLock();

      //! @brief unlock the Mutex
      void Unlock();

      //! @brief check whether the mutex is locked
      //! @return bool; whether the mutex is already locked
      bool TestLock() const;

    //////////////
    // operators//
    //////////////

      //! @brief assignment operator, deletes internal mutex and creates a new one
      //! @param MUTEX parent mutex; both this mutex and MUTEX are asserted to be unlocked
      Mutex &operator=( const Mutex &MUTEX);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief function used internally to delete any existing mutex and create a new one
      //! also asserts that the mutex is not locked, since reinitializing a locked mutex results
      //! in undefined behavior
      void Reinititialize();

    }; // class Mutex

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_MUTEX_H_
