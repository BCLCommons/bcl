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

#ifndef BCL_SCHED_SCOPE_LOCK_H_
#define BCL_SCHED_SCOPE_LOCK_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_mutex.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScopeLock
    //! @brief is a class to lock critical scopes against concurrent thread access
    //!
    //! It can be used very easy in scenarios like
    //! @verbatim
    //! { // begin of critical section
    //!   MutexType m;
    //!   _init_mutex( &m);
    //!   _lock_mutex( &m);
    //!   // do some critical stuff that cannot run concurrent
    //!   _unlock_mutex( &m);
    //!   _destroy_mutex( &m);
    //! } // end of critical section
    //! @endverbatim
    //!
    //! which will easily be transferred into this
    //! @verbatim
    //! { // begin of critical section
    //!   MutexSingleThreaded mt_policy;
    //!   ScopeLock lock( &mt_policy); // constructor locks the mutex
    //!   // do your critical stuff here
    //! } // end of critical section - destructor of ScopeLock will unlock the mutex
    //! @endverbatim
    //!
    //! @see @link example_sched_scope_lock.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScopeLock
    {
    private:

      Mutex m_Mutex; //! the mutex to be used

    public:

      //! @brief construct ScopeLock from MUTEX - locks the mutex
      ScopeLock()
      {
        m_Mutex.Lock();
      }

      //! @brief destruct ScopeLock - unlocks the mutex
      ~ScopeLock()
      {
        m_Mutex.Unlock();
      }

    }; // class ScopeLock

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_SCOPE_LOCK_H_ 
