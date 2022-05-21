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

#ifndef BCL_UTIL_CLEANABLE_INTERFACE_H_
#define BCL_UTIL_CLEANABLE_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "sched/bcl_sched_mutex.h"

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CleanableInterface
    //! @brief derives classes that require a cleanup even on abnormal termination by exit()
    //! @details normally, when objects are declared in a scope, and the scope is passed, those objects are destructed,
    //! and their destructor is called. This is not happening, when a program is terminated with a call to exit(), which
    //! is called in the BCL_Assert or when using external libraries, like boinc that has a graphics loop that never exits
    //!
    //! This cause problems, e.g. their is an open file or a databank connection. Those will never be closed, and can
    //! cause trouble over time.
    //! Even more sever are Semaphores, that have a Lock, and do not unlock, and other processes are waiting for them.
    //!
    //! The CleanableInterface only defines one function, void CleanUp(), that needs to be overwritten in the derived
    //! classes. The implementation takes care of registering each instance of that class with the nested Cleanables
    //! class, that stores pointer to all instances of cleanables. The first time, the Cleanables nested class is
    //! instantiated, the static Cleanup of the Cleanables class is registered with the POSIX atexit( &CleanUp) call.
    //! this CleanUp function, when called, iterates over all remaining instances of Cleanable derived classes, and
    //! and calls the Cleanup function on them.
    //!
    //! Each derived class implements the Cleanup function in a way, that makes sure, that e.g. Databank connections are
    //! closed, file handles are closed, or Semaphores are unlocked and released.
    //!
    //! No o meaningful example can be written, since only the when exit() is called, the effect will be visible
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Feb 25, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CleanableInterface
    {

    private:

    //////////
    // data //
    //////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Cleanables
      //! @brief internal singleton class that stores the objects that have been registered as cleanable
      //!
      //! @remarks example unnecessary
      //! @author woetzen
      //! @date Feb 25, 2010
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class Cleanables
      {
      private:

      //////////
      // data //
      //////////

        sched::Mutex                    m_Mutex;
        std::set< CleanableInterface *> m_Cleanables;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief private constructor
        //! registers the overall cleanup function
        Cleanables();

      public:

      //////////
      // data //
      //////////

        //! @brief access to the only instance of Cleanables
        static Cleanables &GetCleanables();

      ////////////////
      // operations //
      ////////////////

        //! @brief register the CLEANABLE
        void Register( CleanableInterface *CLEANABLE);

        //! @brief unregister the CLEANABLE
        void UnRegister( CleanableInterface *CLEANABLE);

      private:

        //! @brief clean all cleanables that are registered
        void CleanAllCleanables();

      //////////////////////
      // helper functions //
      //////////////////////

        //! @brief Cleanup all cleanables
        static void AtExitClean();

        //! @brief signal handler, for non exit calls
        //! @param SIGNAL
        static void SignalHandler( int SIGNAL);

      }; // class Cleanables

    protected:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! registers this with the cleanables
      CleanableInterface();

      //! @brief virtual destructor
      virtual ~CleanableInterface();

    ////////////////
    // operations //
    ////////////////

      //! @brief virtual function that will be called if cleanable is registered
      virtual void CleanUp() = 0;

    }; // class CleanableInterface

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_CLEANABLE_INTERFACE_H_
