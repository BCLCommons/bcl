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
#include "util/bcl_util_cleanable_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically
#include <map>
#include <signal.h>
#include <stdlib.h>

#ifndef _SIGHANDLER_T_DEFINED
#define _SIGHANDLER_T_DEFINED
//! @brief signal_handler_t
typedef void( *sighandler_t)( int);
#endif

namespace bcl
{
  namespace util
  {

    //! @brief access to singelton map that stores previous signal handlers
    //! @return map that has the signal handler function pointer for the bcl overwritten signal handler
    static std::map< int, sighandler_t> &GetSignalHandlerMap()
    {
      static std::map< int, sighandler_t> *s_signal_handler_map( new std::map< int, sighandler_t>());

      return *s_signal_handler_map;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief private constructor
    //! registers the overall cleanup function
    CleanableInterface::Cleanables::Cleanables()
    {
      atexit( &AtExitClean);

      #ifndef BCL_NO_OS_SIGNAL_HANDLING
      GetSignalHandlerMap()[ SIGILL ] = ::signal( SIGILL , &SignalHandler);
      GetSignalHandlerMap()[ SIGINT ] = ::signal( SIGINT , &SignalHandler);
      GetSignalHandlerMap()[ SIGFPE ] = ::signal( SIGFPE , &SignalHandler);
      GetSignalHandlerMap()[ SIGABRT] = ::signal( SIGABRT, &SignalHandler);
      GetSignalHandlerMap()[ SIGSEGV] = ::signal( SIGSEGV, &SignalHandler);
      GetSignalHandlerMap()[ SIGTERM] = ::signal( SIGTERM, &SignalHandler);
      #endif
    }

  //////////
  // data //
  //////////

    //! @brief access to the only instance of Cleanables
    CleanableInterface::Cleanables &CleanableInterface::Cleanables::GetCleanables()
    {
      // only instance of cleanables - needs to exist to he absolute end, so it is a new pointer, that will never be
      // deleted and only cleaned up by the operating system
      static Cleanables *s_cleanables( new Cleanables());

      // return
      return *s_cleanables;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief register the CLEANABLE
    void CleanableInterface::Cleanables::Register( CleanableInterface *CLEANABLE)
    {
      m_Mutex.Lock();
      m_Cleanables.insert( CLEANABLE);
      m_Mutex.Unlock();
    }

    //! @brief unregister the CLEANABLE
    void CleanableInterface::Cleanables::UnRegister( CleanableInterface *CLEANABLE)
    {
      m_Mutex.Lock();
      m_Cleanables.erase( CLEANABLE);
      m_Mutex.Unlock();
    }

    //! @brief clean all cleanables that are registered
    void CleanableInterface::Cleanables::CleanAllCleanables()
    {
      // iterate over all cleanables and call clean
      while( !m_Cleanables.empty())
      {
        CleanableInterface *current( *m_Cleanables.begin());
        m_Cleanables.erase( m_Cleanables.begin());
        current->CleanUp();
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Cleanup all cleanables
    void CleanableInterface::Cleanables::AtExitClean()
    {
      GetCleanables().CleanAllCleanables();
    }

    //! @brief signal handler, for non exit calls
    //! @param SIGNAL
    void CleanableInterface::Cleanables::SignalHandler( int SIGNAL)
    {
      static volatile sig_atomic_t flag( 0);
      if( flag != 0)
      {
        return;
      }

      // set flag, so that signal is not processed twice
      flag = 1;
      BCL_MessageTop
      (
        "caught signal: " + util::Format()( SIGNAL) + " cleaning! " + util::CallStack().String()
      );

      // clean all the things that needed to be cleaned
      GetCleanables().CleanAllCleanables();
      std::map< int, sighandler_t>::iterator itr( GetSignalHandlerMap().find( SIGNAL));

      // reset the signal handler
      ::signal( SIGNAL, itr->second);

      // remove the itr
      GetSignalHandlerMap().erase( itr);

      // raise the signal again
      raise( SIGNAL);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! registers this with the cleanables
    CleanableInterface::CleanableInterface()
    {
      Cleanables::GetCleanables().Register( this);
    }

    //! @brief virtual destructor
    CleanableInterface::~CleanableInterface()
    {
      Cleanables::GetCleanables().UnRegister( this);
    }

  } // namespace util
} // namespace bcl
