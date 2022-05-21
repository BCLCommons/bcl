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

#ifndef BCL_SIGNAL_SIGNAL_H_
#define BCL_SIGNAL_SIGNAL_H_

// include the namespace header
#include "bcl_signal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_signal_connection.h"
#include "bcl_signal_signal_base.h"

// external includes - sorted alphabetically
#include <list>

namespace bcl
{
  namespace signal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal0
    //! @brief is a class that was connected to Slots derived classes, that have registered there
    //! member functions that will be invoked if Emit or operator is called. It can call void function with 0 arguments
    //!
    //! @see @link example_signal_signal.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Signal0 :
      public SignalBase0
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Signal0()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL the signal to be copied
      Signal0( const Signal0 &SIGNAL) :
        SignalBase0( SIGNAL)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief connect a SLOT derived class and one of its void, 0 argument member functions to be called on that signal
      //! @param DESTINATION pointer to the destination class instance
      //! @param MEMBER_FUNCTION_PTR pointer to void member function with 0 arguments
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of this Signal, when false, they will be not inherited
      //! @param PUSH_FRONT is true, the Connection will be inserted to the front of the connection list, so that it will be called first, when signal is emitted
      template< typename t_DestinationType>
      void Connect
      (
        t_DestinationType *DESTINATION,
        void ( t_DestinationType::*MEMBER_FUNCTION_PTR)(),
        const bool COPY_CONNECTION = false,
        const bool PUSH_FRONT = false
      )
      {
        // lock that critical section
        sched::ScopeLock lock;

        // create new connection
        Connection0< t_DestinationType> *conn =
        new Connection0< t_DestinationType>( DESTINATION, MEMBER_FUNCTION_PTR, COPY_CONNECTION);

        // normal connection is added to end of connection list
        if( !PUSH_FRONT)
        {
          // add connection to signals slots
          SignalBase0::m_ConnectedSlots.push_back( conn);
        }
        // this connection hs to be notified first, since following notifications might depend on it
        else
        {
          SignalBase0::m_ConnectedSlots.push_front( conn);
        }

        // add connection to destinations slots, so that the destination can cleanup properly when necessary
        DESTINATION->SignalConnect( this);
      }

      //! @brief emit the signal to all connected slots
      void Emit()
      {
        // lock that critical section
        sched::ScopeLock lock;

        // iterators to iterate over all connections
        SignalBase0::const_iterator itr_next, itr = SignalBase0::m_ConnectedSlots.begin();
        SignalBase0::const_iterator itr_end = SignalBase0::m_ConnectedSlots.end();

        // iterate over all connected slots
        while( itr != itr_end)
        {
          itr_next = itr;
          ++itr_next;

          ( *itr)->Emit();

          itr = itr_next;
        }
      }

    //////////////
    // operator //
    //////////////

      //! @brief convenience function to emit signal
      //! @see Emit
      void operator()()
      {
        Emit();
      }

    }; // template class Signal0

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal1
    //! @brief is a class that was connected to Slots derived classes, that have registered there
    //! member functions that will be invoked if Emit or operator is called. It can call void function with 1 arguments
    //!
    //! @tparam t_ArgType1 type of first argument of registered member functions that are called on the signal
    //!
    //! @see Signal0
    //! @see @link example_signal_signal.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1>
    class Signal1 :
      public SignalBase1< t_ArgType1>
    {
    public:

      //! @brief default constructor
      Signal1()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL the signal to be copied
      Signal1( const Signal1< t_ArgType1> &SIGNAL) :
        SignalBase1< t_ArgType1>( SIGNAL)
      {
      }

      //! @brief connect a SLOT derived class and one of its void, 0 argument member functions to be called on that signal
      //! @param DESTINATION pointer to the destination class instance
      //! @param MEMBER_FUNCTION_PTR pointer to void member function with 1 arguments
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of this Signal, when false, they will be not inherited
      //! @param PUSH_FRONT is true, the Connection will be inserted to the front of the connection list, so that it will be called first, when signal is emitted
      template< typename t_DestinationType>
      void Connect
      (
        t_DestinationType *DESTINATION,
        void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1),
        const bool COPY_CONNECTION = false,
        const bool PUSH_FRONT = false
      )
      {
        // lock that critical section
        sched::ScopeLock lock;

        // create new connection
        Connection1< t_DestinationType, t_ArgType1> *
          conn = new Connection1< t_DestinationType, t_ArgType1>( DESTINATION, MEMBER_FUNCTION_PTR, COPY_CONNECTION);

        // normal connection is added to end of connection list
        if( !PUSH_FRONT)
        {
          // add connection to signals slots
          SignalBase1< t_ArgType1>::m_ConnectedSlots.push_back( conn);
        }
        // this connection hs to be notified first, since following notifications might depend on it
        else
        {
          SignalBase1< t_ArgType1>::m_ConnectedSlots.push_front( conn);
        }

        DESTINATION->SignalConnect( this);
      }

      //! @brief emit the signal to all connected slots
      //! @param A1 first argument to the function
      void Emit( t_ArgType1 A1)
      {
        sched::ScopeLock lock;
        typename SignalBase1< t_ArgType1>::const_iterator itr_next, itr = SignalBase1< t_ArgType1>::m_ConnectedSlots.begin();
        typename SignalBase1< t_ArgType1>::const_iterator itr_end = SignalBase1< t_ArgType1>::m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          itr_next = itr;
          ++itr_next;

          ( *itr)->Emit( A1);

          itr = itr_next;
        }
      }

      //! @brief convenience function to emit signal
      //! @see Emit
      void operator()( t_ArgType1 A1)
      {
        Emit( A1);
      }

    }; // template class Signal1

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal2
    //! @brief is a class that was connected to Slots derived classes, that have registered there
    //! member functions that will be invoked if Emit or operator is called. It can call void function with 2 arguments
    //!
    //! @see Signal0
    //! @see @link example_signal_signal.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2>
    class Signal2 :
      public SignalBase2< t_ArgType1, t_ArgType2>
    {
    public:

      //! @brief default constructor
      Signal2()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL the signal to be copied
      Signal2( const Signal2< t_ArgType1, t_ArgType2> &SIGNAL) :
        SignalBase2< t_ArgType1, t_ArgType2>( SIGNAL)
      {
      }

      //! @brief connect a SLOT derived class and one of its void, 0 argument member functions to be called on that signal
      //! @param DESTINATION pointer to the destination class instance
      //! @param MEMBER_FUNCTION_PTR pointer to void member function with 2 arguments
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of this Signal, when false, they will be not inherited
      //! @param PUSH_FRONT is true, the Connection will be inserted to the front of the connection list, so that it will be called first, when signal is emitted
      template< typename t_DestinationType>
      void Connect
      (
        t_DestinationType *DESTINATION,
        void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2),
        const bool COPY_CONNECTION = false,
        const bool PUSH_FRONT = false
      )
      {
        // lock that critical section
        sched::ScopeLock lock;

        // create new connection
        Connection2< t_DestinationType, t_ArgType1, t_ArgType2> *
          conn = new Connection2< t_DestinationType, t_ArgType1, t_ArgType2>( DESTINATION, MEMBER_FUNCTION_PTR);

        // normal connection is added to end of connection list
        if( !PUSH_FRONT)
        {
          // add connection to signals slots
          SignalBase2< t_ArgType1, t_ArgType2>::m_ConnectedSlots.push_back( conn);
        }
        // this connection hs to be notified first, since following notifications might depend on it
        else
        {
          SignalBase2< t_ArgType1, t_ArgType2>::m_ConnectedSlots.push_front( conn);
        }

        DESTINATION->SignalConnect( this);
      }

      //! @brief emit the signal to all connected slots
      //! @param A1 first argument to the function
      //! @param A2 second argument to the function
      void Emit( t_ArgType1 A1, t_ArgType2 A2)
      {
        sched::ScopeLock lock;
        typename SignalBase2< t_ArgType1, t_ArgType2>::const_iterator itr_next, itr = SignalBase2< t_ArgType1, t_ArgType2>::m_ConnectedSlots.begin();
        typename SignalBase2< t_ArgType1, t_ArgType2>::const_iterator itr_end = SignalBase2< t_ArgType1, t_ArgType2>::m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          itr_next = itr;
          ++itr_next;

          ( *itr)->Emit( A1, A2);

          itr = itr_next;
        }
      }

      //! @brief convenience function to emit signal
      //! @see Emit
      void operator()( t_ArgType1 A1, t_ArgType2 A2)
      {
        Emit( A1, A2);
      }

    }; // template class Signal2

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal3
    //! @brief is a class that was connected to Slots derived classes, that have registered there
    //! member functions that will be invoked if Emit or operator is called. It can call void function with 3 arguments
    //!
    //! @see Signal0
    //! @see @link example_signal_signal.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class Signal3 :
      public SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>
    {
    public:

      //! @brief default constructor
      Signal3()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL the signal to be copied
      Signal3( const Signal3< t_ArgType1, t_ArgType2, t_ArgType3> &SIGNAL) :
        SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>( SIGNAL)
      {
      }

      //! @brief connect a SLOT derived class and one of its void, 0 argument member functions to be called on that signal
      //! @param DESTINATION pointer to the destination class instance
      //! @param MEMBER_FUNCTION_PTR pointer to void member function with 3 arguments
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of this Signal, when false, they will be not inherited
      //! @param PUSH_FRONT is true, the Connection will be inserted to the front of the connection list, so that it will be called first, when signal is emitted
      template< typename t_DestinationType>
      void Connect
      (
        t_DestinationType *DESTINATION,
        void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2, t_ArgType3),
        const bool COPY_CONNECTION = false,
        const bool PUSH_FRONT = false
      )
      {
        // lock that critical section
        sched::ScopeLock lock;

        // create new connection
        Connection3< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3> *
          conn = new Connection3< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3>
                 (
                   DESTINATION, MEMBER_FUNCTION_PTR, COPY_CONNECTION
                 );

        // normal connection is added to end of connection list
        if( !PUSH_FRONT)
        {
          // add connection to signals slots
          SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::m_ConnectedSlots.push_back( conn);
        }
        // this connection hs to be notified first, since following notifications might depend on it
        else
        {
          SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::m_ConnectedSlots.push_front( conn);
        }

        DESTINATION->SignalConnect( this);
      }

      //! @brief emit the signal to all connected slots
      //! @param A1 first argument to the function
      //! @param A2 second argument to the function
      //! @param A3 third argument to the function
      void Emit( t_ArgType1 A1, t_ArgType2 A2, t_ArgType3 A3)
      {
        sched::ScopeLock lock;
        typename SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::const_iterator itr_next, itr = SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::m_ConnectedSlots.begin();
        typename SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::const_iterator itr_end = SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3>::m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          itr_next = itr;
          ++itr_next;

          ( *itr)->Emit( A1, A2, A3);

          itr = itr_next;
        }
      }

      //! @brief convenience function to emit signal
      //! @see Emit
      void operator()( t_ArgType1 A1, t_ArgType2 A2, t_ArgType3 A3)
      {
        Emit( A1, A2, A3);
      }

    }; // template class Signal3

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal4
    //! @brief is a class that was connected to Slots derived classes, that have registered there
    //! member functions that will be invoked if Emit or operator is called. It can call void function with 4 arguments
    //!
    //! @tparam t_ArgType1 type of first argument of registered member functions that are called on the signal
    //!
    //! @see Signal0
    //! @see @link example_signal_signal.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class Signal4 :
      public SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>
    {
    public:

      //! @brief default constructor
      Signal4()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL the signal to be copied
      Signal4( const Signal4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> &SIGNAL) :
        SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>( SIGNAL)
      {
      }

      //! @brief connect a SLOT derived class and one of its void, 0 argument member functions to be called on that signal
      //! @param DESTINATION pointer to the destination class instance
      //! @param MEMBER_FUNCTION_PTR pointer to void member function with 4 arguments
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of this Signal, when false, they will be not inherited
      //! @param PUSH_FRONT is true, the Connection will be inserted to the front of the connection list, so that it will be called first, when signal is emitted
      template< typename t_DestinationType>
      void Connect
      (
        t_DestinationType *DESTINATION,
        void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4),
        const bool COPY_CONNECTION = false,
        const bool PUSH_FRONT = false
      )
      {
        // lock that critical section
        sched::ScopeLock lock;

        // create new connection
        Connection4< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> *
          conn = new Connection4< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>
                 (
                   DESTINATION, MEMBER_FUNCTION_PTR, COPY_CONNECTION
                 );

        // normal connection is added to end of connection list
        if( !PUSH_FRONT)
        {
          // add connection to signals slots
          SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::m_ConnectedSlots.push_back( conn);
        }
        // this connection hs to be notified first, since following notifications might depend on it
        else
        {
          SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::m_ConnectedSlots.push_front( conn);
        }

        DESTINATION->SignalConnect( this);
      }

      //! @brief emit the signal to all connected slots
      //! @param A1 first argument to the function
      //! @param A2 second argument to the function
      //! @param A3 third argument to the function
      //! @param A4 fourth argument to the function
      void Emit( t_ArgType1 A1, t_ArgType2 A2, t_ArgType3 A3, t_ArgType4 A4)
      {
        sched::ScopeLock lock;
        typename SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::const_iterator itr_next, itr = SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::m_ConnectedSlots.begin();
        typename SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::const_iterator itr_end = SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>::m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          itr_next = itr;
          ++itr_next;

          ( *itr)->Emit( A1, A2, A3, A4);

          itr = itr_next;
        }
      }

      //! @brief convenience function to emit signal
      //! @see Emit
      void operator()( t_ArgType1 A1, t_ArgType2 A2, t_ArgType3 A3, t_ArgType4 A4)
      {
        Emit( A1, A2, A3, A4);
      }

    }; // template class Signal4

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_SIGNAL_H_
