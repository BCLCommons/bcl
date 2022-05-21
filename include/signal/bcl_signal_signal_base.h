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

#ifndef BCL_SIGNAL_SIGNAL_BASE_H_
#define BCL_SIGNAL_SIGNAL_BASE_H_

// include the namespace header
#include "bcl_signal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_signal_connection_interface.h"
#include "bcl_signal_slots.h"
#include "sched/bcl_sched_scope_lock.h"

// external includes - sorted alphabetically
#include <list>

namespace bcl
{
  namespace signal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase
    //! @brief SignalBase provides functionality for destinations, that need to be notified to connect to them
    //! @details derived classes do maintain a list of connections to destinations, that will get notified, once a
    //!          signal is emitted
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SignalBase
    {
    public:

      //! @brief default constructor
      SignalBase()
      {
      }

      virtual void SlotDisconnect( Slots *SLOT_PTR) = 0;
      virtual bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR) = 0;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase0
    //! @brief TODO: add brief class comment
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SignalBase0 :
      public SignalBase
    {
    public:

      typedef std::list< ConnectionInterface0 *> ConnectionsListType;
      typedef ConnectionsListType::const_iterator const_iterator;
      typedef ConnectionsListType::iterator iterator;

    protected:

    //////////
    // data //
    //////////

      ConnectionsListType m_ConnectedSlots;

    public:

      //! @brief default constructor
      SignalBase0()
      {
      }

      //! @brief copy constructor
      //! @param SIGNAL_BASE source to be copied
      SignalBase0( const SignalBase0 &SIGNAL_BASE) :
        SignalBase( SIGNAL_BASE)
      {
        // lock the scope
        sched::ScopeLock lock;
        const_iterator itr = SIGNAL_BASE.m_ConnectedSlots.begin();
        const_iterator itr_end = SIGNAL_BASE.m_ConnectedSlots.end();

        // iterate over all connections
        while( itr != itr_end)
        {
          // if they can be copied
          if( ( *itr)->CanBeCopied())
          {
            // connect this to the destination
            ( *itr)->GetDestination()->SignalConnect( this);

            // insert a clone of the connection
            m_ConnectedSlots.push_back( ( *itr)->Clone());
          }

          ++itr;
        }
      }

      //! @brief destructor
      //! @details dsiconnects all connections properly
      virtual ~SignalBase0()
      {
        DisconnectAll();
      }

      //! @brief disconnect all destinations and remove connections
      void DisconnectAll()
      {
        // lock the scope
        sched::ScopeLock lock;

        const_iterator itr = m_ConnectedSlots.begin();
        const_iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          // call SignalDisconnect on the destination with this as argument
          ( *itr)->GetDestination()->SignalDisconnect( this);

          // delete the connection
          delete *itr;

          ++itr;
        }

        m_ConnectedSlots.erase( m_ConnectedSlots.begin(), m_ConnectedSlots.end());
      }

      //! @brief disconnect connection from a destination, which is a Slots derived class
      //! @param SLOTS_PTR ptr to the slots derived object
      void Disconnect( Slots *SLOTS_PTR)
      {
        // lock the scope
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        // iterate over connections
        while( itr != itr_end)
        {
          // if the connections destination matches the given SLOTS_PTR
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            // delete the connection and remove it
            delete *itr;
            m_ConnectedSlots.erase( itr);

            // disconnect this from the slots
            SLOTS_PTR->SignalDisconnect( this);

            // only one connection to each slot
            return;
          }

          ++itr;
        }
      }

      //! @brief called by a sender to request disconnection
      //! @param SLOTS_PTR Slots derived class that requests the disconnection
      void SlotDisconnect( Slots *SLOTS_PTR)
      {
        // lock the scope
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        // iterate over all connections
        while( itr != itr_end)
        {
          // copy to the next
          iterator itr_next = itr;
          ++itr_next;

          // delete and remove, if the connection point to the given slot
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
          }

          itr = itr_next;
        }
      }

      bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        bool copied( false);
        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == OLD_SLOT_PTR && ( *itr)->CanBeCopied())
          {
            m_ConnectedSlots.push_back( ( *itr)->Duplicate( NEW_SLOT_PTR));
            copied = true;
          }

          ++itr;
        }

        return copied;
      }

    }; // class SignalBase0

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase1
    //! @brief TODO: add brief class comment
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1>
    class SignalBase1 :
      public SignalBase
    {
    public:

      typedef std::list< ConnectionInterface1< t_ArgType1> *> ConnectionsListType;
      typedef typename ConnectionsListType::const_iterator const_iterator;
      typedef typename ConnectionsListType::iterator iterator;

    protected:

    //////////
    // data //
    //////////

      //! @brief list of slots that are connected and will recive the signal on emit
      ConnectionsListType m_ConnectedSlots;

    public:

      SignalBase1()
      {
      }

      SignalBase1( const SignalBase1< t_ArgType1> &SIGNAL_BASE) :
        SignalBase( SIGNAL_BASE)
      {
      }

      bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        bool copied( false);
        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == OLD_SLOT_PTR && ( *itr)->CanBeCopied())
          {
            m_ConnectedSlots.push_back( ( *itr)->Duplicate( NEW_SLOT_PTR));
            copied = true;
          }

          ++itr;
        }

        return copied;
      }

      virtual ~SignalBase1()
      {
        DisconnectAll();
      }

      void DisconnectAll()
      {
        sched::ScopeLock lock;
        const_iterator itr = m_ConnectedSlots.begin();
        const_iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          ( *itr)->GetDestination()->SignalDisconnect( this);
          delete *itr;

          ++itr;
        }

        m_ConnectedSlots.erase( m_ConnectedSlots.begin(), m_ConnectedSlots.end());
      }

      void Disconnect( Slots *SLOTS_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
            SLOTS_PTR->SignalDisconnect( this);
            return;
          }

          ++itr;
        }
      }

      void SlotDisconnect( Slots *SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          iterator itr_next = itr;
          ++itr_next;

          if( ( *itr)->GetDestination() == SLOT_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
          }

          itr = itr_next;
        }
      }

    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase2
    //! @brief TODO: add brief class comment
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2>
    class SignalBase2 :
      public SignalBase
    {
    public:

      typedef std::list< ConnectionInterface2< t_ArgType1, t_ArgType2> *> ConnectionsListType;
      typedef typename ConnectionsListType::const_iterator const_iterator;
      typedef typename ConnectionsListType::iterator iterator;

    protected:

    //////////
    // data //
    //////////

      ConnectionsListType m_ConnectedSlots;

    public:

      SignalBase2()
      {
      }

      SignalBase2( const SignalBase2< t_ArgType1, t_ArgType2> &SIGNAL_BASE) :
        SignalBase( SIGNAL_BASE)
      {
        sched::ScopeLock lock;
        const_iterator itr = SIGNAL_BASE.m_ConnectedSlots.begin();
        const_iterator itr_end = SIGNAL_BASE.m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->CanBeCopied())
          {
            ( *itr)->GetDestination()->SignalConnect( this);
            m_ConnectedSlots.push_back( ( *itr)->Clone());
          }

          ++itr;
        }
      }

      bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        bool copied( false);
        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == OLD_SLOT_PTR && ( *itr)->CanBeCopied())
          {
            m_ConnectedSlots.push_back( ( *itr)->Duplicate( NEW_SLOT_PTR));
            copied = true;
          }

          ++itr;
        }

        return copied;
      }

      virtual ~SignalBase2()
      {
        DisconnectAll();
      }

      void DisconnectAll()
      {
        sched::ScopeLock lock;
        const_iterator itr = m_ConnectedSlots.begin();
        const_iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          ( *itr)->GetDestination()->SignalDisconnect( this);
          delete *itr;

          ++itr;
        }

        m_ConnectedSlots.erase( m_ConnectedSlots.begin(), m_ConnectedSlots.end());
      }

      void Disconnect( Slots *SLOTS_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
            SLOTS_PTR->SignalDisconnect( this);
            return;
          }

          ++itr;
        }
      }

      void SlotDisconnect( Slots *SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          iterator itr_next = itr;
          ++itr_next;

          if( ( *itr)->GetDestination() == SLOT_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
          }

          itr = itr_next;
        }
      }

    }; // template class SignalBase2

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase3
    //! @brief TODO: add brief class comment
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class SignalBase3 :
      public SignalBase
    {
    public:

      typedef std::list< ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3> *> ConnectionsListType;
      typedef typename ConnectionsListType::const_iterator const_iterator;
      typedef typename ConnectionsListType::iterator iterator;

    protected:

    //////////
    // data //
    //////////

      ConnectionsListType m_ConnectedSlots;

    public:
      SignalBase3()
      {
      }

      SignalBase3( const SignalBase3< t_ArgType1, t_ArgType2, t_ArgType3> &SIGNAL_BASE) :
        SignalBase( SIGNAL_BASE)
      {
        sched::ScopeLock lock;
        const_iterator itr = SIGNAL_BASE.m_ConnectedSlots.begin();
        const_iterator itr_end = SIGNAL_BASE.m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->CanBeCopied())
          {
            ( *itr)->GetDestination()->SignalConnect( this);
            m_ConnectedSlots.push_back( ( *itr)->Clone());
          }

          ++itr;
        }
      }

      bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        bool copied( false);
        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == OLD_SLOT_PTR && ( *itr)->CanBeCopied())
          {
            m_ConnectedSlots.push_back( ( *itr)->Duplicate( NEW_SLOT_PTR));
            copied = true;
          }

          ++itr;
        }

        return copied;
      }

      virtual ~SignalBase3()
      {
        DisconnectAll();
      }

      void DisconnectAll()
      {
        sched::ScopeLock lock;
        const_iterator itr = m_ConnectedSlots.begin();
        const_iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          ( *itr)->GetDestination()->SignalDisconnect( this);
          delete *itr;

          ++itr;
        }

        m_ConnectedSlots.erase( m_ConnectedSlots.begin(), m_ConnectedSlots.end());
      }

      void Disconnect( Slots *SLOTS_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
            SLOTS_PTR->SignalDisconnect( this);
            return;
          }

          ++itr;
        }
      }

      void SlotDisconnect( Slots *SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          iterator itr_next = itr;
          ++itr_next;

          if( ( *itr)->GetDestination() == SLOT_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
          }

          itr = itr_next;
        }
      }

    }; // template class SignalBase3

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SignalBase4
    //! @brief TODO: add brief class comment
    //!
    //! @remarks example_unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class SignalBase4 :
      public SignalBase
    {
    public:

      typedef std::list< ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> *> ConnectionsListType;
      typedef typename ConnectionsListType::const_iterator const_iterator;
      typedef typename ConnectionsListType::iterator iterator;

    protected:

    //////////
    // data //
    //////////

      ConnectionsListType m_ConnectedSlots;

    public:

      SignalBase4()
      {
      }

      SignalBase4( const SignalBase4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> &SIGNAL_BASE) :
        SignalBase( SIGNAL_BASE)
      {
        sched::ScopeLock lock;
        const_iterator itr = SIGNAL_BASE.m_ConnectedSlots.begin();
        const_iterator itr_end = SIGNAL_BASE.m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->CanBeCopied())
          {
            ( *itr)->GetDestination()->SignalConnect( this);
            m_ConnectedSlots.push_back( ( *itr)->Clone());
          }

          ++itr;
        }
      }

      bool SlotDuplicate( const Slots *OLD_SLOT_PTR, Slots *NEW_SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        bool copied( false);
        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == OLD_SLOT_PTR && ( *itr)->CanBeCopied())
          {
            m_ConnectedSlots.push_back( ( *itr)->Duplicate( NEW_SLOT_PTR));
            copied = true;
          }

          ++itr;
        }

        return copied;
      }

      virtual ~SignalBase4()
      {
        DisconnectAll();
      }

      void DisconnectAll()
      {
        sched::ScopeLock lock;
        const_iterator itr = m_ConnectedSlots.begin();
        const_iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          ( *itr)->GetDestination()->SignalDisconnect( this);
          delete *itr;

          ++itr;
        }

        m_ConnectedSlots.erase( m_ConnectedSlots.begin(), m_ConnectedSlots.end());
      }

      void Disconnect( Slots *SLOTS_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          if( ( *itr)->GetDestination() == SLOTS_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
            SLOTS_PTR->SignalDisconnect( this);
            return;
          }

          ++itr;
        }
      }

      void SlotDisconnect( Slots *SLOT_PTR)
      {
        sched::ScopeLock lock;
        iterator itr = m_ConnectedSlots.begin();
        iterator itr_end = m_ConnectedSlots.end();

        while( itr != itr_end)
        {
          iterator itr_next = itr;
          ++itr_next;

          if( ( *itr)->GetDestination() == SLOT_PTR)
          {
            delete *itr;
            m_ConnectedSlots.erase( itr);
          }

          itr = itr_next;
        }
      }

    }; // template class SignalBase4

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_SIGNAL_BASE_H_
