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
#include "signal/bcl_signal_slots.h"

// includes from bcl - sorted alphabetically
#include "signal/bcl_signal_signal_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace signal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Slots::Slots()
    {
    }

    //! @brief copy constructor
    //! duplicates all the senders to this Slots
    //! @param SLOTS the argument slots
    Slots::Slots( const Slots &SLOTS)
    {
      sched::ScopeLock lock;
      std::set< SignalBase *>::const_iterator it = SLOTS.m_Senders.begin();
      std::set< SignalBase *>::const_iterator itEnd = SLOTS.m_Senders.end();

      while( it != itEnd)
      {
        if( ( *it)->SlotDuplicate( &SLOTS, this))
        {
          m_Senders.insert( *it);
        }
        ++it;
      }
    }

    //! @brief destructor disconnects all sender
    Slots::~Slots()
    {
      DisconnectAll();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief connect a signal
    //! @brief SENDER the signal source
    void Slots::SignalConnect( SignalBase *SENDER)
    {
      sched::ScopeLock lock;
      m_Senders.insert( SENDER);
    }

    //! @brief disconnect a given signal
    //! @brief SENDER the source of the signal that needs to be disconnected
    void Slots::SignalDisconnect( SignalBase *SENDER)
    {
      sched::ScopeLock lock;
      m_Senders.erase( SENDER);
    }

    //! @brief disconnects all senders
    //! this is necessary so that this slot does not receives signals if it gets destroyed
    void Slots::DisconnectAll()
    {
      sched::ScopeLock lock;
      std::set< SignalBase *>::const_iterator it = m_Senders.begin();
      std::set< SignalBase *>::const_iterator itEnd = m_Senders.end();

      // disconnect signals
      while( it != itEnd)
      {
        ( *it)->SlotDisconnect( this);
        ++it;
      }

      // remove all signals from sender list
      m_Senders.erase( m_Senders.begin(), m_Senders.end());
    }

  } // namespace signal
} // namespace bcl
