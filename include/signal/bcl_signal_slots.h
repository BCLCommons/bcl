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

#ifndef BCL_SIGNAL_SLOTS_H_
#define BCL_SIGNAL_SLOTS_H_

// include the namespace header
#include "bcl_signal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace signal
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Slots
    //! @brief is a base class for classes, that need to be able to handle signals
    //! @details Signals are emitted from class with Signal handlers - e.g. a Light is derived from Slots, and registers
    //! its member functions like "TurnOn()" with a Signal member in a class "Switch". When the Switch emits a signal
    //! through the Signal member, the functions will be called on the Light - the class with the slots.
    //!
    //! @see @link example_signal_slots.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Slots
    {
    private:

    //////////
    // data //
    //////////

      //! @brief a list of all senders that send signal to this Slots class
      std::set< SignalBase *> m_Senders;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Slots();

      //! @brief copy constructor
      //! duplicates all the senders to this Slots
      //! @param SLOTS the argument slots
      Slots( const Slots &SLOTS);

      //! @brief destructor disconnects all sender
      virtual ~Slots();

    ////////////////
    // operations //
    ////////////////

      //! @brief connect a signal
      //! @brief SENDER the signal source
      void SignalConnect( SignalBase *SENDER);

      //! @brief disconnect a given signal
      //! @brief SENDER the source of the signal that needs to be disconnected
      void SignalDisconnect( SignalBase *SENDER);

      //! @brief disconnects all senders
      //! this is necessary so that this slot does not receives signals if it gets destroyed
      void DisconnectAll();

    }; // class Slots

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_SLOTS_H_
