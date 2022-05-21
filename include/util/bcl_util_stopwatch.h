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

#ifndef BCL_UTIL_STOPWATCH_H
#define BCL_UTIL_STOPWATCH_H

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_message.h"
#include "bcl_util_time.h"
#include "bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Stopwatch
    //! @brief tracks the total amount of time spent in functions that are called repeatedly
    //! optionally notifying the user whenever some amount of time since the last output has passed
    //! @details
    //! void foo()
    //! {
    //!   static Stopwatch s_FooTimer( "total time spent in foo()");
    //!   s_FooTimer.Start();
    //!   ... do something that takes a while in here ...
    //!   s_FooTimer.Stop();
    //! };
    //!
    //! @see @link example_util_stopwatch.cpp @endlink
    //! @author woetzen
    //! @date 06/02/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Stopwatch :
      public ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      std::string m_Description;            //!< a string containing information on what is timed
      Time        m_TimerStart;             //!< time when the object was last started - is set to zero, while not running
      Time        m_TotalTime;              //!< total runtime

      Time        m_LastNotificationTime;   //!< the last time when notification was issued
      Time        m_NotificationInterval;   //!< time that needs to elapse before the next notification is issued

      Message::MessageLevelEnum m_MessageLevel; //!< message level used when printing messages

      bool        m_PrintOnDestruction;     //!< true if the stopwatch should print its usual message when destroyed

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! set the Timer Start to the current time
      //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
      Stopwatch( const bool PRINT_ON_DESTRUCTION = true);

      //! @brief constructor with description and notification options
      //! set the Timer Start to the current time
      //! @param DESCRIPTION string that describes what the stopwatch is timing
      //! @param NOTIFICATION_INTERVAL the minimum interval between messages
      //! @param MESSAGE_LEVEL the message level to use when printing messages
      //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
      //! @param AUTO_START whether to start the stopwatch upon construction
      explicit Stopwatch
      (
        const std::string &DESCRIPTION,
        const Time &NOTIFICATION_INTERVAL = Time( std::numeric_limits< size_t>::max(), 0),
        const Message::MessageLevel &MESSAGE_LEVEL = Message::e_Verbose,
        const bool PRINT_ON_DESTRUCTION = true,
        const bool AUTO_START = true
      );

      //! @brief constructor with description and message level
      //! set the Timer Start to the current time
      //! @param DESCRIPTION string that describes what the stopwatch is timing
      //! @param MESSAGE_LEVEL the message level to use when printing messages
      //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
      explicit Stopwatch
      (
        const std::string &DESCRIPTION,
        const Message::MessageLevel &MESSAGE_LEVEL,
        const bool PRINT_ON_DESTRUCTION,
        const bool AUTO_START = true
      );

      //! virtual copy constructor
      Stopwatch *Clone() const;

      //! @brief destructor
      //! outputs the duration of the scope in which it was initialized if m_PrintWhenDestroyed is true
      virtual ~Stopwatch();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return description for that Stopwatch
      //! @return description
      const std::string &GetDescription() const;

      //! @brief return the starting time
      //! @return Time object of the current time, when this Stopwatch was constructed
      const Time &GetLastStartTime() const;

      //! @brief return the total accumulative Time
      Time GetTotalTime() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief is stopwatch running
      //! @return true, is the stopwatch is currently running
      bool IsRunning() const;

      //! @brief start this stopwatch again.
      void Start();

      //! @brief stop this stopwatch and print out a message, if enough time has passed
      void Stop();

      //! @brief stops the timer, and sets total time to zero
      void Reset();

      //! @brief return the time since the stopwatch was last started
      Time GetProcessDuration() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write a message that contains the total time that this stopwatch has accumulated
      void WriteMessage() const;

    protected:

      //! Write StopWatch to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! Read StopWatch fom std::istream
      std::istream &Read( std::istream &ISTREAM);

    }; // class Stopwatch

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_STOPWATCH_H
