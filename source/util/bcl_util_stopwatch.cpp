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
#include "util/bcl_util_stopwatch.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! set the Timer Start to the current time
    //! @param PRINT_ON_DESTRUCTION true if the current time of this stopwatch should be printed when destroyed
    Stopwatch::Stopwatch( const bool PRINT_ON_DESTRUCTION) :
      m_Description(),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( std::numeric_limits< size_t>::max(), 0),
      m_MessageLevel( Message::e_Verbose),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      Start();
    }

    //! @brief constructor with description and notification options
    //! set the Timer Start to the current time
    //! @param DESCRIPTION string that describes what the stopwatch is timing
    //! @param NOTIFICATION_INTERVAL the minimum interval between messages
    //! @param MESSAGE_LEVEL the message level to use when printing messages
    //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
    //! @param AUTO_START whether to start the stopwatch upon construction
    Stopwatch::Stopwatch
    (
      const std::string &DESCRIPTION,
      const Time &NOTIFICATION_INTERVAL,
      const Message::MessageLevel &MESSAGE_LEVEL,
      const bool PRINT_ON_DESTRUCTION,
      const bool AUTO_START
    ) :
      m_Description( DESCRIPTION),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( NOTIFICATION_INTERVAL),
      m_MessageLevel( MESSAGE_LEVEL),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      if( AUTO_START)
      {
        Start();
      }
    }

    //! @brief constructor with description and message level
    //! set the Timer Start to the current time
    //! @param DESCRIPTION string that describes what the stopwatch is timing
    //! @param MESSAGE_LEVEL the message level to use when printing messages
    //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
    Stopwatch::Stopwatch
    (
      const std::string &DESCRIPTION,
      const Message::MessageLevel &MESSAGE_LEVEL,
      const bool PRINT_ON_DESTRUCTION,
      const bool AUTO_START
    ) :
      m_Description( DESCRIPTION),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( std::numeric_limits< size_t>::max(), 0),
      m_MessageLevel( MESSAGE_LEVEL),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      if( AUTO_START)
      {
        Start();
      }
    }

    //! virtual copy constructor
    Stopwatch *Stopwatch::Clone() const
    {
      return new Stopwatch( *this);
    }

    //! virtual destructor - will output the lifetime of this object and though the duration of the scope in which it was initilized
    Stopwatch::~Stopwatch()
    {
      if( m_PrintOnDestruction)
      {
        Stop();
        // only write the message if the stopwatch was actually used.
        if( m_TotalTime.GetSecondsFractional())
        {
          WriteMessage();
        }
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Stopwatch::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return description for that Stopwatch
    //! @return description
    const std::string &Stopwatch::GetDescription() const
    {
      return m_Description;
    }

    //! @brief return the starting time
    //! @return Time object of the current time, when this Stopwatch was constructed
    const Time &Stopwatch::GetLastStartTime() const
    {
      return m_TimerStart;
    }

    //! @brief return the total accumulative Time
    Time Stopwatch::GetTotalTime() const
    {
      // if running,
      return IsRunning() ? m_TotalTime + GetProcessDuration() : m_TotalTime;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is stopwatch running
    //! @return true, is the stopwatch is currently running
    bool Stopwatch::IsRunning() const
    {
      return !m_TimerStart.IsZero();
    }

    //! @brief start this stopwatch again.
    void Stopwatch::Start()
    {
      // check if not already running
      if( IsRunning())
      {
        return;
      }

      m_TimerStart = Time::GetCurrent();
    }

    //! @brief stop this stopwatch and print out a message if enough time has passed
    void Stopwatch::Stop()
    {
      // check if already stopped
      if( !IsRunning())
      {
        return;
      }

      // add difference between timer start and current time to total time
      m_TotalTime += Time::GetCurrent() - m_TimerStart;
      m_TimerStart.SetZero(); // set start time to zero, to indicate, that the stopwatch is not running

      // check for how long no notification was issued
      if( ( m_TotalTime - m_LastNotificationTime) > m_NotificationInterval)
      {
        m_LastNotificationTime = m_TotalTime;
        WriteMessage();
      }
    }

    //! @brief stops the timer, and sets total time to zero
    void Stopwatch::Reset()
    {
      Stop();
      m_TotalTime.SetZero();
      m_LastNotificationTime.SetZero();
    }

    //! @brief return the time since the stop watch was last started
    Time Stopwatch::GetProcessDuration() const
    {
      // calculate the current duration if it is running
      if( IsRunning())
      {
        return Time::GetCurrent() - m_TimerStart;
      }
      // return a zero Time object
      else
      {
        return Time();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write a message that contains the total time that this stopwatch has accumulated
    void Stopwatch::WriteMessage() const
    {
      Time total_time( GetTotalTime());
      BCL_Message
      (
        m_MessageLevel,
        m_Description + " has run for " + util::Format()( total_time.GetSecondsFractional()) + " seconds"
      );
    }

    //! @brief Write StopWatch to std::ostream
    std::ostream &Stopwatch::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Description,          OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TimerStart,           OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TotalTime,            OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LastNotificationTime, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NotificationInterval, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MessageLevel,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PrintOnDestruction,   OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief Read StopWatch fom std::istream
    std::istream &Stopwatch::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Description,          ISTREAM);
      io::Serialize::Read( m_TimerStart,           ISTREAM);
      io::Serialize::Read( m_TotalTime,            ISTREAM);
      io::Serialize::Read( m_LastNotificationTime, ISTREAM);
      io::Serialize::Read( m_NotificationInterval, ISTREAM);
      io::Serialize::Read( m_MessageLevel,         ISTREAM);
      io::Serialize::Read( m_PrintOnDestruction,   ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace util
} // namespace bcl
