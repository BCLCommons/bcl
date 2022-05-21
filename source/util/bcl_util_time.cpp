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
#include "util/bcl_util_time.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

// http://gethelp.devx.com/techtips/cpp_pro/10min/10min0500.asp
// http://www.linuxjournal.com/node/5574/print
// http://synesis.com.au/software/unixem.html

#if defined (__GNUC__) && !defined(__MINGW32__)
  #include <sys/time.h>
  #include <unistd.h>
  #define msecond_sleep( MSECONDS) usleep( MSECONDS * 1000)
#elif defined (_MSC_VER)
  #include <sys/timeb.h>
  #include <sys/types.h>
  #include <time.h>
  #include <windows.h>
  #include <winsock.h>
  void gettimeofday( struct timeval *t, void *timezone)
  {       struct _timeb timebuffer;
          _ftime64_s( &timebuffer);
          t->tv_sec = long( timebuffer.time); // explicit cast to long, otherwise C4244 compiler warning in VS
          t->tv_usec = 1000 * timebuffer.millitm;
  };
  //in windows sleep is written with capital Sleep
  #define msecond_sleep( MSECONDS) Sleep( MSECONDS)
#elif defined(__MINGW32__)
  #include <sys/time.h>
  #include <windows.h>
  //in windows sleep is written with capital Sleep
  #define msecond_sleep( MSECONDS) Sleep( MSECONDS)
#endif

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Time::Time() :
      m_Seconds( 0),
      m_Microseconds( 0)
    {
    }

    //! construct Time from Seconds and Microseconds
    Time::Time( const size_t &SECONDS, const size_t &MICROSECONDS) :
      m_Seconds( SECONDS + MICROSECONDS / s_MicroSecondsPerSecond),
      m_Microseconds( MICROSECONDS % s_MicroSecondsPerSecond)
    {
    }

    //! virtual copy constructor
    Time *Time::Clone() const
    {
      return new Time( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Time::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get seconds
    const size_t &Time::GetSeconds() const
    {
      return m_Seconds;
    }

    //! get microseconds
    const size_t &Time::GetMicroSeconds() const
    {
      return m_Microseconds;
    }

    //! set this to current time
    Time &Time::SetToCurrentTime()
    {
      return *this = GetCurrent();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is the time zero
    //! @return true is time is zero
    bool Time::IsZero() const
    {
      return m_Seconds == 0 && m_Microseconds == 0;
    }

    //! @brief set the time to zero time
    void Time::SetZero()
    {
      m_Seconds = 0;
      m_Microseconds = 0;
    }

    //! return the current system time
    Time Time::GetCurrent()
    {
      struct timeval tv;
      gettimeofday( &tv, NULL);
      return Time( tv.tv_sec, tv.tv_usec);
    }

    //! return time in the format like Tue Nov 18 10:54:53 2008
    std::string Time::GetTimeAsDate() const
    {
      const time_t currenttime( m_Seconds);
      return TrimString( std::string( ctime( &currenttime)));
    }

    //! returns the time in the format : hour:minute:second
    std::string Time::GetTimeAsHourMinuteSecond() const
    {
      static const Format s_format( Format().W( 2).R().Fill( '0'));
      std::string time_string;
      time_string = s_format( long( m_Seconds) / long( 3600)) + ":"
                    + s_format( ( long( m_Seconds) % long( 3600)) / long( 60)) + ":" + s_format( long( m_Seconds) % long( 60));

      return time_string;
    }

    //! return time as days hour:minute:second
    std::string Time::GetTimeAsDayHourMinuteSecond() const
    {
      static const Format s_format( Format().W( 2).R().Fill( '0'));
      std::string time_string;
      time_string = s_format(   long( m_Seconds) / long( 3600 * 24)) + "d " +
                    s_format( ( long( m_Seconds) % long( 3600 * 24)) / long( 3600)) + ":" +
                    s_format( ( long( m_Seconds) % long( 3600)) / long( 60)) + ":" +
                    s_format(   long( m_Seconds) % long( 60));

      return time_string;
    }

    //! returns the time in the format hour:minute:second.milliseconds
    std::string Time::GetTimeAsHourMinuteSecondMilliSeconds() const
    {
      const std::string time_string
      (
        Format().W( 2).R().Fill( '0')( int( m_Seconds) / 3600) +
        ":" + Format().W( 2).R().Fill( '0')( ( int( m_Seconds) % 3600) / 60) +
        ":" + Format().W( 2).R().Fill( '0')( int( m_Seconds) % 60) +
        "." + Format().W( 3).R().Fill( '0')( m_Microseconds / 1000)
      );

      return time_string;
    }

    //! @brief get seconds as fractional seconds.fractional
    //! @return double as seconds, and microseconds as fraction of a second
    double Time::GetSecondsFractional() const
    {
      return double( m_Seconds) + double( m_Microseconds) / double( s_MicroSecondsPerSecond);
    }

    //! @brief get the time in total milliseconds
    //! @return sum of seconds an microsenconds converted to milliseconds
    size_t Time::GetTotalMilliseconds() const
    {
      return m_Seconds * s_MilliSecondsPerSecond + m_Microseconds / s_MicroSecondsPerMiliSecond;
    }

    //! @brief get the time in total microseconds
    //! @return sum of seconds converted to microseconds and microseconds
    size_t Time::GetTotalMicroseconds() const
    {
      return m_Seconds * s_MicroSecondsPerSecond + m_Microseconds;
    }

    //! @brief delays the process
    //! @param DELAY the time to delay by
    void Time::Delay( const Time &DELAY)
    {
      msecond_sleep( DELAY.GetTotalMilliseconds());
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief subtract time from this time
    Time &Time::operator -=( const Time &TIME)
    {
      //if TIME is larger then this, this is set to zero and returned
      if( TIME.m_Seconds > m_Seconds || ( TIME.m_Seconds == m_Seconds && TIME.m_Microseconds > m_Microseconds))
      {
        return *this = Time();
      }

      //Seconds are larger but microseonds are smaller
      if( TIME.m_Microseconds > m_Microseconds)
      {
        return *this = Time( m_Seconds - TIME.m_Seconds - 1, s_MicroSecondsPerSecond - ( TIME.m_Microseconds - m_Microseconds));
      }

      //Seconds and Microseconds are larger
      return *this = Time( m_Seconds - TIME.m_Seconds, m_Microseconds - TIME.m_Microseconds);
    }

    //! @brief add time to this time
    Time &Time::operator +=( const Time &TIME)
    {
      // add microseconds
      m_Microseconds += TIME.m_Microseconds;

      // add seconds
      m_Seconds += TIME.m_Seconds;

      // if Microseconds exceed
      if( m_Microseconds >= s_MicroSecondsPerSecond)
      {
        m_Seconds += m_Microseconds / s_MicroSecondsPerSecond;
        m_Microseconds %= s_MicroSecondsPerSecond;
      }

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! Write Time to std::ostream
    std::ostream &Time::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write time as seconds.xxxxxx
      io::Serialize::Write( m_Seconds, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Microseconds, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! Read Time from std::istream
    std::istream &Time::Read( std::istream &ISTREAM)
    {
      //read member
      io::Serialize::Read( m_Seconds, ISTREAM);
      io::Serialize::Read( m_Microseconds, ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Create time from a fractional minutes string, e.g. 35.75 means 35 minutes and 45 seconds
    //! @param TIME the time object to initialize
    //! @param INIT the string to use to initialize the time object
    //! @param ERR_STREAM stream to output errors to
    //! @return true on success
    bool Time::CreateTimeFromMinutesString( Time &TIME, const std::string &INIT, std::ostream &ERR_STREAM)
    {
      double minutes;
      if( !TryConvertFromString( minutes, INIT, ERR_STREAM))
      {
        return false;
      }
      // convert minute to seconds
      TIME.m_Seconds = size_t( minutes * s_SecondsPerMinute);
      TIME.m_Microseconds = size_t( minutes * s_SecondsPerMinute * s_MicroSecondsPerSecond) % s_MicroSecondsPerSecond;
      return true;
    }

    //! @brief Convert the given time object to a minutes string
    //! @return the converted string
    std::string Time::ConvertTimeToMinutesString( const Time &TIME)
    {
      return
        Format()
        (
          ( double( TIME.m_Seconds) + double( TIME.m_Microseconds) / double( s_MicroSecondsPerSecond))
          / double( s_SecondsPerMinute)
        );
    }

    //! @brief Time object from number of hours
    //! @param HOURS the number of hours
    //! @return a Time object with seconds for that many hours
    Time Time::CreateTimeFromHours( const size_t HOURS)
    {
      return Time( HOURS * s_SecondsPerHour, 0);
    }

    //! @brief create time object from compiler __DATE__ and __TIME__ macro
    //! @param DATE format has to be: Mmm dd yyyy as it comes from the __DATE__ macro
    //! @param TIME format has to be: hh::mm::ss as is come from the __TIME__ macro
    Time Time::CreateTimeFromCompilerMacro( const std::string &DATE, const std::string &TIME)
    {
      // split data and time into its components
      const storage::Vector< std::string> splitted_date( SplitString( DATE));
      const storage::Vector< std::string> splitted_time( SplitString( TIME, ":"));

      // array containing all three letter months
      static const std::string s_months[] =
      {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
      };

      // get current timeinfo and modify it to the user's choice
      time_t resulting_time;
      struct tm *timeinfo;
      time( &resulting_time);
      timeinfo = localtime( &resulting_time);
      timeinfo->tm_year = ConvertStringToNumericalValue< size_t>( splitted_date( 2)) - 1900;
      timeinfo->tm_mon = size_t( std::find( s_months, s_months + 12, splitted_date( 0)) - s_months);
      timeinfo->tm_mday = ConvertStringToNumericalValue< size_t>( splitted_date( 1));
      timeinfo->tm_hour = ConvertStringToNumericalValue< size_t>( splitted_time( 0));
      timeinfo->tm_min  = ConvertStringToNumericalValue< size_t>( splitted_time( 1));
      timeinfo->tm_sec  = ConvertStringToNumericalValue< size_t>( splitted_time( 2));
      BCL_Assert( timeinfo->tm_mon < 12, "unknown month supplied: " + splitted_date( 0));

      // call mktime to convert to real time
      resulting_time = mktime( timeinfo);

      // create time object from seconds
      return Time( size_t( resulting_time), 0);
    }

    //! @brief operator for subtracting times.
    //! if time 1 smaller time 2 difference will be zero
    //! @param TIME_LHS time that is usually larger than the time that is subtracted
    //! @param TIME_RHS smaller time that is subtracted
    //! @return difference of times - if left was smaller than right, it will be time zero
    Time operator -( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return Time( TIME_LHS).operator -=( TIME_RHS);
    }

    //! @brief operator for adding times.
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return summation of times
    Time operator +( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return Time( TIME_LHS).operator +=( TIME_RHS);
    }

    //! @brief equal two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is equal to RHS time, false otherwise
    bool operator ==( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds() && TIME_LHS.GetMicroSeconds() == TIME_RHS.GetMicroSeconds();
    }

    //! @brief not equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time not equal than RHS time, false otherwise
    bool operator !=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return TIME_LHS.GetSeconds() != TIME_RHS.GetSeconds() || TIME_LHS.GetMicroSeconds() != TIME_RHS.GetMicroSeconds();
    }

    //! @brief smaller than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller than RHS time, false otherwise
    bool operator <( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() < TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() < TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief smaller equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller or equal than RHS time, false otherwise
    bool operator <=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() < TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() <= TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief larger than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger than RHS time, false otherwise
    bool operator >( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() > TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() > TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief larger equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger or equal than RHS time, false otherwise
    bool operator >=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() > TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() >= TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

  } // namespace util
} // namespace bcl
