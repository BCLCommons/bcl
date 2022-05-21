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

#ifndef BCL_UTIL_TIME_H_
#define BCL_UTIL_TIME_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Time
    //! @brief class is used to store time in seconds plus microseconds
    //! @details It provides functions to get the time, the current systemtime, read and write operation and operator -
    //! also look at:  http://www.informit.com/guides/content.asp?g=cplusplus&seqNum=64
    //!
    //! @see @link example_util_time.cpp @endlink
    //! @author woetzen
    //! @date 06/02/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Time :
      public ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t m_Seconds;      //!< Seconds
      size_t m_Microseconds; //!< Microseconds

    public:

      static const size_t s_MicroSecondsPerSecond     = 1000000; //!< microseconds per second (10^6)
      static const size_t s_MilliSecondsPerSecond     = 1000;    //!< milliseconds per second (10^3)
      static const size_t s_MicroSecondsPerMiliSecond = 1000;    //!< microseconds per millisecond( 10^3)
      static const size_t s_SecondsPerMinute          = 60;      //!< seconds per minute
      static const size_t s_SecondsPerHour            = 3600;    //!< seconds per hour

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Time();

      //! construct Time from Seconds and Microseconds
      Time( const size_t &SECONDS, const size_t &MICROSECONDS);

      //! copy constructor
      Time *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! get seconds
      const size_t &GetSeconds() const;

      //! get microseconds
      const size_t &GetMicroSeconds() const;

      //! set this to current time
      Time &SetToCurrentTime();

    ////////////////
    // operations //
    ////////////////

      //! @brief is the time zero
      //! @return true is time is zero
      bool IsZero() const;

      //! @brief set the time to zero time
      void SetZero();

      //! return the current system time
      static Time GetCurrent();

      //! return time in the format like Tue Nov 18 10:54:53 2008
      std::string GetTimeAsDate() const;

      //! returns the time in the format : hour:minute:second
      std::string GetTimeAsHourMinuteSecond() const;

      //! return time as days hour:minute:second
      std::string GetTimeAsDayHourMinuteSecond() const;

      //! @brief get time in a readable format with Milliseconds
      //! returns the time in the format : hour:minute:second.milliseconds
      std::string GetTimeAsHourMinuteSecondMilliSeconds() const;

      //! @brief get seconds as fractional seconds.fractional
      //! @return double as seconds, and microseconds as fraction of a second
      double GetSecondsFractional() const;

      //! @brief get the time in total milliseconds
      //! @return sum of seconds an microsenconds converted to milliseconds
      size_t GetTotalMilliseconds() const;

      //! @brief get the time in total microseconds
      //! @return sum of seconds converted to microseconds and microseconds
      size_t GetTotalMicroseconds() const;

      //! @brief delays the process
      //! @param DELAY the time to delay by (maximal milli second precision)
      static void Delay( const Time &DELAY);

    ///////////////
    // operators //
    ///////////////

      //! @brief subtract time from this time
      Time &operator -=( const Time &TIME);

      //! @brief add time to this time
      Time &operator +=( const Time &TIME);

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! Write Time to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! Read Time from std::istream
      std::istream &Read( std::istream &ISTREAM);

    //////////////////////
    // helper functions //
    //////////////////////
    public:

      //! @brief Create time from a fractional minutes string, e.g. 35.75 means 35 minutes and 45 seconds
      //! @param TIME the time object to initialize
      //! @param INIT the string to use to initialize the time object
      //! @param ERR_STREAM stream to output errors to
      //! @return true on success
      static bool CreateTimeFromMinutesString( Time &TIME, const std::string &INIT, std::ostream &ERR_STREAM);

      //! @brief Convert the given time object to a minutes string
      //! @return the converted string
      static std::string ConvertTimeToMinutesString( const Time &TIME);

      //! @brief Time object from number of hours
      //! @param HOURS the number of hours
      //! @return a Time object with seconds for that many hours
      static Time CreateTimeFromHours( const size_t HOURS);

      //! @brief create time object from compiler DATE and TIME macro
      //! @param DATE format has to be: Mmm dd yyyy as it comes from the DATE macro
      //! @param TIME format has to be: hh:mm:ss as is come from the TIME macro
      static Time CreateTimeFromCompilerMacro( const std::string &DATE, const std::string &TIME);

    }; //class Time

    //! @brief operator for subtracting times.
    //! if time 1 smaller time 2 difference will be zero
    //! @param TIME_LHS time that is usually larger than the time that is subtracted
    //! @param TIME_RHS smaller time that is subtracted
    //! @return difference of times - if left was smaller than right, it will be time zero
    BCL_API Time operator -( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief operator for adding times.
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return summation of times
    BCL_API Time operator +( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief equal two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is equal to RHS time, false otherwise
    BCL_API bool operator ==( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief not equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time not equal than RHS time, false otherwise
    BCL_API bool operator !=( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief smaller than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller than RHS time, false otherwise
    BCL_API bool operator <( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief smaller equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller or equal than RHS time, false otherwise
    BCL_API bool operator <=( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief larger than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger than RHS time, false otherwise
    BCL_API bool operator >( const Time &TIME_LHS, const Time &TIME_RHS);

    //! @brief larger equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger or equal than RHS time, false otherwise
    BCL_API bool operator >=( const Time &TIME_LHS, const Time &TIME_RHS);

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_TIME_H_
