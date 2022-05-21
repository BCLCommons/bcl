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

#ifndef BCL_UTIL_LOGGERS_H_
#define BCL_UTIL_LOGGERS_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_enumerate.h"
#include "bcl_util_logger_interface.h"
#include "bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief for Logger
    typedef Enum< ShPtr< LoggerInterface>, Loggers> Logger;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Loggers
    //! @brief manages the system-wide list of enumerated loggers
    //!
    //! @see @link example_util_loggers.cpp @endlink
    //! @author heinzes1
    //! @date Jan 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Loggers :
      public Enumerate< ShPtr< LoggerInterface>, Loggers>
    {
      friend class Enumerate< ShPtr< LoggerInterface>, Loggers>;

    public:

    //////////
    // data //
    //////////

      //! default logger writes to std::cout or std::err
      Logger e_Default;

      //! the current logger
      ShPtr< LoggerInterface> m_CurrentLogger;

      //! the enum for the current logger
      Logger m_CurrentLoggerEnum;

    //////////
    // data //
    //////////

      //! @brief command line flag to be used to set Logger over the command line
      //! @return ShPtr to a FlagInterface which is used to set Logger
      static ShPtr< command::FlagInterface> &GetFlagLogger();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      Loggers();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Initialize the logger from the command line flag
      static void UpdateCurrentLoggerFromCommandLineFlag();

      //! @brief Initialize and return the current logger
      //! @return one and only reference to one of the loggers
      LoggerInterface &GetCurrentLogger();

    }; // class Loggers

    //! @brief get enumerated list of Loggers
    BCL_API
    Loggers &GetLoggers();

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< LoggerInterface>, Loggers>;

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_LOGGERS_H_
