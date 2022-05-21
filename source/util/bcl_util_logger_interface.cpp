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
#include "util/bcl_util_logger_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_loggers.h"

// external includes - sorted alphabetically

#if defined(_WIN32)
  #include <windows.h>
#else
  // posix window size from ioctl
  #include <errno.h>
  #include <sys/ioctl.h>
#endif

namespace bcl
{
  namespace util
  {

    //! @brief get the default maximum line width to be used for a particular logger
    //! this value should be used whenever the output is not to a terminal or the terminal window size is unavailable
    size_t LoggerInterface::GetDefaultMaxLineWidth()
    {
      return 120;
    }

    //! @brief get the line width of the terminal
    //! this value should be used whenever the output is to a terminal
    size_t LoggerInterface::GetTerminalLineWidth()
    {
      if( GetIgnoreTerminalLineWidth())
      {
        return GetDefaultMaxLineWidth();
      }
      size_t actual_cols( 0);
#if defined(_WIN32)
      // Windows
      CONSOLE_SCREEN_BUFFER_INFO buffer_rows_columns;

      if( !GetConsoleScreenBufferInfo( GetStdHandle( STD_OUTPUT_HANDLE), &buffer_rows_columns))
      {
        // probably not running in a terminal
        actual_cols = LoggerInterface::GetDefaultMaxLineWidth();
      }
      else
      {
        // running in a terminal; get the actual width of the window by subtracting the column positions
        // NOTE: when calling this function on an executable run with Wine, be aware that wine does not pass linux console
        //       information to the application.  Likewise, this value will always be 80 unless run on A. a windows box or
        //       B. running it with a proper windows console using wineconsole --backend=user.  If doing B, you may want
        //       to set wineconsole up to not exit at program completion, this setting is in ~/.wine/user.reg under
        //       ExitOnDie, which should be set to 0 to keep the console open
        actual_cols = buffer_rows_columns.srWindow.Right - buffer_rows_columns.srWindow.Left + 1;
      }
#else
      // Linux / Apple
      struct winsize window_rows_columns;

      // cache the old error number and reset it to detect errors in ioctl,
      // since its return value is non-standardized (implementations may return 0 on success, others -1, etc.)
      const int old_errno( errno);
      errno = 0;
      ioctl( 0, TIOCGWINSZ, &window_rows_columns);
      if( !errno)
      {
        // get the column width, if it could be retrieved successfully
        actual_cols = window_rows_columns.ws_col;
      }
      else
      {
        // probably not running in a terminal
        actual_cols = LoggerInterface::GetDefaultMaxLineWidth();
      }
      errno = old_errno;
#endif
      static const size_t s_min_window_size( 60);
      if( actual_cols < s_min_window_size)
      {
        // if the terminal size is really small (say, smaller than 60), than most well formatted output will look bad
        // anyway, so scale up the output
        actual_cols = std::max( actual_cols, size_t( 1));

        // try to still keep actual columns as an integer multiple of its real value
        actual_cols *= ( ( s_min_window_size - 1) / actual_cols + 1);
      }

      return actual_cols;
    }

    //! @brief set a flag to ignore the terminal's line width when reporting the default max line width
    //! @param IGNORE_TERMINAL_LINE_WIDTH true to ignore the terminal's line width when reporting GetTerminalLineWidth
    void LoggerInterface::SetIgnoreTerminalLineWidth( const bool &IGNORE_TERMINAL_LINE_WIDTH)
    {
      GetIgnoreTerminalLineWidth() = IGNORE_TERMINAL_LINE_WIDTH;
    }

    //! @brief get currently used Logger
    LoggerInterface &GetLogger()
    {
      return Loggers::GetEnums().GetCurrentLogger();
    }

  } // namespace util
} // namespace bcl
