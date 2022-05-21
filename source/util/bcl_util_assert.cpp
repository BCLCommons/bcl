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
#include "util/bcl_util_assert.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically
#include <stdlib.h>

namespace bcl
{
  namespace util
  {

    //! an error message is output to the logger and program is terminated
    //! @param ERRORMSG message describing the error
    //! @param ERROR_CODE the error code that should be passed to exit
    //! @param FILE_NAME the filename the assert was raised
    //! @param LINE_NUMBER the line number in the file
    //! @param FUNCTION_NAME the name of the function, the assert was raised
    //! @param DISPLAY_CALLSTACK whether to display the call stack, if possible
    void Assert::Exit
    (
      const std::string &ERRORMSG,
      const int ERROR_CODE,
      const char *FILE_NAME,
      const int LINE_NUMBER,
      const char *FUNCTION_NAME,
      const bool &DISPLAY_CALLSTACK
    )
    {
      // remove path from filename
      std::string file_name( FILE_NAME);
      file_name = file_name.substr( file_name.find_last_of( PATH_SEPARATOR) + 1);

      // if exit code does not indicate an error
      if( !DISPLAY_CALLSTACK)
      {
        if( !ERRORMSG.empty())
        {
          GetLogger() << "===> " << ERRORMSG << std::endl;
        }
      }
      else
      {
        std::cerr << "BCL Error: "   << ERROR_CODE
                  << " | function: " << FUNCTION_NAME
                  << " | file: "     << file_name
                  << " | line: "     << LINE_NUMBER
                  << "\n===> "       << ERRORMSG << std::endl;

        // stack till here
        const CallStack stack( 1);

        // output error
        std::cerr                    << stack.String();

        // if the logger is not the default logger, write the error out to the logger as well
        if( &GetLogger() != &**GetLoggers().e_Default)
        {
          // output error also to logger, in case std::cerr is not being monitored
          GetLogger() << "BCL Error: "   << ERROR_CODE
                      << " | function: " << FUNCTION_NAME
                      << " | file: "     << file_name
                      << " | line: "     << LINE_NUMBER
                      << "\n===> "       << ERRORMSG << std::endl;
        }
      }

      // abort program
      GetRuntimeEnvironment().Finalize( ERROR_CODE);
      exit( ERROR_CODE);
    }

  } // namespace util
} // namespace bcl
