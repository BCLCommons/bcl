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

#ifndef BCL_UTIL_ASSERT_H_
#define BCL_UTIL_ASSERT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Assert
    //! @brief Class used by BCL_Assert and BCL_Exit to write error-related information before exiting
    //!
    //! @see @link example_util_assert.cpp @endlink
    //! @author woetzen
    //! @date Nov 6, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Assert
    {

    public:

      //! an error message is output to the logger and program is terminated
      //! @param ERRORMSG message describing the error
      //! @param ERROR_CODE the error code that should be passed to exit
      //! @param FILE_NAME the filename the assert was raised
      //! @param LINE_NUMBER the line number in the file
      //! @param FUNCTION_NAME the name of the function, the assert was raised
      //! @param DISPLAY_CALLSTACK whether to display the call stack, if possible
      static void Exit
      (
        const std::string &ERRORMSG,
        const int ERROR_CODE,
        const char *FILE_NAME,
        const int LINE_NUMBER,
        const char *FUNCTION_NAME,
        const bool &DISPLAY_CALLSTACK
      );

    }; // class Assert

  } // namespace util

  // Use this version when a violation of the assert condition is not expected and indicates a bug in the code
  // this version prints off a callstack
  #define BCL_Assert( BCL_ASSERT_CONDITION, BCL_ASSERT_DESCRIPTION)  { if( !( BCL_ASSERT_CONDITION)){ util::Assert::Exit( BCL_ASSERT_DESCRIPTION, -1, __FILE__, __LINE__, __FUNCTION__, true);}}

  // Use this version when a violation of the assert condition indicates that the user gave us bad input
  // this version does not try to print off a callstack
  #define BCL_UserAssert( BCL_ASSERT_CONDITION, BCL_ASSERT_DESCRIPTION)  { if( !( BCL_ASSERT_CONDITION)){ util::Assert::Exit( BCL_ASSERT_DESCRIPTION, -1, __FILE__, __LINE__, __FUNCTION__, false);}}

  #define BCL_Exit( BCL_EXIT_DESCRIPTION, EXIT_CODE)                 { util::Assert::Exit( BCL_EXIT_DESCRIPTION, EXIT_CODE, __FILE__, __LINE__, __FUNCTION__, EXIT_CODE);}

  // When exiting, either without error or because of some error that is clearly on the user's part, the call stack
  // should not be printed to avoid confusion that the fault lies in the code
  #define BCL_ExitWithoutCallstack( BCL_EXIT_DESCRIPTION, EXIT_CODE) { util::Assert::Exit( BCL_EXIT_DESCRIPTION, EXIT_CODE, __FILE__, __LINE__, __FUNCTION__, false);}
} // namespace bcl

#endif //BCL_UTIL_ASSERT_H_
