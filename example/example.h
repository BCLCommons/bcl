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

#ifndef EXAMPLE_H_
#define EXAMPLE_H_

// include forward header of this class

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "io/bcl_io.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "example_check_macros.h"
#include "example_interface.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_table.hpp"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @class ExampleClass
  //! @brief This class is the basic Example class
  //! For adding an example to the bcl do:\n
  //! 1) derive your example from ExampleInterface
  //! 2) overwrite all abstract functions (Clone() - virtual copy constructor and Run() containing your example code)
  //! 3) make a static member in your example: static const ExampleClass::EnumType Example{your Example class name}_Instance
  //! 4) initialize this instance by calling the constructor outside the class:<br>
  //!    const ExampleClass::EnumType Example{your Example class name}::Example{your Example class name}_Instance( GetExamples().AddEnum( "{your Example class name}", Example{your Example class name}());
  //!
  //! All examples have to check their results by itself \n
  //! - give a warning for noncritical result changes \n
  //! - Assert if you know the result for certain
  //! - all output has to be a message - no std::cout please
  //! - use util::Format in the util::Message for objects
  //!
  //! @author meilerj, woetzen
  //!
  //! @date 21.08.2005
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClass :
    public util::Enumerate< util::ShPtr< ExampleInterface>, ExampleClass>
  {
      friend class util::Enumerate< util::ShPtr< ExampleInterface>, ExampleClass>;
  public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExampleResult
    //! @brief collect all results for all examples that define checks in that data structure
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class ExampleResult :
      public storage::Table< std::string>
    {

    private:

    //////////
    // data //
    //////////

      //! @brief the header for the result table
      //! @return the header for the result table
      static const storage::Vector< std::string> &GetHeader();

      //! @brief count number of successful checks
      size_t m_NumberSuccesses;

      //! @brief count number of errors in checks
      size_t m_NumberErrors;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExampleResult();

    /////////////////
    // data access //
    /////////////////

      //! @brief return the total number of errors
      //! @return number of errors
      size_t GetNumberErrors() const
      {
        return m_NumberErrors;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief insert an examples test case result into the table
      //! @param EXAMPLE_NAME name of the example the test case is generated in
      //! @param LINE_NUMBER line number where test is generated
      //! @param SUCCESS true for a successful test
      //! @param ERROR_STRING if test is unsuccessful this will be the error message
      void InsertTest
      (
        const std::string &EXAMPLE_NAME,
        const size_t LINE_NUMBER,
        const bool SUCCESS,
        const std::string &ERROR_STRING
      );

      //! insert default row for total result of example
      void InsertDefaultTotalRow( const std::string &EXAMPLE_NAME);

      //! @brief write all results
      //! @param OSTREAM stream to be written to
      std::ostream &WriteResults( std::ostream &OSTREAM);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief generate row name for example
      static const std::string GenerateTotalRowName( const std::string &EXAMPLE_NAME);

    }; // class ExampleResult

  //////////
  // data //
  //////////

    storage::Vector< std::string> m_Namespaces;

  private:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExampleClass();

  public:

  //////////
  // data //
  //////////

    //! flag to change example path
    static util::ShPtr< command::FlagInterface> &GetExamplePathFlag();

    //! collection of all example results
    static ExampleResult &GetResults();

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const;

    //! @brief get example path
    //! @return path where example files are written and read from
    const std::string &GetExamplePath() const;

    //! @brief add a new example to the examples
    //! @param EXAMPLE_I ExampleInterface derived instance of that example
    //! @return reference to the inserted EnumData object
    EnumType &AddEnum( const ExampleInterface &EXAMPLE_I);

    //! @brief return the namespaces examples are available for
    //! @return vector of strings of namespaces
    const storage::Vector< std::string> &GetNamespaces() const
    {
      return m_Namespaces;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief remove "bcl::Example" from the static class name of the example
    //! @param EXAMPLE_CLASS_NAME the class name returned form GetClassIdentifier
    //! @return a shorter name - removed "bcl::Example" from the front
    static std::string RemoveUnnecessaryString( const std::string &EXAMPLE_CLASS_NAME);

    //! @brief get first capital word from a string of the form "WordWord" -> "Word"
    //! @param WORD1_WORD2 string of the Form "WordoneWordtwo"
    //! @return string of the form "WordOne"
    static std::string FirstWord( const std::string &WORD1_WORD2);

    //! @brief concatenate STRINGS onto separate lines, squeeze multiple entries per line until line-width would exceed threshold
    //! @param STRINGS a vector of strings.  These may contain empty strings or strings with new-lines in them
    //! @param THRESHOLD entries will be added to each line so long as the line's width stays under this value
    //! @return the resulting string
    static std::string AddNewLinesWhenLineWidthIsExcessive
    (
      const storage::Vector< std::string> &STRINGS,
      const size_t &THRESHOLD
    );

    //! @brief open a file, report if the file is missing
    //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
    //! @param LINE_NUMBER line number of the example where the example check was called
    //! @param FILENAME file the example check was called on
    //! @param ISTREAM stream that is to be opened
    //! @param OPEN_MODE the manner in which FILENAME should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return the resulting string
    static bool ExampleMustOpenInputFile
    (
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const std::string &FILENAME,
      io::IFStream &ISTREAM,
      const std::ios_base::openmode OPEN_MODE = std::ios::in
    );

    //! @brief open an output file, report if the file could not be opened for writing
    //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
    //! @param LINE_NUMBER line number of the example where the example check was called
    //! @param FILENAME file the example check was called on
    //! @param OSTREAM stream that is to be opened
    //! @param OPEN_MODE the manner in which FILENAME should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return the resulting string
    static bool ExampleMustOpenOutputFile
    (
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const std::string &FILENAME,
      io::OFStream &OSTREAM,
      const std::ios_base::openmode OPEN_MODE = std::ios::out
    );

    //! @brief ExampleCheckReport a function called by macros to report example checks, optionally on success
    //! Example checks are checks that ACTUAL_RESULT == EXPECTED_RESULT
    //! This function is used with macros to automatically inform the user at the point of failure, what was executed, and
    //! what were the results, and that they were not equal
    //! @param EXAMPLE_CHECK_NAME description of what was checked; used only if the example checks something indirectly
    //! for example, checking that clone worked right based on whether the object is defined
    //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
    //! @param LINE_NUMBER passed via preprocessor __LINE__
    //! @param ACTUAL_CODE The code that was executed to get the result on the left hand side of the ==
    //! @param ACTUAL_RESULT The object that resulted from executing ACTUAL_CODE
    //! @param EXPECTED_CODE The code that was executed to get the result on the right hand side of the ==
    //! @param EXPECTED_RESULT The object that resulted from executing EXPECTED_CODE
    //! @return true iff ACTUAL_RESULT == EXPECTED_RESULT
    template< typename t_DataType1, typename t_DataType2>
    static bool ExampleCheckReport
    (
      const std::string &EXAMPLE_CHECK_NAME,
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const char *const &ACTUAL_CODE,
      t_DataType1 const &ACTUAL_RESULT,
      const char *const &EXPECTED_CODE,
      t_DataType2 const &EXPECTED_RESULT
    )
    {
      return HandleResult
             (
               EXAMPLE_CHECK_NAME,
               EXAMPLE_NAME,
               LINE_NUMBER,
               ACTUAL_CODE,
               util::Format()( ACTUAL_RESULT),
               EXPECTED_CODE,
               util::Format()( EXPECTED_RESULT),
               ACTUAL_RESULT == EXPECTED_RESULT
             );
    } // ExampleCheckReport

    //! @brief ExampleCheckReport a function called by macros to report example checks, optionally on success
    //! Example checks are checks that ACTUAL_RESULT == EXPECTED_RESULT
    //! This function is used with macros to automatically inform the user at the point of failure, what was executed, and
    //! what were the results, and that they were not equal
    //! @param EXAMPLE_CHECK_NAME description of what was checked; used only if the example checks something indirectly
    //! for example, checking that clone worked right based on whether the object is defined
    //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
    //! @param LINE_NUMBER passed via preprocessor __LINE__
    //! @param ACTUAL_CODE The code that was executed to get the result on the left hand side of the ==
    //! @param ACTUAL_RESULT The object that resulted from executing ACTUAL_CODE
    //! @param EXPECTED_CODE The code that was executed to get the result on the right hand side of the ==
    //! @param EXPECTED_RESULT The object that resulted from executing EXPECTED_CODE
    //! @param RELATIVE_TOLERANCE relative maximum that the results may differ by
    //! @param ABSOLUTE_TOLERANCE absolute maximum that the results may differ by
    //! @return true iff ACTUAL_RESULT == EXPECTED_RESULT
    template< typename t_DataType1, typename t_DataType2>
    static bool
    ExampleCheckReportWithTolerance
    (
      const std::string &EXAMPLE_CHECK_NAME,
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const char *const &ACTUAL_CODE,
      t_DataType1 const &ACTUAL_RESULT,
      const char *const &EXPECTED_CODE,
      t_DataType2 const &EXPECTED_RESULT,
      const double &ABSOLUTE_TOLERANCE,
      const double &RELATIVE_TOLERANCE
    )
    {
      return HandleResult
             (
               EXAMPLE_CHECK_NAME,
               EXAMPLE_NAME,
               LINE_NUMBER,
               ACTUAL_CODE,
               util::Format()( ACTUAL_RESULT),
               EXPECTED_CODE,
               util::Format()( EXPECTED_RESULT),
               math::EqualWithinTolerance( ACTUAL_RESULT, EXPECTED_RESULT, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE),
               RELATIVE_TOLERANCE,
               ABSOLUTE_TOLERANCE
             );
    } // ExampleCheckReport

    //! @brief overload for comparing size_t's to integers, which the compiler will otherwise complain about
    //! @param LINE_NUMBER, EXAMPLE_CHECK_NAME, EXAMPLE_NAME, ACTUAL_CODE, ACTUAL_RESULT, EXPECTED_CODE, EXPECTED_RESULT see above
    //! @return true iff ACTUAL_RESULT == EXPECTED_RESULT
    template< typename t_Integer>
    static bool ExampleCheckReport
    (
      const std::string &EXAMPLE_CHECK_NAME,
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const char *const &ACTUAL_CODE,
      const size_t &ACTUAL_RESULT,
      const char *const &EXPECTED_CODE,
      const t_Integer &EXPECTED_RESULT
    )
    {
      return ExampleCheckReport< size_t, size_t>
             (
               EXAMPLE_CHECK_NAME,
               EXAMPLE_NAME,
               LINE_NUMBER,
               ACTUAL_CODE,
               ACTUAL_RESULT,
               EXPECTED_CODE,
               size_t( EXPECTED_RESULT)
             );
    } // ExampleCheckReport

    //! @brief overload for comparing size_t's to size_t's, which catches comparisons of size_t's to anything else
    //! @param EXAMPLE_CHECK_NAME, EXAMPLE_NAME, LINE_NUMBER see above
    //! @param ACTUAL_CODE, ACTUAL_RESULT, EXPECTED_CODE, EXPECTED_RESULT see above
    //! @return true iff ACTUAL_RESULT == EXPECTED_RESULT
    static bool ExampleCheckReport
    (
      const std::string &EXAMPLE_CHECK_NAME,
      const std::string &EXAMPLE_NAME,
      const size_t &LINE_NUMBER,
      const char *const &ACTUAL_CODE,
      const size_t &ACTUAL_RESULT,
      const char *const &EXPECTED_CODE,
      const size_t &EXPECTED_RESULT
    );

    //! @brief Try to register a file for this example to read or write
    //! @param EXAMPLE_NAME the example that is trying to access the file
    //! @param FILENAME the file that the example is attempting to access
    //! @param OPEN_MODE the manner in which FILENAME should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @param ERR_STREAM the stream to write error messages out to
    //! @return true if the file is openable by this example
    static bool RequestExampleFileAccess
    (
      const std::string &EXAMPLE_NAME,
      const std::string &FILENAME,
      const std::ios_base::openmode OPEN_MODE,
      std::ostream &ERR_STREAM
    );

  private:

    //! @brief replace line breaks with given string
    //! @param STRING the string with line breaks
    //! @param NEW_LINE_PREFIX the string that is prepends to every new line
    //! @return a string with NEW_LINE_PREFIX inserted after each '\n' character
    static std::string FollowNewlinesWith( const std::string &STRING, const std::string &NEW_LINE_PREFIX);

    //! @brief HandleReportedResult a helper function for Example::ExampleCheckReport
    //! @param EXAMPLE_CHECK_NAME description of what was checked; used only if the example checks something indirectly
    //! for example, checking that clone worked right based on whether the object is defined
    //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
    //! @param LINE_NUMBER passed via preprocessor __LINE__
    //! @param ACTUAL_CODE The code that was executed to get the result on the left hand side of the ==
    //! @param ACTUAL_RESULT serialized object that resulted from executing ACTUAL_CODE
    //! @param EXPECTED_CODE The code that was executed to get the result on the right hand side of the ==
    //! @param EXPECTED_RESULT serialized object that resulted from executing EXPECTED_CODE
    //! @param SUCCESS did the actual result and expected result objects compare equal
    //! @param RELATIVE_TOLERANCE relative maximum that the results may differ by
    //! @param ABSOLUTE_TOLERANCE absolute maximum that the results may differ by (for vectors, the max that any single value may differ by)
    static bool HandleResult
    (
      const std::string &EXAMPLE_CHECK_NAME,
      const std::string &EXAMPLE_NAME,
      const size_t      &LINE_NUMBER,
      const char *const &ACTUAL_CODE,
      const std::string &ACTUAL_RESULT,
      const char *const &EXPECTED_CODE,
      const std::string &EXPECTED_RESULT,
      const bool        &SUCCESS,
      const double      &RELATIVE_TOLERANCE = 0.0,
      const double      &ABSOLUTE_TOLERANCE = 0.0
    );

    //! @brief get the map from file name to example that wrote to it
    //! @return the map from file name to example that wrote to it
    //! @detail to prevent race conditions when examples are run in parallel, examples may not write to the same
    //! @detail output file that another example has read or written to
    //! @detail the base path for all entries is the main bcl path
    static storage::Map< std::string, std::string> &GetExampleOutputFileMap();

    //! @brief get the map from file name to an example that read it
    //! @return the map from file name to an example that read it
    //! @detail to prevent race conditions, an example cannot read a file that a different example writes
    //! @detail the base path for all entries is the main bcl path
    static storage::Map< std::string, std::string> &GetExampleInputFileMap();

  }; // class ExampleClass

  ExampleClass &GetExamples();

} // namespace bcl

#endif // EXAMPLE_H_
