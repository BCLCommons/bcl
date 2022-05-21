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
#include "example.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

//! path definition for example files
#if defined (__GNUC__)
  #define BCL_EXAMPLE_PATH_DEFAULT "example/example_files/"
#elif defined (_MSC_VER)
  #define BCL_EXAMPLE_PATH_DEFAULT "../../../example/example_files/"
#endif

namespace bcl
{
//////////
// data //
//////////

  //! @brief the header for the result table
  //! @return the header for the result table
  const storage::Vector< std::string> &ExampleClass::ExampleResult::GetHeader()
  {
    static const storage::Vector< std::string> s_header
    (
      storage::Vector< std::string>::Create( "success", "error", "error_string")
    );
    return s_header;
  }

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

  //! @brief default constructor
  ExampleClass::ExampleClass() :
    util::Enumerate< util::ShPtr< ExampleInterface>, ExampleClass>( false)
  {
  }

  //! @brief default constructor
  ExampleClass::ExampleResult::ExampleResult() :
    storage::Table< std::string>( storage::TableHeader( GetHeader())),
    m_NumberSuccesses( 0),
    m_NumberErrors( 0)
  {
  }

////////////////
// operations //
////////////////

  //! @brief insert a examples testcase result into the table
  //! @param EXAMPLE_NAME name of the example the testcase is generated in
  //! @param LINE_NUMBER line number where test is generated
  //! @param SUCCESS true for a successful test
  //! @param ERROR_STRING if test is unsuccessful this will be the error message
  void ExampleClass::ExampleResult::InsertTest
  (
    const std::string &EXAMPLE_NAME,
    const size_t LINE_NUMBER,
    const bool SUCCESS,
    const std::string &ERROR_STRING
  )
  {
    // count total up
    const std::string total_name( RemoveUnnecessaryString( EXAMPLE_NAME) + "_total");
    std::string &field( operator[]( total_name)( size_t( !SUCCESS)));
    field = util::Format()( util::ConvertStringToNumericalValue< size_t>( field) + 1);

    if( SUCCESS)
    {
      ++m_NumberSuccesses;
    }
    else
    {
      // add particular test case
      const std::string name( RemoveUnnecessaryString( EXAMPLE_NAME) + '_' + util::Format().W( 5).Fill( '0').R()( LINE_NUMBER));

      ++m_NumberErrors;
      InsertRow( name, storage::Vector< std::string>::Create( "_", "X", ERROR_STRING), true);
    }
  }

  //! @brief generate row name for example
  const std::string ExampleClass::ExampleResult::GenerateTotalRowName( const std::string &EXAMPLE_NAME)
  {
    return EXAMPLE_NAME + "_total";
  }

  //! insert default row for total result of example
  void ExampleClass::ExampleResult::InsertDefaultTotalRow( const std::string &EXAMPLE_NAME)
  {
    InsertRow( GenerateTotalRowName( EXAMPLE_NAME), storage::Vector< std::string>::Create( "0", "0", ""));
  }

  //! @brief write all results
  //! @param OSTREAM stream to be written to
  std::ostream &ExampleClass::ExampleResult::WriteResults( std::ostream &OSTREAM)
  {
    InsertRow
    (
      "total",
      storage::Vector< std::string>::Create( util::Format()( m_NumberSuccesses), util::Format()( m_NumberErrors), "")
    );

    // end
    return WriteFormatted( OSTREAM);
  }

/////////////////
// data access //
/////////////////

  //! @brief returns class name
  //! @return the class name as const ref std::string
  const std::string &ExampleClass::GetClassIdentifier() const
  {
    return GetStaticClassName( *this);
  }

  //! @brief add a new example to the examples
  //! @param EXAMPLE_I ExampleInterface derived instance of that example
  //! @return reference to the inserted EnumData object
  ExampleClass::EnumType &ExampleClass::AddEnum
  (
    const ExampleInterface &EXAMPLE_I
  )
  {
    const std::string name( RemoveUnnecessaryString( EXAMPLE_I.GetClassIdentifier()));
    const std::string namespace_name( FirstWord( name));

    if( std::find( m_Namespaces.Begin(), m_Namespaces.End(), namespace_name) == m_Namespaces.End())
    {
      m_Namespaces.PushBack( namespace_name);
    }

    return util::Enumerate< util::ShPtr< ExampleInterface>, ExampleClass>::AddEnum
    (
      name, util::ShPtr< ExampleInterface>( EXAMPLE_I.Clone())
    );
  }

  //! collection of all example results
  ExampleClass::ExampleResult &ExampleClass::GetResults()
  {
    static ExampleClass::ExampleResult s_example_result;
    return s_example_result;
  }

  //! flag to change example path
  util::ShPtr< command::FlagInterface> &ExampleClass::GetExamplePathFlag()
  {
    // static example path flag
    static util::ShPtr< command::FlagInterface> s_example_path_flag
    (
      new command::FlagStatic
      (
        "example_path",
        "change path for reading and writing example files",
        command::Parameter
        (
          "path",
          "relative or absolute example path",
          command::ParameterCheckFileInSearchPath( "example/example_files", "./", io::Directory::e_Dir),
          ""
        )
      )
    );

    // end
    return s_example_path_flag;
  }

  //! @brief get example path
  //! @return path where example files are written and read from
  const std::string &ExampleClass::GetExamplePath() const
  {
    return GetExamplePathFlag()->GetFirstParameter()->GetValue();
  }

//////////////////////
// helper functions //
//////////////////////

  //! @brief remove "bcl::Example" from the static class name of the example
  //! @param EXAMPLE_CLASS_NAME the class name returned form GetClassIdentifier
  //! @return a shorter name - removed "bcl::Example" from the front
  std::string ExampleClass::RemoveUnnecessaryString( const std::string &EXAMPLE_CLASS_NAME)
  {
    // the string to be shortened by
    static const std::string s_unnecessary_string( "bcl::Example");

    // Return the part of the name after the unneccessary part
    return EXAMPLE_CLASS_NAME.substr( s_unnecessary_string.length());
  }

  //! @brief get first capital word from a string of the form "WordWord" -> "Word"
  //! @param WORD1_WORD2 string of the Form "WordoneWordtwo"
  //! @return string of the form "WordOne"
  std::string ExampleClass::FirstWord( const std::string &WORD1_WORD2)
  {
    // if word is empty
    if( WORD1_WORD2.empty())
    {
      return WORD1_WORD2;
    }

    // search for the first capital letter after the first letter
    for( std::string::const_iterator itr( WORD1_WORD2.begin() + 1), itr_end( WORD1_WORD2.end()); itr != itr_end; ++itr)
    {
      if( std::isupper( *itr))
      {
        // return the first word before the capital letter
        return std::string( WORD1_WORD2.begin(), itr);
      }
    }

    // end the WORD since no first word was found
    return WORD1_WORD2;
  }

  //! @brief concatenate STRINGS onto separate lines, squeeze multiple entries per line until line-width would exceed threshold
  //! @param STRINGS a vector of strings.  These may contain empty strings or strings with new-lines in them
  //! @param THRESHOLD entries will be added to each line so long as the line's width stays under this value
  //! @return the resulting string
  std::string ExampleClass::AddNewLinesWhenLineWidthIsExcessive
  (
    const storage::Vector< std::string> &STRINGS,
    const size_t &THRESHOLD
  )
  {
    std::string result; // string to hold the result of concatenation
    size_t width( 0); // number of chars since the last new-line

    // this variable will track when the last string had several lines but did not end in a new-line
    bool last_str_was_continued_multiline( false);

    for
    (
      storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
      itr != itr_end;
      ++itr
    )
    {
      // check whether to add a new line before appending the next string
      // if a string was spread out over multiple lines, but did not end with a new line, then it usually looks better
      // to separate the previous line from the next
      // otherwise, add a new line if appending the current string would make the width >= the threshold
      if( last_str_was_continued_multiline || ( width > 0 && width + itr->size() >= THRESHOLD))
      {
        result += '\n'; // add a new line
        width = 0;      // now on a new line, so the width is 0
      }
      // if this string and the current line are not empty, add a space
      else if( itr->size() > 0 && width > 0)
      {
        result += ' '; // add a space to separate this string from the previous
        ++width;       // the current line is now one character longer with the space
      }

      const size_t last_newline( itr->rfind( '\n')); // find the first new line from the right
      last_str_was_continued_multiline =   // update whether the string had a new-line, but no new-line at the end
        bool
        (
           last_newline != std::string::npos  // there was a new line
           && last_newline != itr->size() - 1 // the last new-line was not the last char in the string
        );

      result += *itr;          // append the current string
      width += itr->size(); // and add its size to the width of the current line
    }

    return result;
  }

  //! @brief Try to register a file for this example to read or write
  //! @param EXAMPLE_NAME the example that is trying to access the file
  //! @param FILENAME the file that the example is attempting to access
  //! @param OPEN_MODE the manner in which FILENAME should be opened
  //!        for explanation on the types and use of open modes please see
  //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
  //! @param ERR_STREAM the stream to write error messages out to
  //! @return true if the file is openable by this example
  bool ExampleClass::RequestExampleFileAccess
  (
    const std::string &EXAMPLE_NAME,
    const std::string &FILENAME,
    const std::ios_base::openmode OPEN_MODE,
    std::ostream &ERR_STREAM
  )
  {
    // test for whether this file was already written to
    storage::Map< std::string, std::string>::const_iterator itr( GetExampleOutputFileMap().Find( FILENAME));

    // The example/example_files/bcl_objects directory often becomes cluttered, so after major refactorings, the
    // following set of commands can be used to clean it up:
    // BCL_Debug( FILENAME);
    // Then compile, then:
    // ./build/linux64_release/bcl-all-static.exe Examples | grep "DEBUG: FILENAME YIELDED: example/example_files/bcl_objects/" > used_files.txt
    // cat used_files.txt | sort | sed 's/DEBUG: FILENAME YIELDED: //' | uniq > used_files_full.txt
    // find example/example_files/bcl_objects/ -type f | sort > all_files.txt
    // diff -y used_files_full.txt all_files.txt --width 800 | grep '|' | awk '{print "\""$3"\""}' | tr '\n' ' ' > files_to_remove.txt
    // diff -y used_files_full.txt all_files.txt --width 800 | grep '>' | awk '{print "\""$2"\""}' | tr '\n' ' ' >> files_to_remove.txt
    // At this point, look over the files to remove and ensure that it contains only files that really should be removed
    // cat files_to_remove.txt | xargs svn rm
    // rm files_to_remove.txt used_files.txt used_files_full.txt all_files.txt
    // The same procedure could be used for the input and output files as well, but in this case it is necessary to run
    // the application examples.  Currently the app examples cannot be run serially; it is necessary to run them with
    // different calls to each app example, due to interference of the protein app examples on various global flags.
    if( itr != GetExampleOutputFileMap().End() && itr->second != EXAMPLE_NAME)
    {
      ERR_STREAM << EXAMPLE_NAME << " tried to access " << FILENAME << ", which was written to by " << itr->second;
      return false;
    }

    // Examples cannot write to files that were written or read by other examples, otherwise there is a race
    // condition
    if( OPEN_MODE & std::ios::out)
    {
      itr = GetExampleInputFileMap().Find( FILENAME);
      if( itr != GetExampleInputFileMap().End())
      {
        if( itr->second != EXAMPLE_NAME)
        {
          ERR_STREAM << EXAMPLE_NAME << " tried to write to " << FILENAME << ", which was read by " << itr->second;
          return false;
        }
      }
      else
      {
        // a new output file for the given example
        GetExampleOutputFileMap()[ FILENAME] = EXAMPLE_NAME;
      }
      return true;
    }
    else if( itr == GetExampleOutputFileMap().End())
    {
      // a new input file for the given example
      GetExampleInputFileMap()[ FILENAME] = EXAMPLE_NAME;
    }
    return true;
  }

  //! @brief overload for comparing size_t's to size_t's, which catches comparisons of size_t's to anything else
  //! @param LINE_NUMBER, EXAMPLE_CHECK_NAME, EXAMPLE_NAME, ACTUAL_CODE, ACTUAL_RESULT, EXPECTED_CODE, EXPECTED_RESULT see above
  //! @return true iff ACTUAL_RESULT == EXPECTED_RESULT
  bool ExampleClass::ExampleCheckReport
  (
    const std::string &EXAMPLE_CHECK_NAME,
    const std::string &EXAMPLE_NAME,
    const size_t &LINE_NUMBER,
    const char *const &ACTUAL_CODE,
    const size_t &ACTUAL_RESULT,
    const char *const &EXPECTED_CODE,
    const size_t &EXPECTED_RESULT
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

  //! @brief replace line breaks with given string
  //! @param STRING the string with line breaks
  //! @param NEW_LINE_PREFIX the string that is prepends to every new line
  //! @return a string with NEW_LINE_PREFIX inserted after each '\n' character
  std::string ExampleClass::FollowNewlinesWith( const std::string &STRING, const std::string &NEW_LINE_PREFIX)
  {
    if( STRING.find( '\n') == std::string::npos) // no new-lines, just return
    {
      return STRING;
    }
    else
    {
      return util::ReplaceString( STRING, "\n", "\n" + NEW_LINE_PREFIX);
    }
  }

  //! @brief open a file, report if the file is missing
  //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
  //! @param LINE_NUMBER line number of the example where the example check was called
  //! @param FILENAME file the example check was called on
  //! @param ISTREAM stream that is to be opened
  //! @return the resulting string
  //! @param OPEN_MODE the manner in which FILENAME should be opened
  //!        for explanation on the types and use of open modes please see
  //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
  //! @return the resulting string
  bool ExampleClass::ExampleMustOpenInputFile
  (
    const std::string &EXAMPLE_NAME,
    const size_t &LINE_NUMBER,
    const std::string &FILENAME,
    io::IFStream &ISTREAM,
    const std::ios_base::openmode OPEN_MODE
  )
  {
    // whether an input file opens or not does not constitute a real example check, so only insert the test if it failed

    // first, check whether this example is allowed to access the given file, this is used to prevent race conditions
    // when examples are run in parallel
    std::ostringstream message_stream;
    bool had_error( false);
    if( !RequestExampleFileAccess( EXAMPLE_NAME, FILENAME, OPEN_MODE, message_stream))
    {
      had_error = true;
    }
    else if( !io::File::TryOpenIFStream( ISTREAM, FILENAME, OPEN_MODE))
    {
      message_stream << "Example input file " << FILENAME
                     << " could not be opened; remainder of " << EXAMPLE_NAME
                     << " was skipped";
      had_error = true;
    }
    if( had_error)
    {
      // add the test
      ExampleClass::GetResults().InsertTest( EXAMPLE_NAME, LINE_NUMBER, false, message_stream.str());
      BCL_MessageCrt( "ERROR: " + message_stream.str());
      return false;
    }
    return true;
  }

  //! @brief open an output file, report if the file could not be opened for writing
  //! @param EXAMPLE_NAME result of GetClassIdentifier() called within example
  //! @param LINE_NUMBER line number of the example where the example check was called
  //! @param FILENAME file the example check was called on
  //! @param OSTREAM stream that is to be opened
  //! @param OPEN_MODE the manner in which FILENAME should be opened
  //!        for explanation on the types and use of open modes please see
  //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
  //! @return the resulting string
  bool ExampleClass::ExampleMustOpenOutputFile
  (
    const std::string &EXAMPLE_NAME,
    const size_t &LINE_NUMBER,
    const std::string &FILENAME,
    io::OFStream &OSTREAM,
    const std::ios_base::openmode OPEN_MODE
  )
  {
    // whether an output file opens or not does not constitute a real example check, so only insert the test if it failed

    // first, check whether this example is allowed to access the given file, this is used to prevent race conditions
    // when examples are run in parallel
    std::ostringstream message_stream;
    bool had_error( false);
    if( !RequestExampleFileAccess( EXAMPLE_NAME, FILENAME, OPEN_MODE, message_stream))
    {
      had_error = true;
    }
    else if( !io::File::TryOpenOFStream( OSTREAM, FILENAME, OPEN_MODE))
    {
      message_stream << "Example output file " << FILENAME
                     << " could not be opened; remainder of " << EXAMPLE_NAME
                     << " was skipped";
      had_error = true;
    }
    if( had_error)
    {
      // add the test
      ExampleClass::GetResults().InsertTest( EXAMPLE_NAME, LINE_NUMBER, false, message_stream.str());
      BCL_MessageCrt( "ERROR: " + message_stream.str());
      return false;
    }
    return true;
  }

  //! @brief HandleReportedResult a helper function for ExampleClass::ExampleCheckReport
  //! Example checks are checks that ACTUAL_RESULT == EXPECTED_RESULT
  //! This function is used with macros to automatically inform the user at the point of failure, what was executed, and
  //! what were the results, and that they were not equal
  //! To use this function, include example_candidate_macros.h and use one of them
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
  bool ExampleClass::HandleResult
  (
    const std::string &EXAMPLE_CHECK_NAME,
    const std::string &EXAMPLE_NAME,
    const size_t      &LINE_NUMBER,
    const char *const &ACTUAL_CODE,
    const std::string &ACTUAL_RESULT,
    const char *const &EXPECTED_CODE,
    const std::string &EXPECTED_RESULT,
    const bool        &SUCCESS,
    const double      &RELATIVE_TOLERANCE,
    const double      &ABSOLUTE_TOLERANCE
  )
  {
    // maximum line width in the output. Does not count the indent string, which is ~10 characters.
    const size_t split_line_width( 110);

    // maximum length of output actual_result except at debug level
    const size_t max_output_length_non_debug( 440);

    // remove any extraneous space from the codes
    const std::string expected_code( util::TrimString( EXPECTED_CODE));
    const std::string actual_code( util::TrimString( ACTUAL_CODE));

    // store whether EXPECTED_CODE was a quoted string equivalent to EXPECTED_RESULT
    const bool expected_code_was_string
    (
      expected_code.size() == EXPECTED_RESULT.size() + 2 // is it the right size
      && ( expected_code[ 0] == '\'' || expected_code[ 0] == '\"') // did it start and end with single or double quotes
      && expected_code[ expected_code.size() - 1] == expected_code[ 0]
      && expected_code.substr( 1, expected_code.size() - 2) == EXPECTED_RESULT // is everything in between the quotes the same as the result
    );

    // store whether EXPECTED_RESULT is essentially equivalent to expected_code
    // in which case writing out both EXPECTED_RESULT and expected_code would be redundant
    const bool expected_is_a_simple_value
               (
                 EXPECTED_RESULT == expected_code
                 || expected_code_was_string
                 || expected_code == "true"
                 || expected_code == "false"
               );

    const std::string tolerance_string
    (
      RELATIVE_TOLERANCE > 0.0
      ? " (+/- " + util::Format()( RELATIVE_TOLERANCE) + "%)"
      : ABSOLUTE_TOLERANCE > 0.0
        ? " (+/- " + util::Format()( ABSOLUTE_TOLERANCE) + ")"
        : ""
    );

    if( SUCCESS) // Inform the user
    {
      static const std::string s_GoodTestLine( "              ");

      // report the result to the ExampleClass
      ExampleClass::GetResults().InsertTest( EXAMPLE_NAME, LINE_NUMBER, true, "");

      // report the successful test to the user.  Here are a few example macros and what should be output if
      // everything worked correctly
      //
      // BCL_ExampleIndirectCheck carbon.GetAtomicNumber(), 6, "Set atom type");
      // =std=bcl::app=> Set atom type test succeeded: carbon.GetAtomicNumber() == 6
      //
      // BCL_ExampleCheck( carbon.GetAtomicNumber(), 6);
      // =std=bcl::app=> successful example check: carbon.GetAtomicNumber() == 6
      //
      // BCL_ExampleCheck( carbon.GetAtomicNumber(), Atom( carbon).GetAtomicNumber());
      // the message would read:
      // =std=bcl::app=> successful example check: carbon.GetAtomicNumber() == Atom( carbon).GetAtomicNumber() == 6
      //
      // BCL_ExampleCheckWithinTolerance( result, expected, 0.001)
      // =std=bcl::app=> successful example check: result == expected == 6 (within tolerance of: 0.001)
      storage::Vector< std::string> msg_segments;
      msg_segments.AllocateMemory( 6);
      msg_segments.PushBack
                   (
                     EXAMPLE_CHECK_NAME.size() == 0
                     ? std::string( "successful example check:\n")
                     : EXAMPLE_CHECK_NAME + " test succeeded:\n"
                   );

      msg_segments.PushBack( actual_code); // the snippet of code used to compute the actual result
      if( !expected_is_a_simple_value) // show the code that was produced the expected result
      {
        msg_segments.PushBack( "==");
        msg_segments.PushBack( expected_code); // the snippet of code used to compute the expected result
      }

      // only append the actual result if it would be less than 4 lines
      if
      (
        ACTUAL_RESULT.size() > max_output_length_non_debug
        && util::GetMessenger().GetMessageVerbosity() != util::Message::e_Detail
      )
      {
        msg_segments.PushBack
        (
          "output of result suppressed (result was " + util::Format()( ACTUAL_RESULT.size())
          + " bytes), set message verbosity to Detail to override"
        );
      }
      else
      {
        msg_segments.PushBack( "==");
        msg_segments.PushBack( ACTUAL_RESULT); // the actual result string
      }

      // Append snippet with tolerance information if tolerance was used
      if( tolerance_string.size())
      {
        msg_segments.PushBack( tolerance_string);
      }

      std::string msg( AddNewLinesWhenLineWidthIsExcessive( msg_segments, split_line_width));

      // remove a new line if there is only 1 and it would fit cleanly on one line
      if( msg.size() - s_GoodTestLine.size() < split_line_width && std::count( msg.begin(), msg.end(), '\n') == 1)
      {
        msg[ msg.find( '\n')] = ' ';
      }

      BCL_MessageStd( FollowNewlinesWith( msg, s_GoodTestLine));
    }
    else // failure.  write out an error message immediately and put the failed test into the table
    {
      // string to put before every line in the table if the example failed
      static const std::string s_ErrLine( "  XXXX        ");

      // If the macro call were
      // BCL_ExampleCheck( carbon.GetAtomicNumber(), 6);
      // and carbon.GetAtomicNumber() returned 7, the message would read:
      // =std=bcl::app=> FAILED example check: carbon.GetAtomicNumber() returned 7 instead of 6

      // BCL_ExampleCheck( carbon.GetAtomicNumber(), Atom( carbon).GetAtomicNumber());
      // the message would read:
      // =std=bcl::app=> FAILED example check:
      //  XXXX           carbon.GetAtomicNumber() returned 6 but Atom( carbon).GetAtomicNumber() returned 7
      storage::Vector< std::string> msg_segments;
      msg_segments.AllocateMemory( 8);

      // if the user did not pass an example check name, write "FAILED example check"
      // otherwise, the example is indirect, and what's really being tested is EXAMPLE_CHECK_NAME, so write
      // EXAMPLE_CHECK_NAME + " test FAILED:\n"
      msg_segments.PushBack
                   (
                     EXAMPLE_CHECK_NAME.size() == 0
                     ? std::string( "FAILED example check:\n")
                     : EXAMPLE_CHECK_NAME + " test FAILED:\n"
                   );
      msg_segments.PushBack( actual_code); // the snippet of code used to compute the actual result
      msg_segments.PushBack( "returned");
      msg_segments.PushBack( ACTUAL_RESULT); // the actual result (object)

      if( expected_is_a_simple_value) // then we only need to output the expected result
      {
        msg_segments.PushBack( "instead of");
      }
      else
      {
        msg_segments.PushBack( "but");
        msg_segments.PushBack( expected_code); // the snippet of code used to compute the expected result
        msg_segments.PushBack( "returned");
      }
      msg_segments.PushBack( EXPECTED_RESULT); // the expected result

      // Or for the cases with tolerance
      // BCL_ExampleCheckWithinTolerance( result(), expected(), 0.001);
      // the message would read:
      // =std=bcl::app=> FAILED example check:
      //  XXXX           result() returned 6 but expected() returned 7 (with tolerance = 0.001)

      // Append snippet with tolerance information if tolerance was used
      if( tolerance_string.size())
      {
        msg_segments.PushBack( tolerance_string);
      }

      std::string msg( AddNewLinesWhenLineWidthIsExcessive( msg_segments, split_line_width));

      // replace a new line with a space if there is only 1 new line and the rest of the message
      // would fit cleanly on the line right after failed example check
      if( msg.size() - s_ErrLine.size() < split_line_width && std::count( msg.begin(), msg.end(), '\n') == 1)
      {
        msg[ msg.find( '\n')] = ' ';
      }

      // add the test
      ExampleClass::GetResults().InsertTest( EXAMPLE_NAME, LINE_NUMBER, false, FollowNewlinesWith( msg, s_ErrLine));

      // print the error message immediately as well, along with the line number and example as well, so the user
      // can easily locate the failed example in the messages.  Also let them see the complete error message and line
      // number at this point
      BCL_MessageCrt
      (
        FollowNewlinesWith
        (
          "\nExample " + RemoveUnnecessaryString( EXAMPLE_NAME)
          + " failed at line #" + util::Format()( LINE_NUMBER)
          + '\n' + msg,
          s_ErrLine
        )
      );
    }
    return SUCCESS;
  }

  //! @brief get the map from file name to example that wrote to it
  //! @return the map from file name to example that wrote to it
  //! @detail to prevent race conditions when examples are run in parallel, examples may not write to the same
  //! @detail output file that another example has read or written to
  //! @detail the base path for all entries is the main bcl path
  storage::Map< std::string, std::string> &ExampleClass::GetExampleOutputFileMap()
  {
    static storage::Map< std::string, std::string> s_map;
    return s_map;
  }

  //! @brief get the map from file name to an example that read it
  //! @return the map from file name to an example that read it
  //! @detail to prevent race conditions, an example cannot read a file that a different example writes
  //! @detail the base path for all entries is the main bcl path
  storage::Map< std::string, std::string> &ExampleClass::GetExampleInputFileMap()
  {
    static storage::Map< std::string, std::string> s_map;
    return s_map;
  }

  ExampleClass &GetExamples()
  {
    return ExampleClass::GetEnums();
  }

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< ExampleInterface>, ExampleClass>;

  } // namespace util
} // namespace bcl
