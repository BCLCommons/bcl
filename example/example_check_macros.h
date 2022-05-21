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

#ifndef EXAMPLE_CHECK_MACROS_H_
#define EXAMPLE_CHECK_MACROS_H_

/////////////////////////////////////////////////////////////////////////////////////
//!
//! @brief This header contains macros that perform example checks
//!
//! @author mendenjl
//!
//! @date 06/14/2010
//!
/////////////////////////////////////////////////////////////////////////////////////
// Because macros are scopeless, it would be misleading to put them in namespace bcl

// A macro that standardizes and automatically writes appropriate messages for BCL_Example Checks
// This version optionally reports the actual code used in successful tests, and the resulting objects
// Always reports the actual code used in failed tests, and the resulting objects
#define BCL_General_Example_Check( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME) \
  ExampleClass::ExampleCheckReport( std::string( EXAMPLE_CHECK_NAME), GetClassIdentifier(), __LINE__, \
    #ACTUAL_RESULT, ACTUAL_RESULT, \
    #EXPECTED_RESULT, EXPECTED_RESULT)

// return from example with error code 1 if BCL_General_Example_Check failed
#define BCL_General_Example_Assert( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME) \
  if( !BCL_General_Example_Check( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME)) { return 1;}

// A macro that standardizes and automatically writes appropriate messages for BCL_Example Checks
// This version optionally reports the actual code used in successful tests, and the resulting objects
// Always reports the actual code used in failed tests, and the resulting objects
#define BCL_General_Example_Check_With_Tolerance( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE) \
  ExampleClass::ExampleCheckReportWithTolerance( std::string( EXAMPLE_CHECK_NAME), GetClassIdentifier(), __LINE__, \
    #ACTUAL_RESULT, ACTUAL_RESULT, \
    #EXPECTED_RESULT, EXPECTED_RESULT, \
    ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE)

// A macro that standardizes and in many cases automatically writes appropriate messages for BCL_Example Checks

// This version reports to the user what was tested and the result, so they know what functionality can be used
// and can easily see what went wrong (if anything)
#define BCL_ExampleCheck( ACTUAL_RESULT, EXPECTED_RESULT) \
  BCL_General_Example_Check( ACTUAL_RESULT, EXPECTED_RESULT, "")

// Like BCL_ExampleCheck, but has an extra string parameter, which should indicate what this example
// is indirectly checking
#define BCL_ExampleIndirectCheck(ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME) \
  BCL_General_Example_Check( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME)

// This version reports to the user what was tested and the result, so they know what functionality can be used
// and can easily see what went wrong (if anything)
#define BCL_ExampleCheckWithinTolerance( ACTUAL_RESULT, EXPECTED_RESULT, TOLERANCE) \
  BCL_General_Example_Check_With_Tolerance( ACTUAL_RESULT, EXPECTED_RESULT, "", std::numeric_limits< double>::epsilon(), TOLERANCE)

// Like BCL_ExampleCheck, but has an extra string parameter, which should indicate what this example
// is indirectly checking
#define BCL_ExampleIndirectCheckWithinTolerance( ACTUAL_RESULT, EXPECTED_RESULT, TOLERANCE, EXAMPLE_CHECK_NAME) \
  BCL_General_Example_Check_With_Tolerance( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME, std::numeric_limits< double>::epsilon(), TOLERANCE)

// This version reports to the user what was tested and the result, so they know what functionality can be used
// and can easily see what went wrong (if anything)
#define BCL_ExampleCheckWithinAbsTolerance( ACTUAL_RESULT, EXPECTED_RESULT, TOLERANCE) \
  BCL_General_Example_Check_With_Tolerance( ACTUAL_RESULT, EXPECTED_RESULT, "", TOLERANCE, 0.0)

// Like BCL_ExampleCheck, but has an extra string parameter, which should indicate what this example
// is indirectly checking
#define BCL_ExampleIndirectCheckWithinAbsTolerance( ACTUAL_RESULT, EXPECTED_RESULT, TOLERANCE, EXAMPLE_CHECK_NAME) \
  BCL_General_Example_Check_With_Tolerance( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME, TOLERANCE, 0.0)

// Similar to BCL_ExampleCheck, but immediately returns from the example (with 1) to indicate a severe failure
#define BCL_ExampleAssert( ACTUAL_RESULT, EXPECTED_RESULT) \
  BCL_General_Example_Assert( ACTUAL_RESULT, EXPECTED_RESULT, "")

// Similar to BCL_ExampleCheck, but immediately returns from the example (with 1) to indicate a severe failure
#define BCL_ExampleIndirectAssert( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME) \
  BCL_General_Example_Assert( ACTUAL_RESULT, EXPECTED_RESULT, EXAMPLE_CHECK_NAME)

// Try opening an input file; does not count as an example check because the presence of an
// input file does not indicate whether the class works (except e.g. io::IFstream)
// on failure, counts as a failed example check and prints out the filename could not be opened for reading
#define BCL_ExampleMustOpenInputFile( STREAM, FILENAME) \
  if( !ExampleClass::ExampleMustOpenInputFile( GetClassIdentifier(), __LINE__, FILENAME, STREAM)) { return 1;}

// Try opening a binary input file; does not count as an example check because the presence of an
// input file does not indicate whether the class works (except e.g. io::IFstream)
// on failure, counts as a failed example check and prints out the filename could not be opened for reading
#define BCL_ExampleMustOpenBinaryInputFile( STREAM, FILENAME) \
  if( !ExampleClass::ExampleMustOpenInputFile( GetClassIdentifier(), __LINE__, FILENAME, STREAM, std::ios::binary))\
  { return 1;}

// Try opening an output file; does not count as an example check because the ability to open an output file does not
// not test the class (except e.g. io::OFstream)
// on failure, counts as a failed example check and prints out the filename that could not be opened for writing
#define BCL_ExampleMustOpenOutputFile( STREAM, FILENAME) \
  if( !ExampleClass::ExampleMustOpenOutputFile( GetClassIdentifier(), __LINE__, FILENAME, STREAM)) { return 1;}

// Try opening a binary output file; does not count as an example check because the ability to open an output file does
// not test the class (except e.g. io::OFstream)
// on failure, counts as a failed example check and prints out the filename that could not be opened for writing
#define BCL_ExampleMustOpenBinaryOutputFile( STREAM, FILENAME) \
  if( !ExampleClass::ExampleMustOpenOutputFile( GetClassIdentifier(), __LINE__, FILENAME, STREAM, std::ios::binary))\
  { return 1;}

// Example check that which uses a severity level, condition, and complete error string
#define BCL_Example_Check( CONDITION, ERROR_STRING) \
  ExampleClass::GetResults().InsertTest( GetClassIdentifier(), __LINE__, ( CONDITION), ERROR_STRING);

#endif // EXAMPLE_CHECK_MACROS_H_
