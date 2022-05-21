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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "util/bcl_util_string_functions.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically
#include <cstring>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_string_functions.cpp
  //!
  //! @author woetzen, mendenjl
  //! @date Nov 06, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilStringFunctions :
    public ExampleInterface
  {

  public:

    ExampleUtilStringFunctions *Clone() const
    {
      return new ExampleUtilStringFunctions( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      std::string first( " aa bb cc dd ee ");
      std::string second( "zaazbbzcczddzeez");
      std::string third( " \n");
      BCL_MessageStd( "first string: " + first);
      BCL_MessageStd( "split first by spaces: " + util::Format()( util::SplitString( first)));
      BCL_MessageStd( "trim first: " + util::TrimString( first));
      BCL_MessageStd( "remove spaces from first: " + util::RemoveSpacesFromString( first));
      BCL_MessageStd( "trim third: " + util::TrimString( third));

      BCL_MessageStd( "second string: " + second);
      BCL_MessageStd( "split second by \"z\": " + util::Format()( util::SplitString( second, "z")));

      std::string a( "1234.5");
      BCL_MessageStd( "string a: " + a);
      BCL_MessageStd( "is numerical? " + std::string( util::IsNumerical( a) ? "true" : "false"));

      std::string d( "7000.067");
      BCL_MessageStd( "string d: " + d);
      BCL_MessageStd
      (
        "is numerical? " + std::string( util::IsNumerical( d) ? "true" : "false")
      );

      BCL_MessageStd
      (
        "a + d = " +
        util::Format()
        (
          util::ConvertStringToNumericalValue< double>( a)
          + util::ConvertStringToNumericalValue< double>( d)
        )
      );

      std::string b( "ab12cd34ef ");
      BCL_MessageStd( "string b: " + b);
      BCL_MessageStd( "is numerical? " + std::string( util::IsNumerical( b) ? "true" : "false"));

      BCL_ExampleCheck( util::TrimString( " h e l l o "), "h e l l o");
      BCL_ExampleCheck( util::RemoveSpacesFromString( " h e l l o "), "hello");
      BCL_ExampleCheck
      (
        util::ConvertStringToNumericalValue< double>( " 5.5"),
        double( 5.5)
      );

      // StartsWith / EndsWith
      BCL_ExampleCheck( util::StartsWith( "Apple", "App"), true);
      BCL_ExampleCheck( util::StartsWith( "Apple", "Apple"), true);
      BCL_ExampleCheck( util::StartsWith( "Apple", "AppleMac"), false);
      BCL_ExampleCheck( util::StartsWith( "Apple", "ple"), false);
      BCL_ExampleCheck( util::EndsWith( "Apple", "Apple"), true);
      BCL_ExampleCheck( util::EndsWith( "Apple", "App"), false);
      BCL_ExampleCheck( util::EndsWith( "Apple", "ple"), true);
      BCL_ExampleCheck( util::EndsWith( "Apple", ""), true);
      BCL_ExampleCheck( util::StartsWith( "Apple", ""), true);

      BCL_ExampleCheck( util::ToLower( "ApPlEs And ORANGEs"), "apples and oranges");
      BCL_ExampleCheck( util::ToUpper( "appLeS And ORANGEs"), "APPLES AND ORANGES");

      // make a vector containing strings to test util::IsNumerical with.
      // For all of these strings, IsNumerical(STRING) should return true
      std::vector< std::string> should_be_convertable_to_double;
      should_be_convertable_to_double.push_back( "1");
      should_be_convertable_to_double.push_back( "2e3");
      should_be_convertable_to_double.push_back( "4e-5");
      should_be_convertable_to_double.push_back( " -6.04362354 ");
      should_be_convertable_to_double.push_back( "7.e8  ");
      should_be_convertable_to_double.push_back( "  +9.0e+10");
      should_be_convertable_to_double.push_back( "  +9.0e+10");
      should_be_convertable_to_double.push_back( "+5e-07");
      should_be_convertable_to_double.push_back( " +5.5");
      should_be_convertable_to_double.push_back( " -5.5");
      should_be_convertable_to_double.push_back( "5.");
      should_be_convertable_to_double.push_back( "\t5.5\t");

      for( size_t test_number = 0; test_number < should_be_convertable_to_double.size(); test_number++)
      {
        std::string &test_string = should_be_convertable_to_double[ test_number];
        BCL_MessageStd
        (
          "is " + test_string + " a numerical value? : "
          + util::Format()( util::IsNumerical( test_string))
          + ", it converts to "
          + util::Format()( util::ConvertStringToNumericalValue<double>( test_string))
        );
        BCL_Example_Check
        (
          util::IsNumerical( should_be_convertable_to_double[ test_number]),
          test_string + " should have been recognized as a double!"
        );
      }

      // make a vector containing strings to test util::IsNumerical with.
      // For all of these strings, IsNumerical(STRING) should return false
      std::vector< std::string> should_not_be_convertable_to_double;
      should_not_be_convertable_to_double.push_back( "1e");
      should_not_be_convertable_to_double.push_back( "2e3+");
      should_not_be_convertable_to_double.push_back( "4.+e-5");
      should_not_be_convertable_to_double.push_back( " 6.04362354.5 ");
      should_not_be_convertable_to_double.push_back( "7.e8.5  ");
      should_not_be_convertable_to_double.push_back( "  ++9.0e+10");
      should_not_be_convertable_to_double.push_back( "  -+9.0e+10");
      should_not_be_convertable_to_double.push_back( "+5.a");
      should_not_be_convertable_to_double.push_back( "a5");
      should_not_be_convertable_to_double.push_back( "5. 5");
      should_not_be_convertable_to_double.push_back( " +5..5");
      should_not_be_convertable_to_double.push_back( " +-5.5");
      should_not_be_convertable_to_double.push_back( " -+5.5");
      should_not_be_convertable_to_double.push_back( " --5.5");
      should_not_be_convertable_to_double.push_back( " ++5.5");
      should_not_be_convertable_to_double.push_back( " +5.5+");

      for( size_t test_number = 0; test_number < should_not_be_convertable_to_double.size(); test_number++)
      {
        std::string &test_string = should_not_be_convertable_to_double[ test_number];
        BCL_MessageStd
        (
          "was " + test_string
          + " seen as a numerical value? : "
          + util::Format()( util::IsNumerical( test_string))
        );
        BCL_Example_Check
        (
          !util::IsNumerical( test_string),
          test_string + " should not have been recognized as a double!"
        );
      }

      // make a vector containing pairs of strings to test util::LengthOfUnsignedIntegerType with.  the second member
      // is the expected answer
      std::vector< std::pair< std::string, size_t> > size_t_length_test;
      size_t_length_test.push_back( std::pair< std::string, size_t>( " 1 ", std::strlen( " 1 ")));
      size_t_length_test.push_back( std::pair< std::string, size_t>( "1", std::strlen( "1")));
      size_t_length_test.push_back( std::pair< std::string, size_t>( " 100 A B C D ", std::strlen( " 100 ")));
      size_t_length_test.push_back( std::pair< std::string, size_t>( " -100 A B C D ", std::strlen( ""))); // negative number can't be converted to
      // size_t
      size_t_length_test.push_back( std::pair< std::string, size_t>( " A B C D 100 ", std::strlen( ""))); // can't convert " A " to size_t!

      for( size_t test_number = 0; test_number < size_t_length_test.size(); test_number++)
      {
        BCL_ExampleIndirectCheck
        (
          util::LengthOfUnsignedIntegerType( size_t_length_test[ test_number].first),
          size_t_length_test[ test_number].second,
          "util::LengthOfUnsignedIntegerType starting at index 0 in " + size_t_length_test[ test_number].first
        );
      }

      // StringListFromIStream
      std::istringstream iss( "  ");
      BCL_ExampleIndirectCheck( util::StringListFromIStream( iss).GetSize(), 0, "iss.str( \"  \")");
      iss.clear();
      iss.str( "'  '");
      BCL_ExampleIndirectCheck( util::StringListFromIStream( iss).GetSize(), 1, "iss.str( \"'  '\")");
      iss.clear();
      iss.str( "'  \\'' abra \"ca da be\r'\\\"\\\nra\"");
      BCL_ExampleIndirectCheck
      (
        util::StringListFromIStream( iss),
        storage::Vector< std::string>::Create( "  '", "abra", "ca da be\r\'\"\nra"),
        "iss.str( \"'  \\'' abra \"ca da be\r\'\\\"\nra\"\")"
      );

      // make a vector containing pairs of strings to test util::LengthOfIntegerType with.  the second member
      // is the expected answer
      std::vector< std::pair< std::string, size_t> > long_length_test;
      long_length_test.push_back( std::pair< std::string, size_t>( " +1 ", std::strlen( " +1 ")));
      long_length_test.push_back( std::pair< std::string, size_t>( "-1-", std::strlen( "-1")));
      long_length_test.push_back( std::pair< std::string, size_t>( " 100 A B C D ", std::strlen( " 100 ")));

      // negative numbers can be converted to long, so length should be non-zero here
      long_length_test.push_back( std::pair< std::string, size_t>( " -100 A B C D ", std::strlen( " -100 ")));
      long_length_test.push_back( std::pair< std::string, size_t>( " A B C D 100 ", 0UL)); // can't convert " A " to size_t!

      for( size_t test_number = 0; test_number < long_length_test.size(); test_number++)
      {
        BCL_Example_Check
        (
          util::LengthOfIntegerType( long_length_test[ test_number].first) == long_length_test[ test_number].second,
          long_length_test[ test_number].first + " was the wrong size! (was "
          + util::Format()( util::LengthOfIntegerType( long_length_test[ test_number].first))
          + " instead of " + util::Format()( long_length_test[ test_number].second) + ")"
        );
      }

      // try out join
      BCL_ExampleCheck( util::Join( " ", storage::Vector< std::string>::Create( "Hello", "World")), "Hello World");
      BCL_ExampleCheck( util::Join( " ", storage::Vector< std::string>::Create( "Hello World")), "Hello World");
      BCL_ExampleCheck
      (
        util::Join( ", ", storage::Vector< std::string>::Create( "Red", "White", "Blue")),
        "Red, White, Blue"
      );

      first  = "PQRSTALPHAGAM7 MAKKD";
      second = "XYZLAPHAGAM7 MAZ";
      BCL_MessageStd( "first string:  " + first);
      BCL_MessageStd( "second string: " + second);
      storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> > trip
      (
        util::MaximumMatchingSubstring( first, second, 3)
      );
      BCL_MessageStd( "maximum match of second within first: " + trip.First());
      BCL_MessageStd( "test 1: " + first.substr( trip.Second().First(), trip.Second().Second()));
      BCL_MessageStd( "test 2: " + second.substr( trip.Third().First(), trip.Third().Second()));
      trip = util::MaximumMatchingSubstring( second, first, 3);
      BCL_MessageStd( "maximum match of first within second: " + trip.First());
      BCL_MessageStd( "test 1: " + second.substr( trip.Second().First(), trip.Second().Second()));
      BCL_MessageStd( "test 2: " + first.substr( trip.Third().First(), trip.Third().Second()));

      return 0;
    } //end ExampleClass::ExampleUtilStringFunctions

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleUtilStringFunctions

  const ExampleClass::EnumType ExampleUtilStringFunctions::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilStringFunctions())
  );
} // namespace bcl

