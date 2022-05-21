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
#include "io/bcl_io_serialization_builtin.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_builtin.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationBuiltin :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationBuiltin *Clone() const
    {
      return new ExampleIoSerializationBuiltin( *this);
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
      // make a double and a test string for the label handlers to change
      double test_double( 5.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // SerializationBuiltin( const double *DATA)
      io::SerializationBuiltin< double> serialization_handler_dbl( &test_double);

      // label handlers made with this constructor are typically only used inside containers
      // since they don't provide info about what is labeled
      BCL_ExampleCheck( serialization_handler_dbl.GetLabel().ToString(), "5");
      std::ostringstream err_stream_dbl;
      // set serialization_handler_dbl with a legitimate double inside its range
      serialization_handler_dbl.TryRead( util::ObjectDataLabel( "0.5"), err_stream_dbl);
      BCL_ExampleIndirectCheck
      (
        err_stream_dbl.str(),
        "",
        "serialization_handler_dbl.TryRead( util::ObjectDataLabel( 0.5), error_from_serialization_handler_dbl)"
      );
      // check that set parameter succeeded
      BCL_ExampleIndirectCheck( test_double, double( 0.5), "TryRead (valid case)");
      err_stream_dbl.str( "");
      // now try something that isn't even a double
      BCL_ExampleCheck
      (
        serialization_handler_dbl.TryRead( util::ObjectDataLabel( "test"), err_stream_dbl),
        false
      );
      BCL_ExampleIndirectCheck
      (
        err_stream_dbl.str().size() > 0,
        true,
        "serialization_handler_dbl.TryRead( util::ObjectDataLabel( test), error_from_serialization_handler_dbl)"
      );

      // try nan, and +inf
      err_stream_dbl.str( "");
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "NAN"), err_stream_dbl), true);
      BCL_ExampleIndirectCheck( util::IsDefined( test_double), false, "Setting NAN");
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "nan"), err_stream_dbl), true);
      BCL_ExampleIndirectCheck( util::IsDefined( test_double), false, "Setting nan");
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "+inf"), err_stream_dbl), true);
      BCL_ExampleIndirectCheck( test_double > std::numeric_limits< double>::max(), true, "Setting +inf");
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "inf"), err_stream_dbl), true);
      BCL_ExampleIndirectCheck( test_double > std::numeric_limits< double>::max(), true, "Setting inf");
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "-INF"), err_stream_dbl), true);
      BCL_ExampleIndirectCheck( test_double < -std::numeric_limits< double>::max(), true, "Setting -INF");

      // make another double and string for a fully parameterized label handler
      double my_val( 3.14);

      // SerializationBuiltin( const command::ParameterCheckInterface &PARAMETER_CHECK, const t_DataType *DATA)
      io::SerializationBuiltin< double> serialization_handler_dbl_param( &my_val);

      // label handlers made with this constructor are commonly used in parameterized objects
      BCL_ExampleCheck( serialization_handler_dbl_param.GetLabel().ToString(), "3.14");

      // make a couple test labels to try with the object data label version
      util::ObjectDataLabel test_trivial( "test");
      util::ObjectDataLabel test_basic_object( "Alpha(Omega,Delta,Kappa)");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      io::SerializationBuiltin< util::ObjectDataLabel> serialization_handler_test
      (
        &test_trivial
      );
      io::SerializationBuiltin< util::ObjectDataLabel> serialization_handler_test_basic_object
      (
        &test_basic_object
      );
      BCL_ExampleCheck( serialization_handler_test.GetLabel(), test_trivial);
      BCL_ExampleCheck( serialization_handler_test_basic_object.GetLabel(), test_basic_object);

      // TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM) for strings
      std::ostringstream error_from_serialization_handler;
      // try reading a legitimate object
      BCL_ExampleCheck
      (
        serialization_handler_test.TryRead( test_basic_object, error_from_serialization_handler),
        true
      );
      BCL_ExampleIndirectCheck( test_trivial, test_basic_object, "TryRead");

      // specialized tests for string type
      std::string test_string_needs_quotes( "test=should-have-quotes");
      std::string test_string_escape_quotes( "test=should-escape-\"quotes");

      io::SerializationBuiltin< std::string> serialization_handler_should_escape
      (
        &test_string_needs_quotes
      );
      BCL_ExampleCheck( serialization_handler_should_escape.GetLabel().ToString(), "\"" + test_string_needs_quotes + "\"");
      test_string_needs_quotes = test_string_escape_quotes;
      BCL_ExampleCheck( serialization_handler_should_escape.GetLabel().ToString(), "\"test=should-escape-\\\"quotes\"");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationBuiltin

  const ExampleClass::EnumType ExampleIoSerializationBuiltin::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationBuiltin())
  );

} // namespace bcl
