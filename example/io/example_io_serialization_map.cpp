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
#include "io/bcl_io_serialization_map.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_map.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationMap :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationMap *Clone() const
    {
      return new ExampleIoSerializationMap( *this);
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
      // make some vectors for the label to set
      storage::Map< std::string, size_t> map_string_to_int;
      storage::Map< std::string, double> constants_map;

      std::string test_map_string_to_int_string( "(alpha=5, beta=10, gamma=49)");
      std::string test_constants_map_string( "(c = 3.0e8, pi = 3.14159)");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a label handler that will set map_string_to_int
      io::SerializationMap< storage::Map< std::string, size_t> > map_string_to_int_handler
      (
        *io::Serialization::GetAgent( ( std::string *)( NULL)),
        *io::Serialization::GetAgentWithRange( ( size_t *)( NULL), size_t( 0), size_t( 50)),
        &map_string_to_int
      );

      // construct a label handler that will set constants_map
      io::SerializationMap< storage::Map< std::string, double> >
      constants_map_handler
      (
        *io::Serialization::GetAgent( ( std::string *)( NULL)),
        *io::Serialization::GetAgent( ( double *)( NULL)),
        &constants_map
      );

      // Check that () are not shown for empty maps
      BCL_ExampleCheck( map_string_to_int_handler.GetLabel().ToString(), "\"\"");
      BCL_ExampleCheck( constants_map_handler.GetLabel().ToString(), "\"\"");

      std::ostringstream error_from_map_string_to_int_handler_handler;
      // set coefficient_handler with a legitimate vector of size_ts inside its range
      map_string_to_int_handler.TryRead
      (
        util::ObjectDataLabel( test_map_string_to_int_string),
        error_from_map_string_to_int_handler_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_map_string_to_int_handler_handler.str(),
        "",
        "coefficient_handler.TryRead( util::ObjectDataLabel( test_map_string_to_int_string), ...)"
      );
      storage::Map< std::string, size_t> expected_map;
      expected_map[ "alpha"] = 5;
      expected_map[ "beta"] = 10;
      expected_map[ "gamma"] = 49;
      // check that set parameter succeeded
      BCL_ExampleIndirectCheck( map_string_to_int, expected_map, "TryRead (valid case)");
      error_from_map_string_to_int_handler_handler.str( "");
      // now try a double outside the range
      map_string_to_int_handler.TryRead
      (
        util::ObjectDataLabel( "(alpha=51)"),
        error_from_map_string_to_int_handler_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_map_string_to_int_handler_handler.str().size() > 0,
        true,
        "coefficient_handler.TryRead( util::ObjectDataLabel( 51), error_from_coefficient_handler)"
      );
      // check that set parameter failed
      BCL_ExampleIndirectCheck
      (
        map_string_to_int,
        expected_map,
        "TryRead (invalid case, should not change internal map)"
      );
      error_from_map_string_to_int_handler_handler.str( "");
      // now try something that isn't even a size_t
      map_string_to_int_handler.TryRead
      (
        util::ObjectDataLabel( "(alpha)"),
        error_from_map_string_to_int_handler_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_map_string_to_int_handler_handler.str().size() > 0,
        true,
        "map_string_to_int_handler.TryRead( util::ObjectDataLabel( \"alpha\"), ...)"
      );

      error_from_map_string_to_int_handler_handler.str( "");
      // check that set parameter failed
      BCL_ExampleIndirectCheck
      (
        map_string_to_int,
        expected_map,
        "TryRead (invalid case, should not change internal map)"
      );

      // TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM) for 2a different type of map
      std::ostringstream error_from_constants_map_handler;
      // set coefficient_2d_handler with a legitimate 2d vector with #s inside its range
      constants_map_handler.TryRead
      (
        util::ObjectDataLabel( test_constants_map_string),
        error_from_constants_map_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_constants_map_handler.str(),
        "",
        "constants_map_handler.TryRead( "
        + test_constants_map_string + ", error_from_constants_map_handler)"
      );

      // test the get string functionality; which converts the handled object back into a string

      // label handlers made with this constructor are commonly used in parameterized objects
      BCL_ExampleCheck( map_string_to_int_handler.GetLabel().ToString(), "(alpha=5,beta=10,gamma=49)");
      BCL_ExampleCheckWithinTolerance( constants_map[ "c"], 3.0e8, 0.001);
      BCL_ExampleCheckWithinTolerance( constants_map[ "pi"], 3.14159, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationMap

  const ExampleClass::EnumType ExampleIoSerializationMap::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationMap())
  );

} // namespace bcl
