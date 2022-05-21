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
#include "io/bcl_io_serialization_container.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_container.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationContainer :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationContainer *Clone() const
    {
      return new ExampleIoSerializationContainer( *this);
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
      storage::Vector< size_t> coefficients;
      storage::Vector< storage::Vector< size_t> > coefficients_2d;
      storage::Vector< math::Range< size_t> > ranges;

      std::string test_ranges_string( "([0], \"[8,9)\", \"(49,57)\" )");
      std::string test_coefficient_string( "(0, 9, 49)");
      std::string test_coefficients_2d_string( "( (0, 9, 123), (16, 12) )");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a label handler that will set the test coefficients vector
      io::SerializationContainer< storage::Vector< size_t> > coefficient_handler
      (
        *io::Serialization::GetAgentWithRange< size_t>( 0, 50),
        &coefficients
      );

      // SerializationContainer( const command::ParameterCheckInterface &PARAMETER_CHECK, const std::string *DATA)
      io::SerializationContainer< storage::Vector< storage::Vector< size_t> > >
      coefficient_2d_handler
      (
        io::SerializationContainer< storage::Vector< size_t> >
        (
          *io::Serialization::GetAgentWithRange< size_t>( 0, 150)
        ),
        &coefficients_2d
      );

      // Serialization handlers for containers should output a null value
      // if there is nothing in the container
      BCL_ExampleCheck( coefficient_handler.GetLabel().ToString(), "\"\"");
      BCL_ExampleCheck( coefficient_2d_handler.GetLabel().ToString(), "\"\"");

      std::ostringstream error_from_coefficient_handler;
      // set coefficient_handler with a legitimate vector of size_ts inside its range
      coefficient_handler.TryRead
      (
        util::ObjectDataLabel( test_coefficient_string),
        error_from_coefficient_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_coefficient_handler.str(),
        "",
        "coefficient_handler.TryRead( util::ObjectDataLabel( test_coefficient_string), error_from_coefficient_handler)"
      );
      // check that set parameter succeeded
      BCL_ExampleIndirectCheck
      (
        coefficients,
        storage::Vector< size_t>::Create( 0, 9, 49),
        "TryRead (valid case)"
      );
      error_from_coefficient_handler.str( "");
      // now try a double outside the range
      coefficient_handler.TryRead( util::ObjectDataLabel( "(51)"), error_from_coefficient_handler);
      BCL_ExampleIndirectCheck
      (
        error_from_coefficient_handler.str().size() > 0,
        true,
        "coefficient_handler.TryRead( util::ObjectDataLabel( 51), error_from_coefficient_handler)"
      );
      // check that set parameter failed
      BCL_ExampleIndirectCheck
      (
        coefficients,
        storage::Vector< size_t>::Create( 0, 9, 49),
        "TryRead (invalid case)"
      );
      error_from_coefficient_handler.str( "");
      // now try something that isn't even a size_t
      coefficient_handler.TryRead( util::ObjectDataLabel( "(alpha)"), error_from_coefficient_handler);
      BCL_ExampleIndirectCheck
      (
        error_from_coefficient_handler.str().size() > 0,
        true,
        "coefficient_handler.TryRead( util::ObjectDataLabel( \"alpha\"), error_from_coefficient_handler"
      );

      error_from_coefficient_handler.str( "");
      // check that set parameter failed
      BCL_ExampleIndirectCheck
      (
        coefficients,
        storage::Vector< size_t>::Create( 0, 9, 49),
        "TryRead (invalid case)"
      );

      // TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM) for 2d vectors
      // (any depth of vectors is possible; this just illustrates that with 2d vectors for simplicity)
      std::ostringstream error_from_coefficients_2d_handler;
      // set coefficient_2d_handler with a legitimate 2d vector with #s inside its range
      coefficient_2d_handler.TryRead
      (
        util::ObjectDataLabel( test_coefficients_2d_string),
        error_from_coefficients_2d_handler
      );
      BCL_ExampleIndirectCheck
      (
        error_from_coefficients_2d_handler.str(),
        "",
        "coefficient_2d_handler.TryRead( "
        + test_coefficients_2d_string + ", error_from_coefficients_2d_handler)"
      );

      // create the expected vector
      storage::Vector< storage::Vector< size_t> > expected_coefficients_2d
      (
        storage::Vector< storage::Vector< size_t> >::Create
        (
          storage::Vector< size_t>::Create( 0, 9, 123),
          storage::Vector< size_t>::Create( 16, 12)
        )
      );

      // check that set parameter succeeded
      BCL_ExampleIndirectCheck
      (
        coefficients_2d,
        expected_coefficients_2d,
        "TryRead 2d vector (valid case)"
      );

      // test the get string functionality; which converts the handled object back into a string

      // label handlers made with this constructor are commonly used in parameterized objects
      BCL_ExampleCheck( coefficient_handler.GetLabel().ToString(), "(0,9,49)");
      BCL_ExampleCheck( coefficient_2d_handler.GetLabel().ToString(), "((0,9,123),(16,12))");

      // construct a label handler that will set the test ranges vector
      util::OwnPtr< io::SerializationBase< storage::Vector< math::Range< size_t> > > > ranges_handler
      (
        io::Serialization::GetAgent( &ranges)
      );

      // Serialization handlers for containers should output a null value
      // if there is nothing in the container
      BCL_ExampleCheck( ranges_handler->GetLabel().ToString(), "\"\"");

      std::ostringstream error_from_handler;
      // set ranges with a legitimate vector of ranges
      BCL_ExampleCheck
      (
        ranges_handler->TryRead( util::ObjectDataLabel( test_ranges_string), error_from_handler),
        true
      );
      // check that set parameter succeeded
      BCL_ExampleIndirectCheck
      (
        ranges,
        storage::Vector< math::Range< size_t> >::Create
        (
          math::Range< size_t>( 0, 0),
          math::Range< size_t>( math::RangeBorders::e_LeftClosed, 8, 9, math::RangeBorders::e_RightOpen),
          math::Range< size_t>( math::RangeBorders::e_LeftOpen, 49, 57, math::RangeBorders::e_RightOpen)
        ),
        "TryRead (valid case)"
      );

      // try GetLabel
      BCL_ExampleCheck( ranges_handler->GetLabel().ToString(), "(\"[0,0]\",\"[8,9)\",\"(49,57)\")");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationContainer

  const ExampleClass::EnumType ExampleIoSerializationContainer::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationContainer())
  );

} // namespace bcl
