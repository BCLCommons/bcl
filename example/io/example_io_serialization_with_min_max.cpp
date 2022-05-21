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
#include "io/bcl_io_serialization_with_min_max.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_with_min_max.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationWithMinMax :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationWithMinMax *Clone() const
    {
      return new ExampleIoSerializationWithMinMax( *this);
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

      // SerializationWithMinMax( const command::ParameterCheckInterface &PARAMETER_CHECK, const std::string *DATA)
      io::SerializationWithMinMax< double> serialization_handler_dbl( &test_double, 0.0, 1.0);

      // label handlers made with this constructor are typically only used inside containers
      // since they don't provide info about what is labeled
      BCL_ExampleCheck( serialization_handler_dbl.GetLabel().ToString(), "5");

      std::ostringstream error_from_dbl;
      // set serialization_handler_dbl with a legitimate double inside its range
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "0.5"), error_from_dbl), true);
      BCL_ExampleIndirectCheck
      (
        error_from_dbl.str(),
        "",
        "serialization_handler_dbl.TryRead( util::ObjectDataLabel( 0.5), error_from_serialization_handler_dbl)"
      );
      // check that set parameter succeeded
      BCL_ExampleIndirectCheck( test_double, double( 0.5), "TryRead (valid case)");
      error_from_dbl.str( "");
      // now try a double outside the range
      BCL_ExampleCheck( serialization_handler_dbl.TryRead( util::ObjectDataLabel( "1.5"), error_from_dbl), false);
      BCL_ExampleIndirectCheck
      (
        error_from_dbl.str().size() > 0,
        true,
        "serialization_handler_dbl.TryRead( util::ObjectDataLabel( 1.5), error_from_serialization_handler_dbl)"
      );
      error_from_dbl.str( "");
      // now try something that isn't even a double
      BCL_ExampleCheck
      (
        serialization_handler_dbl.TryRead( util::ObjectDataLabel( "test"), error_from_dbl),
        false
      );
      BCL_ExampleIndirectCheck
      (
        error_from_dbl.str().size() > 0,
        true,
        "serialization_handler_dbl.TryRead( util::ObjectDataLabel( test), error_from_serialization_handler_dbl)"
      );
      error_from_dbl.str( "");

      // make another double and string for a fully parameterized label handler
      double my_val( 3.14);

      // SerializationWithMinMax( const command::ParameterCheckInterface &PARAMETER_CHECK, const t_DataType *DATA)
      io::SerializationWithMinMax< double> serialization_handler_dbl_param( &my_val, 0.0, 10.0);

      // label handlers made with this constructor are commonly used in parameterized objects
      BCL_ExampleCheck( serialization_handler_dbl_param.GetLabel().ToString(), "3.14");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationWithMinMax

  const ExampleClass::EnumType ExampleIoSerializationWithMinMax::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationWithMinMax())
  );

} // namespace bcl
