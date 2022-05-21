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
#include "io/bcl_io_serialization_with_check.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_file_existence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_with_check.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationWithCheck :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationWithCheck *Clone() const
    {
      return new ExampleIoSerializationWithCheck( *this);
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
      std::string test_string( "test");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // SerializationWithCheck< std::string>( const command::ParameterCheckInterface &PARAMETER_CHECK, const std::string *DATA)
      io::SerializationWithCheck< std::string> serialization_handler_file
      (
        command::ParameterCheckFileExistence(),
        &test_string
      );

      // label handlers made with this constructor are typically only used inside containers
      // since they don't provide info about what is labeled
      BCL_ExampleCheck( serialization_handler_file.GetLabel().ToString(), test_string);

      // TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM) for strings
      std::ostringstream error_from_serialization_handler_file;
      // set serialization_handler_dbl with a legitimate double inside its range
      BCL_ExampleCheck
      (
        serialization_handler_file.TryRead( util::ObjectDataLabel( "0.5"), error_from_serialization_handler_file),
        false
      );
      BCL_ExampleIndirectCheck
      (
        error_from_serialization_handler_file.str().size() > 0,
        true,
        "serialization_handler_file.TryRead( util::ObjectDataLabel( 0.5), error_from_serialization_handler_file)"
      );

      // now write a file
      WriteBCLObject( serialization_handler_file);
      // get the location of that file
      const std::string valid_filename( GetExampleOutputPathForBCLObject( serialization_handler_file));

      error_from_serialization_handler_file.str( "");
      // now try that file with the label handler
      serialization_handler_file.TryRead( util::ObjectDataLabel( valid_filename), error_from_serialization_handler_file);
      BCL_ExampleIndirectCheck
      (
        error_from_serialization_handler_file.str(),
        "",
        "serialization_handler_file.TryRead( util::ObjectDataLabel( valid_filename), error_from_serialization_handler_file)"
      );
      // check that set parameter failed
      BCL_ExampleIndirectCheck( test_string, valid_filename, "TryRead file (valid case)");

      // make another double and string for a fully parameterized label handler
      std::string filename( "unknown");
      io::SerializationWithCheck< std::string> serialization_handler_file_param
      (
        command::ParameterCheckFileExistence(),
        &filename
      );

      // label handlers made with this constructor are commonly used in parameterized objects
      BCL_ExampleCheck( serialization_handler_file_param.GetLabel().ToString(), "unknown");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationWithCheck

  const ExampleClass::EnumType ExampleIoSerializationWithCheck::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationWithCheck())
  );

} // namespace bcl
