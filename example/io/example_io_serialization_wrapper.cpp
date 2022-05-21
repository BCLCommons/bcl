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
#include "io/bcl_io_serialization_wrapper.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialization_wrapper.cpp
  //!
  //! @author mendenjl
  //! @date Nov 01, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerializationWrapper :
    public ExampleInterface
  {
  public:

    ExampleIoSerializationWrapper *Clone() const
    {
      return new ExampleIoSerializationWrapper( *this);
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
      math::Range< size_t> range;

      std::string test_range_string[ 3] = { "[0]", "[8,9)", "(49,57)"};
      math::Range< size_t> correct_ranges[ 3] =
      {
        math::Range< size_t>( 0, 0),
        math::Range< size_t>( math::RangeBorders::e_LeftClosed, 8, 9, math::RangeBorders::e_RightOpen),
        math::Range< size_t>( math::RangeBorders::e_LeftOpen, 49, 57, math::RangeBorders::e_RightOpen)
      };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a label handler that will set the test ranges vector
      io::SerializationWrapper< math::Range< size_t> > range_handler;

      std::ostringstream error_from_handler;
      for( int i( 0); i < 3; ++i)
      {
        // set ranges with a legitimate range string
        BCL_ExampleCheck
        (
          range_handler.TryReadObject( range, util::ObjectDataLabel( "", test_range_string[ i]), error_from_handler),
          true
        );
        // check that set parameter succeeded
        BCL_ExampleIndirectCheck( range, correct_ranges[ i], "TryRead (valid case)");
        BCL_ExampleIndirectCheck( range.GetLabel(), correct_ranges[ i].GetLabel(), "TryRead (valid case)");
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerializationWrapper

  const ExampleClass::EnumType ExampleIoSerializationWrapper::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerializationWrapper())
  );

} // namespace bcl
