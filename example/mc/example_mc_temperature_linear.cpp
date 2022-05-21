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
#include "mc/bcl_mc_temperature_linear.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_tracker.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_temperature_linear.cpp
  //!
  //! @author mendenjl
  //! @date Nov 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMcTemperatureLinear :
    public ExampleInterface
  {
  public:

    ExampleMcTemperatureLinear *Clone() const
    {
      return new ExampleMcTemperatureLinear( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      mc::TemperatureLinear default_constructor;

      // test parameterized constructor
      const double start_temp( 101.0);
      const double end_temp( 1.0);
      const size_t iterations( 25);
      mc::TemperatureLinear param_constructor( start_temp, end_temp, iterations);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck( param_constructor.GetLastCalculatedTemperature(), start_temp);

      // create a tracker to update get temperature with
      opti::Tracker< double, double> tracker_double( opti::e_SmallerIsBetter);
      tracker_double.Track
      (
        util::ShPtr< storage::Pair< double, double> >( new storage::Pair< double, double>( 50.0, 50.0)),
        opti::e_Improved
      );
      tracker_double.Track
      (
        util::ShPtr< storage::Pair< double, double> >( new storage::Pair< double, double>( 50.0, 30.0)),
        opti::e_Improved
      );

      // test get temperature
      BCL_ExampleCheckWithinAbsTolerance
      (
        param_constructor.GetTemperature( tracker_double),
        93.0,
        0.0001
      );
      BCL_ExampleIndirectCheck
      (
        param_constructor.GetTemperature( tracker_double),
        param_constructor.GetTemperature( tracker_double),
        "GetTemperature should return a constant value for a given tracker"
      );

      param_constructor.Reset();
      BCL_ExampleIndirectCheck( param_constructor.GetLastCalculatedTemperature(), start_temp, "Reset");

      // construct using serializer
      const std::string parameters( param_constructor.GetString( true));
      mc::TemperatureLinear serialized_temp;
      serialized_temp.AssertRead( parameters);
      BCL_ExampleCheck
      (
        serialized_temp.GetTemperature( tracker_double), param_constructor.GetTemperature( tracker_double)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( param_constructor, mc::TemperatureLinear()),
        true,
        "TemperatureLinear I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMcTemperatureLinear

  const ExampleClass::EnumType ExampleMcTemperatureLinear::s_Instance
  (
    GetExamples().AddEnum( ExampleMcTemperatureLinear())
  );

} // namespace bcl
