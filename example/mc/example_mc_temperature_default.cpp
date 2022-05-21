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
#include "mc/bcl_mc_temperature_default.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_tracker.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_temperature_default.cpp
  //!
  //! @author mendenjl
  //! @date Nov 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMcTemperatureDefault :
    public ExampleInterface
  {
  public:

    ExampleMcTemperatureDefault *Clone() const
    {
      return new ExampleMcTemperatureDefault( *this);
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
      mc::TemperatureDefault default_constructor;

      // test parameterized constructor
      const double start_temp( 100.0);
      mc::TemperatureDefault param_constructor( start_temp);

      // construct from parameters
      const std::string parameters( param_constructor.GetLabel().ToString());
      mc::TemperatureDefault param_temp;
      param_temp.AssertRead( parameters);
      BCL_ExampleCheck
      (
        param_constructor.GetLastCalculatedTemperature(), param_temp.GetLastCalculatedTemperature()
      );

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
      BCL_ExampleCheck( param_constructor.GetTemperature( tracker_double), start_temp);
      BCL_ExampleIndirectCheck
      (
        param_constructor.GetTemperature( tracker_double),
        param_constructor.GetTemperature( tracker_double),
        "GetTemperature should return a constant value for a given tracker"
      );

      param_constructor.Reset();
      BCL_ExampleIndirectCheck( param_constructor.GetLastCalculatedTemperature(), start_temp, "Reset");

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( param_constructor, mc::TemperatureDefault()),
        true,
        "TemperatureLinear I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMcTemperatureDefault

  const ExampleClass::EnumType ExampleMcTemperatureDefault::s_Instance
  (
    GetExamples().AddEnum( ExampleMcTemperatureDefault())
  );

} // namespace bcl
