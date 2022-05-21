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
#include "mc/bcl_mc_temperature_accepted.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_tracker.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_temperature_accepted.cpp
  //!
  //! @author mendenjl
  //! @date Nov 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMcTemperatureAccepted :
    public ExampleInterface
  {
  public:

    ExampleMcTemperatureAccepted *Clone() const
    {
      return new ExampleMcTemperatureAccepted( *this);
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
      mc::TemperatureAccepted default_constructor;

      // test parameterized constructor
      const double start_fraction( 0.5);
      const double end_fraction( 0.25);
      const double start_temp( 100.0);
      const size_t iterations( 30);
      mc::TemperatureAccepted param_constructor( start_fraction, end_fraction, iterations, start_temp, 1);

      // construct from serialization
      const std::string parameters( param_constructor.GetString( true));
      mc::TemperatureAccepted serialized_temp;
      serialized_temp.AssertRead( parameters);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      util::ShPtr< storage::Pair< double, double> > arbitrary_result( new storage::Pair< double, double>( 50.0, 30.0));

      BCL_ExampleCheck( param_constructor.GetLastCalculatedTemperature(), start_temp);

      storage::Vector< storage::Triplet< double, opti::StepStatusEnum, double> > test_step_statuses_temps
      (
        storage::Vector< storage::Triplet< double, opti::StepStatusEnum, double> >::Create
        (
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 0.0), opti::e_Improved, 70.0),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 1.0), opti::e_Improved, 48.6358),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 2.0), opti::e_Accepted, 33.6515),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 3.0), opti::e_Rejected, 23.495),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 4.0), opti::e_Rejected, 16.774),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 5.0), opti::e_Rejected, 12.1569),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 6.0), opti::e_Improved, 8.63073),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 7.0), opti::e_Improved, 5.75716),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 8.0), opti::e_Accepted, 3.60181),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 9.0), opti::e_Accepted, 2.3429),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 10.0), opti::e_Accepted, 1.87549),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 11.0), opti::e_Accepted, 1.74565),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 12.0), opti::e_Accepted, 1.45418),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 13.0), opti::e_Accepted, 0.836078),
          storage::Triplet< double, opti::StepStatusEnum, double>( -sin( 14.0), opti::e_Accepted, 0.156511)
        )
      );

      // create a tracker to update get temperature with
      opti::Tracker< double, double> tracker_double( opti::e_SmallerIsBetter);
      for
      (
        auto itr( test_step_statuses_temps.Begin()), itr_end( test_step_statuses_temps.End());
        itr != itr_end;
        ++itr
      )
      {
        tracker_double.Track( arbitrary_result, itr->Second());
        param_constructor.TrackDelta( itr->First());
        BCL_ExampleIndirectCheckWithinAbsTolerance
        (
          param_constructor.GetTemperature( tracker_double),
          itr->Third(),
          0.1,
          "Temperature after "
          + util::Format()( tracker_double.GetCounts()( opti::e_Rejected)) + " rejected, "
          + util::Format()( tracker_double.GetCounts()( opti::e_Improved)) + " improved, "
          + util::Format()( tracker_double.GetCounts()( opti::e_Accepted)) + " accepted steps "
        );
      }

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
        TestBCLObjectIOForSymmetry( param_constructor, mc::TemperatureAccepted()),
        true,
        "TemperatureAccepted I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMcTemperatureAccepted

  const ExampleClass::EnumType ExampleMcTemperatureAccepted::s_Instance
  (
    GetExamples().AddEnum( ExampleMcTemperatureAccepted())
  );

} // namespace bcl
