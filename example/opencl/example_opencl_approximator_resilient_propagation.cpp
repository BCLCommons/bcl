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
#include "opencl/bcl_opencl_approximator_resilient_propagation.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_approximator_resilient_propagation.cpp
  //!
  //! @author mendenjl
  //! @date Sep 04, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclApproximatorResilientPropagation :
    public ExampleInterface
  {
  public:

    ExampleOpenclApproximatorResilientPropagation *Clone() const
    {
      return new ExampleOpenclApproximatorResilientPropagation( *this);
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
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      static const size_t s_number_values( 4);

      //initializing
      storage::VectorND< 2, linal::Vector< float> > features_results_setup[ s_number_values] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0, 0.0), linal::MakeVector< float>( 0.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0, 1.0), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 1.0, 0.0), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 1.0, 1.0), linal::MakeVector< float>( 0.0))
      };

      // set up the data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < 4; ++count)
      {
        features_results->PushBack( features_results_setup[ count]);
      }

      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( *features_results)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      storage::Vector< size_t> architecture( 1, size_t( 3));

      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function
      (
        new model::ObjectiveFunctionWrapper( features_results_frds, util::ObjectDataLabel( "RMSD"))
      );

      // construct from properties
      util::ShPtr< model::ApproximatorBase> sp_iterate
      (
        new opencl::ApproximatorResilientPropagation
        (
          features_results_frds,
          architecture,
          25,                            // steps per call
          opencl::GetTools().GetFirstCommandQueue()
        )
      );

      sp_objective_function->SetData( features_results_frds, sp_iterate->GetRescaleFeatureDataSet());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_ExampleCheck( sp_iterate->GetClassIdentifier(), "bcl::opencl::ApproximatorResilientPropagation");
      sp_iterate->SetObjectiveFunction( sp_objective_function);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // iterate (train) the network one step
      storage::Pair< util::ShPtr< model::Interface>, float>::s_Instance.IsDefined();
      sp_iterate->Next();
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        improved_pair( sp_iterate->GetCurrentApproximation());
      model::FeatureDataSet< float> pred( improved_pair->First()->operator()( *features_results_frds->GetFeaturesPtr()));
      model::FeatureDataSet< float> exp( features_results_frds->GetResultsPtr()->GetMatrix());

      // check the output of the untrained network
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "NN predicts \n"
//        + util::Format()( pred) +
//        " for \n" + util::Format()( exp) + " before training!"
//      );

      util::ShPtr< model::ObjectiveFunctionWrapper> sp_approx_objective_function
      (
        new model::ObjectiveFunctionWrapper( features_results_frds, util::ObjectDataLabel( "RMSD"))
      );

      sp_approx_objective_function->SetData( features_results_frds, sp_iterate->GetRescaleFeatureDataSet());

      // terminate after 1 iterations
      sp_iterate->SetCriterion( opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float>( 4));

      util::Stopwatch timer;
      // call minimizer till criterion met
      sp_iterate->Approximate();
      const storage::Pair< util::ShPtr< model::Interface>, float> result( *sp_iterate->GetCurrentApproximation());
      BCL_MessageStd( "Time to predict using gpu: " + util::Format()( timer.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds()));

      // check the output of the trained network
      linal::Vector< float> pred_result( features_results_frds->GetResultSize());
      pred_result = ( result.First()->operator ()( *features_results_frds->GetFeaturesPtr())( 0));
      linal::Vector< float> actual_result( features_results_frds->GetResultSize());
      actual_result = ( features_results_frds->GetResults().DeScale()( 0));

//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "NN predicts \n"
//        + util::Format()( result.First()->operator ()( *features_results_frds->GetFeaturesPtr()))
//        + " for \n" + util::Format()( *features_results_frds->GetResultsPtr())
//      );

      BCL_Example_Check
      (
        ( pred_result - actual_result).Norm() < 0.05,
        "The deviation between the results and ( 0.0) should be smaller than 0.05 but is: "
        + util::Format()
        (
          ( pred_result - actual_result).Norm()
        )
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclApproximatorResilientPropagation

  const ExampleClass::EnumType ExampleOpenclApproximatorResilientPropagation::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclApproximatorResilientPropagation())
  );

} // namespace bcl
